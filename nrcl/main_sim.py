# Noise-resilient contrastive encoder plus a q-value detector, on the
# simulated data and the v1 CSV format. Also dumps the embeddings used by

import os
import sys
import argparse
import math
import numpy as np
import pandas as pd
from typing import Dict, Tuple, List, Optional

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import lightning as L
from torch.utils.data import DataLoader

import umap
import matplotlib.pyplot as plt

from sklearn.metrics import confusion_matrix
from scipy import stats
from statsmodels.stats.multitest import multipletests, fdrcorrection
from eval import LightGOFAE, get_novelty_info
from myutils import Database
from paths import EXPERIMENT_DATA, SIMULATED_DATA


# Encoder + Triplet training

class TSEncoder(nn.Module):
    def __init__(self, input_dim: int = 2, latent_dim: int = 64):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv1d(input_dim, 64, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Conv1d(64, 128, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.AdaptiveAvgPool1d(1),
        )
        self.proj = nn.Linear(128, latent_dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = self.conv(x).squeeze(-1)  # [B, 128]
        z = self.proj(h)              # [B, D]
        return z

class NRTriplet(L.LightningModule):
    """
    Per batch: take the centroid of the normal set, treat the clean_ratio
    fraction of noisy samples farthest from it as confident outliers, and
    train with a triplet hinge loss.
    """
    def __init__(self, latent_dim: int = 64, margin: float = 1.0, clean_ratio: float = 0.05, lr: float = 1e-3):
        super().__init__()
        self.save_hyperparameters()
        self.encoder = TSEncoder(latent_dim=latent_dim)
        self.margin = margin
        self.clean_ratio = clean_ratio
        self.lr = lr

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def training_step(self, batch, batch_idx):
        norm_x, noisy_x = batch
        z_norm = self(norm_x)
        z_noisy = self(noisy_x)

        center = z_norm.mean(dim=0)  # [D]
        dist = torch.norm(z_noisy - center, dim=1)

        Bn = z_noisy.size(0)
        K = max(1, int(Bn * self.clean_ratio))
        _, idx = torch.topk(dist, K)

        anchors = z_norm[torch.randperm(z_norm.size(0), device=z_norm.device)[:K]]
        positives = z_norm[torch.randperm(z_norm.size(0), device=z_norm.device)[:K]]
        negatives = z_noisy[idx]

        pos_d = F.pairwise_distance(anchors, positives)
        neg_d = F.pairwise_distance(anchors, negatives)

        loss = F.relu(pos_d - neg_d + self.margin).mean()

        self.log("train_loss", loss, prog_bar=True)
        self.log("avg_noisy_dist", dist.mean(), prog_bar=False)
        return loss

    def configure_optimizers(self):
        return optim.Adam(self.parameters(), lr=self.lr)

class ELRNoiseResilientContrastive(L.LightningModule):
    def __init__(self, num_noisy_samples, latent_dim=64, margin=1.0, clean_ratio=0.5, beta=0.7, lambd_elr=3.0):
        super().__init__()
        self.save_hyperparameters()
        self.encoder = TSEncoder(latent_dim=latent_dim)
        
        # ELR Hyperparameters
        self.beta = beta           # Momentum for memory update
        self.lambd_elr = lambd_elr # Strength of ELR regularization
        
        # ELR Memory Buffer: Stores historical similarity to the normal center
        # num_noisy_samples gives each noisy sequence its own slot.
        self.register_buffer("normality_memory", torch.zeros(num_noisy_samples))
        
        self.margin = margin
        self.clean_ratio = clean_ratio
        self.l1_lambda = 1e-4
        self.umap_weight = 1
        

    def forward(self, x):
        return self.encoder(x)

    def training_step(self, batch, batch_idx):
        # The DataLoader yields (norm_x, noisy_x, noisy_indices).
        norm_x, noisy_x, noisy_idx = batch
        
        norm_emb = self(norm_x)
        noisy_emb = self(noisy_x)

        norm_center = norm_emb.mean(dim=0)
        distances = torch.norm(noisy_emb - norm_center, dim=1)# eu distance
        num_to_keep = int(len(noisy_x) * self.clean_ratio)
        _, true_outlier_indices = torch.topk(distances, num_to_keep)

        # Triplets for the contrastive loss.
        anchors = norm_emb[torch.randperm(len(norm_emb))[:num_to_keep]]
        positives = norm_emb[torch.randperm(len(norm_emb))[:num_to_keep]]
        negatives = noisy_emb[true_outlier_indices]

        pos_dist = F.pairwise_distance(anchors, positives)
        neg_dist = F.pairwise_distance(anchors, negatives)
        contrastive_loss = F.relu(pos_dist - neg_dist + self.margin).mean()

        # Asymmetric ELR regularization.
        # Calculate current normality (similarity to the current normal centroid)
        norm_center = norm_emb.mean(dim=0, keepdim=True)
        # Cosine similarity: 1.0 means it looks perfectly normal
        current_sim = F.cosine_similarity(noisy_emb, norm_center)
        
        # Fetch historical normality from memory
        past_sim = self.normality_memory[noisy_idx]
        
        # Squared error between the current and remembered normality: pushing a
        # sample the memory still considers normal is penalized, which protects
        # mislabelled normals from the triplet loss.
        elr_reg = F.mse_loss(current_sim, past_sim)
        
        # Update the memory buffer (moving average).
        with torch.no_grad():
            new_sim = self.beta * past_sim + (1 - self.beta) * current_sim
            self.normality_memory[noisy_idx] = new_sim

        l1_norm = sum(p.abs().sum() for name, p in self.named_parameters() if 'weight' in name)
        
        # Only apply ELR after a small warm-up (e.g., epoch 5)
        lambd = self.lambd_elr if self.current_epoch > 5 else 0.1
        
        total_loss = contrastive_loss + (lambd * elr_reg) + (self.l1_lambda * l1_norm)

        self.log("train_loss", total_loss, prog_bar=True)
        self.log("elr_reg", elr_reg, prog_bar=False)
        
        return total_loss

    def configure_optimizers(self):
        return optim.Adam(self.parameters(), lr=1e-3)



@torch.no_grad
def get_embeddings(model: nn.Module, loader: DataLoader, device: str = "cuda") -> np.ndarray:
    model.eval().to(device)
    outs = []
    for batch in loader:
        x = batch[0].to(device)
        z = model(x).detach().cpu().numpy()
        outs.append(z)
    return np.concatenate(outs, axis=0)


# UMAP plotting

def umap_plot(z_all: np.ndarray, is_outlier: np.ndarray, save_path: str = "umap.png"):
    reducer = umap.UMAP(
        n_neighbors=15, min_dist=0.1, n_components=2,
        metric="euclidean", random_state=42
    )
    emb2 = reducer.fit_transform(z_all)

    plt.figure(figsize=(10, 8))
    plt.scatter(emb2[~is_outlier, 0], emb2[~is_outlier, 1], s=5, alpha=0.25, label="normal")
    plt.scatter(emb2[is_outlier, 0], emb2[is_outlier, 1], s=10, alpha=0.8, label="G4_nonb")
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"[UMAP] saved to {save_path}")


# Statistical detector: distance-feature fusion

def compute_distance_features(
    z_ref: np.ndarray, 
    z_query: np.ndarray, 
    knn_k: int = 50,
    ref_stats: Optional[Dict] = None
) -> Tuple[np.ndarray, Dict]:
    """
    Three distance features per query point:
    1. Centroid distance (Euclidean to the reference mean)
    2. Mahalanobis distance (in the embedding space)
    3. kNN mean distance (local density)

    Returns (features, ref_stats), where features is (N_query, 3) with columns
    [d_cent, d_mahal, d_knn] and ref_stats holds mu, inv_cov and the kNN index.
    """
    # Fit the reference statistics on z_ref unless they were passed in.
    if ref_stats is None:
        mu = z_ref.mean(axis=0, keepdims=True)
        
        cov = np.cov(z_ref.T) + 1e-6 * np.eye(z_ref.shape[1])
        inv_cov = np.linalg.pinv(cov)
        
        from sklearn.neighbors import NearestNeighbors
        nnbr = NearestNeighbors(n_neighbors=min(knn_k + 1, len(z_ref)), algorithm="auto").fit(z_ref)
        
        ref_stats = {
            "mu": mu,
            "inv_cov": inv_cov,
            "nnbr": nnbr
        }

    # Feature 1: Centroid Distance
    diff = z_query - ref_stats["mu"]
    feat_cent = np.linalg.norm(diff, axis=1)


    feat_mahal2 = np.einsum("ni,ij,nj->n", diff, ref_stats["inv_cov"], diff)
    feat_mahal = np.sqrt(np.maximum(feat_mahal2, 0))

    # Feature 3: KNN Mean Distance
    dists, _ = ref_stats["nnbr"].kneighbors(z_query)
    # On a self-query the first neighbour is the point itself (distance 0); at
    # inference z_query is disjoint from the index, so the mean is taken over
    # all k returned neighbours either way.
    feat_knn = dists.mean(axis=1)

    features = np.stack([feat_cent, feat_mahal, feat_knn], axis=1)
    
    return features, ref_stats

def fit_null_stats(z_bdna: np.ndarray, knn_k: int = 50) -> Dict[str, np.ndarray]:
    """
    Build empirical null distributions for each statistic on B-DNA only.

    Returns dict: stat_name -> null_values (array of length N_bdna)
    """
    # centroid distance
    mu = z_bdna.mean(axis=0, keepdims=True)
    centroid_dist = np.linalg.norm(z_bdna - mu, axis=1)

    # covariance, jittered so pinv stays well behaved
    cov = np.cov(z_bdna.T)
    cov = cov + 1e-6 * np.eye(cov.shape[0])
    inv_cov = np.linalg.pinv(cov)
    diff = z_bdna - mu
    mahal2 = np.einsum("ni,ij,nj->n", diff, inv_cov, diff)

    # knn mean distance (in embedding space)
    from sklearn.neighbors import NearestNeighbors
    nnbr = NearestNeighbors(n_neighbors=min(knn_k + 1, len(z_bdna)), algorithm="auto").fit(z_bdna)
    dists, _ = nnbr.kneighbors(z_bdna)
    # exclude self at index 0
    knn_mean = dists[:, 1:].mean(axis=1)

    return {
        "centroid_dist": centroid_dist,
        "mahalanobis2": mahal2,
        "knn_mean_dist": knn_mean,
    }

def fit_fusion_model(features_null: np.ndarray) -> Dict:
    """
    Fuse the three distance features into one statistic by taking a second
    Mahalanobis distance, this time in the 3-d distance-feature space.

    Takes features_null (N_bdna, 3); returns the fusion mean and inverse
    covariance together with the resulting null distribution.
    """
    mu_feat = features_null.mean(axis=0, keepdims=True)
    cov_feat = np.cov(features_null.T) + 1e-8 * np.eye(features_null.shape[1])  # jitter guards against a singular covariance
    inv_cov_feat = np.linalg.pinv(cov_feat)
    
    diff = features_null - mu_feat
    dist_fusion_null = np.einsum("ni,ij,nj->n", diff, inv_cov_feat, diff)
    
    return {
        "mu_feat": mu_feat,
        "inv_cov_feat": inv_cov_feat,
        "null_dist_distribution": dist_fusion_null
    }


def empirical_upper_tail_p(null_vals: np.ndarray, vals: np.ndarray) -> np.ndarray:
    """
    Upper-tail empirical p-values:
      p = (1 + #{null >= v}) / (n+1)
    """
    null_sorted = np.sort(null_vals)
    n = len(null_sorted)
    # count(null >= v) = n - left insertion point of v
    left_idx = np.searchsorted(null_sorted, vals, side="left")
    ge = n - left_idx
    p = (ge + 1.0) / (n + 1.0)
    return p


def detect_outliers_qvalue(
    z_bdna_null: np.ndarray,
    z_mixed: np.ndarray,
    alpha: float = 0.20,
    knn_k: int = 50,
) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
      1) Use multiple distances as a feature vector in a low-dimensional
         'statistic space' (centroid, Mahalanobis^2, kNN mean).
      2) Fit a multivariate null on B-DNA in this statistic space.
      3) Compute a single Mahalanobis distance D(x) for every sample.
      4) Use empirical upper-tail p-values of D(x) and apply BH-FDR.

    Returns:
      pred_outlier (0/1) for mixed (bdna+nonb),
      debug dict containing per-stat p-values, D_null, D_mixed,
      final p-values, q-values, threshold used.
    """
    # 1. Per-statistics for B-DNA (null)
    null_stats = fit_null_stats(z_bdna_null, knn_k=knn_k)
    # shape [N_null, 3]
    stat_null_mat = np.stack(
        [
            null_stats["centroid_dist"],
            null_stats["mahalanobis2"],
            null_stats["knn_mean_dist"],
        ],
        axis=1,
    )

    # 2. Per-statistics for mixed set (using BDNA as reference)
    # centroid distance in embedding space
    mu_embed = z_bdna_null.mean(axis=0, keepdims=True)
    centroid_dist_m = np.linalg.norm(z_mixed - mu_embed, axis=1)

    # Mahalanobis^2 in embedding space
    cov_embed = np.cov(z_bdna_null.T) + 1e-6 * np.eye(z_bdna_null.shape[1])
    inv_cov_embed = np.linalg.pinv(cov_embed)
    diff_embed = z_mixed - mu_embed
    mahal2_m = np.einsum("ni,ij,nj->n", diff_embed, inv_cov_embed, diff_embed)

    # kNN mean distance in embedding space (neighbors built on BDNA)
    from sklearn.neighbors import NearestNeighbors
    nnbr = NearestNeighbors(
        n_neighbors=min(knn_k + 1, len(z_bdna_null)), algorithm="auto"
    ).fit(z_bdna_null)
    dists_m, _ = nnbr.kneighbors(z_mixed)
    knn_mean_m = dists_m[:, 1:].mean(axis=1)

    # shape [N_mixed, 3]
    stat_mixed_mat = np.stack(
        [centroid_dist_m, mahal2_m, knn_mean_m],
        axis=1,
    )

    # 3. Optional: per-stat empirical p-values (for analysis only)
    # Diagnostics and ablation only; the FDR decision does not use these.
    p_each = {
        "centroid_dist": empirical_upper_tail_p(
            null_stats["centroid_dist"], centroid_dist_m
        ),
        "mahalanobis2": empirical_upper_tail_p(
            null_stats["mahalanobis2"], mahal2_m
        ),
        "knn_mean_dist": empirical_upper_tail_p(
            null_stats["knn_mean_dist"], knn_mean_m
        ),
    }

    # 4. Build a single unified statistic in 'stat space'
    # Fit multivariate normal on stat_null_mat
    mu_stat = stat_null_mat.mean(axis=0, keepdims=True)           # [1, 3]
    cov_stat = np.cov(stat_null_mat.T)
    cov_stat = cov_stat + 1e-6 * np.eye(cov_stat.shape[0])
    inv_cov_stat = np.linalg.pinv(cov_stat)

    # Mahalanobis^2 in statistic space
    diff_null_stat = stat_null_mat - mu_stat                      # [N_null, 3]
    D_null = np.einsum("ni,ij,nj->n", diff_null_stat, inv_cov_stat, diff_null_stat)

    diff_mixed_stat = stat_mixed_mat - mu_stat                    # [N_mixed, 3]
    D_mixed = np.einsum("ni,ij,nj->n", diff_mixed_stat, inv_cov_stat, diff_mixed_stat)

    # 5. Empirical p-values from D_null
    # Upper-tail: larger D(x) => more outlier-like
    p_final = empirical_upper_tail_p(D_null, D_mixed)

    # 6. BH-FDR on this single family of p-values
    reject, q_values, _, _ = multipletests(p_final, alpha=alpha, method="fdr_bh")

    # q < alpha
    thr = alpha
    pred = (q_values < thr).astype(int)

    debug = {
        # per-stat p-values (for ablation / plots)
        "p_each": p_each,
        # unified-statistic Mahalanobis^2 in stat-space
        "D_null": D_null,
        "D_mixed": D_mixed,
        # final p-values & q-values actually used for FDR decision
        "p_final": p_final,
        "q": q_values,
        "threshold": np.array([thr], dtype=float),
    }
    return pred, debug


def get_raw_embeddings(x):
    x = x.reshape(-1, 50 * 2)
    return x

def get_gofae_embeddings(loader):
    gofae = LightGOFAE()
    x, _, _ = get_novelty_info(gofae, loader)
    print(x.shape)
    return x
# Metrics helper

def stat_binary(y_true_outlier01: np.ndarray, pred_outlier01: np.ndarray):
    # confusion_matrix wants labels [0,1]
    tn, fp, fn, tp = confusion_matrix(y_true_outlier01, pred_outlier01, labels=[0,1]).ravel().tolist()
    # tp = outlier detected
    print(f"TP(outlier): {tp}, FP: {fp}, FN: {fn}, TN(normal): {tn}")
    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-12)
    print(f"Precision: {precision:.4f}, Recall: {recall:.4f}, F1: {f1:.4f}")
    return precision, recall, f1


# Main

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--nonb_type', type=str, default='G_Quadruplex_Motif',
                        help='Non-B DNA structure type')
    parser.add_argument('--alpha', type=float, default=0.2,
                        help='FDR significance level')
    args = parser.parse_args()

    database = Database(
                 data_path=SIMULATED_DATA,
                 data_type='simulation',
                 win_size=50,
                 nonb_ratio=0.10,
                 nonb_type=args.nonb_type,
                 poison_train=False,
                 poison_val=True,
                 use_label=True)

    # Train the encoder
    model = ELRNoiseResilientContrastive(
                num_noisy_samples=len(database.train_nonb),
                latent_dim=64,
                clean_ratio=0.05,
                beta=0.8,        # momentum for the ELR memory buffer
                lambd_elr=10.0,  # strength of the noise protection
            )
    trainer = L.Trainer(max_epochs=10, accelerator="auto", strategy="ddp_find_unused_parameters_true")
    trainer.fit(model, database.get_contrastive_dataloader_elr())

    # Embeddings
    z_bdna_test = get_embeddings(model, database.test_loader_bdna)   # bdna-only set
    z_mixed = get_embeddings(model, database.test_loader_nonb)       # mixed outlier test loader

    # test raw features
    z_bdna_test_raw = get_raw_embeddings(database.test_bdna)
    z_mixed_raw = get_raw_embeddings(database.test_nonb)

    # test gofae embed
    z_bdna_test_gofae = get_gofae_embeddings(database.test_loader_bdna)
    z_mixed_gofae = get_gofae_embeddings(database.test_loader_nonb)
    # Build true labels for simulation (optional)
    # nonb dna: 1 bdna:0
    testset = database.test_sim # total dataset
    nonb_type = database.nonb_type # pick up the nonb data
    testset = testset[testset['label'].isin(nonb_type)]
    true_labels = testset['true_label'].tolist() 
    y_true = np.array([1 if (_ != 'bdna') else 0 for _ in true_labels])   # nonb: 1 bdna:0
    
    np.savez('embedding_plot_paper_CL.npz', embed=z_mixed, label=y_true)
    np.savez('embedding_plot_paper_raw.npz', embed=z_mixed_raw, label=y_true)
    np.savez('embedding_plot_paper_gofae.npz', embed=z_mixed_gofae, label=y_true)


    # Unified-statistic q-value detector
    fdr_alpha = args.alpha
    pred, dbg = detect_outliers_qvalue(
        z_bdna_null=z_bdna_test,
        z_mixed=z_mixed,
        alpha=fdr_alpha,
        knn_k=50,
    )

    print("q-value detector (Mahalanobis in stat space)")
    med_each = {k: float(np.median(v)) for k, v in dbg["p_each"].items()}
    print("median per-stat p:", med_each)
    print("median p_final:", float(np.median(dbg["p_final"])),
          "median q:", float(np.median(dbg["q"])))
    print("threshold used:", dbg["threshold"])

    # Save FDR detection results
    nonb_type_str = database.nonb_type[0] if isinstance(database.nonb_type, list) else database.nonb_type
    output_dir = os.path.join('pu_results', nonb_type_str)
    os.makedirs(output_dir, exist_ok=True)

    # Rejected sample indices (indices into test_nonb where pred==1)
    rejected_indices = np.where(pred == 1)[0]

    # Save rejected sample indices
    pd.DataFrame(rejected_indices, columns=['rejected_idx']).to_csv(
        os.path.join(output_dir, nonb_type_str + '_test_rejected_idx.csv'), index=False)

    # For experimental data: filter BED file to keep only FDR-rejected regions
    if database.data_type != 'simulation':
        bed_path = os.path.join(database.data_path, nonb_type_str + '_test.bed')
        if os.path.exists(bed_path):
            test_bed = pd.read_csv(bed_path, sep='\t')
            test_bed_filtered = test_bed.iloc[rejected_indices]
            test_bed_filtered.to_csv(
                os.path.join(output_dir, nonb_type_str + '_test_filtered.bed'),
                sep='\t', index=False)
            print(f"[SAVE] Filtered BED: {len(test_bed_filtered)} / {len(test_bed)} regions kept")

    # Save test summary
    total = len(database.test_nonb)
    discoveries = int(pred.sum())
    test_summary = pd.DataFrame({
        'total': [total],
        'discoveries': [discoveries],
        'fdr_level': [fdr_alpha]
    })
    test_summary.to_csv(os.path.join(output_dir, 'test_summary_upper.csv'), index=False)

    print(f"[SAVE] Summary: total={total}, discoveries={discoveries}, fdr_level={fdr_alpha}")
    print(f"[SAVE] All results saved to {output_dir}")