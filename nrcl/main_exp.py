"""
pu_contrastive_exp_v2.py
Adapted from pu_contrastive_exp.py to support the new G_Quadruplex_Motif
data format (with metadata columns: chr, strand, win_start50, win_end50,
motif_start, motif_end, features).

Reads either CSV format through data_loader.Database and writes the detections
with their motif coordinates, so no separate .bed file is needed downstream.
"""

import os
import sys
import argparse
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
from statsmodels.stats.multitest import multipletests

from data_loader import Database
from paths import EXPERIMENT_DATA


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
        h = self.conv(x).squeeze(-1)
        z = self.proj(h)
        return z


class NRTriplet(L.LightningModule):
    def __init__(self, latent_dim: int = 64, margin: float = 1.0,
                 clean_ratio: float = 0.05, lr: float = 1e-3):
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

        center = z_norm.mean(dim=0)
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
        return loss

    def configure_optimizers(self):
        return optim.Adam(self.parameters(), lr=self.lr)


class ELRNoiseResilientContrastive(L.LightningModule):
    def __init__(self, num_noisy_samples, latent_dim=64, margin=1.0,
                 clean_ratio=0.5, beta=0.7, lambd_elr=3.0):
        super().__init__()
        self.save_hyperparameters()
        self.encoder = TSEncoder(latent_dim=latent_dim)

        self.beta = beta
        self.lambd_elr = lambd_elr
        self.register_buffer("normality_memory", torch.zeros(num_noisy_samples))

        self.margin = margin
        self.clean_ratio = clean_ratio
        self.l1_lambda = 1e-4

    def forward(self, x):
        return self.encoder(x)

    def training_step(self, batch, batch_idx):
        norm_x, noisy_x, noisy_idx = batch
        norm_emb = self(norm_x)
        noisy_emb = self(noisy_x)

        norm_center = norm_emb.mean(dim=0)
        distances = torch.norm(noisy_emb - norm_center, dim=1)
        num_to_keep = int(len(noisy_x) * self.clean_ratio)
        _, true_outlier_indices = torch.topk(distances, num_to_keep)

        anchors = norm_emb[torch.randperm(len(norm_emb))[:num_to_keep]]
        positives = norm_emb[torch.randperm(len(norm_emb))[:num_to_keep]]
        negatives = noisy_emb[true_outlier_indices]

        pos_dist = F.pairwise_distance(anchors, positives)
        neg_dist = F.pairwise_distance(anchors, negatives)
        contrastive_loss = F.relu(pos_dist - neg_dist + self.margin).mean()

        # ELR regularization
        norm_center2 = norm_emb.mean(dim=0, keepdim=True)
        current_sim = F.cosine_similarity(noisy_emb, norm_center2)
        past_sim = self.normality_memory[noisy_idx]
        elr_reg = F.mse_loss(current_sim, past_sim)

        with torch.no_grad():
            new_sim = self.beta * past_sim + (1 - self.beta) * current_sim
            self.normality_memory[noisy_idx] = new_sim

        l1_norm = sum(p.abs().sum() for name, p in self.named_parameters() if 'weight' in name)
        lambd = self.lambd_elr if self.current_epoch > 5 else 0.1
        total_loss = contrastive_loss + (lambd * elr_reg) + (self.l1_lambda * l1_norm)

        self.log("train_loss", total_loss, prog_bar=True)
        self.log("elr_reg", elr_reg, prog_bar=False)
        return total_loss

    def configure_optimizers(self):
        return optim.Adam(self.parameters(), lr=1e-3)


# Embedding extraction

@torch.no_grad()
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
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2,
                         metric="euclidean", random_state=42)
    emb2 = reducer.fit_transform(z_all)
    plt.figure(figsize=(10, 8))
    plt.scatter(emb2[~is_outlier, 0], emb2[~is_outlier, 1], s=5, alpha=0.25, label="normal")
    plt.scatter(emb2[is_outlier, 0], emb2[is_outlier, 1], s=10, alpha=0.8, label="G4_nonb")
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"[UMAP] saved to {save_path}")


# Statistical Detector

def fit_null_stats(z_bdna: np.ndarray, knn_k: int = 50) -> Dict[str, np.ndarray]:
    mu = z_bdna.mean(axis=0, keepdims=True)
    centroid_dist = np.linalg.norm(z_bdna - mu, axis=1)

    cov = np.cov(z_bdna.T) + 1e-6 * np.eye(z_bdna.shape[1])
    inv_cov = np.linalg.pinv(cov)
    diff = z_bdna - mu
    mahal2 = np.einsum("ni,ij,nj->n", diff, inv_cov, diff)

    from sklearn.neighbors import NearestNeighbors
    nnbr = NearestNeighbors(n_neighbors=min(knn_k + 1, len(z_bdna)),
                            algorithm="auto").fit(z_bdna)
    dists, _ = nnbr.kneighbors(z_bdna)
    knn_mean = dists[:, 1:].mean(axis=1)

    return {
        "centroid_dist": centroid_dist,
        "mahalanobis2": mahal2,
        "knn_mean_dist": knn_mean,
    }


def empirical_upper_tail_p(null_vals: np.ndarray, vals: np.ndarray) -> np.ndarray:
    null_sorted = np.sort(null_vals)
    n = len(null_sorted)
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
    Feature-fusion detector:
      1) Compute centroid-dist, Mahalanobis^2, kNN-mean in embedding space
      2) Fit multivariate null in this 3D statistic space
      3) Mahalanobis distance in stat space -> empirical p-values -> BH-FDR
    """
    # null stats
    null_stats = fit_null_stats(z_bdna_null, knn_k=knn_k)
    stat_null_mat = np.stack([
        null_stats["centroid_dist"],
        null_stats["mahalanobis2"],
        null_stats["knn_mean_dist"],
    ], axis=1)

    # mixed stats
    mu_embed = z_bdna_null.mean(axis=0, keepdims=True)
    centroid_dist_m = np.linalg.norm(z_mixed - mu_embed, axis=1)

    cov_embed = np.cov(z_bdna_null.T) + 1e-6 * np.eye(z_bdna_null.shape[1])
    inv_cov_embed = np.linalg.pinv(cov_embed)
    diff_embed = z_mixed - mu_embed
    mahal2_m = np.einsum("ni,ij,nj->n", diff_embed, inv_cov_embed, diff_embed)

    from sklearn.neighbors import NearestNeighbors
    nnbr = NearestNeighbors(
        n_neighbors=min(knn_k + 1, len(z_bdna_null)), algorithm="auto"
    ).fit(z_bdna_null)
    dists_m, _ = nnbr.kneighbors(z_mixed)
    knn_mean_m = dists_m[:, 1:].mean(axis=1)

    stat_mixed_mat = np.stack([centroid_dist_m, mahal2_m, knn_mean_m], axis=1)

    # per-stat p (diagnostic)
    p_each = {
        "centroid_dist": empirical_upper_tail_p(null_stats["centroid_dist"], centroid_dist_m),
        "mahalanobis2": empirical_upper_tail_p(null_stats["mahalanobis2"], mahal2_m),
        "knn_mean_dist": empirical_upper_tail_p(null_stats["knn_mean_dist"], knn_mean_m),
    }

    # unified stat: Mahalanobis in stat-space
    mu_stat = stat_null_mat.mean(axis=0, keepdims=True)
    cov_stat = np.cov(stat_null_mat.T) + 1e-6 * np.eye(stat_null_mat.shape[1])
    inv_cov_stat = np.linalg.pinv(cov_stat)

    diff_null_stat = stat_null_mat - mu_stat
    D_null = np.einsum("ni,ij,nj->n", diff_null_stat, inv_cov_stat, diff_null_stat)

    diff_mixed_stat = stat_mixed_mat - mu_stat
    D_mixed = np.einsum("ni,ij,nj->n", diff_mixed_stat, inv_cov_stat, diff_mixed_stat)

    # empirical p -> BH-FDR
    p_final = empirical_upper_tail_p(D_null, D_mixed)
    reject, q_values, _, _ = multipletests(p_final, alpha=alpha, method="fdr_bh")

    thr = alpha
    pred = (q_values < thr).astype(int)

    debug = {
        "p_each": p_each,
        "D_null": D_null,
        "D_mixed": D_mixed,
        "p_final": p_final,
        "q": q_values,
        "threshold": np.array([thr], dtype=float),
    }
    return pred, debug


# Metrics helper

def stat_binary(y_true_outlier01: np.ndarray, pred_outlier01: np.ndarray):
    tn, fp, fn, tp = confusion_matrix(
        y_true_outlier01, pred_outlier01, labels=[0, 1]).ravel().tolist()
    print(f"TP(outlier): {tp}, FP: {fp}, FN: {fn}, TN(normal): {tn}")
    precision = tp / max(tp + fp, 1)
    recall = tp / max(tp + fn, 1)
    f1 = 2 * precision * recall / max(precision + recall, 1e-12)
    print(f"Precision: {precision:.4f}, Recall: {recall:.4f}, F1: {f1:.4f}")
    return precision, recall, f1


# embed
from eval import LightGOFAE, get_novelty_info
def get_raw_embeddings(x):
    x = x.reshape(-1, 50 * 2)
    return x

def get_gofae_embeddings(loader):
    gofae = LightGOFAE()
    x, _, _ = get_novelty_info(gofae, loader)
    print(x.shape)
    return x



# Main

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--nonb_type', type=str, default='G_Quadruplex_Motif')
    parser.add_argument('--alpha', type=float, default=0.2, help='FDR level')
    parser.add_argument('--data_path', type=str,
                        default=EXPERIMENT_DATA)
    parser.add_argument('--epochs', type=int, default=10)
    parser.add_argument('--clean_ratio', type=float, default=0.05)
    parser.add_argument('--beta', type=float, default=0.8)
    parser.add_argument('--lambd_elr', type=float, default=10.0)
    parser.add_argument('--smooth', type=str, default=None, choices=[None, 'median', 'mean'],
                        help='Smoothing method for translocation time signals (None/median/mean)')
    parser.add_argument('--smooth_kernel', type=int, default=5,
                        help='Kernel size for smoothing (must be odd)')
    args = parser.parse_args()

    # Load data (auto-detects old vs new format)
    database = Database(
        data_path=args.data_path,
        data_type='experimental',
        nonb_type=args.nonb_type,
        poison_train=False,
        poison_val=True,
        smooth=args.smooth,
        smooth_kernel=args.smooth_kernel,
    )

    # Train encoder
    model = ELRNoiseResilientContrastive(
        num_noisy_samples=len(database.train_nonb),
        latent_dim=64,
        clean_ratio=args.clean_ratio,
        beta=args.beta,
        lambd_elr=args.lambd_elr,
    )
    trainer = L.Trainer(
        max_epochs=args.epochs,
        accelerator="auto",
        strategy="ddp_find_unused_parameters_true",
    )
    trainer.fit(model, database.get_contrastive_dataloader_elr())

    # Get Embeddings: raw, gofae, ours
    z_bdna_test = get_embeddings(model, database.test_loader_bdna)
    z_mixed = get_embeddings(model, database.test_loader_nonb)


    # Detect
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

    # Save results as BED file with metadata
    nonb_type_str = database.nonb_type[0]
    output_dir = os.path.join('pu_results', nonb_type_str)
    os.makedirs(output_dir, exist_ok=True)

    rejected_indices = np.where(pred == 1)[0]
    total = len(pred)
    discoveries = int(pred.sum())

    meta_nonb = database.get_test_nonb_metadata()
    bed_path = os.path.join(output_dir, f'{nonb_type_str}_detected.bed')

    if meta_nonb is not None:
        # New format: metadata available — build BED from it
        meta_nonb = meta_nonb.reset_index(drop=True)
        rejected_meta = meta_nonb.iloc[rejected_indices].copy()

        # Parse sequence from features field:
        #   e.g. "...;sequence=gggttgggtggggaggg"
        def extract_sequence(feat_str):
            if pd.isna(feat_str):
                return ''
            for part in str(feat_str).split(';'):
                if part.startswith('sequence='):
                    return part.split('=', 1)[1]
            return ''

        bed_df = pd.DataFrame({
            'chr':          rejected_meta['chr'].values,
            'win_start':    rejected_meta['win_start50'].values,
            'win_end':      rejected_meta['win_end50'].values,
            'strand':       rejected_meta['strand'].values,
            'motif_start':  rejected_meta['motif_start'].values,
            'motif_end':    rejected_meta['motif_end'].values,
            'sequence':     rejected_meta['features'].apply(extract_sequence).values,
            'p_value':      dbg["p_final"][rejected_indices],
            'q_value':      dbg["q"][rejected_indices],
        })
        bed_df.to_csv(bed_path, sep='\t', index=False)
    else:
        # Old format: no metadata, save indices only
        bed_df = pd.DataFrame({
            'rejected_idx': rejected_indices,
            'p_value':      dbg["p_final"][rejected_indices],
            'q_value':      dbg["q"][rejected_indices],
        })
        bed_df.to_csv(bed_path, sep='\t', index=False)

    print(f"\nDetection summary: {nonb_type_str}")
    print(f"  Total test windows:  {total}")
    print(f"  Discoveries (FDR<{fdr_alpha}): {discoveries}")
    print(f"  Discovery rate:      {discoveries/max(total,1):.4f}")
    print(f"  Output BED:          {bed_path}")
