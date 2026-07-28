import numpy as np
import os
import pandas as pd
import sys
sys.path.append('../../')
import time
import torch
import torch.nn.functional as F
import torch.optim as optim
from pathlib import Path
from torch.utils.data import DataLoader
from torch.utils.data import TensorDataset
from myutils import Database

from sklearn.covariance import MinCovDet
from statsmodels.stats.multitest import fdrcorrection


import lightning as L
from architecture.Encoders import Encoder_P1, Encoder_P2
from architecture.Decoders import Decoder
from utilities.IOUtilities import restore_model
from core.Evaluation import get_novelty_info, NoveltyStats, get_union, get_intersection, get_symdif, \
    get_jaccard, get_idx, disc_ratio, compute_empirical, compute_true_sim_counts, compute_accuracy_metrics



class MODEL_CONFIG:
    path = './results_newdata_with_stiefel_tr/sim_GOFAEDND/sw/nonb_ratio_0.1/config_1/use_min/'
    n_z = 64
    n_proj = 64
    h_dim = 32
    stiefel = True
    loss_fn = 'mae'


class LightGOFAE(L.LightningModule):
    def __init__(self):
        super().__init__()
        self.encoder_p1 = Encoder_P1(n_inputs=2, out_dim=MODEL_CONFIG.n_z, h_dim=MODEL_CONFIG.h_dim, Stiefel=MODEL_CONFIG.stiefel)
        if MODEL_CONFIG.stiefel:
            self.encoder_p2 = Encoder_P2(in_dim=MODEL_CONFIG.h_dim * 24, z_dim=MODEL_CONFIG.n_z)
        self.decoder = Decoder(n_inputs=MODEL_CONFIG.n_z, n_outputs=2, h_dim=MODEL_CONFIG.h_dim)
        restore_model(MODEL_CONFIG.path, self.encoder_p1, self.encoder_p2, self.decoder, MODEL_CONFIG.stiefel, device='cuda')
        
    def forward(self, x):
        y = self.encoder_p1(x)
        if MODEL_CONFIG.stiefel:
            z = self.encoder_p2(y)
        y = self.decoder(z)
        return z, y 
    

def get_novelty_info(model, loader):
    loss_fn_eval = MODEL_CONFIG.loss_fn

    with torch.no_grad():
        store_recon_loss = []
        store_orig = []
        store_recon = []

        store_z_codes = np.empty((0, MODEL_CONFIG.n_z), float)
        for batch_idx, batch in enumerate(loader):
            x = batch[0]
            x = x.to('cuda')
            x_embed, x_recon = model(x)

            # calculate reconstruction loss
            if loss_fn_eval == 'mse':
                recon_resid = F.mse_loss(x_recon, x, reduction='none')
                recon_loss = torch.sum(recon_resid, dim=1)
            elif loss_fn_eval == 'mae':
                recon_resid = F.l1_loss(x_recon, x, reduction='none')
                recon_loss = torch.sum(recon_resid, dim=(1, 2))
            elif loss_fn_eval == 'smae':
                recon_resid = F.smooth_l1_loss(x_recon, x, reduction='none')
                recon_loss = torch.sum(recon_resid, dim=1)
            elif loss_fn_eval == 'rmse':
                recon_loss = torch.sqrt(torch.sum(F.mse_loss(x_recon, x, reduction='none'), dim=1))
            elif loss_fn_eval == 'max':
                recon_loss = torch.amax(torch.abs(x_recon - x), dim=1)
            else:
                sys.exit("Reconstruction Loss Not Defined")

            # Store the original input
            store_orig.extend(x.detach().cpu().numpy())

            # Store the reconstruction
            store_recon.extend(x_recon.detach().cpu().numpy())

            # Store the total reconstruction loss per observation
            store_recon_loss.extend(recon_loss.detach().cpu().numpy())

            # Store the encodings
            store_z_codes = np.append(store_z_codes, x_embed.detach().cpu().numpy(), axis=0)

    return store_z_codes, np.array(store_recon_loss), (np.array(store_orig), np.array(store_recon))


def fit_test_novelty_stats(dataloader, encoder_p1, encoder_p2, decoder):
    latent_code, rec, full = get_novelty_info((encoder_p1, encoder_p2, decoder))
    start_time_MCD1 = time.time()

    robust_cov_latent = MinCovDet().fit(latent_code)
    end_time_MCD1 = time.time()
    total_time_MCD1 = end_time_MCD1 - start_time_MCD1
    print('Total time for MCD1 computation on Test BDNA is {:.2f} minutes'.format(total_time_MCD1 / 60.),
            flush=True)

    test_MH = robust_cov_latent.mahalanobis(latent_code) ** (0.5)
    test_bivdat = np.hstack((test_MH.reshape(-1, 1), rec.reshape(-1, 1)))
    start_time_MCD2 = time.time()
    MH2 = MinCovDet().fit(test_bivdat)
    end_time_MCD2 = time.time()
    total_time_MCD2 = end_time_MCD2 - start_time_MCD2
    print('Total time for MCD2 computation on Test BDNA is {:.2f} minutes'.format(total_time_MCD2 / 60.),
            flush=True)
    MH2_test = MH2.mahalanobis(test_bivdat) ** (0.5)

    return latent_code, rec, full, test_MH, MH2_test

def transform_with_novelty_stats(dataloader, model):
    latent_code, rec, full = get_novelty_info(model, dataloader)
    robust_cov_latent = MinCovDet().fit(latent_code)

    MH = robust_cov_latent.mahalanobis(latent_code) ** (0.5)
    bivdat = np.hstack((MH.reshape(-1, 1), rec.reshape(-1, 1)))
    MH2 = MinCovDet().fit(bivdat)
    MH2_data = MH2.mahalanobis(bivdat) ** (0.5)
    return latent_code, rec, full, MH, MH2_data


def evaluate_acc(database, model):
    latent_code_nonb, rec_nonb, full_nonb, MH_nonb, MH2_data_nonb = transform_with_novelty_stats(database.test_loader_nonb, model)
    latent_code_b, rec_b, full_b, MH_b, MH2_data_b = transform_with_novelty_stats(database.test_loader_bdna, model)
    emp_dist_test_U, sorted_idx_test_U = compute_empirical(MH2_data_b, MH2_data_nonb, tail='upper')
    rejected_test_U, _ = fdrcorrection(emp_dist_test_U, alpha=0.03, method='indep', is_sorted=True)
    tn_idxs, fp_idxs, fn_idxs, tp_idxs, real_idx = compute_true_sim_counts(test_sim=database.test_sim, nonb_type='G_Quadruplex_Motif', \
                                                         idxs=sorted_idx_test_U[rejected_test_U])
    precision, recall, fscore = compute_accuracy_metrics(len(tn_idxs), len(fp_idxs), len(fn_idxs), len(tp_idxs))
    print(precision, recall, fscore)





if __name__ == "__main__":

    gofae = LightGOFAE()  
    # a = torch.rand([64, 2, 50]).to('cuda')
    # print(gofae(a)[0].shape, gofae(a)[1].shape)
    database = Database()

    evaluate_acc(database, gofae)