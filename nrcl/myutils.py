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
from torch.utils.data import TensorDataset, Dataset
from torch.optim.lr_scheduler import ReduceLROnPlateau
from core.Losses import discriminative_loss


from core.Dataset import import_data, reprocess_data2
from utils import load_data2
from paths import EXPERIMENT_DATA, SIMULATED_DATA

import lightning as L
from architecture.Encoders import Encoder_P1, Encoder_P2
from architecture.Decoders import Decoder

import warnings
warnings.filterwarnings("ignore")



# class LightAutoEncoder(L.LightningModule):
#     """
#     Pytorch imp for conv1d is (batch, n_features, seq_len)
#     """
#     def __init__(self, encoder, decoder, stiefel=False, encoder_p2=None):
#         super().__init__()
#         self.save_hyperparameters()
#         self.encoder_p1 = encoder
#         self.stiefel = stiefel
#         if stiefel:
#             self.encoder_p2 = encoder_p2
#         self.decoder = decoder
#         self.loss_fn = torch.nn.MSELoss()
        
#     def forward(self, x):
#         y = self.encoder_p1(x)
#         if self.stiefel:
#             y = self.encoder_p2(y)
#         y = self.decoder(y)
#         return y

#     def training_step(self, batch, batch_idx):
#         x = batch[0] if isinstance(batch, (list, tuple)) else batch
#         y = self(x)
#         loss = self.loss_fn(y, x)
#         self.log(f"loss", loss, on_step=True, on_epoch=True, prog_bar=True)
#         return loss 
    
#     def validation_step(self, batch, batch_idx):
#         pass

#     def test_step(self, batch, batch_idx):
#         pass


#     def configure_optimizers(self):
#         optimizer = torch.optim.AdamW(self.parameters(), lr=1e-3)
#         scheduler = torch.optim.lr_scheduler.OneCycleLR(
#             optimizer,
#             max_lr=1e-3,
#             total_steps=self.trainer.estimated_stepping_batches,
#             pct_start=0.3, 
#             div_factor=10, 
#             final_div_factor=100 
#         )
        
#         lr_scheduler_config = {
#             "scheduler": scheduler,
#             "interval": "step", 
#         }
#         return {"optimizer": optimizer, "lr_scheduler": lr_scheduler_config}


def reprocess_data(nonb_type, train, val, test, poison_train=False, poison_val=True):
    if not isinstance(nonb_type, list):
        nonb_type = [nonb_type]

    train_bdna = train[train["label"] == 'bdna'].drop(['label', 'true_label'], axis=1).to_numpy()
    val_bdna = val[val["label"] == 'bdna'].drop(['label', 'true_label'], axis=1).to_numpy()
    test_bdna = test[test["label"] == 'bdna'].drop(['label', 'true_label'], axis=1).to_numpy()
    train_nonb = train[train["label"].isin(nonb_type)].drop(['label', 'true_label'], axis=1).to_numpy()
    val_nonb = val[val["label"].isin(nonb_type)].drop(['label', 'true_label'], axis=1).to_numpy()
    test_nonb = test[test["label"].isin(nonb_type)].drop(['label', 'true_label'], axis=1).to_numpy()

    train_bdna_poison = []
    val_bdna_poison = []

    if poison_train:
        poison_train_count = len(train_nonb)
        train_bdna_poison = train_bdna[:poison_train_count]
        train_bdna = train_bdna[poison_train_count:]

        train_bdna_poison = np.array([train_bdna_poison[:, :50], train_bdna_poison[:, 50:]]).transpose(1, 0, 2)

    if poison_val:
        poison_val_count = len(val_nonb)
        val_bdna_poison = val_bdna[:poison_val_count]
        val_bdna = val_bdna[poison_val_count:]
        # print(train_bdna.shape)
        train_bdna = np.array([train_bdna[:, :50], train_bdna[:, 50:]]).transpose(1, 0, 2)
        # print(train_bdna.shape)
        val_bdna = np.array([val_bdna[:, :50], val_bdna[:, 50:]]).transpose(1, 0, 2)
        test_bdna = np.array([test_bdna[:, :50], test_bdna[:, 50:]]).transpose(1, 0, 2)

        train_nonb = np.array([train_nonb[:, :50], train_nonb[:, 50:]]).transpose(1, 0, 2)
        val_nonb = np.array([val_nonb[:, :50], val_nonb[:, 50:]]).transpose(1, 0, 2)
        test_nonb = np.array([test_nonb[:, :50], test_nonb[:, 50:]]).transpose(1, 0, 2)

        val_bdna_poison = np.array([val_bdna_poison[:, :50], val_bdna_poison[:, 50:]]).transpose(1, 0, 2)


    return train_bdna, val_bdna, test_bdna, train_nonb, val_nonb, test_nonb, train_bdna_poison, val_bdna_poison


class PUContrastiveDataset(Dataset):
    def __init__(self, normal_ds, noisy_ds):
        self.normal_ds = normal_ds 
        self.noisy_ds = noisy_ds 
        
    def __len__(self):
        return len(self.normal_ds)

    def __getitem__(self, idx):
        x_anchor = self.normal_ds[idx]
        pos_idx = torch.randint(0, len(self.normal_ds), (1,)).item()
        x_pos = self.normal_ds[pos_idx]
        noisy_idx = idx % len(self.noisy_ds)
        x_noisy = self.noisy_ds[noisy_idx]
        return (
            x_anchor.astype(np.float32), 
            x_pos.astype(np.float32), 
            x_noisy.astype(np.float32)
        )

class NRContrastiveDataset(Dataset):
    def __init__(self, normal_ds, noisy_ds):
        self.normal_ds = normal_ds
        self.noisy_ds = noisy_ds

    def __len__(self):
        return len(self.normal_ds)

    def __getitem__(self, idx):
        x_anchor = self.normal_ds[idx]
        noisy_idx = idx % len(self.noisy_ds)
        x_noisy = self.noisy_ds[noisy_idx]
        return (
            x_anchor.astype(np.float32), # bdna data
            x_noisy.astype(np.float32)   # nonb data
        )


class IndexedNoisyDataset(torch.utils.data.Dataset):
    def __init__(self, normal_ds, noisy_ds):
        self.normal_ds = normal_ds
        self.noisy_ds = noisy_ds

    def __len__(self):
        return len(self.normal_ds)

    def __getitem__(self, idx):
        x_anchor = self.normal_ds[idx]
        noisy_idx = idx % len(self.noisy_ds)
        x_noisy = self.noisy_ds[noisy_idx]
        
        return (
            x_anchor.astype(np.float32), # bdna data
            x_noisy.astype(np.float32),  # nonb data
            noisy_idx
        )

class Database():
    # now fix to simulation
    def __init__(self, data_path=None,
                 data_type='simulation',
                 win_size=50,
                 nonb_ratio=0.10,
                #  nonb_type=['G_Quadruplex_Motif', 'Short_Tandem_Repeat'],
                #nonb_type='Short_Tandem_Repeat',
                 nonb_type='G_Quadruplex_Motif',
                 poison_train=False,
                 poison_val=True,
                 use_label=True):
        self.data_path = data_path
        self.data_type = data_type
        self.win_size = win_size
        self.nonb_ratio = nonb_ratio
        if not isinstance(nonb_type, list):
            self.nonb_type = [nonb_type]
        self.poison_train = poison_train
        self.poison_val = poison_val
        self.use_label = use_label
        self.batch_size = 256
        self.num_workers = 4
        if self.data_path is None:
            self.data_path = SIMULATED_DATA if data_type == 'simulation' else EXPERIMENT_DATA
        if data_type == 'simulation':
            self.train_sim, self.val_sim, self.test_sim = load_data2(folder=self.data_path, 
                                                                     win_size=self.win_size, 
                                                                     n_bdna=200_000, 
                                                                     n_nonb=20_000,
                                                                     nonb_ratio=self.nonb_ratio)
            self.train_bdna, self.val_bdna, self.test_bdna, self.train_nonb, self.val_nonb, \
            self.test_nonb, self.train_bdna_poison, self.val_bdna_poison = reprocess_data(nonb_type=self.nonb_type, 
                                                                                           train=self.train_sim,
                                                                                            val=self.val_sim, 
                                                                                            test=self.test_sim,
                                                                                            poison_train=self.poison_train,
                                                                                            poison_val=self.poison_val)
        else:
            self.nonb_type = nonb_type
            self.train_bdna, self.val_bdna, self.test_bdna, self.train_nonb, self.val_nonb, \
            self.test_nonb, self.train_bdna_poison, self.val_bdna_poison = import_data(self.data_path, nonb_type=self.nonb_type,
                                                                        poison_train=self.poison_train,
                                                                        poison_val=self.poison_val)
            
        
        if use_label:     # don't mixup b and nonb in training set.
            self.X_train_outlier = self.train_nonb
            self.outlier_label_train = -np.ones((self.X_train_outlier.shape[0]))
            self.real_label_train = np.ones((self.train_bdna.shape[0]))
            self.outlier_dataset_train = TensorDataset(torch.tensor(self.X_train_outlier.astype(np.float32)), \
                                                    torch.tensor(self.outlier_label_train.astype(np.float32)))
            self.train_dataset = TensorDataset(torch.tensor(self.train_bdna.astype(np.float32)), \
                                            torch.tensor(self.real_label_train.astype(np.float32)))
            
            self.train_loader = DataLoader(dataset=self.train_dataset, batch_size=self.batch_size, pin_memory=True, \
                                      num_workers=self.num_workers, shuffle=True, drop_last=True)

            self.val_dataset = TensorDataset(torch.tensor(self.val_bdna.astype(np.float32)))
            self.val_loader = DataLoader(dataset=self.val_dataset, batch_size=self.batch_size, pin_memory=True, \
                                num_workers=self.num_workers, shuffle=False, drop_last=True)

            self.test_dataset_bdna = TensorDataset(torch.tensor(self.test_bdna.astype(np.float32)))
            self.test_loader_bdna = DataLoader(dataset=self.test_dataset_bdna, batch_size=self.batch_size, pin_memory=True,
                                         num_workers=self.num_workers, shuffle=False, drop_last=False)
            self.test_dataset_nonb = TensorDataset(torch.tensor(self.test_nonb.astype(np.float32)))
            self.test_loader_nonb = DataLoader(dataset=self.test_dataset_nonb, batch_size=self.batch_size, pin_memory=True,
                                         num_workers=self.num_workers, shuffle=False, drop_last=False)
            
    def get_contrastive_dataloader(self):
        dataset = PUContrastiveDataset(normal_ds=self.train_bdna, noisy_ds=self.train_nonb)
        dataloader = DataLoader(dataset=dataset, batch_size=self.batch_size, pin_memory=True, shuffle=True, drop_last=True)
        return dataloader
    

    def get_contrastive_dataloader_nr(self):
        dataset = NRContrastiveDataset(normal_ds=self.train_bdna, noisy_ds=self.train_nonb)
        dataloader = DataLoader(dataset=dataset, batch_size=512, pin_memory=True, shuffle=True, drop_last=True)
        return dataloader
    
    def get_contrastive_dataloader_elr(self):
        dataset = IndexedNoisyDataset(normal_ds=self.train_bdna, noisy_ds=self.train_nonb)
        dataloader = DataLoader(dataset=dataset, batch_size=512, pin_memory=True, shuffle=True, drop_last=True)
        return dataloader


if __name__ == "__main__":
    database = Database()
    # print(database.train_sim.columns)
    # print(database.train_sim['label'].unique())
    # print(database.train_sim['true_label'].unique())

    # enc = Encoder_P1()
    # dec = Decoder()
    # model = LightAutoEncoder(encoder=enc, decoder=dec)
    # trainer = L.Trainer(max_epochs=20, accelerator="auto")
    # trainer.fit(model, train_dataloaders=database.train_loader, val_dataloaders=database.val_loader)


    for batch in database.get_contrastive_dataloader():
        print(batch[0].shape, batch[1].shape, batch[2].shape)


    print(database.test_sim)


