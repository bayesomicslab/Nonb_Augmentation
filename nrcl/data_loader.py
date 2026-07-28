"""
data_loader.py
Unified data loader that supports both OLD format and NEW format.
"""

import os
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, TensorDataset, Dataset
from typing import Dict, List, Optional, Tuple, Union
from scipy.ndimage import median_filter as _median_filter
from paths import EXPERIMENT_DATA


# metadata columns in new format
NEW_FORMAT_META_COLS = ['chr', 'strand', 'win_start50', 'win_end50',
                        'motif_start', 'motif_end', 'features']
FEATURE_COLS = [f'forward_{i}' for i in range(50)] + [f'reverse_{i}' for i in range(50)]


# Format detection & CSV reading

def _detect_format(csv_path: str) -> str:
    """Peek at the header to decide old vs new format."""
    df_head = pd.read_csv(csv_path, nrows=0)
    if 'chr' in df_head.columns:
        return 'new'
    return 'old'


def _resolve_val_path(data_path: str, nonb_type: str) -> str:
    """Try _val.csv first (new), fall back to _validation.csv (old)."""
    val_path = os.path.join(data_path, f'{nonb_type}_val.csv')
    if os.path.exists(val_path):
        return val_path
    val_path2 = os.path.join(data_path, f'{nonb_type}_validation.csv')
    if os.path.exists(val_path2):
        return val_path2
    raise FileNotFoundError(
        f'Validation file not found. Tried:\n  {val_path}\n  {val_path2}')


def _read_csv_auto(csv_path: str) -> pd.DataFrame:
    """Read CSV with auto-detection of index column."""
    fmt = _detect_format(csv_path)
    if fmt == 'new':
        return pd.read_csv(csv_path)          # no index_col
    else:
        return pd.read_csv(csv_path, index_col=0)  # old format has numeric index


def _extract_features_and_meta(
    df: pd.DataFrame,
    fmt: str,
) -> Tuple[np.ndarray, Optional[pd.DataFrame]]:
    """
    Extract (N, 100) feature array and optionally the metadata DataFrame.

    Returns:
        features: np.ndarray of shape (N, 100)
        metadata: pd.DataFrame with NEW_FORMAT_META_COLS + 'label', or None for old format
    """
    if fmt == 'new':
        features = df[FEATURE_COLS].to_numpy(dtype=np.float64)
        meta = df[NEW_FORMAT_META_COLS + ['label']].copy()
        return features, meta
    else:
        # old format: everything except 'label' is feature
        features = df.drop(columns=['label']).to_numpy(dtype=np.float64)
        return features, None


# Core data import

def import_data(
    data_path: str,
    nonb_type: Union[str, List[str]],
    poison_train: bool = False,
    poison_val: bool = True,
    smooth: Optional[str] = None,
    smooth_kernel: int = 5,
) -> Dict:
    """
    Load train/val/test CSVs with auto format detection.

    Returns a dict:
        train_bdna, val_bdna, test_bdna         : np.ndarray (N, 2, 50)
        train_nonb, val_nonb, test_nonb         : np.ndarray (N, 2, 50)
        train_bdna_poison, val_bdna_poison       : np.ndarray or []
        metadata: {
            'train_nonb': pd.DataFrame or None,
            'val_nonb':   pd.DataFrame or None,
            'test_nonb':  pd.DataFrame or None,
            'test_bdna':  pd.DataFrame or None,
        }
    """
    if isinstance(nonb_type, list):
        nonb_type_str = nonb_type[0]
        nonb_type_set = set(nonb_type)
    else:
        nonb_type_str = nonb_type
        nonb_type_set = {nonb_type}

    # resolve file paths
    train_path = os.path.join(data_path, f'{nonb_type_str}_train.csv')
    val_path = _resolve_val_path(data_path, nonb_type_str)
    test_path = os.path.join(data_path, f'{nonb_type_str}_test.csv')

    # detect format from train file
    fmt = _detect_format(train_path)
    print(f'[data_loader] Detected format: {fmt}  (nonb_type={nonb_type_str})')

    # read CSVs
    train_df = _read_csv_auto(train_path)
    val_df = _read_csv_auto(val_path)
    test_df = _read_csv_auto(test_path)

    # split by label
    def split_and_extract(df):
        bdna_df = df[df['label'] == 'bdna']
        nonb_df = df[df['label'].isin(nonb_type_set)]
        bdna_feat, bdna_meta = _extract_features_and_meta(bdna_df, fmt)
        nonb_feat, nonb_meta = _extract_features_and_meta(nonb_df, fmt)
        return bdna_feat, nonb_feat, bdna_meta, nonb_meta

    train_bdna_raw, train_nonb_raw, _, meta_train_nonb = split_and_extract(train_df)
    val_bdna_raw, val_nonb_raw, _, meta_val_nonb = split_and_extract(val_df)
    test_bdna_raw, test_nonb_raw, meta_test_bdna, meta_test_nonb = split_and_extract(test_df)

    # poison logic
    train_bdna_poison = []
    val_bdna_poison = []

    if poison_train:
        n_poison = len(train_nonb_raw)
        train_bdna_poison = train_bdna_raw[:n_poison]
        train_bdna_raw = train_bdna_raw[n_poison:]
        train_bdna_poison = _reshape_to_2ch(train_bdna_poison)

    if poison_val:
        n_poison = len(val_nonb_raw)
        val_bdna_poison = val_bdna_raw[:n_poison]
        val_bdna_raw = val_bdna_raw[n_poison:]
        val_bdna_poison = _reshape_to_2ch(val_bdna_poison)

    # optional smoothing to suppress single-position outliers
    if smooth is not None:
        print(f'[data_loader] Smoothing: method={smooth}, kernel={smooth_kernel}')
        train_bdna_raw = _smooth_signals(train_bdna_raw, smooth, smooth_kernel)
        val_bdna_raw = _smooth_signals(val_bdna_raw, smooth, smooth_kernel)
        test_bdna_raw = _smooth_signals(test_bdna_raw, smooth, smooth_kernel)
        train_nonb_raw = _smooth_signals(train_nonb_raw, smooth, smooth_kernel)
        val_nonb_raw = _smooth_signals(val_nonb_raw, smooth, smooth_kernel)
        test_nonb_raw = _smooth_signals(test_nonb_raw, smooth, smooth_kernel)

    # reshape all to (N, 2, 50)
    train_bdna = _reshape_to_2ch(train_bdna_raw)
    val_bdna = _reshape_to_2ch(val_bdna_raw)
    test_bdna = _reshape_to_2ch(test_bdna_raw)
    train_nonb = _reshape_to_2ch(train_nonb_raw)
    val_nonb = _reshape_to_2ch(val_nonb_raw)
    test_nonb = _reshape_to_2ch(test_nonb_raw)

    return {
        'train_bdna': train_bdna,
        'val_bdna': val_bdna,
        'test_bdna': test_bdna,
        'train_nonb': train_nonb,
        'val_nonb': val_nonb,
        'test_nonb': test_nonb,
        'train_bdna_poison': train_bdna_poison,
        'val_bdna_poison': val_bdna_poison,
        'metadata': {
            'train_nonb': meta_train_nonb,
            'val_nonb': meta_val_nonb,
            'test_nonb': meta_test_nonb,
            'test_bdna': meta_test_bdna,
        },
    }


def _smooth_signals(arr: np.ndarray, method: str = 'median', kernel: int = 5) -> np.ndarray:
    """Smooth translocation-time signals to suppress single-position outliers.

    arr is (N, 100): the first 50 columns are forward, the last 50 reverse.
    method is 'median' or 'mean'; kernel is the (odd) window size.
    """
    if arr.ndim != 2 or arr.shape[1] != 100:
        return arr
    smoothed = np.empty_like(arr)
    for ch_start in (0, 50):  # forward, then reverse
        block = arr[:, ch_start:ch_start + 50]
        if method == 'median':
            smoothed[:, ch_start:ch_start + 50] = _median_filter(block, size=(1, kernel))
        else:  # mean / moving average
            from scipy.ndimage import uniform_filter1d
            smoothed[:, ch_start:ch_start + 50] = uniform_filter1d(block, size=kernel, axis=1)
    return smoothed


def _reshape_to_2ch(arr: np.ndarray) -> np.ndarray:
    """(N, 100) -> (N, 2, 50): channel 0 = forward, channel 1 = reverse."""
    if arr.ndim == 2 and arr.shape[1] == 100:
        return np.stack([arr[:, :50], arr[:, 50:]], axis=1)
    return arr


# Dataset wrappers

class NRContrastiveDataset(Dataset):
    def __init__(self, normal_ds: np.ndarray, noisy_ds: np.ndarray):
        self.normal_ds = normal_ds
        self.noisy_ds = noisy_ds

    def __len__(self):
        return len(self.normal_ds)

    def __getitem__(self, idx):
        x_anchor = self.normal_ds[idx]
        noisy_idx = idx % len(self.noisy_ds)
        x_noisy = self.noisy_ds[noisy_idx]
        return x_anchor.astype(np.float32), x_noisy.astype(np.float32)


class IndexedNoisyDataset(Dataset):
    def __init__(self, normal_ds: np.ndarray, noisy_ds: np.ndarray):
        self.normal_ds = normal_ds
        self.noisy_ds = noisy_ds

    def __len__(self):
        return len(self.normal_ds)

    def __getitem__(self, idx):
        x_anchor = self.normal_ds[idx]
        noisy_idx = idx % len(self.noisy_ds)
        x_noisy = self.noisy_ds[noisy_idx]
        return x_anchor.astype(np.float32), x_noisy.astype(np.float32), noisy_idx


# Database: unified high-level interface

class Database:
    def __init__(
        self,
        data_path: str = EXPERIMENT_DATA,
        data_type: str = 'experimental',
        nonb_type: Union[str, List[str]] = 'G_Quadruplex_Motif',
        poison_train: bool = False,
        poison_val: bool = True,
        batch_size: int = 256,
        num_workers: int = 4,
        smooth: Optional[str] = None,
        smooth_kernel: int = 5,
    ):
        self.data_path = data_path
        self.data_type = data_type
        self.batch_size = batch_size
        self.num_workers = num_workers

        if not isinstance(nonb_type, list):
            self.nonb_type = [nonb_type]
        else:
            self.nonb_type = nonb_type

        # load data
        result = import_data(
            data_path,
            nonb_type=self.nonb_type,
            poison_train=poison_train,
            poison_val=poison_val,
            smooth=smooth,
            smooth_kernel=smooth_kernel,
        )

        self.train_bdna = result['train_bdna']
        self.val_bdna = result['val_bdna']
        self.test_bdna = result['test_bdna']
        self.train_nonb = result['train_nonb']
        self.val_nonb = result['val_nonb']
        self.test_nonb = result['test_nonb']
        self.train_bdna_poison = result['train_bdna_poison']
        self.val_bdna_poison = result['val_bdna_poison']
        self.metadata = result['metadata']   # dict of DataFrames

        print(f'[Database] train_bdna={self.train_bdna.shape}, train_nonb={self.train_nonb.shape}')
        print(f'[Database] val_bdna={self.val_bdna.shape}, val_nonb={self.val_nonb.shape}')
        print(f'[Database] test_bdna={self.test_bdna.shape}, test_nonb={self.test_nonb.shape}')
        if self.metadata['test_nonb'] is not None:
            print(f'[Database] metadata preserved for test_nonb ({len(self.metadata["test_nonb"])} rows)')

        # build dataloaders
        self._build_loaders()

    def _build_loaders(self):
        def to_tensor(arr):
            return torch.tensor(arr.astype(np.float32))

        self.train_dataset = TensorDataset(to_tensor(self.train_bdna))
        self.train_loader = DataLoader(
            self.train_dataset, batch_size=self.batch_size,
            pin_memory=True, num_workers=self.num_workers,
            shuffle=True, drop_last=True)

        self.val_dataset = TensorDataset(to_tensor(self.val_bdna))
        self.val_loader = DataLoader(
            self.val_dataset, batch_size=self.batch_size,
            pin_memory=True, num_workers=self.num_workers,
            shuffle=False, drop_last=True)

        self.test_dataset_bdna = TensorDataset(to_tensor(self.test_bdna))
        self.test_loader_bdna = DataLoader(
            self.test_dataset_bdna, batch_size=self.batch_size,
            pin_memory=True, num_workers=self.num_workers,
            shuffle=False, drop_last=False)

        self.test_dataset_nonb = TensorDataset(to_tensor(self.test_nonb))
        self.test_loader_nonb = DataLoader(
            self.test_dataset_nonb, batch_size=self.batch_size,
            pin_memory=True, num_workers=self.num_workers,
            shuffle=False, drop_last=False)

    def get_contrastive_dataloader_nr(self):
        dataset = NRContrastiveDataset(self.train_bdna, self.train_nonb)
        return DataLoader(dataset, batch_size=512, pin_memory=True,
                          shuffle=True, drop_last=True)

    def get_contrastive_dataloader_elr(self):
        dataset = IndexedNoisyDataset(self.train_bdna, self.train_nonb)
        return DataLoader(dataset, batch_size=512, pin_memory=True,
                          shuffle=True, drop_last=True)

    # helpers for result saving

    def get_test_nonb_metadata(self) -> Optional[pd.DataFrame]:
        return self.metadata.get('test_nonb', None)

    def get_test_bdna_metadata(self) -> Optional[pd.DataFrame]:
        return self.metadata.get('test_bdna', None)
