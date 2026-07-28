import os

import numpy as np
import pandas as pd


# Columns in the new (ONT) experimental CSVs that are *not* model features.
# Everything else (forward_0..forward_49, reverse_0..reverse_49) is a feature.
_META_COLS = ['chr', 'strand', 'win_start50', 'win_end50',
              'motif_start', 'motif_end', 'features', 'label']
# Genomic-coordinate columns used to (re)build the BED of flagged windows.
_BED_COLS = ['chr', 'strand', 'win_start50', 'win_end50', 'motif_start', 'motif_end']


def _feature_columns(df):
    """Return the ordered list of feature columns (forward_*/reverse_*).

    Falls back to "every column that is not metadata" so the loader keeps
    working if the feature naming convention ever changes.
    """
    fwd = sorted([c for c in df.columns if str(c).startswith('forward_')],
                 key=lambda c: int(str(c).split('_')[1]))
    rev = sorted([c for c in df.columns if str(c).startswith('reverse_')],
                 key=lambda c: int(str(c).split('_')[1]))
    if fwd and rev:
        return fwd + rev
    return [c for c in df.columns if c not in _META_COLS]


def _maybe_write_test_bed(path, nonb_type, test, mask):
    """Write a BED-like CSV of the test non-B windows (in the same row order
    as ``test_nonb``) so the downstream filtering step has coordinates.

    The new experimental data ships coordinates inside the main CSV rather than
    as a separate ``*_test_bed.csv`` file, so we materialize it here (cached).
    """
    bed_path = os.path.join(path, str(nonb_type) + '_test_bed.csv')
    if os.path.exists(bed_path):
        return
    bed_cols = [c for c in _BED_COLS if c in test.columns]
    if not bed_cols:
        return
    test.loc[mask, bed_cols].reset_index(drop=True).to_csv(bed_path, index=False)


def import_data(path, nonb_type, poison_train=False, poison_val=True):
    train = pd.read_csv(path + str(nonb_type) + '_train.csv')
    val = pd.read_csv(path + str(nonb_type) + '_val.csv')
    test = pd.read_csv(path + str(nonb_type) + '_test.csv')

    feat_cols = _feature_columns(train)

    # float32 keeps the (large) feature arrays at half the memory of the pandas
    # float64 default; the model casts to float32 downstream anyway.
    train_bdna = train.loc[train["label"] == 'bdna', feat_cols].to_numpy(dtype=np.float32)
    val_bdna = val.loc[val["label"] == 'bdna', feat_cols].to_numpy(dtype=np.float32)
    test_bdna = test.loc[test["label"] == 'bdna', feat_cols].to_numpy(dtype=np.float32)
    train_nonb = train.loc[train["label"] == nonb_type, feat_cols].to_numpy(dtype=np.float32)
    val_nonb = val.loc[val["label"] == nonb_type, feat_cols].to_numpy(dtype=np.float32)
    test_nonb_mask = test["label"] == nonb_type
    test_nonb = test.loc[test_nonb_mask, feat_cols].to_numpy(dtype=np.float32)

    # Materialize the BED coordinates for the flagged-window export step.
    _maybe_write_test_bed(path, nonb_type, test, test_nonb_mask)

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

        train_bdna = np.array([train_bdna[:, :50], train_bdna[:, 50:]]).transpose(1, 0, 2)
        val_bdna = np.array([val_bdna[:, :50], val_bdna[:, 50:]]).transpose(1, 0, 2)
        test_bdna = np.array([test_bdna[:, :50], test_bdna[:, 50:]]).transpose(1, 0, 2)

        train_nonb = np.array([train_nonb[:, :50], train_nonb[:, 50:]]).transpose(1, 0, 2)
        val_nonb = np.array([val_nonb[:, :50], val_nonb[:, 50:]]).transpose(1, 0, 2)
        test_nonb = np.array([test_nonb[:, :50], test_nonb[:, 50:]]).transpose(1, 0, 2)

        val_bdna_poison = np.array([val_bdna_poison[:, :50], val_bdna_poison[:, 50:]]).transpose(1, 0, 2)

    return train_bdna, val_bdna, test_bdna, train_nonb, val_nonb, test_nonb, train_bdna_poison, val_bdna_poison


def reprocess_data2(nonb_type, train, val, test, poison_train=False, poison_val=True):
    train_bdna = train[train["label"] == 'bdna'].drop(['label', 'true_label'], axis=1).to_numpy()
    val_bdna = val[val["label"] == 'bdna'].drop(['label', 'true_label'], axis=1).to_numpy()
    test_bdna = test[test["label"] == 'bdna'].drop(['label', 'true_label'], axis=1).to_numpy()
    train_nonb = train[train["label"] == nonb_type].drop(['label', 'true_label'], axis=1).to_numpy()
    val_nonb = val[val["label"] == nonb_type].drop(['label', 'true_label'], axis=1).to_numpy()
    test_nonb = test[test["label"] == nonb_type].drop(['label', 'true_label'], axis=1).to_numpy()

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

        train_bdna = np.array([train_bdna[:, :50], train_bdna[:, 50:]]).transpose(1, 0, 2)
        val_bdna = np.array([val_bdna[:, :50], val_bdna[:, 50:]]).transpose(1, 0, 2)
        test_bdna = np.array([test_bdna[:, :50], test_bdna[:, 50:]]).transpose(1, 0, 2)

        train_nonb = np.array([train_nonb[:, :50], train_nonb[:, 50:]]).transpose(1, 0, 2)
        val_nonb = np.array([val_nonb[:, :50], val_nonb[:, 50:]]).transpose(1, 0, 2)
        test_nonb = np.array([test_nonb[:, :50], test_nonb[:, 50:]]).transpose(1, 0, 2)

        val_bdna_poison = np.array([val_bdna_poison[:, :50], val_bdna_poison[:, 50:]]).transpose(1, 0, 2)

    return train_bdna, val_bdna, test_bdna, train_nonb, val_nonb, test_nonb, train_bdna_poison, val_bdna_poison



