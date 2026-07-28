"""Data locations, resolved once so nothing else has to hard-code an absolute path.

"""

import os

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(_HERE)

#: Root under which all inputs are expected.
DATA_ROOT = os.environ.get("NONB_DATA_ROOT", os.path.join(REPO_ROOT, "data"))

#: ONT experimental windows: ``<nonb_type>_{train,val,test}.csv`` (+ ``_test_bed.csv``).
EXPERIMENT_DATA = os.environ.get(
    "NONB_EXPERIMENT_DATA", os.path.join(DATA_ROOT, "Experiment_Data")
)

#: Simulated translocation-time windows (``forward_*.npy`` / ``reverse_*.npy`` + centered CSVs).
SIMULATED_DATA = os.environ.get(
    "NONB_SIMULATED_DATA", os.path.join(DATA_ROOT, "simulated_data")
)

#: Directory holding the pretrained GoFAE-DND ``last_model.pt`` used for the baseline panel.
GOFAE_CHECKPOINT = os.environ.get(
    "NONB_GOFAE_CHECKPOINT", os.path.join(DATA_ROOT, "gofae_checkpoint")
)
