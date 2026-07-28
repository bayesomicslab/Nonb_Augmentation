#!/bin/bash
# Run the core GoFAE-DND model (gofae.py) on all 7 non-B DNA types using the new
# ONT experimental data (no BED files; coordinates are embedded in the CSVs).
set -e

EXP_DATA_PATH="${NONB_EXPERIMENT_DATA:-/path/to/Experiment_Data}"
EPOCHS="${EPOCHS:-50}"
FDR="${FDR:-0.2}"
GPU="${GPU:-0}"
LOGDIR="run_logs_gofae"
mkdir -p "${LOGDIR}"

NONB_TYPES=(
    "Z_DNA_Motif"
    "G_Quadruplex_Motif"
    "A_Phased_Repeat"
    "Direct_Repeat"
    "Mirror_Repeat"
    "Short_Tandem_Repeat"
    "Inverted_Repeat"
)

for nonb in "${NONB_TYPES[@]}"; do
    echo "GoFAE-DND on ${nonb} (epochs=${EPOCHS}, fdr=${FDR})"
    CUDA_VISIBLE_DEVICES="${GPU}" python -u gofae.py \
        --data_type experimental \
        --nonb_type "${nonb}" \
        --exp_data_path "${EXP_DATA_PATH}" \
        --epochs "${EPOCHS}" \
        --fdr_level "${FDR}" \
        > "${LOGDIR}/${nonb}.log" 2>&1
    echo "Finished ${nonb}; tail of log:"
    tail -5 "${LOGDIR}/${nonb}.log"
    echo ""
done

echo "All 7 nonb_types finished."
