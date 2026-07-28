#!/bin/bash
set -e

NONB_TYPES=(
    "G_Quadruplex_Motif"
    "Short_Tandem_Repeat"
    "Z_DNA_Motif"
    "A_Phased_Repeat"
    "Direct_Repeat"
    "Inverted_Repeat"
    "Mirror_Repeat"
)

ALPHA=0.2

for nonb in "${NONB_TYPES[@]}"; do
    echo "Running ${nonb}"
    python main_sim.py --nonb_type "${nonb}" --alpha "${ALPHA}"
    echo ""
done

echo "All nonb_types finished."
