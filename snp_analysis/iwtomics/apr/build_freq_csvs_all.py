"""Build the per-bin SNP frequency CSVs that plot_apr_iwtomics_all.R expects.

Same logic as the APR notebook's plotting section, pointed at apr_notebook_run/.
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd

WORKDIR = os.path.join(os.environ.get("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "apr_notebook_run")
os.chdir(WORKDIR)

MAX_FLANK   = 100
TRACT_BINS  = 40
SPACER_BINS = 30
PART_LEN = {"Tract1": TRACT_BINS, "Spacer1": SPACER_BINS,
            "Tract2": TRACT_BINS, "Spacer2": SPACER_BINS,
            "Tract3": TRACT_BINS}
scaledposTract  = [round(i, 4) for i in np.arange(1, TRACT_BINS  * 2, 2.0) / (TRACT_BINS  * 2)]
scaledposSpacer = [round(i, 4) for i in np.arange(1, SPACER_BINS * 2, 2.0) / (SPACER_BINS * 2)]

_counts = {k: int(v) for k, v in
           (line.split("\t") for line in Path("split_counts.txt").read_text().strip().splitlines())}
N_MOTIF, N_CTRL = _counts["N_MOTIF"], _counts["N_CTRL"]
print(f"N_MOTIF={N_MOTIF}  N_CTRL={N_CTRL}")


def _build_matrices(intersect_path, n_loci):
    out = {
        "left":    np.zeros((n_loci, MAX_FLANK),   dtype=np.float32),
        "right":   np.zeros((n_loci, MAX_FLANK),   dtype=np.float32),
        "Tract1":  np.zeros((n_loci, TRACT_BINS),  dtype=np.float32),
        "Spacer1": np.zeros((n_loci, SPACER_BINS), dtype=np.float32),
        "Tract2":  np.zeros((n_loci, TRACT_BINS),  dtype=np.float32),
        "Spacer2": np.zeros((n_loci, SPACER_BINS), dtype=np.float32),
        "Tract3":  np.zeros((n_loci, TRACT_BINS),  dtype=np.float32),
    }
    seen = set()
    with open(intersect_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 8: continue
            if len(f) >= 6 and f[5] == ".": continue
            lid, ps, pe, ptype, snp_start = int(f[3]), int(f[1]), int(f[2]), f[4], int(f[6])
            key = (lid, ptype, snp_start)
            if key in seen: continue
            seen.add(key)
            plen = pe - ps
            if ptype == "LeftFlank":
                idx = -(pe - snp_start) + MAX_FLANK
                if 0 <= idx < MAX_FLANK: out["left"][lid, idx] += 1
            elif ptype == "RightFlank":
                idx = snp_start - ps
                if 0 <= idx < MAX_FLANK: out["right"][lid, idx] += 1
            elif ptype in PART_LEN and plen > 0:
                grid = scaledposTract if ptype.startswith("Tract") else scaledposSpacer
                fs = (snp_start     - ps) / plen
                fe = (snp_start + 1 - ps) / plen
                for i, sp in enumerate(grid):
                    if sp >= fs and sp < fe:
                        out[ptype][lid, i] += 1
    return out


print("building matrices (motif) ...")
mot = _build_matrices("NCNRAPhasedRepeats.gnomAD.intersect",      N_MOTIF)
print("building matrices (control) ...")
ctl = _build_matrices("NCNRAPhasedRepeats.gnomAD.ctrl.intersect", N_CTRL)


def _freq(m, n): return m.sum(0) / n

pd.DataFrame({"LeftFlank":  _freq(mot["left"],  N_MOTIF),
              "ctrl_LeftFlank":  _freq(ctl["left"],  N_CTRL)}
             ).to_csv("iwt_freq_LeftFlank.csv",  index=False)
pd.DataFrame({"RightFlank": _freq(mot["right"], N_MOTIF),
              "ctrl_RightFlank": _freq(ctl["right"], N_CTRL)}
             ).to_csv("iwt_freq_RightFlank.csv", index=False)

center_df = {}
for nm in ("Tract1", "Spacer1", "Tract2", "Spacer2", "Tract3"):
    center_df[nm]            = _freq(mot[nm], N_MOTIF)
    center_df[f"ctrl_{nm}"]  = _freq(ctl[nm], N_CTRL)
maxL = max(PART_LEN.values())
center_padded = {k: np.concatenate([v, np.full(maxL - len(v), np.nan)]) for k, v in center_df.items()}
pd.DataFrame(center_padded).to_csv("iwt_freq_center.csv", index=False)
print("wrote iwt_freq_{LeftFlank,RightFlank,center}.csv")
