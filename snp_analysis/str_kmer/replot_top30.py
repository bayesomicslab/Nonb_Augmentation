#!/usr/bin/env python3
"""Plot the top-N STR units from str_unit_mutrate.csv. Plotting only; the rates
are computed by str_unit_mutrate.py."""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

WORK=os.path.join(os.environ.get("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "str_unit_mutrate")
agg=pd.read_csv(f"{WORK}/str_unit_mutrate.csv")
MINN=10; TOPN=12

def period(u): return len(u)
tot=(agg[agg.group.isin(["detected","undetected"])]
     .groupby("unit")["n"].sum().sort_values(ascending=False))
top=list(tot.index[:TOPN])
print("top units:",top)

COL={"detected":"#d62728","undetected":"#1f77b4","control":"0.45"}
ncol=4; nrow=int(np.ceil(len(top)/ncol))
fig,axes=plt.subplots(nrow,ncol,figsize=(3.2*ncol,2.6*nrow))
axes=np.array(axes).reshape(-1)
for ax,unit in zip(axes,top):
    xs=set()
    for grp in ("detected","undetected","control"):
        sub=agg[(agg.group==grp)&(agg.unit==unit)&(agg.n>=MINN)].sort_values("nrep")
        if not len(sub): continue
        xs.update(int(v) for v in sub["nrep"])
        ax.errorbar(sub["nrep"],sub["rate"],yerr=sub["sem"],
                    fmt="o"+("--" if grp=="control" else "-"),
                    color=COL[grp],ms=3,lw=1.3,capsize=2)
    ax.set_title(f"{unit}  ({period(unit)}-mer, n={int(tot[unit]):,})",fontsize=9)
    ax.set_xlabel("number of repeats",fontsize=8); ax.tick_params(labelsize=7)
    ax.set_ylabel("SNPs / bp",fontsize=8)
    if xs:
        xs=sorted(xs)
        ax.set_xticks(xs)
        if len(xs)==1: ax.set_xlim(xs[0]-0.5, xs[0]+0.5)
    else:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
for ax in axes[len(top):]: ax.axis("off")
h=[Line2D([0],[0],color=COL[g],ls=("--" if g=="control" else "-"),marker="o",ms=4,label=g)
   for g in ("detected","undetected","control")]
fig.legend(handles=h,fontsize=9,loc="upper right",ncol=3,frameon=False,bbox_to_anchor=(0.99,1.0))
fig.suptitle(f"STR mutation rate vs number of repeats  (top {TOPN} units, MINN={MINN})  detected / undetected / control",
             fontsize=11,y=1.0)
fig.tight_layout(rect=[0,0,1,0.97])
fig.savefig(f"{WORK}/str_unit_mutrate.png",dpi=180,bbox_inches="tight")
print("saved str_unit_mutrate.png (top 12)")
