#!/usr/bin/env python3
"""One figure per k-mer length, one sub-panel per repeat unit; each panel shows
detected / undetected / control against repeat number, on linear and log10 y.
Per-locus SNP rate is aggregated from counts.dedup.bed into (group, unit, nrep).
"""
import os, sys
import numpy as np
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

WORK=os.path.join(os.environ.get("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "str_unit_mutrate")
COUNTS=f"{WORK}/counts.dedup.bed"
MINN=10
UNIT_MIN={1:0, 2:0, 3:1000, 4:500}   # min total (det+und) loci for a unit to be drawn
COL={"detected":"#d62728","undetected":"#1f77b4","control":"0.45"}
KNAME={1:"mononucleotide",2:"dinucleotide",3:"trinucleotide",4:"tetranucleotide"}

rec=defaultdict(list)          # (grp,unit,nrep) -> rates
tot=defaultdict(int)           # unit -> det+und count
for line in open(COUNTS):
    a=line.rstrip("\n").split("\t")
    if len(a)<5: continue
    grp,unit,nrep=a[3].split("|"); L=int(a[2])-int(a[1])
    if L<=0: continue
    rec[(grp,unit,int(nrep))].append(int(a[4])/L)
    if grp in ("detected","undetected"): tot[unit]+=1
agg={}
for (grp,unit,nrep),v in rec.items():
    v=np.array(v); n=len(v)
    agg[(grp,unit,nrep)]=(n, v.mean(), v.std(ddof=1)/np.sqrt(n) if n>1 else 0.0)

def units_for(k):
    us=sorted([u for u in tot if len(u)==k], key=lambda u:-tot[u])
    keep=[u for u in us if tot[u]>=UNIT_MIN[k]]
    drop=[u for u in us if tot[u]<UNIT_MIN[k]]
    return keep, drop

def plot_k(k, logy):
    keep, drop = units_for(k)
    if not keep: return None
    ncol=min(4,len(keep)); nrow=int(np.ceil(len(keep)/ncol))
    fig,axes=plt.subplots(nrow,ncol,figsize=(3.3*ncol,2.7*nrow),squeeze=False)
    axes=axes.reshape(-1)
    for ax,unit in zip(axes,keep):
        xs=set()
        for grp in ("detected","undetected","control"):
            pts=sorted((nrep,*agg[(grp,unit,nrep)][1:]) for (g,u,nrep) in agg
                       if g==grp and u==unit and agg[(g,u,nrep)][0]>=MINN)
            if not pts: continue
            x=np.array([p[0] for p in pts]); m=np.array([p[1] for p in pts]); se=np.array([p[2] for p in pts])
            xs.update(int(v) for v in x)
            y=np.log10(m) if logy else m
            yerr=(np.log10(m+se)-np.log10(np.clip(m-se,1e-9,None)))/2 if logy else se
            ax.errorbar(x,y,yerr=yerr,fmt="o"+("--" if grp=="control" else "-"),
                        color=COL[grp],ms=3,lw=1.3,capsize=2)
        ax.set_title(f"{unit}  (n={tot[unit]:,})",fontsize=9)
        ax.set_xlabel("repeat number",fontsize=8); ax.tick_params(labelsize=7)
        ax.set_ylabel(r"$\log_{10}$(SNPs/bp)" if logy else "SNPs/bp",fontsize=8)
        if xs:
            xs=sorted(xs); ax.set_xticks(xs)
            if len(xs)==1: ax.set_xlim(xs[0]-0.5,xs[0]+0.5)
        else: ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    for ax in axes[len(keep):]: ax.axis("off")
    h=[Line2D([0],[0],color=COL[g],ls=("--" if g=="control" else "-"),marker="o",ms=4,label=g)
       for g in ("detected","undetected","control")]
    fig.legend(handles=h,fontsize=9,loc="upper right",ncol=3,frameon=False,bbox_to_anchor=(0.99,1.0))
    scale="log10" if logy else "linear"
    fig.suptitle(f"{KNAME[k]} ({k}-mer): mutation rate vs repeat number, per unit  [{scale} y]  "
                 f"detected / undetected / control",fontsize=11,y=1.0)
    fig.tight_layout(rect=[0,0,1,0.96])
    out=f"{WORK}/str_perunit_dedup_k{k}_{scale}.png"
    fig.savefig(out,dpi=170,bbox_inches="tight"); plt.close(fig)
    if drop: print(f"  k{k}: drew {len(keep)} units; dropped {len(drop)} below {UNIT_MIN[k]} loci: {drop}")
    return out

ks=[int(x) for x in sys.argv[1:]] if len(sys.argv)>1 else [1,2,3,4]
for k in ks:
    for logy in (False,True):
        o=plot_k(k,logy)
        if o: print("saved",o)
