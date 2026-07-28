#!/usr/bin/env python3
"""Mutation rate vs repeat number, one line per k-mer length, for detected /
undetected / control. Reads the deduplicated per-locus counts.
"""

import os
import numpy as np
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

WORK=os.path.join(os.environ.get("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "str_unit_mutrate")
COUNTS=f"{WORK}/counts.dedup.bed"
KMERS=[1,2,3]; MINN=10

buck=defaultdict(list)
for line in open(COUNTS):
    a=line.rstrip("\n").split("\t")
    if len(a)<5: continue
    grp,unit,nrep=a[3].split("|"); L=int(a[2])-int(a[1])
    if L<=0: continue
    buck[(grp,len(unit),int(nrep))].append(int(a[4])/L)

agg={}
for (grp,k,nrep),v in buck.items():
    v=np.array(v); n=len(v)
    if n<MINN: continue
    agg[(grp,k,nrep)]=(n, v.mean(), v.std(ddof=1)/np.sqrt(n) if n>1 else 0.0)

KCOL={1:"#1f77b4", 2:"#d62728", 3:"#e8a800"}            # mono / di / tri
KLAB={1:"mononucleotide",2:"dinucleotide",3:"trinucleotide"}
LS  ={"detected":"-","undetected":"--","control":":"}
MK  ={"detected":"v","undetected":"^","control":"o"}

fig,ax=plt.subplots(figsize=(8.2,5.6))
for k in KMERS:
    for grp in ("detected","undetected","control"):
        pts=sorted((nrep,*agg[(grp,k,nrep)][1:]) for (g,kk,nrep) in agg if g==grp and kk==k)
        if len(pts)<2: continue
        x=np.array([p[0] for p in pts]); m=np.array([p[1] for p in pts]); se=np.array([p[2] for p in pts])
        ly=np.log10(m)
        lo=np.log10(np.clip(m-se,1e-9,None)); hi=np.log10(m+se)
        ax.fill_between(x,lo,hi,color=KCOL[k],alpha=0.18,lw=0)
        ax.plot(x,ly,LS[grp],color=KCOL[k],lw=1.6,marker=MK[grp],ms=4,
                markerfacecolor=KCOL[k] if grp!="control" else "white",
                markeredgecolor=KCOL[k])
ax.set_xlabel("Repeat number",fontsize=13)
ax.set_ylabel(r"$\log_{10}$(mutation rate, SNPs/bp)",fontsize=13)
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.tick_params(labelsize=11)

kmer_h=[Line2D([0],[0],color=KCOL[k],lw=2.4,label=KLAB[k]) for k in KMERS]
sig_h =[Line2D([0],[0],color="0.3",ls=LS[g],marker=MK[g],ms=5,mfc=("white" if g=="control" else "0.3"),
               label=g) for g in ("detected","undetected","control")]
# Legends sit outside the axes so they never cover the curves.
leg1=ax.legend(handles=kmer_h,fontsize=10,frameon=True,title="k-mer length",
               loc="upper left",bbox_to_anchor=(1.02,1.0),borderaxespad=0)
ax.add_artist(leg1)
ax.legend(handles=sig_h,fontsize=10,frameon=True,title="signal",
          loc="upper left",bbox_to_anchor=(1.02,0.66),borderaxespad=0)
ax.set_title("STR mutation rate vs repeat number, by k-mer length\n(detected / undetected / control)",fontsize=12)
fig.tight_layout()
fig.savefig(f"{WORK}/str_mutrate_bykmer_line_dedup.png",dpi=180,bbox_inches="tight")
fig.savefig(f"{WORK}/str_mutrate_bykmer_line_dedup.pdf",bbox_inches="tight")
print("saved str_mutrate_bykmer_line_dedup.{png,pdf}")
for k in KMERS:
    for grp in ("detected","undetected","control"):
        xs=sorted(nrep for (g,kk,nrep) in agg if g==grp and kk==k)
        if xs: print(f"  k{k} {grp:11s}: nrep {xs[0]}..{xs[-1]} ({len(xs)} pts)")
