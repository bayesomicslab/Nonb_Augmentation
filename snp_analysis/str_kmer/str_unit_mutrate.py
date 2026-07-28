#!/usr/bin/env python3
"""Mutation rate (SNPs/bp, length-normalised) against number of repeats,
stratified by repeat unit (rotation-merged, strands kept separate), for
detected / undetected / control.

Per locus: rate = (#gnomAD SNPs in motif) / (motif length bp).
Per (group, unit, n_repeats): mean of per-locus rates (+ SEM, n).
Control = MR-style 1:1 length-matched bedtools shuffle (inherits unit/nrep).
"""
import os, re, subprocess, csv as _csv
from collections import defaultdict
import numpy as np, pandas as pd

BASE=os.environ.get("SNP_AF_BASE", "/path/to/snp_af_analysis")
WORK=f"{BASE}/results/str_unit_mutrate"; os.makedirs(WORK,exist_ok=True); os.chdir(WORK)
SNPDIR=f"{BASE}/snp_data/gnomad"
NCNR=f"{BASE}/results/str_notebook_detected_run/NCNR.hg38.bed"
GAPS=f"{BASE}/annotations/hg38_gaps.bed"; CHROMS=f"{BASE}/reference/hg38.chrom.sizes"
BT=os.environ.get("BEDTOOLS", "bedtools")
TESTCSV=os.environ.get("NONB_EXPERIMENT_DATA", "/path/to/Experiment_Data") + "/Short_Tandem_Repeat_test.csv"
DETBED=os.environ.get("NONB_PU_RESULTS", "/path/to/pu_results") + "/Short_Tandem_Repeat/Short_Tandem_Repeat_detected.bed"
MAXLEN=100; MINN=30; TOPN=12
_csv.field_size_limit(10**7)

def period(comp): return sum(int(re.match(r'(\d+)([ACGT])',t).group(1)) for t in comp.split('/'))
def canon_rot(u): return min(u[i:]+u[:i] for i in range(len(u)))  # rotation-merged, strand kept

# 1) feature map: (chr,mstart1,mend) -> (unit, nrep)
print("parse features ...",flush=True)
feat={}
with open(TESTCSV) as f:
    r=_csv.reader(f); h=next(r); ci=h.index('chr'); li=h.index('label'); fi=h.index('features')
    si=h.index('motif_start'); ei=h.index('motif_end')
    for row in r:
        if row[li]!='Short_Tandem_Repeat' or not row[si] or not row[ei]: continue
        d={}; xc=None
        for p in row[fi].split(';'):
            if '=' in p: k,v=p.split('=',1); d[k]=v
            elif re.match(r'x\d+\+\d+',p): xc=p
        if 'composition' not in d or 'sequence' not in d or xc is None: continue
        k=period(d['composition']); seq=d['sequence'].upper()
        if len(seq)<k: continue
        feat[(row[ci],int(float(row[si])),int(float(row[ei])))]=(canon_rot(seq[:k]),int(re.match(r'x(\d+)\+',xc).group(1)))
det=set()
dd=pd.read_csv(DETBED,sep='\t').dropna(subset=['motif_start','motif_end'])
for c,s,e in zip(dd['chr'],dd['motif_start'].astype(int),dd['motif_end'].astype(int)): det.add((c,s,e))
print(f"  loci {len(feat)}, detected {len(det)}",flush=True)

# 2) labelled motif beds (0-based), name=group|unit|nrep, <=100bp
with open("loci.bed","w") as fo, open("all_coords.bed","w") as fc:
    for (c,s1,e),(unit,nrep) in feat.items():
        st=s1-1; L=e-st
        if L<=0 or L>MAXLEN: continue
        grp='detected' if (c,s1,e) in det else 'undetected'
        fo.write(f"{c}\t{st}\t{e}\t{grp}|{unit}|{nrep}\n")
        fc.write(f"{c}\t{st}\t{e}\n")
subprocess.run(f"{BT} intersect -a loci.bed -b {NCNR} -v | sort -k1,1 -k2,2n > loci.ncnr.bed",shell=True,check=True)

# 3) matched control (shuffle, name preserved -> relabel control)
subprocess.run(f"cat all_coords.bed {GAPS} {NCNR} | cut -f1-3 | sort -k1,1 -k2,2n | {BT} merge -i stdin > excl.bed",shell=True,check=True)
subprocess.run(f"{BT} shuffle -i loci.ncnr.bed -excl excl.bed -g {CHROMS} -chrom -f 0 -noOverlapping -seed 42 "
   "| awk 'BEGIN{OFS=\"\\t\"}{split($4,a,\"|\"); print $1,$2,$3,\"control|\"a[2]\"|\"a[3]}' "
   "| sort -k1,1 -k2,2n > ctrl.ncnr.bed",shell=True,check=True)
subprocess.run("cat loci.ncnr.bed ctrl.ncnr.bed | sort -k1,1 -k2,2n > all.sorted.bed",shell=True,check=True)

# 4) per-chrom intersect -c
print("intersect -c vs gnomAD ...",flush=True)
chroms=sorted(set(l.split('\t')[0] for l in open("all.sorted.bed")))
with open("counts.bed","w") as fout:
    for c in chroms:
        snp=f"{SNPDIR}/{c}_snps.bed"
        if not os.path.exists(snp): continue
        cmd=f"awk '$1==\"{c}\"' all.sorted.bed | {BT} intersect -a stdin -b {snp} -c -sorted"
        fout.write(subprocess.run(cmd,shell=True,capture_output=True,text=True).stdout)
        print("  ",c,flush=True)

# 5) per-locus rate -> aggregate per (group,unit,nrep)
print("aggregate ...",flush=True)
rates=defaultdict(list)
for line in open("counts.bed"):
    a=line.rstrip("\n").split("\t")
    if len(a)<5: continue
    grp,unit,nrep=a[3].split("|"); L=int(a[2])-int(a[1]); cnt=int(a[4])
    if L>0: rates[(grp,unit,int(nrep))].append(cnt/L)
rows=[]
for (grp,unit,nrep),v in rates.items():
    v=np.array(v); rows.append(dict(group=grp,unit=unit,nrep=nrep,n=len(v),
        rate=v.mean(), sem=v.std(ddof=1)/np.sqrt(len(v)) if len(v)>1 else 0.0))
agg=pd.DataFrame(rows).sort_values(["unit","group","nrep"])
agg.to_csv("str_unit_mutrate.csv",index=False)

# 6) plot top-N units by total (det+und) loci
tot=defaultdict(int)
for (grp,unit,nrep),v in rates.items():
    if grp in ("detected","undetected"): tot[unit]+=len(v)
top=[u for u,_ in sorted(tot.items(),key=lambda x:-x[1])[:TOPN]]
print("top units:",top,flush=True)

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
COL={"detected":"#d62728","undetected":"#1f77b4","control":"0.45"}
ncol=4; nrow=int(np.ceil(len(top)/ncol))
fig,axes=plt.subplots(nrow,ncol,figsize=(3.4*ncol,2.7*nrow))
axes=np.array(axes).reshape(-1)
for ax,unit in zip(axes,top):
    for grp in ("detected","undetected","control"):
        sub=agg[(agg.group==grp)&(agg.unit==unit)&(agg.n>=MINN)].sort_values("nrep")
        if not len(sub): continue
        ls="--" if grp=="control" else "-"
        ax.errorbar(sub["nrep"],sub["rate"],yerr=sub["sem"],fmt="o"+("--" if grp=="control" else "-"),
                    color=COL[grp],ms=3,lw=1.3,capsize=2)
    ax.set_title(f"{unit}  (n={tot[unit]:,})",fontsize=9)
    ax.set_xlabel("number of repeats",fontsize=8); ax.tick_params(labelsize=7)
    ax.set_ylabel("SNPs / bp",fontsize=8)
for ax in axes[len(top):]: ax.axis("off")
h=[Line2D([0],[0],color=COL[g],ls=("--" if g=="control" else "-"),marker="o",ms=4,label=g) for g in ("detected","undetected","control")]
fig.legend(handles=h,fontsize=9,loc="upper right",ncol=3,frameon=False,bbox_to_anchor=(0.99,1.0))
fig.suptitle("STR mutation rate vs number of repeats detected / undetected / control",fontsize=11,y=1.0)
fig.tight_layout(rect=[0,0,1,0.96])
fig.savefig("str_unit_mutrate.png",dpi=180,bbox_inches="tight")
fig.savefig("str_unit_mutrate.pdf",bbox_inches="tight")
print("saved str_unit_mutrate.{png,pdf} in",WORK,flush=True)
