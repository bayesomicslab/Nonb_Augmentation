#!/usr/bin/env python3
"""STR SNP frequency by k-mer (paper-style IWTomics positional plot).

Mirrors the APR/STR notebook processing (split LeftFlank|Locus|RightFlank ->
flanks -> matched control -> locuschoice -> per-chrom gnomAD intersect ->
per-bin matrices -> IWTomics), but stratified by repeat-unit length (k-mer).

Layout: one row per k-mer (color = k-mer); each row = upstream(bp) | scaled Locus |
downstream(bp); detected = solid, undetected = dashed (same k-mer color),
control = grey dotted; grey shading = IWTomics(detected vs undetected) significant bins.
"""
import os, re, subprocess, random, csv as _csv
from collections import defaultdict
import numpy as np, pandas as pd

# config
BASE   = os.environ.get("SNP_AF_BASE", "/path/to/snp_af_analysis")
WORK   = f"{BASE}/results/str_bykmer_run"
SNPDIR = f"{BASE}/snp_data/gnomad"
NCNR   = f"{BASE}/results/str_notebook_detected_run/NCNR.hg38.bed"
GAPS   = f"{BASE}/annotations/hg38_gaps.bed"
CHROMS = f"{BASE}/reference/hg38.chrom.sizes"
BT     = os.environ.get("BEDTOOLS", "bedtools")
TESTCSV= os.environ.get("NONB_EXPERIMENT_DATA", "/path/to/Experiment_Data") + "/Short_Tandem_Repeat_test.csv"
DETBED = os.environ.get("NONB_PU_RESULTS", "/path/to/pu_results") + "/Short_Tandem_Repeat/Short_Tandem_Repeat_detected.bed"

MAX_FLANK  = 100
LOCUS_BINS = 180
FLANK      = 100          # bp of flank kept
MAXLEN     = 100          # motif length threshold
N_PERM     = 1000
ALPHA      = 0.01
KMERS      = [1, 2, 3, 4, 5, 6]
os.makedirs(WORK, exist_ok=True); os.chdir(WORK)
_csv.field_size_limit(10**7); random.seed(42)

scaledposLocus = [round(i,4) for i in np.arange(1, LOCUS_BINS*2, 2.0)/(LOCUS_BINS*2)]

def sh(cmd): subprocess.run(cmd, shell=True, check=True)

def locuschoice(in_bed, out_bed, seed=42):
    rng=random.Random(seed); last_chrom="chr0"; last_end=0; tochoose=[]; nkept=0
    with open(in_bed) as fi, open(out_bed,"w") as fo:
        for line in fi:
            arr=line.rstrip("\n").split("\t"); chrom=arr[0]; ws=int(arr[1]); we=int(arr[2])
            if chrom==last_chrom:
                if ws<last_end: tochoose.append(arr); last_end=we
                else:
                    fo.write("\t".join(tochoose[rng.randint(0,len(tochoose)-1)])+"\n"); nkept+=1
                    tochoose=[arr]; last_end=we
            else:
                if tochoose: fo.write("\t".join(tochoose[rng.randint(0,len(tochoose)-1)])+"\n"); nkept+=1
                tochoose=[arr]; last_chrom=chrom; last_end=we
        if tochoose: fo.write("\t".join(tochoose[rng.randint(0,len(tochoose)-1)])+"\n"); nkept+=1
    print(f"  locuschoice {in_bed} -> {out_bed}: kept {nkept}"); return nkept

def kmer_of(comp): return sum(int(re.match(r'(\d+)([ACGT])',t).group(1)) for t in comp.split('/'))

# 1. feature map + detected set
print("parsing features ...", flush=True)
feat={}
with open(TESTCSV) as fh:
    r=_csv.reader(fh); h=next(r)
    ci=h.index('chr'); li=h.index('label'); fi=h.index('features'); si=h.index('motif_start'); ei=h.index('motif_end')
    for row in r:
        if row[li]!='Short_Tandem_Repeat' or not row[si] or not row[ei]: continue
        d={}
        for p in row[fi].split(';'):
            if '=' in p: k,v=p.split('=',1); d[k]=v
        if 'composition' not in d: continue
        feat[(row[ci],int(float(row[si])),int(float(row[ei])))]=kmer_of(d['composition'])
det=set()
dd=pd.read_csv(DETBED,sep='\t').dropna(subset=['motif_start','motif_end'])
for c,s,e in zip(dd['chr'],dd['motif_start'].astype(int),dd['motif_end'].astype(int)): det.add((c,s,e))
print(f"  test STR loci {len(feat)}, detected {len(det)}", flush=True)

# 2. labeled motif beds (0-based, name=group|kmer), NCNR, <=100bp
with open("motif.det.bed","w") as fd, open("motif.und.bed","w") as fu, open("motif.all.bed","w") as fa:
    for (c,s1,e),k in feat.items():
        st=s1-1; L=e-st
        if L<=0 or L>MAXLEN: continue
        grp='detected' if (c,s1,e) in det else 'undetected'
        line=f"{c}\t{st}\t{e}\t{grp}|{k}\n"
        (fd if grp=='detected' else fu).write(line)
        fa.write(f"{c}\t{st}\t{e}\tall|{k}\n")
for nm in ("det","und","all"):
    sh(f"{BT} intersect -a motif.{nm}.bed -b {NCNR} -v | sort -k1,1 -k2,2n > motif.{nm}.ncnr.bed")

# 3. matched control = shuffle of all test STR (inherits kmer)
sh(f"cut -f1-3 motif.all.bed | sort -k1,1 -k2,2n > str_all_coords.bed")
sh(f"cat str_all_coords.bed {GAPS} {NCNR} | cut -f1-3 | sort -k1,1 -k2,2n | {BT} merge -i stdin > exclude.bed")
sh(f"{BT} shuffle -i motif.all.ncnr.bed -excl exclude.bed -g {CHROMS} -chrom -f 0 -noOverlapping -seed 42 "
   "| awk 'BEGIN{OFS=\"\\t\"}{split($4,a,\"|\"); print $1,$2,$3,\"control|\"a[2]}' "
   "| sort -k1,1 -k2,2n > motif.ctrl.ncnr.bed")

# 4. flank + locuschoice + split (record locus_id -> kmer)
def flank_choose_split(motif_ncnr, tag, seed):
    # build flanked: chrom f_start f_end name l_start l_end   (clamp f_start at 0)
    with open(motif_ncnr) as fi, open(f"{tag}.flanked.bed","w") as fo:
        for line in fi:
            c,s,e,name=line.rstrip("\n").split("\t")
            s=int(s); e=int(e); fs=max(0,s-FLANK); fe=e+FLANK
            fo.write(f"{c}\t{fs}\t{fe}\t{name}\t{s}\t{e}\n")
    sh(f"sort -k1,1 -k2,2n {tag}.flanked.bed > {tag}.flanked.sorted.bed")
    locuschoice(f"{tag}.flanked.sorted.bed", f"{tag}.nooverlap.bed", seed=seed)
    kmap=[]
    with open(f"{tag}.nooverlap.bed") as fi, open(f"{tag}.split.bed","w") as fo:
        for line in fi:
            c,fs,fe,name,ls,le=line.rstrip("\n").split("\t")
            lid=len(kmap); kmap.append(int(name.split("|")[1]))
            fo.write(f"{c}\t{fs}\t{ls}\t{lid}\tLeftFlank\n")
            fo.write(f"{c}\t{ls}\t{le}\t{lid}\tLocus\n")
            fo.write(f"{c}\t{le}\t{fe}\t{lid}\tRightFlank\n")
    sh(f"sort -k1,1 -k2,2n {tag}.split.bed > {tag}.split.sorted.bed && mv {tag}.split.sorted.bed {tag}.split.bed")
    return np.array(kmap, dtype=int)

print("flank/choose/split ...", flush=True)
kmap_det = flank_choose_split("motif.det.ncnr.bed",  "det",  seed=42)
kmap_und = flank_choose_split("motif.und.ncnr.bed",  "und",  seed=42)
kmap_ctl = flank_choose_split("motif.ctrl.ncnr.bed", "ctrl", seed=43)
print(f"  loci det/und/ctrl = {len(kmap_det)}/{len(kmap_und)}/{len(kmap_ctl)}", flush=True)

# 5. per-chrom gnomAD intersect (-loj)
def intersect(tag):
    chroms=sorted(set(l.split("\t")[0] for l in open(f"{tag}.split.bed")))
    with open(f"{tag}.intersect","w") as fo:
        for c in chroms:
            snp=f"{SNPDIR}/{c}_snps.bed"
            if not os.path.exists(snp): continue
            cmd=f"awk '$1==\"{c}\"' {tag}.split.bed | {BT} intersect -wa -wb -a stdin -b {snp} -loj"
            fo.write(subprocess.run(cmd,shell=True,capture_output=True,text=True).stdout)
    print(f"  intersect {tag} done", flush=True)
intersect("det"); intersect("und"); intersect("ctrl")

# 6. matrices
def build_matrices(intersect_path, n_loci):
    left =np.zeros((n_loci,MAX_FLANK),np.float32)
    locus=np.zeros((n_loci,LOCUS_BINS),np.float32)
    right=np.zeros((n_loci,MAX_FLANK),np.float32)
    seen=set()
    with open(intersect_path) as fh:
        for line in fh:
            f=line.rstrip("\n").split("\t")
            if len(f)<8: continue
            lid=int(f[3]); ps=int(f[1]); pe=int(f[2]); ptype=f[4]
            if len(f)>=6 and f[5]==".": continue
            snp=int(f[6]); plen=pe-ps; key=(lid,ptype,snp)
            if key in seen: continue
            seen.add(key)
            if ptype=="LeftFlank":
                idx=-(pe-snp)+MAX_FLANK
                if 0<=idx<MAX_FLANK: left[lid,idx]+=1
            elif ptype=="RightFlank":
                co=snp-ps
                if 0<=co<MAX_FLANK: right[lid,co]+=1
            elif ptype=="Locus":
                if plen==0: continue
                fsf=(snp-ps)/plen; fef=(snp+1-ps)/plen
                for i,sp in enumerate(scaledposLocus):
                    if fsf<=sp<fef: locus[lid,i]+=1
    return left,locus,right

print("build matrices ...", flush=True)
det_L,det_C,det_R = build_matrices("det.intersect",  len(kmap_det))
und_L,und_C,und_R = build_matrices("und.intersect",  len(kmap_und))
ctl_L,ctl_C,ctl_R = build_matrices("ctrl.intersect", len(kmap_ctl))

# 7. IWTomics (det vs undet)
def iwtomics(d1,d2,n_perm=N_PERM,alpha=ALPHA,seed=42):
    n1,p=d1.shape; n2=d2.shape[0]; n=n1+n2; allm=np.vstack([d1,d2])
    obs=(d1.mean(0)-d2.mean(0))**2; perm=np.zeros((n_perm,p))
    rng=np.random.default_rng(seed)
    for b in range(n_perm):
        idx=rng.permutation(n); perm[b]=(allm[idx[:n1]].mean(0)-allm[idx[n1:]].mean(0))**2
    unadj=(np.sum(perm>=obs[None,:],0)+1)/(n_perm+1)
    oc=np.concatenate([[0],np.cumsum(obs)]); pc=np.hstack([np.zeros((n_perm,1)),np.cumsum(perm,1)])
    adj=unadj.copy()
    for k in range(2,p+1):
        ow=(oc[k:]-oc[:-k])/k; pw=(pc[:,k:]-pc[:,:-k])/k
        wp=(np.sum(pw>=ow[None,:],0)+1)/(n_perm+1)
        for j in range(len(wp)):
            for x in range(j,j+k):
                if wp[j]>adj[x]: adj[x]=wp[j]
    return adj<alpha

def freq(mat, kmap, k):
    m=mat[kmap==k]
    return m.sum(0)/len(m) if len(m) else np.full(mat.shape[1],np.nan)

# 8. plot: rows = k-mer
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec, matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
CMAP=plt.get_cmap("tab10")
COL_CTRL=(0,0,0,0.45); COL_SIG=(0.88,0.88,0.88)
left_x=np.arange(-MAX_FLANK,0); right_x=np.arange(1,MAX_FLANK+1)

def filt_runs(mask,minrun=5):
    out=mask.copy(); n=len(out); i=0
    while i<n:
        if out[i]:
            j=i
            while j<n and out[j]: j+=1
            if j-i<minrun: out[i:j]=False
            i=j
        else: i+=1
    return out

nrow=len(KMERS)
fig=plt.figure(figsize=(11,2.4*nrow))
outer=gridspec.GridSpec(nrow,1,hspace=0.55)
rows=[]
for ri,k in enumerate(KMERS):
    dL,dC,dR=freq(det_L,kmap_det,k),freq(det_C,kmap_det,k),freq(det_R,kmap_det,k)
    uL,uC,uR=freq(und_L,kmap_und,k),freq(und_C,kmap_und,k),freq(und_R,kmap_und,k)
    cL,cC,cR=freq(ctl_L,kmap_ctl,k),freq(ctl_C,kmap_ctl,k),freq(ctl_R,kmap_ctl,k)
    # significance det vs undet (per section)
    dmaskL=det_L[kmap_det==k]; umaskL=und_L[kmap_und==k]
    dmaskC=det_C[kmap_det==k]; umaskC=und_C[kmap_und==k]
    dmaskR=det_R[kmap_det==k]; umaskR=und_R[kmap_und==k]
    okL=len(dmaskL)>1 and len(umaskL)>1
    sigL=filt_runs(iwtomics(dmaskL,umaskL)) if okL else np.zeros(MAX_FLANK,bool)
    sigC=filt_runs(iwtomics(dmaskC,umaskC)) if okL else np.zeros(LOCUS_BINS,bool)
    sigR=filt_runs(iwtomics(dmaskR,umaskR)) if okL else np.zeros(MAX_FLANK,bool)
    allv=np.concatenate([v for v in [dL,dC,dR,uL,uC,uR,cL,cC,cR]])
    allv=allv[(allv>0)&np.isfinite(allv)]
    ylo=allv.min()*0.6 if len(allv) else 1e-3; yhi=allv.max()*1.6 if len(allv) else 1.0
    color=CMAP(ri)
    gs=gridspec.GridSpecFromSubplotSpec(1,3,subplot_spec=outer[ri],width_ratios=[100,180,100],wspace=0.05)
    ax0=fig.add_subplot(gs[0]); ax1=fig.add_subplot(gs[1]); ax2=fig.add_subplot(gs[2])
    def style(ax,sig,xs,lo,hi):
        cl=sig
        if cl.any(): ax.fill_between(xs,lo,hi,where=cl,color=COL_SIG,step="mid",zorder=0,lw=0)
        ax.set_yscale("log"); ax.set_ylim(lo,hi); ax.tick_params(labelsize=7)
    style(ax0,sigL,left_x,ylo,yhi)
    ax0.scatter(left_x,cL,s=5,color=COL_CTRL,zorder=2)
    ax0.plot(left_x,uL,color=color,ls="--",lw=1.1,zorder=3)
    ax0.plot(left_x,dL,color=color,ls="-",lw=1.3,zorder=4)
    ax0.set_xlim(-MAX_FLANK-0.5,-0.5); ax0.set_xticks([-100,-50])
    ax0.set_ylabel(f"{k}-mer\nSNP freq",fontsize=8); ax0.spines["top"].set_visible(False); ax0.spines["right"].set_visible(False)
    xb=np.arange(1,LOCUS_BINS+1); style(ax1,sigC,xb,ylo,yhi)
    ax1.scatter(xb,cC,s=4,color=COL_CTRL,zorder=2)
    ax1.plot(xb,uC,color=color,ls="--",lw=1.1,zorder=3)
    ax1.plot(xb,dC,color=color,ls="-",lw=1.3,zorder=4)
    ax1.set_xlim(0.5,LOCUS_BINS+0.5); ax1.set_xticks([]); ax1.yaxis.set_visible(False)
    for sp in ("left","right","top"): ax1.spines[sp].set_visible(False)
    if ri==0: ax1.set_title("upstream(bp) | scaled STR locus | downstream(bp)",fontsize=9)
    style(ax2,sigR,right_x,ylo,yhi)
    ax2.scatter(right_x,cR,s=5,color=COL_CTRL,zorder=2)
    ax2.plot(right_x,uR,color=color,ls="--",lw=1.1,zorder=3)
    ax2.plot(right_x,dR,color=color,ls="-",lw=1.3,zorder=4)
    ax2.set_xlim(0.5,MAX_FLANK+0.5); ax2.set_xticks([50,100]); ax2.yaxis.set_visible(False)
    ax2.spines["top"].set_visible(False); ax2.spines["left"].set_visible(False)
    rows.append((k,color))

handles=[Line2D([0],[0],color="0.3",ls="-",lw=1.4,label="detected (k-mer color)"),
         Line2D([0],[0],color="0.3",ls="--",lw=1.1,label="undetected"),
         Line2D([0],[0],marker="o",color="w",markerfacecolor=COL_CTRL,markersize=5,label="control"),
         Patch(facecolor=COL_SIG,label="sig (det vs undet)")]
fig.legend(handles=handles,fontsize=8,loc="upper right",frameon=False,bbox_to_anchor=(0.99,1.0))
fig.suptitle("STR SNP frequency along motif + flanks, by k-mer (detected / undetected / control)",fontsize=11,y=1.0)
fig.savefig("str_snp_frequency_by_kmer.png",dpi=200,bbox_inches="tight")
fig.savefig("str_snp_frequency_by_kmer.pdf",bbox_inches="tight")
print("saved str_snp_frequency_by_kmer.{png,pdf} in",WORK, flush=True)
