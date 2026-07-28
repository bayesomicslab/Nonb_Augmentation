#!/usr/bin/env python3
"""STR-by-kmer preparation, stratified by k-mer; mirrors the Python half of the
STR notebook's IWTomics section.
Reuses the per-kmer intersects already in results/str_bykmer_run/.
Outputs (in str_bykmer_run/):
  str_bykmer_freq.csv               tidy per-bin freq: region,kmer,x,detected,undetected,control
  mat_{det,und}_k{K}_{L,C,R}.txt.gz integer count matrices (loci x bins) for official R IWTomics
"""
import os, numpy as np, pandas as pd

WORK=os.path.join(os.environ.get("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "str_bykmer_run")
os.chdir(WORK)
MAX_FLANK=100; LOCUS_BINS=180; KMERS=[1,2,3,4,5,6]
scaledposLocus=[round(i,4) for i in np.arange(1,LOCUS_BINS*2,2.0)/(LOCUS_BINS*2)]

def kmap_from_bed(bed):
    return np.array([int(l.split("\t")[3].split("|")[1]) for l in open(bed)], dtype=int)

def build_matrices(intersect_path, n_loci):
    L=np.zeros((n_loci,MAX_FLANK),np.int32); C=np.zeros((n_loci,LOCUS_BINS),np.int32); R=np.zeros((n_loci,MAX_FLANK),np.int32)
    seen=set()
    with open(intersect_path) as fh:
        for line in fh:
            f=line.rstrip("\n").split("\t")
            if len(f)<8: continue
            if len(f)>=6 and f[5]==".": continue
            lid=int(f[3]); ps=int(f[1]); pe=int(f[2]); ptype=f[4]; snp=int(f[6])
            key=(lid,ptype,snp)
            if key in seen: continue
            seen.add(key); plen=pe-ps
            if ptype=="LeftFlank":
                idx=-(pe-snp)+MAX_FLANK
                if 0<=idx<MAX_FLANK: L[lid,idx]+=1
            elif ptype=="RightFlank":
                idx=snp-ps
                if 0<=idx<MAX_FLANK: R[lid,idx]+=1
            elif ptype=="Locus" and plen>0:
                fsf=(snp-ps)/plen; fef=(snp+1-ps)/plen
                for i,sp in enumerate(scaledposLocus):
                    if fsf<=sp<fef: C[lid,i]+=1
    return {"L":L,"C":C,"R":R}

print("kmer maps ...",flush=True)
km={"det":kmap_from_bed("det.nooverlap.bed"),
    "und":kmap_from_bed("und.nooverlap.bed"),
    "ctrl":kmap_from_bed("ctrl.nooverlap.bed")}
print("build matrices ...",flush=True)
mats={"det":build_matrices("det.intersect",len(km["det"])),
      "und":build_matrices("und.intersect",len(km["und"])),
      "ctrl":build_matrices("ctrl.intersect",len(km["ctrl"]))}

REG_X={"L":np.arange(-MAX_FLANK,0),"C":np.arange(1,LOCUS_BINS+1),"R":np.arange(1,MAX_FLANK+1)}
rows=[]
for k in KMERS:
    for reg in ("L","C","R"):
        d=mats["det"][reg][km["det"]==k]; u=mats["und"][reg][km["und"]==k]; c=mats["ctrl"][reg][km["ctrl"]==k]
        df=mats["det"][reg].shape[1]
        det_f=d.sum(0)/len(d) if len(d) else np.full(df,np.nan)
        und_f=u.sum(0)/len(u) if len(u) else np.full(df,np.nan)
        ctl_f=c.sum(0)/len(c) if len(c) else np.full(df,np.nan)
        for x,a,b,cc in zip(REG_X[reg],det_f,und_f,ctl_f):
            rows.append((reg,k,int(x),a,b,cc))
        np.savetxt(f"mat_det_k{k}_{reg}.txt.gz", d, fmt="%d")
        np.savetxt(f"mat_und_k{k}_{reg}.txt.gz", u, fmt="%d")
        print(f"  k{k} {reg}: det {len(d)}, und {len(u)}, ctrl {len(c)} loci",flush=True)

pd.DataFrame(rows,columns=["region","kmer","x","detected","undetected","control"]).to_csv("str_bykmer_freq.csv",index=False)
print("wrote str_bykmer_freq.csv + mat_*.txt.gz",flush=True)
