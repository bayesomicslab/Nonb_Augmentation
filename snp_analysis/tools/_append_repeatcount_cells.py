#!/usr/bin/env python3
"""Append the Section 7 repeat-count cells to the STR notebook."""
import os
import json

NB = os.environ.get("STR_NOTEBOOK",
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                 "notebooks", "str_snp_frequency_notebook_detected.ipynb"))

md = """## 7. STR SNP frequency vs repeat count (faceted by k-mer)

Non-positional view of STR variation. Each locus is collapsed to a scalar SNP rate and
grouped by **k-mer length** (= sum of `composition` counts in `features`) and **number of
repeats** (= `N` in `xN+M`). Because k-mers occupy non-overlapping repeat-count ranges
(1-mer 10–19, 4-mer 3–4, …), a shared x-axis is meaningless, so the panels are **faceted by k-mer** with
a **per-facet y-scale**.

- x = number of repeats; y = pooled SNP frequency (ΣSNP/Σbp); color = detected / undetected;
  grey dashed = MR-style matched control.
- control = 1:1 length-matched `bedtools shuffle` (excl = STR ∪ gaps ∪ NCNR), inherits each
  locus's k-mer/repeat-count; empirically flat at ~0.23 SNPs/bp.
- significance: per-(k-mer, repeat-count) permutation test of pooled SNP rate, real vs matched
  control (N=1000, α=0.01); `*` marks significant points. Groups with < `RC_NMIN` loci hidden.

Reuses `locuschoice` (5.4), `STR_BED` (5.1), `TEST_CSV` (5.8), and env from the setup cell."""

build = r"""# 7.1-7.4  Build labeled motif beds (detected/undetected) + MR-style matched control
import csv as _csv
_csv.field_size_limit(10**7)

RC_MAXLEN = 100      # motif length threshold (bp)
RC_NMIN   = 50       # min loci per (group,kmer,nrep) to plot / test
RC_NPERM  = 1000
RC_ALPHA  = 0.01

def _rc_period(comp):
    return sum(int(re.match(r'(\d+)([ACGT])', t).group(1)) for t in comp.split('/'))

# 7.1 feature map from test CSV: (chr, mstart_1based, mend) -> (kmer, nrep)
rc_feat = {}
with open(TEST_CSV) as fh:
    r = _csv.reader(fh); h = next(r)
    ci=h.index('chr'); li=h.index('label'); fi=h.index('features')
    si=h.index('motif_start'); ei=h.index('motif_end')
    for row in r:
        if row[li] != 'Short_Tandem_Repeat' or not row[si] or not row[ei]:
            continue
        d={}; xc=None
        for p in row[fi].split(';'):
            if '=' in p: k,v=p.split('=',1); d[k]=v
            elif re.match(r'x\d+\+\d+', p): xc=p
        if 'composition' not in d or xc is None: continue
        rc_feat[(row[ci], int(float(row[si])), int(float(row[ei])))] = (
            _rc_period(d['composition']), int(re.match(r'x(\d+)\+', xc).group(1)))
print(f'test STR loci with features: {len(rc_feat)}')

# 7.2 detected set (1-based) from detected.bed
rc_det = set()
_dd = pd.read_csv(STR_BED, sep='\t').dropna(subset=['motif_start','motif_end'])
for c,s,e in zip(_dd['chr'], _dd['motif_start'].astype(int), _dd['motif_end'].astype(int)):
    rc_det.add((c,s,e))
print(f'detected loci: {len(rc_det)}')

# 7.3 labeled motif bed (0-based, name=group|kmer|nrep, <=100bp) + all-STR exclusion bed
with open('rc_loci.bed','w') as f:
    for (c,s1,e),(kmer,nrep) in rc_feat.items():
        start=s1-1; L=e-start
        if L<=0 or L>RC_MAXLEN: continue
        grp = 'detected' if (c,s1,e) in rc_det else 'undetected'
        f.write(f'{c}\t{start}\t{e}\t{grp}|{kmer}|{nrep}\n')
with open('rc_str_all.bed','w') as f:
    for (c,s1,e) in rc_feat:
        f.write(f'{c}\t{s1-1}\t{e}\n')
subprocess.run('bedtools intersect -a rc_loci.bed -b NCNR.hg38.bed -v '
               '| sort -k1,1 -k2,2n > rc_loci.ncnr.bed', shell=True, check=True)
locuschoice('rc_loci.ncnr.bed', 'rc_loci.nooverlap.bed', seed=42)

# 7.4 MR-style matched control: shuffle (name preserved), excl = STR ∪ gaps ∪ NCNR
subprocess.run('cat rc_str_all.bed "$GAPS_BED" NCNR.hg38.bed | cut -f1-3 '
               '| sort -k1,1 -k2,2n | bedtools merge -i stdin > rc_exclude.bed',
               shell=True, check=True)
subprocess.run('bedtools shuffle -i rc_loci.ncnr.bed -excl rc_exclude.bed '
               '-g "$CHROM_SIZES" -chrom -f 0 -noOverlapping -seed 42 '
               "| awk 'BEGIN{OFS=\"\\t\"}{split($4,a,\"|\"); print $1,$2,$3,\"control|\"a[2]\"|\"a[3]}' "
               '| sort -k1,1 -k2,2n > rc_ctrl.bed', shell=True, check=True)
locuschoice('rc_ctrl.bed', 'rc_ctrl.nooverlap.bed', seed=43)"""

count = r"""# 7.5  Per-locus SNP counts via bedtools intersect -c (per chrom, -sorted)
def rc_intersect_count(in_bed, out_path):
    snp_dir = os.environ['SNP_DIR']
    chroms = sorted(set(l.split('\t')[0] for l in open(in_bed)))
    with open(out_path,'w') as fout:
        for c in chroms:
            snp = f'{snp_dir}/{c}_snps.bed'
            if not os.path.exists(snp):
                print('  skip', c); continue
            cmd = f"awk '$1==\"{c}\"' {in_bed} | bedtools intersect -a stdin -b {snp} -c -sorted"
            fout.write(subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout)
    print(f'{in_bed} -> {out_path}')

rc_intersect_count('rc_loci.nooverlap.bed', 'rc_loci.counts.bed')
rc_intersect_count('rc_ctrl.nooverlap.bed', 'rc_ctrl.counts.bed')"""

agg = r"""# 7.6  Aggregate pooled SNPs/bp + per-group permutation test vs matched control
from collections import defaultdict
rc_loci = defaultdict(list)   # (group,kmer,nrep) -> list of (snp, bp)
def _rc_load(path):
    for line in open(path):
        a=line.rstrip('\n').split('\t')
        if len(a)<5: continue
        grp,kmer,nrep=a[3].split('|')
        rc_loci[(grp,int(kmer),int(nrep))].append((int(a[4]), int(a[2])-int(a[1])))
_rc_load('rc_loci.counts.bed'); _rc_load('rc_ctrl.counts.bed')

def _rc_pooled(pairs):
    s=sum(x[0] for x in pairs); b=sum(x[1] for x in pairs)
    return s/b if b else np.nan

def _rc_perm(real, ctrl, nperm=RC_NPERM, seed=42):
    # two-sample permutation on the pooled SNPs/bp difference (the plotted metric)
    rng=np.random.default_rng(seed)
    allp=np.array(real+ctrl, dtype=float); n1=len(real); n=len(allp)
    obs=abs(_rc_pooled(real)-_rc_pooled(ctrl)); snp=allp[:,0]; bp=allp[:,1]; cnt=0
    for _ in range(nperm):
        idx=rng.permutation(n); a=idx[:n1]; b=idx[n1:]
        if abs(snp[a].sum()/bp[a].sum() - snp[b].sum()/bp[b].sum()) >= obs: cnt+=1
    return (cnt+1)/(nperm+1)

rc_df = pd.DataFrame([dict(group=g,kmer=k,nrep=n,n=len(p),
                           snp=sum(x[0] for x in p), bp=sum(x[1] for x in p),
                           freq=_rc_pooled(p)) for (g,k,n),p in rc_loci.items()])
rc_df['sig']=False; rc_df['pval']=np.nan
for kmer in sorted(rc_df.kmer.unique()):
    for nrep in sorted(rc_df[rc_df.kmer==kmer].nrep.unique()):
        ctrl=rc_loci.get(('control',kmer,nrep),[])
        if len(ctrl)<RC_NMIN: continue
        for grp in ('detected','undetected'):
            real=rc_loci.get((grp,kmer,nrep),[])
            if len(real)<RC_NMIN: continue
            p=_rc_perm(real,ctrl); m=(rc_df.group==grp)&(rc_df.kmer==kmer)&(rc_df.nrep==nrep)
            rc_df.loc[m,'pval']=p; rc_df.loc[m,'sig']=p<RC_ALPHA
rc_df=rc_df.sort_values(['group','kmer','nrep']).reset_index(drop=True)
rc_df.to_csv('str_repeatcount_freq.csv', index=False)
rc_df"""

plot = r"""# 7.7  Faceted plot: per-k-mer panel, per-facet y-scale, significance marks
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
RC_COL={'detected':'#d62728','undetected':'#1f77b4','control':'0.45'}
plot_df=rc_df[rc_df.n>=RC_NMIN].copy()
kmers=sorted(plot_df[plot_df.group!='control'].kmer.unique())
fig,axes=plt.subplots(1,len(kmers),figsize=(3.4*len(kmers),3.8),sharey=False)
if len(kmers)==1: axes=[axes]
for ax,k in zip(axes,kmers):
    vals=[]
    for grp in ('detected','undetected','control'):
        sub=plot_df[(plot_df.group==grp)&(plot_df.kmer==k)].sort_values('nrep')
        if not len(sub): continue
        ls='--' if grp=='control' else '-'
        ax.plot(sub.nrep,sub.freq,ls,color=RC_COL[grp],marker='o',ms=4,lw=1.6)
        vals+=list(sub.freq.values)
        for _,rr in sub[sub.sig].iterrows():
            ax.annotate('*',(rr.nrep,rr.freq),textcoords='offset points',
                        xytext=(0,4),ha='center',fontsize=11,color=RC_COL[grp])
    if vals:
        ax.set_yscale('log'); ax.set_ylim(min(vals)*0.9, max(vals)*1.15)
    ax.set_title(f'{k}-mer',fontsize=11); ax.set_xlabel('number of repeats'); ax.tick_params(labelsize=8)
    xs=sorted(plot_df[plot_df.kmer==k].nrep.unique())
    if xs: ax.set_xticks(xs)
axes[0].set_ylabel('SNP frequency (SNPs/bp)')
handles=[Line2D([0],[0],color=RC_COL[g],ls=('--' if g=='control' else '-'),marker='o',ms=4,lw=1.6,label=g)
         for g in ('detected','undetected','control')]
handles.append(Line2D([0],[0],marker='*',color='0.3',ls='',ms=9,label=f'sig vs ctrl (p<{RC_ALPHA})'))
axes[-1].legend(handles=handles,fontsize=8,frameon=False,loc='best')
fig.suptitle('STR SNP frequency vs repeat count, by k-mer (detected/undetected vs matched control)',fontsize=12)
fig.tight_layout(rect=[0,0,1,0.94])
fig.savefig('str_repeatcount_snp_frequency.pdf', bbox_inches='tight')
fig.savefig('str_repeatcount_snp_frequency.png', dpi=200, bbox_inches='tight')
print('saved str_repeatcount_snp_frequency.{pdf,png}')
plt.show()"""

def code(src, cid):
    return {"cell_type":"code","execution_count":None,"id":cid,"metadata":{},
            "outputs":[],"source":[l+"\n" for l in src.split("\n")[:-1]]+[src.split("\n")[-1]]}
def mdc(src, cid):
    return {"cell_type":"markdown","id":cid,"metadata":{},
            "source":[l+"\n" for l in src.split("\n")[:-1]]+[src.split("\n")[-1]]}

nb=json.load(open(NB))
nb["cells"] += [
    mdc(md,   "rc_section_md"),
    code(build,"rc_build"),
    code(count,"rc_count"),
    code(agg,  "rc_agg"),
    code(plot, "rc_plot"),
]
json.dump(nb, open(NB,"w"), indent=1)
print("appended 5 cells; total now", len(nb["cells"]))
