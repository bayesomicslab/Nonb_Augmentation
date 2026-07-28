# Data


## Stage 1 — `NONB_DATA_ROOT`

```
Experiment_Data/                  ONT windows, per non-B type
  <Motif>_train.csv               columns: chr, strand, win_start50, win_end50,
  <Motif>_val.csv                 motif_start, motif_end, features, label,
  <Motif>_test.csv                forward_0..49, reverse_0..49
simulated_data/                   simulated translocation-time windows
  forward_{bdna,train_100,val_100,test_100}.npy, reverse_*.npy
  <Motif>_centered{,_train,_validation,_test}.csv
```

`<Motif>` is one of `A_Phased_Repeat`, `Direct_Repeat`, `G_Quadruplex_Motif`,
`Inverted_Repeat`, `Mirror_Repeat`, `Short_Tandem_Repeat`, `Z_DNA_Motif`. The `test.csv` files
are 3–4 GB each. Both feature blocks are the per-base translocation times of the 50 bp window on
the forward and reverse strand; `features` is a `key=value;…` string that carries the motif
`sequence` and, for STRs, the repeat `composition`.

Everything above comes from the ONT non-B DNA dataset of
[ONT-nonb-GoFAE-DND](https://github.com/bayesomicslab/ONT-nonb-GoFAE-DND).

## Stage 1 output, consumed by Stage 2 — `NONB_PU_RESULTS`

```
pu_results/<Motif>/<Motif>_detected.bed
```

Written by `nrcl/main_exp.py`. Columns: `chr, win_start, win_end, strand,
motif_start, motif_end, sequence, …`. Coordinates are 1-based inclusive; the notebooks convert
them to 0-based half-open on load.

## Stage 2 — `SNP_AF_BASE`

```
reference/
  hg38.fa, hg38.fa.fai         3.3 GB   UCSC hg38; only the G4 notebook reads it
                                        (bedtools getfasta for undetected motif sequences)
  hg38.chrom.sizes              11 KB   UCSC
annotations/
  hg38_gaps.bed                 28 KB   UCSC gap.txt.gz — assembly gaps, excluded from controls
  coding_regions.bed             5 MB   GENCODE v44 CDS, merged (212,645 intervals)
  repetitive_elements.bed      100 MB   UCSC RepeatMasker hg38.fa.out.gz, merged (4,128,792)
snp_data/gnomad/
  chr1_snps.bed … chrY_snps.bed 57 GB   gnomAD; columns chr, start, end, ref, alt, AF
results/                                created by the notebooks; one directory per run
```

`coding_regions.bed ∪ repetitive_elements.bed` is the **NCNR** mask (non-coding, non-repetitive):
motifs overlapping it are dropped, and it is part of the `bedtools shuffle -excl` union together
with all motifs of the type being analysed and `hg38_gaps.bed`.

This NCNR definition is deliberately narrower than the seven-component `CRgenome` used by
Guiblet et al. (which additionally masks exons ±100 bp, 5 kb upstream/downstream, phastCons
±100 bp and cCREs ±100 bp, covering ~80% of hg38 against ~52% here). The wider mask was
piloted on direct repeats, mirror repeats and Z-DNA and gave a noisier signal on a much smaller
locus set; those pilot notebooks are not included in this repo.

### Rebuilding the annotation tracks

```bash
# assembly gaps
zcat gap.txt.gz | cut -f2-4 | sort -k1,1 -k2,2n | bedtools merge > hg38_gaps.bed

# CDS
zcat gencode.v44.annotation.gtf.gz | awk '$3=="CDS"' \
  | awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5}' | sort -k1,1 -k2,2n | bedtools merge > coding_regions.bed

# repeats
zcat hg38.fa.out.gz | tail -n +4 \
  | awk 'BEGIN{OFS="\t"}{print $5,$6-1,$7}' | sort -k1,1 -k2,2n | bedtools merge > repetitive_elements.bed
```
