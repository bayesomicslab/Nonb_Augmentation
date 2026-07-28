#!/usr/bin/env python3
"""Per-locus variant rate, defined as in the other non-B notebooks:

    rate = (distinct positions carrying >=1 gnomAD variant) / locus length

A plain `bedtools intersect -c` counts multiallelic sites once per ALT record,
which pushes the rate above 1 (up to 118/bp). Variant type is not filtered: any
gnomAD record at a position counts once, matching the notebooks' dedup on
(locus, part, position). gnomAD *_snps.bed is position-sorted, so a single
streaming awk pass collapses duplicate positions.

Reuses the existing all.sorted.bed (loci plus matched control), leaving the
control unchanged.
Output: counts.dedup.bed (chr start end group|unit|nrep distinct_position_count)
"""
import os, subprocess
BASE=os.environ.get("SNP_AF_BASE", "/path/to/snp_af_analysis")
WORK=f"{BASE}/results/str_unit_mutrate"
SNPDIR=f"{BASE}/snp_data/gnomad"
BT=os.environ.get("BEDTOOLS", "bedtools")
ALL=f"{WORK}/all.sorted.bed"
OUT=f"{WORK}/counts.dedup.bed"

chroms=sorted(set(l.split('\t')[0] for l in open(ALL)))
# One record per distinct position; no variant-type filter.
AWK_SNPUNIQ=(r"""awk 'BEGIN{OFS="\t"} """
             r"""!($1==pc && $2==ps){print $1,$2,$3; pc=$1; ps=$2}'""")

open(OUT,"w").close()
for c in chroms:
    snp=f"{SNPDIR}/{c}_snps.bed"
    if not os.path.exists(snp):
        print("  missing:",snp,flush=True); continue
    # -a: this chrom's loci (already sorted in all.sorted.bed); -b: SNP-only unique positions
    cmd=(f"awk '$1==\"{c}\"' {ALL} | "
         f"{BT} intersect -a stdin -b <({AWK_SNPUNIQ} {snp}) -c -sorted >> {OUT}")
    subprocess.run(cmd,shell=True,check=True,executable="/bin/bash")
    print("  done",c,flush=True)
print("wrote",OUT,flush=True)
