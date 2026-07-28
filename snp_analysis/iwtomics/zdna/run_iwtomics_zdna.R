# IWTomics on Z-DNA (model-detected)
# Run with Rscript, or block by block in an interactive session.
# Reference: IWTomics 1.36.0 manual
#   https://bioconductor.org/packages/release/bioc/manuals/IWTomics/man/IWTomics.pdf
#
# Inputs (already produced by the notebook through cell §3.6):
#   $SNP_AF_BASE/results/zdna_notebook_detected_run/
#     NCNRZDNA.gnomAD.intersect           (motif x SNP, bedtools -loj output)
#     NCNRZDNA.gnomAD.shuffled.intersect  (control x SNP)
#     split_counts.txt                    (N_MOTIF, N_CTRL)
#
# Z-DNA has 3 sub-parts per locus: LeftFlank | Locus | RightFlank.
# (No spacer / repeat partition like DR / IR / MR.)
#
# Outputs (saved into the same workdir so a later plot cell can `load()` them):
#   ZDNA_motif_vs_control_IWTomics_left_0.01.RData
#   ZDNA_motif_vs_control_IWTomics_locus_0.01.RData
#   ZDNA_motif_vs_control_IWTomics_right_0.01.RData


# load packages, verify versions
library(IWTomics)
library(GenomicRanges)
cat("IWTomics version:        ", as.character(packageVersion("IWTomics")), "\n")
cat("GenomicRanges version:   ", as.character(packageVersion("GenomicRanges")), "\n")
cat("R version:               ", R.version.string, "\n")


# paths + IWTomics parameters
WORKDIR <- file.path(Sys.getenv("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "zdna_notebook_detected_run")
setwd(WORKDIR)
cat("workdir:", getwd(), "\n")

ALPHA      <- 0.01      # significance threshold (matches DR pipeline _0.01)
N_PERM     <- 1000      # IWTomicsTest default B (manual p.12)
ALIGN      <- "center"  # all loci have identical width => center alignment.
                        # Not "scale": manual p.13 says "scale" needs smooth() first.
MAX_FLANK  <- 100       # bp on each side
LOCUS_BINS <- 180       # locus normalized to 180 fractional bins

# Locus uses the same scaled grid as the Python pipeline (cell 11):
#   scaledposLocus = round((1, 3, 5, ..., 359) / 360, 4)
SCALED_LOCUS <- round(seq(1, LOCUS_BINS * 2 - 1, by = 2) / (LOCUS_BINS * 2), 4)


# helpers

# Read split_counts.txt -> named integer vector c(N_MOTIF=..., N_CTRL=...)
read_counts <- function(path = "split_counts.txt") {
  d <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  out <- as.integer(d[[2]]); names(out) <- d[[1]]
  out
}

# Build per-locus matrices from a bedtools -loj intersect file.
# Returns list(left = n_loci x MAX_FLANK,
#              locus = n_loci x LOCUS_BINS,
#              right = n_loci x MAX_FLANK).
# Mirrors _build_matrices() in zdna_snp_frequency_notebook_detected.ipynb cell 11.
build_matrices <- function(intersect_path, n_loci) {
  left  <- matrix(0, nrow = n_loci, ncol = MAX_FLANK)
  locus <- matrix(0, nrow = n_loci, ncol = LOCUS_BINS)
  right <- matrix(0, nrow = n_loci, ncol = MAX_FLANK)
  cat("reading", intersect_path, "...\n")

  # Read whole file (60 MB) — fine for these sizes.
  d <- read.table(intersect_path, sep = "\t", header = FALSE,
                  stringsAsFactors = FALSE, comment.char = "",
                  quote = "", colClasses = c(
                    V1 = "character", V2 = "integer", V3 = "integer",
                    V4 = "integer",   V5 = "character",
                    V6 = "character", V7 = "integer", V8 = "integer"))
  cat("  ", nrow(d), "rows total\n")

  # Drop -loj no-overlap rows (V6 == ".")
  d <- d[d$V6 != ".", ]
  cat("  ", nrow(d), "rows after dropping no-overlap\n")

  # Dedup on (locus_id, part_type, snp_start) — same as Python `seen` set
  key <- paste(d$V4, d$V5, d$V7, sep = "|")
  d <- d[!duplicated(key), ]
  cat("  ", nrow(d), "rows after dedup\n")

  ps        <- d$V2
  pe        <- d$V3
  lid       <- d$V4 + 1L           # 0-indexed -> 1-indexed for R
  ptype     <- d$V5
  snp_start <- d$V7
  plen      <- pe - ps

  # LeftFlank: coord = -(pe - snp_start), idx = coord + MAX_FLANK in [0, MAX_FLANK)
  is_l <- ptype == "LeftFlank"
  if (any(is_l)) {
    coord <- -(pe[is_l] - snp_start[is_l])
    idx   <- coord + MAX_FLANK + 1L              # 1-indexed column
    keep  <- idx >= 1L & idx <= MAX_FLANK
    rows  <- lid[is_l][keep]; cols <- idx[keep]
    for (k in seq_along(rows)) left[rows[k], cols[k]] <- left[rows[k], cols[k]] + 1
  }

  # RightFlank: coord = snp_start - ps, in [0, MAX_FLANK)
  is_r <- ptype == "RightFlank"
  if (any(is_r)) {
    coord <- snp_start[is_r] - ps[is_r]
    idx   <- coord + 1L                          # 1-indexed column
    keep  <- idx >= 1L & idx <= MAX_FLANK
    rows  <- lid[is_r][keep]; cols <- idx[keep]
    for (k in seq_along(rows)) right[rows[k], cols[k]] <- right[rows[k], cols[k]] + 1
  }

  # Locus: each SNP at fractional [fs, fe) bumps every SCALED_LOCUS bin in that range
  is_c <- ptype == "Locus" & plen > 0
  if (any(is_c)) {
    sub_lid  <- lid[is_c]
    sub_ps   <- ps[is_c]; sub_pe <- pe[is_c]; sub_pl <- plen[is_c]; sub_ss <- snp_start[is_c]
    fs <- (sub_ss     - sub_ps) / sub_pl
    fe <- (sub_ss + 1 - sub_ps) / sub_pl
    for (k in seq_along(sub_lid)) {
      hits <- which(SCALED_LOCUS >= fs[k] & SCALED_LOCUS < fe[k])
      if (length(hits)) locus[sub_lid[k], hits] <- locus[sub_lid[k], hits] + 1
    }
  }
  list(left = left, locus = locus, right = right)
}


# load counts and build matrices for motif & control
counts  <- read_counts("split_counts.txt")
N_MOTIF <- counts["N_MOTIF"]
N_CTRL  <- counts["N_CTRL"]
cat("N_MOTIF =", N_MOTIF, "  N_CTRL =", N_CTRL, "\n")

mot <- build_matrices("NCNRZDNA.gnomAD.intersect",          N_MOTIF)
ctl <- build_matrices("NCNRZDNA.gnomAD.shuffled.intersect", N_CTRL)

cat("\nmatrix shapes (motif):\n")
cat("  left  :", dim(mot$left),  "\n")
cat("  locus :", dim(mot$locus), "\n")
cat("  right :", dim(mot$right), "\n")
cat("matrix shapes (control):\n")
cat("  left  :", dim(ctl$left),  "\n")
cat("  locus :", dim(ctl$locus), "\n")
cat("  right :", dim(ctl$right), "\n")


# helper to run IWTomics on one region
# IWTomics features matrix wants positions x loci (manual p.8: "columns of
# aligned feature measurements"). build_matrices() returns loci x bins, so
# transpose before handing to IWTomics.
run_iwt <- function(m_loci_x_bins, c_loci_x_bins, label) {
  m <- t(m_loci_x_bins)   # bins x loci
  c <- t(c_loci_x_bins)
  nbins  <- nrow(m)
  n_mot  <- ncol(m); n_ctrl <- ncol(c)

  # Dummy GRanges so each locus has identical width = nbins (center alignment OK).
  gr_mot  <- GRanges("chr1",
                     IRanges(start = seq(1, by = nbins, length.out = n_mot),
                             width = nbins))
  gr_ctrl <- GRanges("chr1",
                     IRanges(start = seq(1, by = nbins, length.out = n_ctrl),
                             width = nbins))
  rf <- IWTomicsData(
    GRangesList(motif = gr_mot, control = gr_ctrl),
    list(snp = list(motif = m, control = c)),
    alignment     = ALIGN,
    id_regions    = c("motif", "control"),
    name_regions  = c("motif", "control"),
    id_features   = "snp",
    name_features = "snp"
  )
  res <- IWTomicsTest(rf,
                      id_region1 = "motif", id_region2 = "control",
                      id_features_subset = "snp",
                      statistics = "mean", B = N_PERM)
  adj <- adjusted_pval(res)[[1]][[1]]
  cat(sprintf("%-12s  %4d bins, %4d motif, %4d control, %4d significant\n",
              label, nbins, n_mot, n_ctrl, sum(adj < ALPHA)))
  list(adj = as.numeric(adj), nbins = nbins)
}


# run all three regions and save .RData
left  <- run_iwt(mot$left,  ctl$left,  "left flank")
locus <- run_iwt(mot$locus, ctl$locus, "locus")
right <- run_iwt(mot$right, ctl$right, "right flank")

adj_pvalue_curve_left <- data.frame(
  snp_coord   = seq(-left$nbins, -1),
  adj_pvalue  = left$adj,
  significant = left$adj < ALPHA
)
save(adj_pvalue_curve_left,
     file = "ZDNA_motif_vs_control_IWTomics_left_0.01.RData")
cat("wrote ZDNA_motif_vs_control_IWTomics_left_0.01.RData\n")

adj_pvalue_curve_locus <- data.frame(
  bin         = seq_len(locus$nbins),
  adj_pvalue  = locus$adj,
  significant = locus$adj < ALPHA
)
save(adj_pvalue_curve_locus,
     file = "ZDNA_motif_vs_control_IWTomics_locus_0.01.RData")
cat("wrote ZDNA_motif_vs_control_IWTomics_locus_0.01.RData\n")

adj_pvalue_curve_right <- data.frame(
  snp_coord   = seq_len(right$nbins),
  adj_pvalue  = right$adj,
  significant = right$adj < ALPHA
)
save(adj_pvalue_curve_right,
     file = "ZDNA_motif_vs_control_IWTomics_right_0.01.RData")
cat("wrote ZDNA_motif_vs_control_IWTomics_right_0.01.RData\n")


# verify outputs + summary
files <- list.files(WORKDIR, pattern = "ZDNA.*IWTomics.*\\.RData$",
                    full.names = TRUE)
cat("\nIWTomics outputs in workdir:\n")
print(file.info(files)[, c("size", "mtime")])

cat("\nSummary of significant bins (alpha =", ALPHA, "):\n")
cat(sprintf("  left  : %d/%d\n",
            sum(adj_pvalue_curve_left$significant),  nrow(adj_pvalue_curve_left)))
cat(sprintf("  locus : %d/%d\n",
            sum(adj_pvalue_curve_locus$significant), nrow(adj_pvalue_curve_locus)))
cat(sprintf("  right : %d/%d\n",
            sum(adj_pvalue_curve_right$significant), nrow(adj_pvalue_curve_right)))
cat("
Done.
")
