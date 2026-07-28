# IWTomics on A-Phased Repeats (model-detected)
# Run with Rscript, or block by block in an interactive session.
# Reference: IWTomics 1.36.0 manual
#   https://bioconductor.org/packages/release/bioc/manuals/IWTomics/man/IWTomics.pdf
#
# Inputs (already produced by the notebook through cell §6.6):
#   $SNP_AF_BASE/results/apr_notebook_detected_run/
#     NCNRAPhasedRepeats.gnomAD.intersect           (motif x SNP, bedtools -loj)
#     NCNRAPhasedRepeats.gnomAD.ctrl.intersect      (control x SNP)
#     split_counts.txt                              (N_MOTIF, N_CTRL)
#
# A-Phased Repeats have 7 sub-parts per locus:
#   LeftFlank | Tract1 | Spacer1 | Tract2 | Spacer2 | Tract3 | RightFlank
# Center = Tract1 + Spacer1 + Tract2 + Spacer2 + Tract3 (5 sub-parts).
#
# Outputs (saved into the same workdir so a later plot cell can `load()` them):
#   APhasedRepeats_motif_vs_control_IWTomics_left_0.01.RData
#   APhasedRepeats_motif_vs_control_IWTomics_center_0.01.RData
#   APhasedRepeats_motif_vs_control_IWTomics_right_0.01.RData


# load packages, verify versions
library(IWTomics)
library(GenomicRanges)
cat("IWTomics version:        ", as.character(packageVersion("IWTomics")), "\n")
cat("GenomicRanges version:   ", as.character(packageVersion("GenomicRanges")), "\n")
cat("R version:               ", R.version.string, "\n")


# paths + IWTomics parameters
WORKDIR <- file.path(Sys.getenv("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "apr_notebook_detected_run")
setwd(WORKDIR)
cat("workdir:", getwd(), "\n")

ALPHA      <- 0.01      # significance threshold
N_PERM     <- 1000      # IWTomicsTest default B (manual p.12)
ALIGN      <- "center"  # all loci within a sub-part have identical width.

MAX_FLANK   <- 100
TRACT_BINS  <- 40
SPACER_BINS <- 30

# scaled fractional grids — match Python pipeline (notebook cell 21)
SCALED_TRACT  <- round(seq(1, TRACT_BINS  * 2 - 1, by = 2) / (TRACT_BINS  * 2), 4)
SCALED_SPACER <- round(seq(1, SPACER_BINS * 2 - 1, by = 2) / (SPACER_BINS * 2), 4)


# helpers

# Read split_counts.txt -> named integer vector c(N_MOTIF=..., N_CTRL=...)
read_counts <- function(path = "split_counts.txt") {
  d <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  out <- as.integer(d[[2]]); names(out) <- d[[1]]
  out
}

# Build per-locus matrices from a bedtools -loj intersect file.
# Returns a named list:
#   left, right  : n_loci x MAX_FLANK
#   Tract1/2/3   : n_loci x TRACT_BINS
#   Spacer1/2    : n_loci x SPACER_BINS
# Mirrors _build_matrices() in cell 21 of apr_snp_frequency_notebook_detected.ipynb.
build_matrices <- function(intersect_path, n_loci) {
  out <- list(
    left    = matrix(0, nrow = n_loci, ncol = MAX_FLANK),
    Tract1  = matrix(0, nrow = n_loci, ncol = TRACT_BINS),
    Spacer1 = matrix(0, nrow = n_loci, ncol = SPACER_BINS),
    Tract2  = matrix(0, nrow = n_loci, ncol = TRACT_BINS),
    Spacer2 = matrix(0, nrow = n_loci, ncol = SPACER_BINS),
    Tract3  = matrix(0, nrow = n_loci, ncol = TRACT_BINS),
    right   = matrix(0, nrow = n_loci, ncol = MAX_FLANK)
  )
  cat("reading", intersect_path, "...\n")

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

  # LeftFlank: coord = -(pe - snp_start), 1-indexed col in [1, MAX_FLANK]
  is_l <- ptype == "LeftFlank"
  if (any(is_l)) {
    coord <- -(pe[is_l] - snp_start[is_l])
    idx   <- coord + MAX_FLANK + 1L
    keep  <- idx >= 1L & idx <= MAX_FLANK
    rows  <- lid[is_l][keep]; cols <- idx[keep]
    for (k in seq_along(rows)) out$left[rows[k], cols[k]] <- out$left[rows[k], cols[k]] + 1
  }

  # RightFlank: coord = snp_start - ps, in [0, MAX_FLANK)
  is_r <- ptype == "RightFlank"
  if (any(is_r)) {
    coord <- snp_start[is_r] - ps[is_r]
    idx   <- coord + 1L
    keep  <- idx >= 1L & idx <= MAX_FLANK
    rows  <- lid[is_r][keep]; cols <- idx[keep]
    for (k in seq_along(rows)) out$right[rows[k], cols[k]] <- out$right[rows[k], cols[k]] + 1
  }

  # Tracts and spacers: each SNP at fractional [fs, fe) bumps every grid bin
  for (pname in c("Tract1", "Tract2", "Tract3", "Spacer1", "Spacer2")) {
    is_p <- ptype == pname & plen > 0
    if (!any(is_p)) next
    grid <- if (startsWith(pname, "Tract")) SCALED_TRACT else SCALED_SPACER
    sub_lid <- lid[is_p]
    sub_ps  <- ps[is_p]; sub_pe <- pe[is_p]; sub_pl <- plen[is_p]; sub_ss <- snp_start[is_p]
    fs <- (sub_ss     - sub_ps) / sub_pl
    fe <- (sub_ss + 1 - sub_ps) / sub_pl
    for (k in seq_along(sub_lid)) {
      hits <- which(grid >= fs[k] & grid < fe[k])
      if (length(hits)) out[[pname]][sub_lid[k], hits] <- out[[pname]][sub_lid[k], hits] + 1
    }
  }
  out
}


# load counts and build matrices for motif & control
counts  <- read_counts("split_counts.txt")
N_MOTIF <- counts["N_MOTIF"]
N_CTRL  <- counts["N_CTRL"]
cat("N_MOTIF =", N_MOTIF, "  N_CTRL =", N_CTRL, "\n")

mot <- build_matrices("NCNRAPhasedRepeats.gnomAD.intersect",      N_MOTIF)
ctl <- build_matrices("NCNRAPhasedRepeats.gnomAD.ctrl.intersect", N_CTRL)

cat("\nmatrix shapes (motif):\n")
for (nm in names(mot)) cat(sprintf("  %-7s : %d x %d\n", nm, nrow(mot[[nm]]), ncol(mot[[nm]])))


# helper to run IWTomics on one region
# IWTomics features matrix wants positions x loci (manual p.8). build_matrices()
# returns loci x bins, so transpose before handing to IWTomics.
run_iwt <- function(m_loci_x_bins, c_loci_x_bins, label) {
  m <- t(m_loci_x_bins)
  c <- t(c_loci_x_bins)
  nbins  <- nrow(m)
  n_mot  <- ncol(m); n_ctrl <- ncol(c)

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


# run all 7 regions (left + 5 center sub-parts + right)
cat("\nrunning IWTomicsTest on left flank ...\n")
left  <- run_iwt(mot$left,    ctl$left,    "left flank")

cat("\nrunning IWTomicsTest on right flank ...\n")
right <- run_iwt(mot$right,   ctl$right,   "right flank")

cat("\nrunning IWTomicsTest on Tract1 / Spacer1 / Tract2 / Spacer2 / Tract3 ...\n")
T1 <- run_iwt(mot$Tract1,  ctl$Tract1,  "Tract1")
S1 <- run_iwt(mot$Spacer1, ctl$Spacer1, "Spacer1")
T2 <- run_iwt(mot$Tract2,  ctl$Tract2,  "Tract2")
S2 <- run_iwt(mot$Spacer2, ctl$Spacer2, "Spacer2")
T3 <- run_iwt(mot$Tract3,  ctl$Tract3,  "Tract3")


# save .RData files
adj_pvalue_curve_left <- data.frame(
  snp_coord   = seq(-left$nbins, -1),
  adj_pvalue  = left$adj,
  significant = left$adj < ALPHA
)
save(adj_pvalue_curve_left,
     file = "APhasedRepeats_motif_vs_control_IWTomics_left_0.01.RData")
cat("wrote APhasedRepeats_motif_vs_control_IWTomics_left_0.01.RData\n")

adj_pvalue_curve_right <- data.frame(
  snp_coord   = seq_len(right$nbins),
  adj_pvalue  = right$adj,
  significant = right$adj < ALPHA
)
save(adj_pvalue_curve_right,
     file = "APhasedRepeats_motif_vs_control_IWTomics_right_0.01.RData")
cat("wrote APhasedRepeats_motif_vs_control_IWTomics_right_0.01.RData\n")

# Combine the 5 center sub-parts into one data.frame with a "subpart" column
adj_pvalue_curve_center <- rbind(
  data.frame(subpart = "Tract1",  bin = seq_len(T1$nbins),
             adj_pvalue = T1$adj, significant = T1$adj < ALPHA),
  data.frame(subpart = "Spacer1", bin = seq_len(S1$nbins),
             adj_pvalue = S1$adj, significant = S1$adj < ALPHA),
  data.frame(subpart = "Tract2",  bin = seq_len(T2$nbins),
             adj_pvalue = T2$adj, significant = T2$adj < ALPHA),
  data.frame(subpart = "Spacer2", bin = seq_len(S2$nbins),
             adj_pvalue = S2$adj, significant = S2$adj < ALPHA),
  data.frame(subpart = "Tract3",  bin = seq_len(T3$nbins),
             adj_pvalue = T3$adj, significant = T3$adj < ALPHA)
)
save(adj_pvalue_curve_center,
     file = "APhasedRepeats_motif_vs_control_IWTomics_center_0.01.RData")
cat("wrote APhasedRepeats_motif_vs_control_IWTomics_center_0.01.RData\n")


# verify outputs + summary
files <- list.files(WORKDIR,
                    pattern = "APhasedRepeats.*IWTomics.*\\.RData$",
                    full.names = TRUE)
cat("\nIWTomics outputs in workdir:\n")
print(file.info(files)[, c("size", "mtime")])

cat("\nSummary of significant bins (alpha =", ALPHA, "):\n")
cat(sprintf("  left    : %d/%d\n",
            sum(adj_pvalue_curve_left$significant),  nrow(adj_pvalue_curve_left)))
for (sp in c("Tract1", "Spacer1", "Tract2", "Spacer2", "Tract3")) {
  rows <- adj_pvalue_curve_center[adj_pvalue_curve_center$subpart == sp, ]
  cat(sprintf("  %-7s : %d/%d\n", sp, sum(rows$significant), nrow(rows)))
}
cat(sprintf("  right   : %d/%d\n",
            sum(adj_pvalue_curve_right$significant), nrow(adj_pvalue_curve_right)))
cat("
Done.
")
