# IWTomics on Inverted Repeats (model-detected)
# Run with Rscript, or block by block in an interactive session.
# Reference: IWTomics 1.36.0 manual
#   https://bioconductor.org/packages/release/bioc/manuals/IWTomics/man/IWTomics.pdf
#
# Inputs (already produced by the notebook through cell 23):
#   $SNP_AF_BASE/results/ir_notebook_detected_run/
#     NCNR.InvertedRepeats.gnomAD.intersect              (motif x SNP, bedtools -loj)
#     NCNR.InvertedRepeats.gnomAD.shuffled.intersect     (control x SNP)
#     split_counts.txt                                   (N_MOTIF, N_MOTIF_SPACER,
#                                                         N_CTRL,  N_CTRL_SPACER)
#
# Inverted Repeats have 5 sub-parts per locus:
#   LeftFlank | LeftRepeat | Spacer | RightRepeat | RightFlank
# Center = LeftRepeat + Spacer + RightRepeat (3 sub-parts).
#
# Spacer: zero-spacer loci have empty Spacer rows (plen=0 is skipped), which
# become all-zero rows in the Spacer matrix. The Python pipeline passes the
# full matrix to iwtomics() and divides the plotted frequency by
# N_MOTIF_SPACER; this port keeps that behaviour, so IWTomicsTest also sees
# the full matrix.
#
# Outputs:
#   InvertedRepeats_motif_vs_control_IWTomics_left_0.01.RData
#   InvertedRepeats_motif_vs_control_IWTomics_center_0.01.RData
#   InvertedRepeats_motif_vs_control_IWTomics_right_0.01.RData


# load packages
library(IWTomics)
library(GenomicRanges)
cat("IWTomics version:        ", as.character(packageVersion("IWTomics")), "\n")
cat("GenomicRanges version:   ", as.character(packageVersion("GenomicRanges")), "\n")
cat("R version:               ", R.version.string, "\n")


# paths + parameters
WORKDIR <- file.path(Sys.getenv("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "ir_notebook_detected_run")
setwd(WORKDIR)
cat("workdir:", getwd(), "\n")

ALPHA      <- 0.01
N_PERM     <- 1000
ALIGN      <- "center"

MAX_FLANK   <- 100
REPEAT_BINS <- 50
SPACER_BINS <- 80

SCALED_REPEAT <- round(seq(1, REPEAT_BINS * 2 - 1, by = 2) / (REPEAT_BINS * 2), 4)
SCALED_SPACER <- round(seq(1, SPACER_BINS * 2 - 1, by = 2) / (SPACER_BINS * 2), 4)


# helpers

read_counts <- function(path = "split_counts.txt") {
  d <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  out <- as.integer(d[[2]]); names(out) <- d[[1]]
  out
}

# Build per-locus matrices from a bedtools -loj intersect file.
# Returns a named list:
#   left, right     : n_loci x MAX_FLANK
#   LeftRepeat,
#   RightRepeat     : n_loci x REPEAT_BINS
#   Spacer          : n_loci x SPACER_BINS  (zero-spacer loci -> all-zero rows)
build_matrices <- function(intersect_path, n_loci) {
  out <- list(
    left        = matrix(0, nrow = n_loci, ncol = MAX_FLANK),
    LeftRepeat  = matrix(0, nrow = n_loci, ncol = REPEAT_BINS),
    Spacer      = matrix(0, nrow = n_loci, ncol = SPACER_BINS),
    RightRepeat = matrix(0, nrow = n_loci, ncol = REPEAT_BINS),
    right       = matrix(0, nrow = n_loci, ncol = MAX_FLANK)
  )
  cat("reading", intersect_path, "...\n")

  d <- read.table(intersect_path, sep = "\t", header = FALSE,
                  stringsAsFactors = FALSE, comment.char = "",
                  quote = "", colClasses = c(
                    V1 = "character", V2 = "integer", V3 = "integer",
                    V4 = "integer",   V5 = "character",
                    V6 = "character", V7 = "integer", V8 = "integer"))
  cat("  ", nrow(d), "rows total\n")

  d <- d[d$V6 != ".", ]
  cat("  ", nrow(d), "rows after dropping no-overlap\n")

  key <- paste(d$V4, d$V5, d$V7, sep = "|")
  d <- d[!duplicated(key), ]
  cat("  ", nrow(d), "rows after dedup\n")

  ps        <- d$V2
  pe        <- d$V3
  lid       <- d$V4 + 1L
  ptype     <- d$V5
  snp_start <- d$V7
  plen      <- pe - ps

  is_l <- ptype == "LeftFlank"
  if (any(is_l)) {
    coord <- -(pe[is_l] - snp_start[is_l])
    idx   <- coord + MAX_FLANK + 1L
    keep  <- idx >= 1L & idx <= MAX_FLANK
    rows  <- lid[is_l][keep]; cols <- idx[keep]
    for (k in seq_along(rows)) out$left[rows[k], cols[k]] <- out$left[rows[k], cols[k]] + 1
  }

  is_r <- ptype == "RightFlank"
  if (any(is_r)) {
    coord <- snp_start[is_r] - ps[is_r]
    idx   <- coord + 1L
    keep  <- idx >= 1L & idx <= MAX_FLANK
    rows  <- lid[is_r][keep]; cols <- idx[keep]
    for (k in seq_along(rows)) out$right[rows[k], cols[k]] <- out$right[rows[k], cols[k]] + 1
  }

  for (pname in c("LeftRepeat", "RightRepeat", "Spacer")) {
    is_p <- ptype == pname & plen > 0
    if (!any(is_p)) next
    grid <- if (endsWith(pname, "Repeat")) SCALED_REPEAT else SCALED_SPACER
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


# counts + matrices for motif & control
counts  <- read_counts("split_counts.txt")
N_MOTIF <- counts["N_MOTIF"]
N_CTRL  <- counts["N_CTRL"]
cat("N_MOTIF =", N_MOTIF, "  N_CTRL =", N_CTRL, "\n")

mot <- build_matrices("NCNR.InvertedRepeats.gnomAD.intersect",          N_MOTIF)
ctl <- build_matrices("NCNR.InvertedRepeats.gnomAD.shuffled.intersect", N_CTRL)

cat("\nmatrix shapes (motif):\n")
for (nm in names(mot)) cat(sprintf("  %-12s : %d x %d\n", nm, nrow(mot[[nm]]), ncol(mot[[nm]])))


# IWTomicsTest wrapper
# IWTomics features matrix wants positions x loci (manual p.8). build_matrices
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
  cat(sprintf("%-12s  %4d bins, %5d motif, %5d control, %4d significant\n",
              label, nbins, n_mot, n_ctrl, sum(adj < ALPHA)))
  list(adj = as.numeric(adj), nbins = nbins)
}


# run all 5 regions (left + 3 center sub-parts + right)
cat("\nrunning IWTomicsTest on left flank ...\n")
left  <- run_iwt(mot$left,  ctl$left,  "left flank")

cat("\nrunning IWTomicsTest on right flank ...\n")
right <- run_iwt(mot$right, ctl$right, "right flank")

cat("\nrunning IWTomicsTest on LeftRepeat / Spacer / RightRepeat ...\n")
LR <- run_iwt(mot$LeftRepeat,  ctl$LeftRepeat,  "LeftRepeat")
SP <- run_iwt(mot$Spacer,      ctl$Spacer,      "Spacer")
RR <- run_iwt(mot$RightRepeat, ctl$RightRepeat, "RightRepeat")


# save .RData files
adj_pvalue_curve_left <- data.frame(
  snp_coord   = seq(-left$nbins, -1),
  adj_pvalue  = left$adj,
  significant = left$adj < ALPHA
)
save(adj_pvalue_curve_left,
     file = "InvertedRepeats_motif_vs_control_IWTomics_left_0.01.RData")

adj_pvalue_curve_right <- data.frame(
  snp_coord   = seq_len(right$nbins),
  adj_pvalue  = right$adj,
  significant = right$adj < ALPHA
)
save(adj_pvalue_curve_right,
     file = "InvertedRepeats_motif_vs_control_IWTomics_right_0.01.RData")

adj_pvalue_curve_center <- rbind(
  data.frame(subpart = "LeftRepeat",  bin = seq_len(LR$nbins),
             adj_pvalue = LR$adj, significant = LR$adj < ALPHA),
  data.frame(subpart = "Spacer",      bin = seq_len(SP$nbins),
             adj_pvalue = SP$adj, significant = SP$adj < ALPHA),
  data.frame(subpart = "RightRepeat", bin = seq_len(RR$nbins),
             adj_pvalue = RR$adj, significant = RR$adj < ALPHA)
)
save(adj_pvalue_curve_center,
     file = "InvertedRepeats_motif_vs_control_IWTomics_center_0.01.RData")


# summary
files <- list.files(WORKDIR,
                    pattern = "InvertedRepeats_motif_vs_control_IWTomics.*\\.RData$",
                    full.names = TRUE)
cat("\nIWTomics outputs in workdir:\n")
print(file.info(files)[, c("size", "mtime")])

cat("\nSummary of significant bins (alpha =", ALPHA, "):\n")
cat(sprintf("  left        : %d/%d\n",
            sum(adj_pvalue_curve_left$significant),  nrow(adj_pvalue_curve_left)))
for (sp in c("LeftRepeat", "Spacer", "RightRepeat")) {
  rows <- adj_pvalue_curve_center[adj_pvalue_curve_center$subpart == sp, ]
  cat(sprintf("  %-11s : %d/%d\n", sp, sum(rows$significant), nrow(rows)))
}
cat(sprintf("  right       : %d/%d\n",
            sum(adj_pvalue_curve_right$significant), nrow(adj_pvalue_curve_right)))
cat("
Done.
")
