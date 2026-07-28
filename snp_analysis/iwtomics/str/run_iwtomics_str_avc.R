# IWTomics on Short Tandem Repeats -- all-test vs all-test control
# Run with Rscript, or block by block in an interactive session.
# Reference: IWTomics 1.36.0 manual
#   https://bioconductor.org/packages/release/bioc/manuals/IWTomics/man/IWTomics.pdf
#
# Inputs (already produced by the notebook through cell 35):
#   $SNP_AF_BASE/results/str_notebook_detected_run/
#     NCNR.STR.all.gnomAD.intersect              (all-test STR x SNP)
#     NCNR.STR.all.gnomAD.shuffled.intersect     (all-test control x SNP)
#
# N_ALL and N_ALL_CTRL are auto-detected from V4 of the intersect files.
#
# STR has 3 sub-parts per locus:
#   LeftFlank | Locus | RightFlank
#
# Outputs:
#   STR_alltest_vs_control_IWTomics_left_0.01.RData
#   STR_alltest_vs_control_IWTomics_locus_0.01.RData
#   STR_alltest_vs_control_IWTomics_right_0.01.RData


# load packages
library(IWTomics)
library(GenomicRanges)
cat("IWTomics version:        ", as.character(packageVersion("IWTomics")), "\n")
cat("GenomicRanges version:   ", as.character(packageVersion("GenomicRanges")), "\n")
cat("R version:               ", R.version.string, "\n")


# paths + parameters
WORKDIR <- file.path(Sys.getenv("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "str_notebook_detected_run")
setwd(WORKDIR)
cat("workdir:", getwd(), "\n")

ALPHA      <- 0.01
N_PERM     <- 1000
ALIGN      <- "center"

MAX_FLANK  <- 100
LOCUS_BINS <- 180

SCALED_LOCUS <- round(seq(1, LOCUS_BINS * 2 - 1, by = 2) / (LOCUS_BINS * 2), 4)


# helpers

detect_n_loci <- function(intersect_path) {
  ncol_file <- length(read.table(intersect_path, sep = "\t", header = FALSE,
                                 nrows = 1, comment.char = "", quote = ""))
  cc <- rep("NULL", ncol_file); cc[4] <- "integer"
  d <- read.table(intersect_path, sep = "\t", header = FALSE,
                  stringsAsFactors = FALSE, comment.char = "",
                  quote = "", colClasses = cc)
  n <- max(d$V4) + 1L
  cat(sprintf("  detected n_loci=%d from %s\n", n, intersect_path))
  n
}

build_matrices <- function(intersect_path, n_loci) {
  out <- list(
    left  = matrix(0, nrow = n_loci, ncol = MAX_FLANK),
    Locus = matrix(0, nrow = n_loci, ncol = LOCUS_BINS),
    right = matrix(0, nrow = n_loci, ncol = MAX_FLANK)
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

  is_p <- ptype == "Locus" & plen > 0
  if (any(is_p)) {
    sub_lid <- lid[is_p]
    sub_ps  <- ps[is_p]; sub_pe <- pe[is_p]; sub_pl <- plen[is_p]; sub_ss <- snp_start[is_p]
    fs <- (sub_ss     - sub_ps) / sub_pl
    fe <- (sub_ss + 1 - sub_ps) / sub_pl
    for (k in seq_along(sub_lid)) {
      hits <- which(SCALED_LOCUS >= fs[k] & SCALED_LOCUS < fe[k])
      if (length(hits)) out$Locus[sub_lid[k], hits] <- out$Locus[sub_lid[k], hits] + 1
    }
  }
  out
}


# counts + matrices
N_ALL      <- detect_n_loci("NCNR.STR.all.gnomAD.intersect")
N_ALL_CTRL <- detect_n_loci("NCNR.STR.all.gnomAD.shuffled.intersect")
cat("N_ALL =", N_ALL, "  N_ALL_CTRL =", N_ALL_CTRL, "\n")

# Append N_ALL and N_ALL_CTRL to split_counts.txt if not already present.
if (file.exists("split_counts.txt")) {
  sc <- readLines("split_counts.txt")
  if (!any(grepl("^N_ALL\\b", sc))) {
    cat(paste0("N_ALL\t", N_ALL, "\n"),
        file = "split_counts.txt", append = TRUE)
    cat("appended N_ALL to split_counts.txt\n")
  }
  if (!any(grepl("^N_ALL_CTRL\\b", sc))) {
    cat(paste0("N_ALL_CTRL\t", N_ALL_CTRL, "\n"),
        file = "split_counts.txt", append = TRUE)
    cat("appended N_ALL_CTRL to split_counts.txt\n")
  }
}

all_ <- build_matrices("NCNR.STR.all.gnomAD.intersect",          N_ALL)
ctl  <- build_matrices("NCNR.STR.all.gnomAD.shuffled.intersect", N_ALL_CTRL)


# IWTomicsTest wrapper
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


# run all 3 regions
cat("\nrunning IWTomicsTest on left flank ...\n")
left  <- run_iwt(all_$left,  ctl$left,  "left flank")

cat("\nrunning IWTomicsTest on Locus ...\n")
loc   <- run_iwt(all_$Locus, ctl$Locus, "Locus")

cat("\nrunning IWTomicsTest on right flank ...\n")
right <- run_iwt(all_$right, ctl$right, "right flank")


# save .RData files
adj_pvalue_curve_left <- data.frame(
  snp_coord   = seq(-left$nbins, -1),
  adj_pvalue  = left$adj,
  significant = left$adj < ALPHA
)
save(adj_pvalue_curve_left,
     file = "STR_alltest_vs_control_IWTomics_left_0.01.RData")

adj_pvalue_curve_locus <- data.frame(
  bin         = seq_len(loc$nbins),
  adj_pvalue  = loc$adj,
  significant = loc$adj < ALPHA
)
save(adj_pvalue_curve_locus,
     file = "STR_alltest_vs_control_IWTomics_locus_0.01.RData")

adj_pvalue_curve_right <- data.frame(
  snp_coord   = seq_len(right$nbins),
  adj_pvalue  = right$adj,
  significant = right$adj < ALPHA
)
save(adj_pvalue_curve_right,
     file = "STR_alltest_vs_control_IWTomics_right_0.01.RData")


# summary
files <- list.files(WORKDIR,
                    pattern = "STR_alltest_vs_control_IWTomics.*\\.RData$",
                    full.names = TRUE)
cat("\nIWTomics outputs in workdir:\n")
print(file.info(files)[, c("size", "mtime")])

cat("\nSummary of significant bins (alpha =", ALPHA, "):\n")
cat(sprintf("  left   : %d/%d\n",
            sum(adj_pvalue_curve_left$significant),  nrow(adj_pvalue_curve_left)))
cat(sprintf("  Locus  : %d/%d\n",
            sum(adj_pvalue_curve_locus$significant), nrow(adj_pvalue_curve_locus)))
cat(sprintf("  right  : %d/%d\n",
            sum(adj_pvalue_curve_right$significant), nrow(adj_pvalue_curve_right)))
cat("
Done.
")
