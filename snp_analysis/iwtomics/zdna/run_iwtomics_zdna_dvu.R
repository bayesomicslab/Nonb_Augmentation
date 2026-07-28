# IWTomics on Z-DNA — detected vs undetected
# Run with Rscript, or block by block in an interactive session.
# Reference: IWTomics 1.36.0 manual
#
# Inputs (produced by the notebook through §3.10.3):
#   $SNP_AF_BASE/results/zdna_notebook_detected_run/
#     NCNRZDNA.gnomAD.intersect           (detected motif x SNP, -loj output)
#     NCNRZDNA.gnomAD.undet.intersect     (undetected motif x SNP, -loj output)
#     split_counts.txt                    (N_MOTIF = N_DETECTED)
#
# Z-DNA has 3 sub-parts per locus: LeftFlank | Locus | RightFlank
# (no spacer / repeat partition like DR / IR / MR).
#
# Outputs (saved into the workdir so the plot cell can `load()` them):
#   ZDNA_detected_vs_undetected_IWTomics_left_0.01.RData   -> adj_pvalue_curve_left
#   ZDNA_detected_vs_undetected_IWTomics_locus_0.01.RData  -> adj_pvalue_curve_locus
#   ZDNA_detected_vs_undetected_IWTomics_right_0.01.RData  -> adj_pvalue_curve_right
#
# Side effect: appends N_UNDET to split_counts.txt so the plot cell can read it.


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

ALPHA      <- 0.01
N_PERM     <- 1000
ALIGN      <- "center"
MAX_FLANK  <- 100
LOCUS_BINS <- 180

SCALED_LOCUS <- round(seq(1, LOCUS_BINS * 2 - 1, by = 2) / (LOCUS_BINS * 2), 4)


# helpers

read_counts <- function(path = "split_counts.txt") {
  d <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  out <- as.integer(d[[2]]); names(out) <- d[[1]]
  out
}

# Detect n_loci = max(V4) + 1 by streaming the (large) intersect file.
scan_n_loci <- function(intersect_path) {
  cat("scanning", intersect_path, "for max locus id ...\n")
  con <- file(intersect_path, "r"); on.exit(close(con))
  CHUNK <- 1e6; m <- -1L
  repeat {
    lines <- readLines(con, n = CHUNK, warn = FALSE)
    if (!length(lines)) break
    parts <- strsplit(lines, "\t", fixed = TRUE)
    v4 <- as.integer(vapply(parts, `[`, character(1), 4L))
    chunk_max <- suppressWarnings(max(v4, na.rm = TRUE))
    if (is.finite(chunk_max) && chunk_max > m) m <- chunk_max
  }
  if (m < 0) stop("could not parse V4 from ", intersect_path)
  m + 1L
}

# Mirrors _build_matrices() in the notebook (cell §3.7).
build_matrices <- function(intersect_path, n_loci) {
  left  <- matrix(0, nrow = n_loci, ncol = MAX_FLANK)
  locus <- matrix(0, nrow = n_loci, ncol = LOCUS_BINS)
  right <- matrix(0, nrow = n_loci, ncol = MAX_FLANK)
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
    for (k in seq_along(rows)) left[rows[k], cols[k]] <- left[rows[k], cols[k]] + 1
  }

  is_r <- ptype == "RightFlank"
  if (any(is_r)) {
    coord <- snp_start[is_r] - ps[is_r]
    idx   <- coord + 1L
    keep  <- idx >= 1L & idx <= MAX_FLANK
    rows  <- lid[is_r][keep]; cols <- idx[keep]
    for (k in seq_along(rows)) right[rows[k], cols[k]] <- right[rows[k], cols[k]] + 1
  }

  is_c <- ptype == "Locus" & plen > 0
  if (any(is_c)) {
    sub_lid <- lid[is_c]
    sub_ps  <- ps[is_c]; sub_pe <- pe[is_c]; sub_pl <- plen[is_c]; sub_ss <- snp_start[is_c]
    fs <- (sub_ss     - sub_ps) / sub_pl
    fe <- (sub_ss + 1 - sub_ps) / sub_pl
    for (k in seq_along(sub_lid)) {
      hits <- which(SCALED_LOCUS >= fs[k] & SCALED_LOCUS < fe[k])
      if (length(hits)) locus[sub_lid[k], hits] <- locus[sub_lid[k], hits] + 1
    }
  }
  list(left = left, locus = locus, right = right)
}


# counts + matrices for detected & undetected
counts <- read_counts("split_counts.txt")
N_DET  <- counts["N_MOTIF"]
N_UNDET <- scan_n_loci("NCNRZDNA.gnomAD.undet.intersect")
cat(sprintf("N_DETECTED = %d   N_UNDETECTED = %d\n", N_DET, N_UNDET))

# persist N_UNDET to split_counts.txt so the plot cell can read it back
counts_path  <- "split_counts.txt"
counts_lines <- readLines(counts_path)
if (!any(grepl("^N_UNDET\\b", counts_lines))) {
  writeLines(c(counts_lines, sprintf("N_UNDET\t%d", N_UNDET)), counts_path)
  cat("appended N_UNDET to split_counts.txt\n")
}

det <- build_matrices("NCNRZDNA.gnomAD.intersect",        N_DET)
und <- build_matrices("NCNRZDNA.gnomAD.undet.intersect",  N_UNDET)

cat("\nmatrix shapes (detected):\n")
cat("  left  :", dim(det$left),  "\n"); cat("  locus :", dim(det$locus), "\n"); cat("  right :", dim(det$right), "\n")
cat("matrix shapes (undetected):\n")
cat("  left  :", dim(und$left),  "\n"); cat("  locus :", dim(und$locus), "\n"); cat("  right :", dim(und$right), "\n")


# IWTomicsTest wrapper (detected vs undetected)
run_iwt <- function(m_loci_x_bins, c_loci_x_bins, label) {
  m <- t(m_loci_x_bins); c <- t(c_loci_x_bins)
  nbins  <- nrow(m)
  n_det  <- ncol(m); n_und <- ncol(c)
  gr_det <- GRanges("chr1", IRanges(start = seq(1, by = nbins, length.out = n_det), width = nbins))
  gr_und <- GRanges("chr1", IRanges(start = seq(1, by = nbins, length.out = n_und), width = nbins))
  rf <- IWTomicsData(
    GRangesList(detected = gr_det, undetected = gr_und),
    list(snp = list(detected = m, undetected = c)),
    alignment     = ALIGN,
    id_regions    = c("detected", "undetected"),
    name_regions  = c("detected", "undetected"),
    id_features   = "snp",
    name_features = "snp"
  )
  res <- IWTomicsTest(rf,
                      id_region1 = "detected", id_region2 = "undetected",
                      id_features_subset = "snp",
                      statistics = "mean", B = N_PERM)
  adj <- adjusted_pval(res)[[1]][[1]]
  cat(sprintf("%-12s  %4d bins, %5d det, %6d undet, %4d significant\n",
              label, nbins, n_det, n_und, sum(adj < ALPHA)))
  list(adj = as.numeric(adj), nbins = nbins)
}


# run all three regions and save .RData
left  <- run_iwt(det$left,  und$left,  "left flank")
locus <- run_iwt(det$locus, und$locus, "locus")
right <- run_iwt(det$right, und$right, "right flank")

adj_pvalue_curve_left <- data.frame(
  snp_coord   = seq(-left$nbins, -1),
  adj_pvalue  = left$adj,
  significant = left$adj < ALPHA
)
save(adj_pvalue_curve_left,
     file = "ZDNA_detected_vs_undetected_IWTomics_left_0.01.RData")
cat("wrote ZDNA_detected_vs_undetected_IWTomics_left_0.01.RData\n")

adj_pvalue_curve_locus <- data.frame(
  bin         = seq_len(locus$nbins),
  adj_pvalue  = locus$adj,
  significant = locus$adj < ALPHA
)
save(adj_pvalue_curve_locus,
     file = "ZDNA_detected_vs_undetected_IWTomics_locus_0.01.RData")
cat("wrote ZDNA_detected_vs_undetected_IWTomics_locus_0.01.RData\n")

adj_pvalue_curve_right <- data.frame(
  snp_coord   = seq_len(right$nbins),
  adj_pvalue  = right$adj,
  significant = right$adj < ALPHA
)
save(adj_pvalue_curve_right,
     file = "ZDNA_detected_vs_undetected_IWTomics_right_0.01.RData")
cat("wrote ZDNA_detected_vs_undetected_IWTomics_right_0.01.RData\n")


# verify outputs + summary
files <- list.files(WORKDIR, pattern = "ZDNA_detected_vs_undetected_IWTomics.*\\.RData$",
                    full.names = TRUE)
cat("\nIWTomics outputs in workdir:\n")
print(file.info(files)[, c("size", "mtime")])

cat("\nSummary of significant bins (alpha =", ALPHA, "), detected vs undetected:\n")
cat(sprintf("  left  : %d/%d\n", sum(adj_pvalue_curve_left$significant),  nrow(adj_pvalue_curve_left)))
cat(sprintf("  locus : %d/%d\n", sum(adj_pvalue_curve_locus$significant), nrow(adj_pvalue_curve_locus)))
cat(sprintf("  right : %d/%d\n", sum(adj_pvalue_curve_right$significant), nrow(adj_pvalue_curve_right)))
cat("
Done.
")
