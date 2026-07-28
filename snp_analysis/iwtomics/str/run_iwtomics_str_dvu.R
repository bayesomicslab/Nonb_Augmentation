# IWTomics on Short Tandem Repeats -- detected vs undetected
# Run with Rscript, or block by block in an interactive session.
# Reference: IWTomics 1.36.0 manual
#   https://bioconductor.org/packages/release/bioc/manuals/IWTomics/man/IWTomics.pdf
#
# Inputs (already produced by the notebook through cell 24):
#   $SNP_AF_BASE/results/str_notebook_detected_run/
#     NCNRSTR.gnomAD.intersect              (detected x SNP, bedtools -loj)
#     NCNR.STR.gnomAD.undet.intersect       (undetected x SNP)
#     split_counts.txt                      (N_MOTIF, N_CTRL)  -- N_UNDET is
#                                           auto-detected from V4 of the
#                                           undetected intersect file.
#
# STR has 3 sub-parts per locus:
#   LeftFlank | Locus | RightFlank
#
# Outputs:
#   STR_detected_vs_undetected_IWTomics_left_0.01.RData
#   STR_detected_vs_undetected_IWTomics_locus_0.01.RData
#   STR_detected_vs_undetected_IWTomics_right_0.01.RData


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

read_counts <- function(path = "split_counts.txt") {
  d <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  out <- as.integer(d[[2]]); names(out) <- d[[1]]
  out
}

# Scan V4 of intersect file to find n_loci (max(V4)+1, since lid is 0-indexed).
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

# Build per-locus matrices from a bedtools -loj intersect file.
# Returns a named list with: left, Locus, right.
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
counts  <- read_counts("split_counts.txt")
N_MOTIF <- counts["N_MOTIF"]
N_UNDET <- detect_n_loci("NCNR.STR.gnomAD.undet.intersect")
cat("N_MOTIF =", N_MOTIF, "  N_UNDET =", N_UNDET, "\n")

# Append N_UNDET to split_counts.txt if not already present.
sc <- readLines("split_counts.txt")
if (!any(grepl("^N_UNDET\\b", sc))) {
  cat(paste0("N_UNDET\t", N_UNDET, "\n"),
      file = "split_counts.txt", append = TRUE)
  cat("appended N_UNDET to split_counts.txt\n")
}

det <- build_matrices("NCNRSTR.gnomAD.intersect",            N_MOTIF)
und <- build_matrices("NCNR.STR.gnomAD.undet.intersect",     N_UNDET)

cat("\nmatrix shapes (detected):\n")
for (nm in names(det)) cat(sprintf("  %-8s : %d x %d\n", nm, nrow(det[[nm]]), ncol(det[[nm]])))


# IWTomicsTest wrapper
# IWTomics features matrix wants positions x loci.
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
    GRangesList(detected = gr_mot, undetected = gr_ctrl),
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
  cat(sprintf("%-12s  %4d bins, %5d det, %5d undet, %4d significant\n",
              label, nbins, n_mot, n_ctrl, sum(adj < ALPHA)))
  list(adj = as.numeric(adj), nbins = nbins)
}


# run all 3 regions
cat("\nrunning IWTomicsTest on left flank ...\n")
left  <- run_iwt(det$left,  und$left,  "left flank")

cat("\nrunning IWTomicsTest on Locus ...\n")
loc   <- run_iwt(det$Locus, und$Locus, "Locus")

cat("\nrunning IWTomicsTest on right flank ...\n")
right <- run_iwt(det$right, und$right, "right flank")


# save .RData files
adj_pvalue_curve_left <- data.frame(
  snp_coord   = seq(-left$nbins, -1),
  adj_pvalue  = left$adj,
  significant = left$adj < ALPHA
)
save(adj_pvalue_curve_left,
     file = "STR_detected_vs_undetected_IWTomics_left_0.01.RData")

adj_pvalue_curve_locus <- data.frame(
  bin         = seq_len(loc$nbins),
  adj_pvalue  = loc$adj,
  significant = loc$adj < ALPHA
)
save(adj_pvalue_curve_locus,
     file = "STR_detected_vs_undetected_IWTomics_locus_0.01.RData")

adj_pvalue_curve_right <- data.frame(
  snp_coord   = seq_len(right$nbins),
  adj_pvalue  = right$adj,
  significant = right$adj < ALPHA
)
save(adj_pvalue_curve_right,
     file = "STR_detected_vs_undetected_IWTomics_right_0.01.RData")


# summary
files <- list.files(WORKDIR,
                    pattern = "STR_detected_vs_undetected_IWTomics.*\\.RData$",
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
