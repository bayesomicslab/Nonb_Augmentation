# IWTomics on A-Phased Repeats — detected vs undetected
# Run with Rscript, or block by block in an interactive session.
# Reference: IWTomics 1.36.0 manual
#   https://bioconductor.org/packages/release/bioc/manuals/IWTomics/man/IWTomics.pdf
#
# Inputs (already produced by the notebook through cell §6.8.2):
#   $SNP_AF_BASE/results/apr_notebook_detected_run/
#     NCNRAPhasedRepeats.gnomAD.intersect           (detected motif x SNP)
#     NCNRAPhasedRepeats.gnomAD.undet.intersect     (undetected motif x SNP)
#     split_counts.txt                              (N_MOTIF = N_DETECTED)
#
# A-Phased Repeats have 7 sub-parts per locus:
#   LeftFlank | Tract1 | Spacer1 | Tract2 | Spacer2 | Tract3 | RightFlank
#
# Outputs:
#   APhasedRepeats_detected_vs_undetected_IWTomics_left_0.01.RData
#   APhasedRepeats_detected_vs_undetected_IWTomics_center_0.01.RData
#   APhasedRepeats_detected_vs_undetected_IWTomics_right_0.01.RData


# load packages
library(IWTomics)
library(GenomicRanges)
cat("IWTomics version:        ", as.character(packageVersion("IWTomics")), "\n")
cat("GenomicRanges version:   ", as.character(packageVersion("GenomicRanges")), "\n")
cat("R version:               ", R.version.string, "\n")


# paths + parameters
WORKDIR <- file.path(Sys.getenv("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "apr_notebook_detected_run")
setwd(WORKDIR)
cat("workdir:", getwd(), "\n")

ALPHA      <- 0.01
N_PERM     <- 1000
ALIGN      <- "center"

MAX_FLANK   <- 100
TRACT_BINS  <- 40
SPACER_BINS <- 30

SCALED_TRACT  <- round(seq(1, TRACT_BINS  * 2 - 1, by = 2) / (TRACT_BINS  * 2), 4)
SCALED_SPACER <- round(seq(1, SPACER_BINS * 2 - 1, by = 2) / (SPACER_BINS * 2), 4)


# helpers

read_counts <- function(path = "split_counts.txt") {
  d <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  out <- as.integer(d[[2]]); names(out) <- d[[1]]
  out
}

# Detect n_loci automatically from V4 in a bedtools -loj intersect
# (loci are 0-indexed, so n_loci = max(V4) + 1).
detect_n_loci <- function(intersect_path) {
  cat("scanning", intersect_path, "for max locus id ...\n")
  con <- file(intersect_path, "r")
  on.exit(close(con))
  CHUNK <- 1e6
  m <- -1L
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


# counts + matrices for detected & undetected
counts    <- read_counts("split_counts.txt")
N_DET     <- counts["N_MOTIF"]
N_UNDET   <- detect_n_loci("NCNRAPhasedRepeats.gnomAD.undet.intersect")
cat("N_DETECTED =", N_DET, "  N_UNDETECTED =", N_UNDET, "\n")

# also persist N_UNDET to split_counts.txt so the plotting cell can reuse it
counts_path <- "split_counts.txt"
counts_lines <- readLines(counts_path)
if (!any(grepl("^N_UNDET\\b", counts_lines))) {
  writeLines(c(counts_lines, sprintf("N_UNDET\t%d", N_UNDET)), counts_path)
  cat("appended N_UNDET to split_counts.txt\n")
}

det <- build_matrices("NCNRAPhasedRepeats.gnomAD.intersect",        N_DET)
und <- build_matrices("NCNRAPhasedRepeats.gnomAD.undet.intersect",  N_UNDET)

cat("\nmatrix shapes (detected):\n")
for (nm in names(det)) cat(sprintf("  %-7s : %d x %d\n", nm, nrow(det[[nm]]), ncol(det[[nm]])))
cat("matrix shapes (undetected):\n")
for (nm in names(und)) cat(sprintf("  %-7s : %d x %d\n", nm, nrow(und[[nm]]), ncol(und[[nm]])))


# IWTomicsTest wrapper
# IWTomics features matrix wants positions x loci (manual p.8). build_matrices
# returns loci x bins, so transpose before handing to IWTomics.
run_iwt <- function(m_loci_x_bins, c_loci_x_bins, label) {
  m <- t(m_loci_x_bins)
  c <- t(c_loci_x_bins)
  nbins  <- nrow(m)
  n_det  <- ncol(m); n_und <- ncol(c)

  gr_det <- GRanges("chr1",
                    IRanges(start = seq(1, by = nbins, length.out = n_det),
                            width = nbins))
  gr_und <- GRanges("chr1",
                    IRanges(start = seq(1, by = nbins, length.out = n_und),
                            width = nbins))
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
  cat(sprintf("%-12s  %4d bins, %4d det, %5d undet, %4d significant\n",
              label, nbins, n_det, n_und, sum(adj < ALPHA)))
  list(adj = as.numeric(adj), nbins = nbins)
}


# run all 7 regions
cat("\nrunning IWTomicsTest on left flank ...\n")
left  <- run_iwt(det$left,  und$left,  "left flank")

cat("\nrunning IWTomicsTest on right flank ...\n")
right <- run_iwt(det$right, und$right, "right flank")

cat("\nrunning IWTomicsTest on Tract1 / Spacer1 / Tract2 / Spacer2 / Tract3 ...\n")
T1 <- run_iwt(det$Tract1,  und$Tract1,  "Tract1")
S1 <- run_iwt(det$Spacer1, und$Spacer1, "Spacer1")
T2 <- run_iwt(det$Tract2,  und$Tract2,  "Tract2")
S2 <- run_iwt(det$Spacer2, und$Spacer2, "Spacer2")
T3 <- run_iwt(det$Tract3,  und$Tract3,  "Tract3")


# save .RData files
adj_pvalue_curve_left <- data.frame(
  snp_coord   = seq(-left$nbins, -1),
  adj_pvalue  = left$adj,
  significant = left$adj < ALPHA
)
save(adj_pvalue_curve_left,
     file = "APhasedRepeats_detected_vs_undetected_IWTomics_left_0.01.RData")

adj_pvalue_curve_right <- data.frame(
  snp_coord   = seq_len(right$nbins),
  adj_pvalue  = right$adj,
  significant = right$adj < ALPHA
)
save(adj_pvalue_curve_right,
     file = "APhasedRepeats_detected_vs_undetected_IWTomics_right_0.01.RData")

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
     file = "APhasedRepeats_detected_vs_undetected_IWTomics_center_0.01.RData")


# summary
files <- list.files(WORKDIR,
                    pattern = "APhasedRepeats_detected_vs_undetected_IWTomics.*\\.RData$",
                    full.names = TRUE)
cat("\nIWTomics outputs in workdir:\n")
print(file.info(files)[, c("size", "mtime")])

cat("\nSummary of significant bins (alpha =", ALPHA, "), detected vs undetected:\n")
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
