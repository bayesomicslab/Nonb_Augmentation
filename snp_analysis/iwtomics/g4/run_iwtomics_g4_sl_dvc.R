# IWTomics on G-Quadruplexes (stem/loop sub-annotation) -- detected vs control
# Reference: Supplementary code 2.ipynb section 1 (G-quadruplexes)
#   locus split into 9 parts: LeftFlank | stem1 | loop1 | stem2 | loop2 |
#                             stem3 | loop3 | stem4 | RightFlank
#   stems scaled to 30 bins, loops to 20 bins, flanks 100 bp (unscaled).
#
# Inputs (produced by the notebook through the intersect cell):
#   results/g4_stemloop_run/
#     G4sl.gnomAD.intersect                 (detected x SNP, bedtools -loj)
#     G4sl.gnomAD.shuffled.intersect        (control  x SNP, bedtools -loj)
#     split_counts.txt                      (N_MOTIF, N_CTRL)
#
# Outputs:
#   G4sl_detected_vs_control_IWTomics_left_0.01.RData    -> adj_pvalue_curve_left
#   G4sl_detected_vs_control_IWTomics_right_0.01.RData   -> adj_pvalue_curve_right
#   G4sl_detected_vs_control_IWTomics_center_0.01.RData  -> adj_pvalue_curve_center
#                                                              (rbind of 7 subparts,
#                                                               with a `subpart` column)


# load packages
library(IWTomics)
library(GenomicRanges)
cat("IWTomics version:        ", as.character(packageVersion("IWTomics")), "\n")
cat("GenomicRanges version:   ", as.character(packageVersion("GenomicRanges")), "\n")
cat("R version:               ", R.version.string, "\n")


# paths + parameters
WORKDIR <- file.path(Sys.getenv("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "g4_stemloop_run")
setwd(WORKDIR)
cat("workdir:", getwd(), "\n")

ALPHA     <- 0.01
N_PERM    <- 1000
ALIGN     <- "center"

MAX_FLANK <- 100
STEM_BINS <- 30
LOOP_BINS <- 20

CENTER_PARTS <- c("stem1", "loop1", "stem2", "loop2", "stem3", "loop3", "stem4")
PART_BINS    <- c(stem1 = STEM_BINS, loop1 = LOOP_BINS, stem2 = STEM_BINS,
                  loop2 = LOOP_BINS, stem3 = STEM_BINS, loop3 = LOOP_BINS,
                  stem4 = STEM_BINS)

SCALED_STEM <- round(seq(1, STEM_BINS * 2 - 1, by = 2) / (STEM_BINS * 2), 4)
SCALED_LOOP <- round(seq(1, LOOP_BINS * 2 - 1, by = 2) / (LOOP_BINS * 2), 4)


# helpers

read_counts <- function(path = "split_counts.txt") {
  d <- read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  out <- as.integer(d[[2]]); names(out) <- d[[1]]
  out
}

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
# Columns of the split bed (V1..V5): chrom, start, end, lid, part.
# Returns a named list: left, right, and one matrix per CENTER_PARTS.
build_matrices <- function(intersect_path, n_loci) {
  out <- list(left  = matrix(0, nrow = n_loci, ncol = MAX_FLANK),
              right = matrix(0, nrow = n_loci, ncol = MAX_FLANK))
  for (p in CENTER_PARTS) out[[p]] <- matrix(0, nrow = n_loci, ncol = PART_BINS[[p]])
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
    idx  <- -(pe[is_l] - snp_start[is_l]) + MAX_FLANK + 1L
    keep <- idx >= 1L & idx <= MAX_FLANK
    rows <- lid[is_l][keep]; cols <- idx[keep]
    for (k in seq_along(rows)) out$left[rows[k], cols[k]] <- out$left[rows[k], cols[k]] + 1
  }
  is_r <- ptype == "RightFlank"
  if (any(is_r)) {
    idx  <- (snp_start[is_r] - ps[is_r]) + 1L
    keep <- idx >= 1L & idx <= MAX_FLANK
    rows <- lid[is_r][keep]; cols <- idx[keep]
    for (k in seq_along(rows)) out$right[rows[k], cols[k]] <- out$right[rows[k], cols[k]] + 1
  }

  for (p in CENTER_PARTS) {
    scaled <- if (grepl("stem", p)) SCALED_STEM else SCALED_LOOP
    is_p <- ptype == p & plen > 0
    if (!any(is_p)) next
    sl <- lid[is_p]; sps <- ps[is_p]; spl <- plen[is_p]; sss <- snp_start[is_p]
    fs <- (sss     - sps) / spl
    fe <- (sss + 1 - sps) / spl
    for (k in seq_along(sl)) {
      hits <- which(scaled >= fs[k] & scaled < fe[k])
      if (length(hits)) out[[p]][sl[k], hits] <- out[[p]][sl[k], hits] + 1
    }
  }
  out
}


# counts + matrices
counts  <- read_counts("split_counts.txt")
N_MOTIF <- counts["N_MOTIF"]
N_CTRL  <- counts["N_CTRL"]
cat("N_MOTIF =", N_MOTIF, "  N_CTRL =", N_CTRL, "\n")

det <- build_matrices("G4sl.gnomAD.intersect",          N_MOTIF)
und <- build_matrices("G4sl.gnomAD.shuffled.intersect", N_CTRL)


# IWTomicsTest wrapper
run_iwt <- function(m_loci_x_bins, c_loci_x_bins, label) {
  m <- t(m_loci_x_bins); c <- t(c_loci_x_bins)
  nbins <- nrow(m); n_mot <- ncol(m); n_ctrl <- ncol(c)
  gr_mot  <- GRanges("chr1", IRanges(start = seq(1, by = nbins, length.out = n_mot),  width = nbins))
  gr_ctrl <- GRanges("chr1", IRanges(start = seq(1, by = nbins, length.out = n_ctrl), width = nbins))
  rf <- IWTomicsData(
    GRangesList(detected = gr_mot, control = gr_ctrl),
    list(snp = list(detected = m, control = c)),
    alignment = ALIGN,
    id_regions = c("detected", "control"),
    name_regions = c("detected", "control"),
    id_features = "snp", name_features = "snp")
  res <- IWTomicsTest(rf, id_region1 = "detected", id_region2 = "control",
                      id_features_subset = "snp", statistics = "mean", B = N_PERM)
  adj <- as.numeric(adjusted_pval(res)[[1]][[1]])
  cat(sprintf("%-10s  %3d bins, %5d det, %5d ctrl, %4d significant\n",
              label, nbins, n_mot, n_ctrl, sum(adj < ALPHA)))
  list(adj = adj, nbins = nbins)
}


# run left, right, and 7 center subparts
cat("\nleft flank ...\n");  left  <- run_iwt(det$left,  und$left,  "left")
cat("right flank ...\n");   right <- run_iwt(det$right, und$right, "right")

center_list <- list()
for (p in CENTER_PARTS) {
  cat(p, "...\n")
  center_list[[p]] <- run_iwt(det[[p]], und[[p]], p)
}


# save .RData files
adj_pvalue_curve_left <- data.frame(
  snp_coord = seq(-left$nbins, -1), adj_pvalue = left$adj,
  significant = left$adj < ALPHA)
save(adj_pvalue_curve_left,
     file = "G4sl_detected_vs_control_IWTomics_left_0.01.RData")

adj_pvalue_curve_right <- data.frame(
  snp_coord = seq_len(right$nbins), adj_pvalue = right$adj,
  significant = right$adj < ALPHA)
save(adj_pvalue_curve_right,
     file = "G4sl_detected_vs_control_IWTomics_right_0.01.RData")

adj_pvalue_curve_center <- do.call(rbind, lapply(CENTER_PARTS, function(p) {
  data.frame(subpart = p, bin = seq_len(center_list[[p]]$nbins),
             adj_pvalue = center_list[[p]]$adj,
             significant = center_list[[p]]$adj < ALPHA)
}))
save(adj_pvalue_curve_center,
     file = "G4sl_detected_vs_control_IWTomics_center_0.01.RData")


# summary
cat("\nSummary of significant bins (alpha =", ALPHA, "):\n")
cat(sprintf("  left   : %d/%d\n", sum(adj_pvalue_curve_left$significant),  nrow(adj_pvalue_curve_left)))
cat(sprintf("  right  : %d/%d\n", sum(adj_pvalue_curve_right$significant), nrow(adj_pvalue_curve_right)))
for (p in CENTER_PARTS) {
  s <- adj_pvalue_curve_center[adj_pvalue_curve_center$subpart == p, ]
  cat(sprintf("  %-6s : %d/%d\n", p, sum(s$significant), nrow(s)))
}
cat("
Done.
")
