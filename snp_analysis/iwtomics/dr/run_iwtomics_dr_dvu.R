# IWTomics on Direct Repeats -- detected vs undetected
# Run with Rscript, or block by block in an interactive session.
#
# Inputs (produced by the notebook through Sec 4.10.3):
#   $SNP_AF_BASE/results/dr_notebook_detected_run/
#     NCNR.DirectRepeats.gnomAD.intersect              (detected motif x SNP)
#     NCNR.DirectRepeats.gnomAD.undet.intersect        (undetected motif x SNP)
#     split_counts.txt                                 (N_MOTIF = N_DETECTED)
#
# Direct Repeats have 5 sub-parts per locus:
#   LeftFlank | LeftRepeat | Spacer | RightRepeat | RightFlank
#
# Outputs:
#   DirectRepeats_detected_vs_undetected_IWTomics_left_0.01.RData
#   DirectRepeats_detected_vs_undetected_IWTomics_center_0.01.RData
#   DirectRepeats_detected_vs_undetected_IWTomics_right_0.01.RData
#
# Side effect: appends N_UNDET (and N_UNDET_SPACER) to split_counts.txt so the
# plotting cell can read it back.


# load packages
library(IWTomics)
library(GenomicRanges)
cat("IWTomics version:        ", as.character(packageVersion("IWTomics")), "\n")
cat("GenomicRanges version:   ", as.character(packageVersion("GenomicRanges")), "\n")
cat("R version:               ", R.version.string, "\n")


# paths + parameters
WORKDIR <- file.path(Sys.getenv("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "dr_notebook_detected_run")
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

# Detect (n_loci, n_loci_spacer) from V4/V5 of an intersect.
# n_loci         = max(V4) + 1
# n_loci_spacer  = #distinct V4 with at least one Spacer row of plen>0
scan_loci_counts <- function(intersect_path) {
  cat("scanning", intersect_path, "for max locus id and spacer-loci ...\n")
  con <- file(intersect_path, "r")
  on.exit(close(con))
  CHUNK <- 1e6
  m <- -1L
  spacer_lids <- integer(0)
  repeat {
    lines <- readLines(con, n = CHUNK, warn = FALSE)
    if (!length(lines)) break
    parts <- strsplit(lines, "\t", fixed = TRUE)
    v4 <- as.integer(vapply(parts, `[`, character(1), 4L))
    v5 <- vapply(parts, `[`, character(1), 5L)
    v2 <- as.integer(vapply(parts, `[`, character(1), 2L))
    v3 <- as.integer(vapply(parts, `[`, character(1), 3L))
    chunk_max <- suppressWarnings(max(v4, na.rm = TRUE))
    if (is.finite(chunk_max) && chunk_max > m) m <- chunk_max
    is_sp <- v5 == "Spacer" & (v3 - v2) > 0L
    if (any(is_sp)) spacer_lids <- unique(c(spacer_lids, v4[is_sp]))
  }
  if (m < 0) stop("could not parse V4 from ", intersect_path)
  list(n_loci = m + 1L, n_loci_spacer = length(unique(spacer_lids)))
}

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

# Per-locus mask: TRUE where the Spacer row has start < end. Read from the
# split BED rather than the intersect, so that real spacer loci with zero SNP
# overlaps stay in the denominator. Mirrors `_spacer_loci_mask` on the
# detected-vs-control path.
# Loci with id >= n_loci (e.g. a trailing locus absent from the intersect) are
# dropped so the mask length matches the matrix row count.
spacer_loci_mask <- function(split_bed_path, n_loci) {
  s <- read.table(split_bed_path, sep = "\t", header = FALSE,
                  stringsAsFactors = FALSE,
                  colClasses = c(V1 = "character", V2 = "integer",
                                 V3 = "integer", V4 = "integer", V5 = "character"))
  s <- s[s$V5 == "Spacer" & (s$V3 - s$V2) > 0, ]
  mask <- rep(FALSE, n_loci)
  ids  <- s$V4[s$V4 < n_loci] + 1L
  mask[ids] <- TRUE
  mask
}


# counts + matrices for detected & undetected
counts    <- read_counts("split_counts.txt")
N_DET     <- counts["N_MOTIF"]
und_info  <- scan_loci_counts("NCNR.DirectRepeats.gnomAD.undet.intersect")
N_UNDET         <- und_info$n_loci
N_UNDET_SPACER  <- und_info$n_loci_spacer
cat(sprintf("N_DETECTED = %d   N_UNDETECTED = %d   N_UNDET_SPACER = %d\n",
            N_DET, N_UNDET, N_UNDET_SPACER))

# Persist N_UNDET / N_UNDET_SPACER to split_counts.txt.
counts_path <- "split_counts.txt"
counts_lines <- readLines(counts_path)
new_lines <- counts_lines
if (!any(grepl("^N_UNDET\\b",        new_lines)))
  new_lines <- c(new_lines, sprintf("N_UNDET\t%d",        N_UNDET))
if (!any(grepl("^N_UNDET_SPACER\\b", new_lines)))
  new_lines <- c(new_lines, sprintf("N_UNDET_SPACER\t%d", N_UNDET_SPACER))
if (!identical(new_lines, counts_lines)) {
  writeLines(new_lines, counts_path)
  cat("appended N_UNDET / N_UNDET_SPACER to split_counts.txt\n")
}

det <- build_matrices("NCNR.DirectRepeats.gnomAD.intersect",        N_DET)
und <- build_matrices("NCNR.DirectRepeats.gnomAD.undet.intersect",  N_UNDET)

cat("\nmatrix shapes (detected):\n")
for (nm in names(det)) cat(sprintf("  %-12s : %d x %d\n", nm, nrow(det[[nm]]), ncol(det[[nm]])))
cat("matrix shapes (undetected):\n")
for (nm in names(und)) cat(sprintf("  %-12s : %d x %d\n", nm, nrow(und[[nm]]), ncol(und[[nm]])))


# IWTomicsTest wrapper
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
  cat(sprintf("%-12s  %4d bins, %5d det, %6d undet, %4d significant\n",
              label, nbins, n_det, n_und, sum(adj < ALPHA)))
  list(adj = as.numeric(adj), nbins = nbins)
}


# run all 5 regions
cat("\nrunning IWTomicsTest on left flank ...\n")
left  <- run_iwt(det$left,  und$left,  "left flank")

cat("\nrunning IWTomicsTest on right flank ...\n")
right <- run_iwt(det$right, und$right, "right flank")

cat("\nrunning IWTomicsTest on LeftRepeat / Spacer / RightRepeat ...\n")
LR <- run_iwt(det$LeftRepeat,  und$LeftRepeat,  "LeftRepeat")

# Restrict to loci with a non-zero spacer before testing, so the tested curves
# match the plotted spacer frequency (divided by N_*_SPACER, not N_*). Left in,
# the zero-spacer loci enter as all-zero rows and dilute the means: the spacer
# panel then shows a large detected-vs-undetected gap with almost no
# significant bins.
det_sp_mask <- spacer_loci_mask("NCNR.DirectRepeats.split.bed",       nrow(det$Spacer))
und_sp_mask <- spacer_loci_mask("NCNR.DirectRepeats.undet.split.bed", nrow(und$Spacer))
cat(sprintf("Spacer loci kept: det %d/%d, undet %d/%d\n",
            sum(det_sp_mask), length(det_sp_mask),
            sum(und_sp_mask), length(und_sp_mask)))
SP <- run_iwt(det$Spacer[det_sp_mask, , drop = FALSE],
              und$Spacer[und_sp_mask, , drop = FALSE], "Spacer")
RR <- run_iwt(det$RightRepeat, und$RightRepeat, "RightRepeat")


# save .RData files
adj_pvalue_curve_left <- data.frame(
  snp_coord   = seq(-left$nbins, -1),
  adj_pvalue  = left$adj,
  significant = left$adj < ALPHA
)
save(adj_pvalue_curve_left,
     file = "DirectRepeats_detected_vs_undetected_IWTomics_left_0.01.RData")

adj_pvalue_curve_right <- data.frame(
  snp_coord   = seq_len(right$nbins),
  adj_pvalue  = right$adj,
  significant = right$adj < ALPHA
)
save(adj_pvalue_curve_right,
     file = "DirectRepeats_detected_vs_undetected_IWTomics_right_0.01.RData")

adj_pvalue_curve_center <- rbind(
  data.frame(subpart = "LeftRepeat",  bin = seq_len(LR$nbins),
             adj_pvalue = LR$adj, significant = LR$adj < ALPHA),
  data.frame(subpart = "Spacer",      bin = seq_len(SP$nbins),
             adj_pvalue = SP$adj, significant = SP$adj < ALPHA),
  data.frame(subpart = "RightRepeat", bin = seq_len(RR$nbins),
             adj_pvalue = RR$adj, significant = RR$adj < ALPHA)
)
save(adj_pvalue_curve_center,
     file = "DirectRepeats_detected_vs_undetected_IWTomics_center_0.01.RData")


# summary
files <- list.files(WORKDIR,
                    pattern = "DirectRepeats_detected_vs_undetected_IWTomics.*\\.RData$",
                    full.names = TRUE)
cat("\nIWTomics outputs in workdir:\n")
print(file.info(files)[, c("size", "mtime")])

cat("\nSummary of significant bins (alpha =", ALPHA, "), detected vs undetected:\n")
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
