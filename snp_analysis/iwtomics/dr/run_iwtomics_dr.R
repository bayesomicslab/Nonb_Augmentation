# IWTomics on Direct Repeats (model-detected)
# Run with Rscript, or block by block in an interactive session.
# Reference: IWTomics 1.36.0 manual
#   https://bioconductor.org/packages/release/bioc/manuals/IWTomics/man/IWTomics.pdf
#
# Inputs (written by the notebook):
#   $SNP_AF_BASE/results/dr_notebook_detected_run/
#     iwt_motif_left.tsv  iwt_ctrl_left.tsv     (loci x 100)
#     iwt_motif_right.tsv iwt_ctrl_right.tsv    (loci x 100)
#     iwt_motif_LR.tsv    iwt_ctrl_LR.tsv       (loci x 50)
#     iwt_motif_Sp.tsv    iwt_ctrl_Sp.tsv       (loci x 80; only non-zero-spacer loci)
#     iwt_motif_RR.tsv    iwt_ctrl_RR.tsv       (loci x 50)
#
# Outputs (saved into the same workdir so the notebook plot cell finds them):
#   DirectRepeats_motif_vs_control_IWTomics_left_0.01.RData
#   DirectRepeats_motif_vs_control_IWTomics_right_0.01.RData
#   DirectRepeats_motif_vs_control_IWTomics_center_0.01.RData


# load packages, verify versions
# IWTomics depends on GenomicRanges (manual p.1: "Depends: GenomicRanges").
library(IWTomics)
library(GenomicRanges)
cat("IWTomics version:        ", as.character(packageVersion("IWTomics")), "\n")
cat("GenomicRanges version:   ", as.character(packageVersion("GenomicRanges")), "\n")
cat("R version:               ", R.version.string, "\n")


# paths + IWTomics parameters
# WORKDIR = where the notebook wrote the TSVs and where the .RData files will land
# (so the notebook plot cell can `load()` them without changing paths).
WORKDIR <- file.path(Sys.getenv("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "dr_notebook_detected_run")
setwd(WORKDIR)
cat("workdir:", getwd(), "\n")

ALPHA  <- 0.01      # significance threshold (matches reference notebook _0.01)
N_PERM <- 1000      # IWTomicsTest default B (manual p.12)
ALIGN  <- "center"  # all regions have identical width => center alignment.
                    # Not "scale": manual p.13 says "scale" needs smooth() first.


# helper to read a loci-by-bins TSV into a positions-by-loci matrix
# IWTomics features matrix: "columns of aligned feature measurements" (manual p.8).
# So each column = one locus's curve, each row = one position. The notebook saves them as
# loci x bins (one row per locus), so the matrix is transposed.
read_loci_x_bins <- function(path) {
  m <- as.matrix(read.table(path, sep = "\t", header = FALSE,
                            comment.char = "", check.names = FALSE))
  t(m)   # -> bins x loci  (positions x loci)
}


# Left flank, worked through explicitly; the other four regions reuse the
# function assembled from these steps below.

# read motif & control matrices for left flank
m_left <- read_loci_x_bins("iwt_motif_left.tsv")
c_left <- read_loci_x_bins("iwt_ctrl_left.tsv")
cat("left flank  motif:  ",   nrow(m_left), "positions x", ncol(m_left), "loci\n")
cat("left flank  control:",   nrow(c_left), "positions x", ncol(c_left), "loci\n")

# build dummy GRanges for each region dataset
# IWTomicsData needs a GRangesList of region locations (manual p.6 "regions" slot).
# The coordinates are arbitrary; only the equal width per locus matters,
# width (= number of positions) so center-alignment is unambiguous.
nbins  <- nrow(m_left)
n_mot  <- ncol(m_left)
n_ctrl <- ncol(c_left)

gr_mot  <- GRanges("chr1",
                   IRanges(start = seq(1, by = nbins, length.out = n_mot),
                           width = nbins))
gr_ctrl <- GRanges("chr1",
                   IRanges(start = seq(1, by = nbins, length.out = n_ctrl),
                           width = nbins))
gr_list <- GRangesList(motif = gr_mot, control = gr_ctrl)
cat("GRangesList lengths:", lengths(gr_list), "\n")

# construct IWTomicsData object
# Manual p.7 alternative constructor:
#   IWTomicsData(x, y, alignment='center', id_regions=NULL, name_regions=NULL,
#                id_features=NULL, name_features=NULL, length_features=NULL)
# x = GRangesList, y = list-of-matrix-lists named by feature -> region.
rf_left <- IWTomicsData(
  gr_list,                                       # x
  list(snp = list(motif = m_left, control = c_left)),  # y
  alignment     = ALIGN,
  id_regions    = c("motif", "control"),
  name_regions  = c("motif", "control"),
  id_features   = "snp",
  name_features = "snp"
)

# run the two-sample Interval-Wise Test
# Manual p.12: IWTomicsTest(regionsFeatures, id_region1, id_region2,
#                          id_features_subset, mu=0, statistics="mean",
#                          probs=0.5, max_scale=NULL, paired=FALSE, B=1000)
# motif vs control on the "snp" feature, mean statistic, B=1000.
res_left <- IWTomicsTest(
  rf_left,
  id_region1         = "motif",
  id_region2         = "control",
  id_features_subset = "snp",
  statistics         = "mean",
  B                  = N_PERM
)

# extract adjusted p-values
# Manual p.9: adjusted_pval(x, test, id_features_subset, scale_threshold).
# Default scale_threshold = NULL -> use max possible interval length (= nbins).
# Return shape: list(test1 = list(snp = numeric_vector)).
adj_left <- adjusted_pval(res_left)[[1]][[1]]
cat("left adj p-values length:", length(adj_left), " (expected", nbins, ")\n")
cat("left significant bins (alpha =", ALPHA, "):",
    sum(adj_left < ALPHA), "/", nbins, "\n")

# assemble adj_pvalue_curve_left (matches reference notebook structure)
adj_pvalue_curve_left <- data.frame(
  snp_coord   = seq(-nbins, -1),
  adj_pvalue  = as.numeric(adj_left),
  significant = as.logical(adj_left < ALPHA)
)

# save to .RData (filename matches the reference notebook convention)
save(adj_pvalue_curve_left,
     file = "DirectRepeats_motif_vs_control_IWTomics_left_0.01.RData")
cat("wrote DirectRepeats_motif_vs_control_IWTomics_left_0.01.RData\n")


# Apply the same procedure to the other four regions.
# (right flank, LeftRepeat, Spacer, RightRepeat)
run_iwt <- function(motif_path, ctrl_path, label) {
  m <- read_loci_x_bins(motif_path)
  c <- read_loci_x_bins(ctrl_path)
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

right <- run_iwt("iwt_motif_right.tsv", "iwt_ctrl_right.tsv", "right flank")
LR    <- run_iwt("iwt_motif_LR.tsv",    "iwt_ctrl_LR.tsv",    "LeftRepeat")
Sp    <- run_iwt("iwt_motif_Sp.tsv",    "iwt_ctrl_Sp.tsv",    "Spacer")
RR    <- run_iwt("iwt_motif_RR.tsv",    "iwt_ctrl_RR.tsv",    "RightRepeat")

# save right flank
adj_pvalue_curve_right <- data.frame(
  snp_coord   = seq(1, right$nbins),
  adj_pvalue  = right$adj,
  significant = right$adj < ALPHA
)
save(adj_pvalue_curve_right,
     file = "DirectRepeats_motif_vs_control_IWTomics_right_0.01.RData")
cat("wrote DirectRepeats_motif_vs_control_IWTomics_right_0.01.RData\n")


# Combine the three center sub-parts into adj_pvalue_curve_center.
# (matches reference notebook variable layout)
adj_pvalue_curve_center <- rbind(
  data.frame(subpart = "LeftRepeat",  bin = seq_len(LR$nbins),
             adj_pvalue = LR$adj, significant = LR$adj < ALPHA),
  data.frame(subpart = "Spacer",      bin = seq_len(Sp$nbins),
             adj_pvalue = Sp$adj, significant = Sp$adj < ALPHA),
  data.frame(subpart = "RightRepeat", bin = seq_len(RR$nbins),
             adj_pvalue = RR$adj, significant = RR$adj < ALPHA)
)
save(adj_pvalue_curve_center,
     file = "DirectRepeats_motif_vs_control_IWTomics_center_0.01.RData")
cat("wrote DirectRepeats_motif_vs_control_IWTomics_center_0.01.RData\n")


# Verify outputs.
files <- list.files(WORKDIR, pattern = "DirectRepeats.*IWTomics.*\\.RData$",
                    full.names = TRUE)
cat("\nIWTomics outputs in workdir:\n")
print(file.info(files)[, c("size", "mtime")])

cat("\nSummary of significant bins (alpha =", ALPHA, "):\n")
cat(sprintf("  left  : %d/%d\n",
            sum(adj_pvalue_curve_left$significant),  nrow(adj_pvalue_curve_left)))
cat(sprintf("  right : %d/%d\n",
            sum(adj_pvalue_curve_right$significant), nrow(adj_pvalue_curve_right)))
for (sp in c("LeftRepeat", "Spacer", "RightRepeat")) {
  s <- adj_pvalue_curve_center[adj_pvalue_curve_center$subpart == sp, ]
  cat(sprintf("  %-11s: %d/%d\n", sp, sum(s$significant), nrow(s)))
}
cat("
Done.
")
