# Official IWTomics on STR, detected vs undetected, stratified by k-mer.
# Mirrors run_iwtomics_str_dvu.R but loops over k-mers and reads the per-kmer
# count matrices written by str_bykmer_prep.py.
# Output: str_bykmer_iwt.csv  (region,kmer,x,adj_pvalue,significant)
suppressMessages({ library(IWTomics); library(GenomicRanges) })

WORK <- file.path(Sys.getenv("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "str_bykmer_run")
setwd(WORK)
ALPHA <- 0.01; N_PERM <- 1000; KMERS <- c(1,2,3,4,5,6)
REGIONS <- c("L","C","R")

run_iwt <- function(det_mat, und_mat) {
  # input: loci x bins ; IWTomics wants bins x loci
  m <- t(det_mat); c <- t(und_mat)
  nbins <- nrow(m)
  gr_m <- GRanges("chr1", IRanges(seq(1, by=nbins, length.out=ncol(m)), width=nbins))
  gr_c <- GRanges("chr1", IRanges(seq(1, by=nbins, length.out=ncol(c)), width=nbins))
  rf <- IWTomicsData(GRangesList(detected=gr_m, undetected=gr_c),
                     list(snp=list(detected=m, undetected=c)),
                     alignment="center",
                     id_regions=c("detected","undetected"),
                     name_regions=c("detected","undetected"),
                     id_features="snp", name_features="snp")
  res <- IWTomicsTest(rf, id_region1="detected", id_region2="undetected",
                      id_features_subset="snp", statistics="mean", B=N_PERM)
  as.numeric(adjusted_pval(res)[[1]][[1]])
}

x_of <- function(reg, nbins) {
  if (reg=="L") return(seq(-nbins, -1))   # -100..-1
  return(seq_len(nbins))                   # C: 1..180, R: 1..100
}

out <- list()
for (k in KMERS) for (reg in REGIONS) {
  fd <- sprintf("mat_det_k%d_%s.txt.gz", k, reg)
  fu <- sprintf("mat_und_k%d_%s.txt.gz", k, reg)
  if (!file.exists(fd) || !file.exists(fu)) next
  det_mat <- as.matrix(read.table(gzfile(fd)))
  und_mat <- as.matrix(read.table(gzfile(fu)))
  if (nrow(det_mat) < 2 || nrow(und_mat) < 2) next
  adj <- run_iwt(det_mat, und_mat)
  nbins <- length(adj)
  out[[paste(k,reg)]] <- data.frame(region=reg, kmer=k, x=x_of(reg,nbins),
                                    adj_pvalue=adj, significant=adj<ALPHA)
  cat(sprintf("k%d %s: %d bins, %d significant\n", k, reg, nbins, sum(adj<ALPHA)))
}
res <- do.call(rbind, out)
write.csv(res, "str_bykmer_iwt.csv", row.names=FALSE)
cat("wrote str_bykmer_iwt.csv\n")
