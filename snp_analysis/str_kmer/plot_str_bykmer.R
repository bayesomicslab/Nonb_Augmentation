# ggplot: STR SNP frequency along motif+flanks, faceted by k-mer (rows) x region (cols).
# detected = solid, undetected = dashed (colour = k-mer), control = grey points,
# grey vlines = IWTomics(det vs undet) significant bins.
suppressMessages({ library(ggplot2); library(scales) })

WORK <- file.path(Sys.getenv("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "str_bykmer_run")
setwd(WORK)

freq <- read.csv("str_bykmer_freq.csv")
iwt  <- read.csv("str_bykmer_iwt.csv")

reg_lab <- c(L="upstream (bp)", C="scaled STR locus", R="downstream (bp)")
freq$region <- factor(freq$region, levels=c("L","C","R"), labels=reg_lab)
iwt$region  <- factor(iwt$region,  levels=c("L","C","R"), labels=reg_lab)
freq$kmer_f <- factor(paste0(freq$kmer,"-mer"), levels=paste0(sort(unique(freq$kmer)),"-mer"))
iwt$kmer_f  <- factor(paste0(iwt$kmer,"-mer"),  levels=levels(freq$kmer_f))

# long form for detected/undetected (colour=k-mer, linetype=signal)
long <- rbind(
  data.frame(region=freq$region, kmer_f=freq$kmer_f, kmer=freq$kmer, x=freq$x,
             signal="detected",   value=freq$detected),
  data.frame(region=freq$region, kmer_f=freq$kmer_f, kmer=freq$kmer, x=freq$x,
             signal="undetected", value=freq$undetected))
long <- long[is.finite(long$value) & long$value>0, ]
ctrl <- freq[is.finite(freq$control) & freq$control>0, ]
sigdf <- iwt[iwt$significant, ]

p <- ggplot() +
  geom_vline(data=sigdf, aes(xintercept=x), colour="grey88", linewidth=0.3) +
  geom_point(data=ctrl, aes(x=x, y=control), colour="grey40", size=0.35) +
  geom_line(data=long, aes(x=x, y=value, colour=factor(kmer), linetype=signal),
            linewidth=0.6) +
  facet_grid(kmer_f ~ region, scales="free") +
  scale_y_continuous(trans="log10", labels=label_number()) +
  scale_linetype_manual(name="", values=c(detected="solid", undetected="dashed")) +
  scale_colour_brewer(name="k-mer", palette="Set1", guide="none") +
  labs(x="position", y="SNP frequency",
       title="STR SNP frequency along motif + flanks, by k-mer (detected / undetected / control)") +
  theme_bw(base_size=9) +
  theme(panel.grid.minor=element_blank(),
        strip.background=element_rect(fill="grey95"),
        legend.position="top")

nk <- length(unique(freq$kmer)); H <- 1.6*nk + 1
ggsave("str_snp_frequency_by_kmer_iwtomics.pdf", p, width=11, height=H, bg="white")
ggsave("str_snp_frequency_by_kmer_iwtomics.png", p, width=11, height=H, dpi=200, bg="white")
cat("saved str_snp_frequency_by_kmer_iwtomics.{pdf,png}\n")
