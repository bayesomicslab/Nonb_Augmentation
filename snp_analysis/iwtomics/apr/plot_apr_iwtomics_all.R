# 7-panel R + ggplot plot for ALL labelled A-Phased Repeats (mirror of paper §6.6)
# Reads the .RData from run_iwtomics_apr_all.R and the iwt_freq_*.csv written
# by build_freq_csvs_all.py.

suppressMessages({
  library(ggplot2); library(grid); library(gridExtra)
})

setwd(file.path(Sys.getenv("SNP_AF_BASE", "/path/to/snp_af_analysis"), "results", "apr_notebook_run"))

load("APhasedRepeats_motif_vs_control_IWTomics_left_0.01.RData")
load("APhasedRepeats_motif_vs_control_IWTomics_right_0.01.RData")
load("APhasedRepeats_motif_vs_control_IWTomics_center_0.01.RData")

freqL <- read.csv("iwt_freq_LeftFlank.csv")
freqR <- read.csv("iwt_freq_RightFlank.csv")
freqC <- read.csv("iwt_freq_center.csv")

LeftFlankDf  <- data.frame(Levels = seq(-100, -1, 1),
                           Freq   = freqL$LeftFlank,
                           C_Freq = freqL$ctrl_LeftFlank,
                           Signif = adj_pvalue_curve_left$significant)
RightFlankDf <- data.frame(Levels = seq(1, 100, 1),
                           Freq   = freqR$RightFlank,
                           C_Freq = freqR$ctrl_RightFlank,
                           Signif = adj_pvalue_curve_right$significant)

mk_part <- function(name, n_bins) {
  sig_rows <- adj_pvalue_curve_center[adj_pvalue_curve_center$subpart == name, ]
  data.frame(Levels = seq_len(n_bins),
             Freq   = freqC[[name]][1:n_bins],
             C_Freq = freqC[[paste0("ctrl_", name)]][1:n_bins],
             Signif = sig_rows$significant)
}
T1Df <- mk_part("Tract1",  40)
S1Df <- mk_part("Spacer1", 30)
T2Df <- mk_part("Tract2",  40)
S2Df <- mk_part("Spacer2", 30)
T3Df <- mk_part("Tract3",  40)

allv <- c(LeftFlankDf$Freq, LeftFlankDf$C_Freq,
          RightFlankDf$Freq, RightFlankDf$C_Freq,
          T1Df$Freq, T1Df$C_Freq, S1Df$Freq, S1Df$C_Freq,
          T2Df$Freq, T2Df$C_Freq, S2Df$Freq, S2Df$C_Freq,
          T3Df$Freq, T3Df$C_Freq)
allv <- allv[is.finite(allv) & allv > 0]
y_lo <- min(allv) * 0.5
y_hi <- max(allv) * 2

Y_BREAKS <- c(0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.5, 1)
Y_BREAKS <- Y_BREAKS[Y_BREAKS >= y_lo & Y_BREAKS <= y_hi]

cols <- c("non-B DNA"               = rgb(t(col2rgb("blue")),  alpha = 100, maxColorValue = 255),
          "significant\ndifference" = "grey88",
          "control"                 = rgb(t(col2rgb("black")), alpha = 100, maxColorValue = 255))

p1 <- ggplot(LeftFlankDf, aes(x = Levels)) +
  geom_vline(xintercept = LeftFlankDf$Levels[LeftFlankDf$Signif], col = "gray88") +
  geom_point(aes(y = C_Freq, col = "control"),   size = 1) +
  geom_point(aes(y = Freq,   col = "non-B DNA"), size = 1) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title  = element_text(size = 8, family = "sans"),
        axis.text.x = element_text(size = 8, family = "sans"),
        axis.text.y = element_text(size = 8, family = "sans")) +
  labs(y = "SNP frequency", x = "upstream sequence (bp)") +
  scale_colour_manual(name = "", values = cols) +
  scale_x_continuous(breaks = c(-100, -50, -1), limits = c(-100, -1)) +
  theme(plot.margin = unit(c(20, 5, 1, 1), "pt")) +
  scale_y_continuous(trans = "log", limits = c(y_lo, y_hi), breaks = Y_BREAKS)

base_panel <- function(df, xlab, n_bins, vline_size = 1.1) {
  ggplot(df, aes(x = Levels)) +
    geom_vline(xintercept = df$Levels[df$Signif], col = "gray88", size = vline_size) +
    geom_point(aes(y = C_Freq, col = "control"),   size = 1) +
    geom_point(aes(y = Freq,   col = "non-B DNA"), size = 1) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x  = element_text(size = 8, family = "sans"),
          axis.title.x = element_text(size = 8, family = "sans"),
          axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank()) +
    labs(x = xlab) +
    scale_colour_manual(name = "", values = cols) +
    scale_x_continuous(breaks = c(1, n_bins), limits = c(1, n_bins)) +
    theme(plot.margin = unit(c(20, 5, 1, 5), "pt")) +
    scale_y_continuous(trans = "log", limits = c(y_lo, y_hi), breaks = Y_BREAKS)
}
# Tracts have 40 bins and spacers 30, so sizes 1.1 / 1.0 render each
# significant bin as a distinct grey block.
p2 <- base_panel(T1Df, "tract",  40, 1.1)
p3 <- base_panel(S1Df, "spacer", 30, 1.0)
p4 <- base_panel(T2Df, "tract",  40, 1.1)
p5 <- base_panel(S2Df, "spacer", 30, 1.0)
p6 <- base_panel(T3Df, "tract",  40, 1.1)

p7 <- ggplot(RightFlankDf, aes(x = Levels)) +
  geom_line(aes(x = 0, y = 0, colour = "significant\ndifference")) +
  geom_vline(xintercept = RightFlankDf$Levels[RightFlankDf$Signif], col = "gray88") +
  geom_point(aes(y = C_Freq, col = "control"),   size = 1) +
  geom_point(aes(y = Freq,   col = "non-B DNA"), size = 1) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(),
        axis.text.x  = element_text(size = 8, family = "sans"),
        axis.title.x = element_text(size = 8, family = "sans"),
        legend.text  = element_text(size = 8, family = "sans"),
        legend.position      = c(0.99, 0.98),
        legend.justification  = c("right", "top"),
        legend.background     = element_rect(fill = "white", colour = NA)) +
  labs(x = "downstream sequence (bp)") +
  scale_colour_manual(name = "", values = cols,
                      breaks = c("non-B DNA", "control", "significant\ndifference")) +
  scale_x_continuous(breaks = c(1, 50, 100), limits = c(1, 100)) +
  theme(plot.margin = unit(c(20, 5, 1, 5), "pt")) +
  scale_y_continuous(trans = "log", limits = c(y_lo, y_hi), breaks = Y_BREAKS) +
  guides(shape = guide_legend(order = 3))

g <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7, ncol = 7,
                 widths = unit(c(60, 22, 16, 22, 16, 22, 60),
                               c("mm","mm","mm","mm","mm","mm","mm")),
                 top = textGrob("A-phased repeats (all labelled) vs control",
                                gp = gpar(fontsize = 12, family = "sans")))

ggsave("apr_snp_frequency_notebook_iwtomics.pdf", g, width = 11.5, height = 3.5, dpi = 300)
ggsave("apr_snp_frequency_notebook_iwtomics.png", g, width = 11.5, height = 3.5, dpi = 200)
cat("saved apr_snp_frequency_notebook_iwtomics.{pdf,png}\n")

cat(sprintf("\nLEFT  flank: %d/100 sig (-1 adj_p=%.4f)\n",
            sum(adj_pvalue_curve_left$significant),
            adj_pvalue_curve_left$adj_pvalue[100]))
cat(sprintf("RIGHT flank: %d/100 sig (+1 adj_p=%.4f)\n",
            sum(adj_pvalue_curve_right$significant),
            adj_pvalue_curve_right$adj_pvalue[1]))
for (sp in c("Tract1","Spacer1","Tract2","Spacer2","Tract3")) {
  rows <- adj_pvalue_curve_center[adj_pvalue_curve_center$subpart == sp, ]
  cat(sprintf("%-7s : %d/%d sig\n", sp, sum(rows$significant), nrow(rows)))
}
