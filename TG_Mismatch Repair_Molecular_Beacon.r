# =============================================================================
#  T:G Mismatch Repair Molecular Beacon Assay
# =============================================================================
#
#  Purpose
#  -------
#  Plot real-time fluorescence traces from a molecular-beacon T:G mismatch
#  repair assay in three HLRCC cell-line models, NCCFH1, UOK262 and UOK350
#  with and without doxycycline induction of FH (Ctrl vs Dox), against methylated
#  and unmethylated DNA substrates. The molecular-beacon probe contains a 5'
#  flurophore located in a hair pin structure next to a 3' quencher. Upon
#  succesful repair of the TG mismatch the 5' flurophore is released giving a
#  real time fluorescence trace of TG repair kinetics.
#
#  Output is a 2 x 3 panel figure: rows = substrate (methylated / unmethylated),
#  columns = cell line (NCCFH1 / UOK262 / UOK350). Each panel overlays Ctrl (red)
#  and Dox (blue) traces with LOESS smoothing and a +/- SD ribbon.
#
#  Inputs (CSV files in the working directory)
#  -------------------------------------------
#  All input data are normalized to a maximum fluorescent value derived from a
#  melt curve analysis to normalize for pippeting error in probe loading.
#
#    <CellLine>_Ctrl_Dox_Meth.csv     methylated substrate, Ctrl + Dox columns
#    <CellLine>_Ctrl_Dox_Unmeth.csv   unmethylated substrate, Ctrl + Dox columns
#
#    where <CellLine> is one of: NCCFH1, UOK262, UOK350.
#
#  Each CSV has one row per timepoint and columns containing "Ctrl" or "Dox"
#  in their name for the replicate measurements (e.g. Ctrl_1, Ctrl_2, ...).
#
#  Outputs
#  -------
#  NCCFH1_UOK262_UOK350_Ctrl_Dox_Meth_Unmeth.pdf   composite 2 x 3 panel figure
#
#  Dependencies
#  ------------
#  ggplot2, ggpubr
#
#  Paper
#  -----
#  Metabolite-Driven Succination of MBD4 Links FH Loss to CpG Mutational Bias
#  in HLRCC.
#
#  Author: Ben Stansfield
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Environment
# -----------------------------------------------------------------------------
rm(list = ls())

library(ggplot2)
library(ggpubr)

# Working directory containing the raw beacon CSVs.
# Update this path to match your local environment.
setwd("")

# -----------------------------------------------------------------------------
# 1. Helpers
# -----------------------------------------------------------------------------

# Add per-timepoint mean and SD across Ctrl replicates and Dox replicates,
# identified by column-name prefix, plus a Measurement index column.
summarise_traces <- function(df) {
  ctrl_cols <- grep("Ctrl", colnames(df))
  dox_cols  <- grep("Dox",  colnames(df))

  df$ctrl_mean   <- rowMeans(df[, ctrl_cols])
  df$dox_mean    <- rowMeans(df[, dox_cols])
  df$ctrl_sd     <- apply(df[, ctrl_cols], 1, sd)
  df$dox_sd      <- apply(df[, dox_cols],  1, sd)
  df$Measurement <- seq_len(nrow(df))
  df
}

# Build a single-panel ggplot: LOESS-smoothed Ctrl (red) and Dox (blue) traces
# with +/- SD ribbons. Legend is suppressed for the composite figure.
plot_traces <- function(df) {
  ggplot(df, aes(x = Measurement)) +
    geom_smooth(aes(y = ctrl_mean), color = "red",
                method = "loess", se = FALSE, lwd = 1) +
    geom_ribbon(aes(ymin = ctrl_mean - ctrl_sd,
                    ymax = ctrl_mean + ctrl_sd),
                fill = "red", alpha = 0.2) +
    geom_smooth(aes(y = dox_mean), color = "blue",
                method = "loess", se = FALSE, lwd = 1) +
    geom_ribbon(aes(ymin = dox_mean - dox_sd,
                    ymax = dox_mean + dox_sd),
                fill = "blue", alpha = 0.2) +
    theme_classic() +
    theme(legend.position = "none") +
    scale_y_continuous(expand = c(0, 0), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
    scale_x_continuous(expand = c(0, 0), limits = c(1, 87)) +
    labs(x = "Measurement", y = "Relative Fluorescence")
}

# -----------------------------------------------------------------------------
# 2. Read fluorescence data for every cell line x substrate combination
# -----------------------------------------------------------------------------
cell_lines <- c("NCCFH1", "UOK262", "UOK350")
substrates <- c("Meth", "Unmeth")

# Named list keyed by "<CellLine>_<Substrate>" (e.g. "NCCFH1_Meth")
traces <- list()
for (cell in cell_lines) {
  for (sub in substrates) {
    key <- paste(cell, sub, sep = "_")
    traces[[key]] <- read.csv(paste0(cell, "_Ctrl_Dox_", sub, ".csv"),
                              header = TRUE)
  }
}

# -----------------------------------------------------------------------------
# 3. Summarise replicates and build per-panel ggplots
# -----------------------------------------------------------------------------
traces <- lapply(traces, summarise_traces)
plots  <- lapply(traces, plot_traces)

# -----------------------------------------------------------------------------
# 4. Compose 2 x 3 panel figure and save
# -----------------------------------------------------------------------------
# Row 1: methylated substrate    (NCCFH1, UOK262, UOK350)
# Row 2: unmethylated substrate  (NCCFH1, UOK262, UOK350)
p1 <- ggarrange(plots[["NCCFH1_Meth"]],   plots[["UOK262_Meth"]],   plots[["UOK350_Meth"]],
                plots[["NCCFH1_Unmeth"]], plots[["UOK262_Unmeth"]], plots[["UOK350_Unmeth"]],
                ncol = 3, nrow = 2, align = "hv")

ggsave(p1, file = "NCCFH1_UOK262_UOK350_Ctrl_Dox_Meth_Unmeth.pdf",
       width = 9, height = 5)