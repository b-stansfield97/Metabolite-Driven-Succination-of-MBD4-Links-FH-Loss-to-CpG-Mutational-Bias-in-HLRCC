# =============================================================================
#  TCGA UCEC TET-Mutant SBS Mutational Signature Analysis
# =============================================================================
#
#  Purpose
#  -------
#  Comparator analysis for the MBD4-mutant cohorts: extract de novo SBS
#  signatures from TCGA UCEC tumours carrying somatic mutations in any of
#  TET1, TET2 or TET3. TET-mutant tumours are expected to enrich for 5mC
#  deamination–driven mutational processes (e.g. COSMIC SBS1) and serve as a
#  reference for interpreting the CpG bias seen with MBD4 loss.
#
#  Inputs
#  ------
#  TCGA UCEC ensemble-masked MAF files in the working directory matching
#  "*ensemble_masked.maf.gz". Each MAF is expected to contain a `CONTEXT`
#  column (11 nt centred on the variant).
#
#  Outputs
#  -------
#    Tet_Signature_Estimate.pdf          goodness-of-fit sweep
#    tet_enriched_signatures.pdf         extracted signature spectra
#    tet_exposures.pdf                   refit exposures
#    Tet_mut_sig.rda                     workspace image 
#
#  Dependencies
#  ------------
#  sigfit (>= 2.0), pheatmap
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

# -----------------------------------------------------------------------------
# 1. Helper: read and concatenate MAF files
# -----------------------------------------------------------------------------
read_data <- function(file_paths) {
  data_list <- list()

  for (i in seq_along(file_paths)) {
    df <- read.delim2(file_paths[i],
                      sep          = "\t",
                      comment.char = "#",
                      colClasses   = "character")
    data_list[[i]] <- df
    print(i)
  }

  data <- do.call(rbind, data_list)
  return(data)
}

# -----------------------------------------------------------------------------
# 2. Load UCEC MAFs and filter to SNVs
# -----------------------------------------------------------------------------
UCEC_path <- dir(pattern = "*ensemble_masked.maf.gz")
UCEC_data <- read_data(UCEC_path)

UCEC_data_filtered <- UCEC_data[which(UCEC_data$VARIANT_CLASS == "SNV"), ]

# -----------------------------------------------------------------------------
# 3. Retain TET1/TET2/TET3-mutant tumours and build trinuc context
# -----------------------------------------------------------------------------
genes_of_interest <- c("TET1", "TET2", "TET3")

# Identify tumours with any somatic mutation in a TET gene, then keep
# *all* variants from those tumours (not just the TET variant itself).
samples_of_interest <- unique(UCEC_data_filtered$Tumor_Sample_Barcode[
  UCEC_data_filtered$Hugo_Symbol %in% genes_of_interest])

tet_filtered_data <- UCEC_data_filtered[
  UCEC_data_filtered$Tumor_Sample_Barcode %in% samples_of_interest, ]

# GDC `CONTEXT` is 11 nt centred on the variant; chars 5-7 give the trinuc
tet_filtered_data$trinuc_context <- substr(tet_filtered_data$CONTEXT, 5, 7)

# -----------------------------------------------------------------------------
# 4. Build sigfit input matrix
# -----------------------------------------------------------------------------
tet_matrix <- matrix(0, ncol = 4, nrow = nrow(tet_filtered_data))
colnames(tet_matrix) <- c("Sample", "Ref", "Alt", "Trinuc")

tet_matrix[, "Sample"] <- tet_filtered_data$Tumor_Sample_Barcode
tet_matrix[, "Ref"]    <- tet_filtered_data$Tumor_Seq_Allele1
tet_matrix[, "Alt"]    <- tet_filtered_data$Tumor_Seq_Allele2
tet_matrix[, "Trinuc"] <- tet_filtered_data$trinuc_context

# -----------------------------------------------------------------------------
# 5. De novo signature extraction (sweep nsignatures = 2:10)
# -----------------------------------------------------------------------------
library(sigfit)
data("cosmic_signatures_v2")

counts_tet <- build_catalogues(tet_matrix)

tet_MCMC_signatures <- extract_signatures(counts      = counts_tet,
                                          nsignatures = 2:10,
                                          iter        = 5000,
                                          warmup      = 1000,
                                          seed        = 1756)

pdf(file = "Tet_Signature_Estimate.pdf", height = 10, width = 10)
plot_gof(tet_MCMC_signatures)
dev.off()

## Goodness-of-fit indicates 3 signatures is optimal for this cohort.

# -----------------------------------------------------------------------------
# 6. Refit with optimal nsignatures = 3
# -----------------------------------------------------------------------------
tet_3_MCMC_signatures <- extract_signatures(counts      = counts_tet,
                                            nsignatures = 3,
                                            iter        = 10000,
                                            warmup      = 1000,
                                            seed        = 1756)

tet_extracted_sigs <- retrieve_pars(tet_3_MCMC_signatures, par = "signatures")

# -----------------------------------------------------------------------------
# 7. Match extracted signatures to COSMIC v2
# -----------------------------------------------------------------------------
tet_matched_signatures <- match_signatures(tet_extracted_sigs,
                                           cosmic_signatures_v2,
                                           stat = "cosine")
# Matches in this cohort: sig 1 -> COSMIC 10, sig 2 -> COSMIC 6, sig 3 -> COSMIC 20.

matched_cosmic <- c("Signature 10", "Signature 6", "Signature 20")
cosmic_matches <- cosmic_signatures_v2[matched_cosmic, ]

# -----------------------------------------------------------------------------
# 8. Refit matched COSMIC signatures for accurate exposures
# -----------------------------------------------------------------------------
tet_mcmc_samples_refit <- fit_signatures(counts     = counts_tet,
                                         signatures = cosmic_matches,
                                         iter       = 5000,
                                         warmup     = 1000,
                                         seed       = 1756)

tet_exposures  <- retrieve_pars(tet_mcmc_samples_refit, par = "exposures",  hpd_prob = 0.95)
tet_signatures <- retrieve_pars(tet_3_MCMC_signatures,  par = "signatures")

# -----------------------------------------------------------------------------
# 9. Plot signature spectra and exposures
# -----------------------------------------------------------------------------
pdf(file = "tet_enriched_signatures.pdf", height = 3, width = 9)
plot_spectrum(tet_signatures$mean)
dev.off()

pdf(file = "tet_exposures.pdf")
plot_exposures(tet_mcmc_samples_refit)
dev.off()

save.image("Tet_mut_sig.rda")