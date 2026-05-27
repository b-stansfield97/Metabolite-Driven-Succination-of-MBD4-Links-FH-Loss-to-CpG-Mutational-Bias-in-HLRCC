# =============================================================================
#  HLRCC Whole-Exome SBS Mutational Signature Analysis
# =============================================================================
#
#  Purpose
#  -------
#  Extract de novo single-base-substitution (SBS) mutational signatures from
#  HLRCC tumour whole-exome MAF files, match the extracted signatures to the
#  COSMIC v2 reference catalogue, and visualise per-tumour exposures and
#  cosine similarity to COSMIC.
#
#  Inputs
#  ------
#  HLRCC tumour MAF files in the working directory matching "*.maf".
#  Each MAF is expected to contain a `ref_context` column holding the
#  reference allele +/- 10 nt.
#
#  Outputs
#  -------
#    HLRCC_Signature_Estimate.pdf         goodness-of-fit across nsignatures
#    HLRCC_enriched_signatures.pdf        extracted signature spectra
#    HLRCC_signature_exposure.pdf         per-tumour exposures (refit)
#    HLRCC_cosine_sim_heatmap.pdf         cosine similarity vs COSMIC v2
#    HLRCC_mut_sig.rda                    saved workspace image
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
# Reads every MAF in `file_paths`, tags each with a sequential `sample` id,
# and row-binds them into a single data frame.
read_data <- function(file_paths) {
  data_list <- list()

  for (i in seq_along(file_paths)) {
    df <- read.delim2(file_paths[i],
                      sep          = "\t",
                      comment.char = "#",
                      colClasses   = "character")
    df$sample <- paste0("Tumor_sample_", i)
    data_list[[i]] <- df
    print(i)
  }

  data <- do.call(rbind, data_list)
  return(data)
}

# -----------------------------------------------------------------------------
# 2. Load HLRCC MAFs and filter to SNPs
# -----------------------------------------------------------------------------
HLRCC_path <- dir(pattern = "*.maf")
HLRCC_data <- read_data(HLRCC_path)

# Retain only the columns needed for downstream signature extraction
cols <- c("Hugo_Symbol", "Variant_Classification", "Variant_Type",
          "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "ref_context", "sample")
HLRCC_data <- HLRCC_data[, cols]

# SBS signatures are defined for single-base substitutions only
HLRCC_data_filtered <- HLRCC_data[which(HLRCC_data$Variant_Type == "SNP"), ]

# -----------------------------------------------------------------------------
# 3. Build trinucleotide context for each variant
# -----------------------------------------------------------------------------
# Sanity check: the 11th character of ref_context should equal the reference
# allele (ref_context = reference allele +/- 10 nt).
check <- substr(HLRCC_data_filtered$ref_context, 11, 11)
check == HLRCC_data_filtered$Tumor_Seq_Allele1

# Trinucleotide context = reference allele +/- 1 nt
HLRCC_data_filtered$trinuc_context <- substr(HLRCC_data_filtered$ref_context, 10, 12)

# -----------------------------------------------------------------------------
# 4. Build sigfit input matrix
# -----------------------------------------------------------------------------
# sigfit expects a 4-column matrix: Sample, Ref, Alt, Trinuc.
HLRCC_matrix <- matrix(0, ncol = 4, nrow = nrow(HLRCC_data_filtered))
colnames(HLRCC_matrix) <- c("Sample", "Ref", "Alt", "Trinuc")

HLRCC_matrix[, "Sample"] <- HLRCC_data_filtered$sample
HLRCC_matrix[, "Ref"]    <- HLRCC_data_filtered$Tumor_Seq_Allele1
HLRCC_matrix[, "Alt"]    <- HLRCC_data_filtered$Tumor_Seq_Allele2
HLRCC_matrix[, "Trinuc"] <- HLRCC_data_filtered$trinuc_context

# -----------------------------------------------------------------------------
# 5. De novo signature extraction (sweep nsignatures = 2:10)
# -----------------------------------------------------------------------------
library(sigfit)
data("cosmic_signatures_v2")

# Build the 96-channel trinucleotide-context catalogue
counts_HLRCC <- build_catalogues(HLRCC_matrix)

# Initial sweep to identify the optimal number of signatures via GoF
HLRCC_MCMC_signatures <- extract_signatures(counts      = counts_HLRCC,
                                            nsignatures = 2:10,
                                            iter        = 5000,
                                            warmup      = 1000,
                                            seed        = 1756)

pdf(file = "HLRCC_Signature_Estimate.pdf", height = 10, width = 10)
plot_gof(HLRCC_MCMC_signatures)
dev.off()

## Goodness-of-fit indicates 3 signatures is optimal for this cohort.

# -----------------------------------------------------------------------------
# 6. Refit with optimal nsignatures = 3 (longer MCMC chain)
# -----------------------------------------------------------------------------
HLRCC_3_MCMC_signatures <- extract_signatures(counts      = counts_HLRCC,
                                              nsignatures = 3,
                                              iter        = 50000,
                                              warmup      = 1000,
                                              seed        = 1756)

# Pull out the extracted signature spectra
HLRCC_extracted_sigs <- retrieve_pars(HLRCC_3_MCMC_signatures, par = "signatures")

# -----------------------------------------------------------------------------
# 7. Match extracted signatures to COSMIC v2
# -----------------------------------------------------------------------------
HLRCC_matched_signatures <- match_signatures(HLRCC_extracted_sigs,
                                             cosmic_signatures_v2,
                                             stat = "cosine")
# Matches in this cohort: sig 1 -> COSMIC 1, sig 2 -> COSMIC 19, sig 3 -> COSMIC 3.

matched_cosmic <- c("Signature 1", "Signature 19", "Signature 3")
cosmic_matches <- cosmic_signatures_v2[matched_cosmic, ]

# -----------------------------------------------------------------------------
# 8. Refit matched COSMIC signatures to obtain accurate exposures
# -----------------------------------------------------------------------------
HLRCC_mcmc_samples_refit <- fit_signatures(counts     = counts_HLRCC,
                                           signatures = cosmic_matches,
                                           iter       = 50000,
                                           warmup     = 1000,
                                           seed       = 1756)

HLRCC_exposures  <- retrieve_pars(HLRCC_mcmc_samples_refit, par = "exposures",  hpd_prob = 0.95)
HLRCC_signatures <- retrieve_pars(HLRCC_3_MCMC_signatures,  par = "signatures")

# -----------------------------------------------------------------------------
# 9. Plot signature spectra and per-tumour exposures
# -----------------------------------------------------------------------------
pdf(file = "HLRCC_enriched_signatures.pdf", height = 3, width = 9)
plot_spectrum(HLRCC_signatures$mean)
dev.off()

pdf(file = "HLRCC_signature_exposure.pdf")
plot_exposures(HLRCC_mcmc_samples_refit)
dev.off()

# -----------------------------------------------------------------------------
# 10. Cosine similarity heatmap vs COSMIC v2
# -----------------------------------------------------------------------------
library(pheatmap)

# Row-normalised cosine similarity between two signature matrices
compute_cosine_similarity <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)

  # Row-wise L2 norm
  norm_A <- sqrt(rowSums(A^2))
  norm_B <- sqrt(rowSums(B^2))

  A_unit <- A / norm_A
  B_unit <- B / norm_B

  # Cosine similarity: A x t(B)
  sim <- A_unit %*% t(B_unit)
  return(sim)
}

cos_sim <- compute_cosine_similarity(HLRCC_extracted_sigs$mean, cosmic_signatures_v2)

# Heatmap: extracted signatures (columns) vs COSMIC v2 references (rows)
cols <- colorRampPalette(c("black", "white", "red"))(74)
pheatmap(t(cos_sim),
         color        = cols,
         breaks       = seq(0, 1, length.out = 75),
         cluster_rows = FALSE,
         cluster_col  = FALSE,
         filename     = "HLRCC_cosine_sim_heatmap.pdf",
         height       = 7,
         width        = 3)


# -----------------------------------------------------------------------------
# 11. Save workspace
# -----------------------------------------------------------------------------
save.image("HLRCC_mut_sig.rda")