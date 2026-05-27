# =============================================================================
#  TCGA MBD4-Mutant SBS Mutational Signature Analysis (COAD + UCEC)
# =============================================================================
#
#  Purpose
#  -------
#  Identify the SBS mutational signatures operating in TCGA tumours harbouring
#  somatic MBD4 mutations, as a genetic phenocopy of fumarate-driven MBD4
#  succination / loss-of-function. MBD4-mutant cases are stratified by tumour
#  type (COAD, UCEC) and MLH1 co-mutant cases are excluded to avoid confounding
#  by mismatch-repair deficiency from a distinct cause.
#
#  Workflow
#  --------
#   1. Aggregate ensemble-masked MAFs across COAD and UCEC.
#   2. Retain SNVs in MBD4-mutant tumours; exclude MLH1 co-mutants.
#   3. Build trinucleotide-context catalogues per tumour type.
#   4. Extract de novo SBS signatures (sigfit MCMC) and choose the optimum.
#   5. Match extracted signatures to COSMIC v2 and refit for exposures.
#
#  Inputs
#  ------
#  Directory `main_dir` containing one sub-directory per TCGA project
#  (e.g. TCGA-COAD, TCGA-UCEC), each holding `*ensemble_masked.maf.gz` files.
#
#  Outputs (PDFs in the working directory)
#  ---------------------------------------
#    MBD4_MT_<COHORT>__Signature_Estimate.pdf         GoF sweep
#    MBD4_<COHORT>_optimal_signature_spectrum.pdf     extracted spectra
#    <COHORT>_refit_signature_exposure.pdf            refit exposures
#    MBD4_Mutant_Tumors_cosine_matrix.pdf             cosine vs COSMIC v2
#    MBD4_mutation_Signature_analysis.rda             full workspace image
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
# 1. Aggregate TCGA MAFs by project (COAD, UCEC)
# -----------------------------------------------------------------------------
# Each project directory contains one or more ensemble-masked MAF files.
# We row-bind all MAFs within a project, then retain only SNVs.
main_dir <- "/blue/aikseng.ooi/bstansfield/TCGA_Analysis/MAF/MBD4"
subdirs  <- list.dirs(main_dir, recursive = FALSE)

tcga_data <- list()
for (i in seq_along(subdirs)) {

  maf_files <- list.files(subdirs[i],
                          pattern    = "*ensemble_masked.maf.gz",
                          full.names = TRUE)

  if (length(maf_files) == 0) {
    message("No MAF files in: ", subdirs[i])
    next
  }

  # Read the first MAF in the project (use tryCatch to skip malformed files)
  subdir_maf <- tryCatch({
    read.delim2(maf_files[1],
                sep          = "\t",
                comment.char = "#",
                colClasses   = "character")
  }, error = function(e) {
    message("Error reading: ", maf_files[1], " - skipping this file")
    return(NULL)
  })

  if (is.null(subdir_maf)) next

  # Append the remaining MAFs (if any) to the first
  if (length(maf_files) > 1) {
    for (y in 2:length(maf_files)) {
      df <- tryCatch({
        read.delim2(maf_files[y],
                    sep          = "\t",
                    comment.char = "#",
                    colClasses   = "character")
      }, error = function(e) {
        message("Error reading: ", maf_files[y], " - skipping this file")
        return(NULL)
      })

      if (!is.null(df)) {
        subdir_maf <- rbind(subdir_maf, df)
      }
      print(y)
    }
  }

  # Retain SNVs only and store keyed by project name
  data_filtered <- subdir_maf[which(subdir_maf$VARIANT_CLASS == "SNV"), ]
  tcga_data[[basename(subdirs[i])]] <- data_filtered
  print(i)
}

# -----------------------------------------------------------------------------
# 2. Retain MBD4-mutant tumours; exclude MLH1 co-mutants
# -----------------------------------------------------------------------------
# MLH1 inactivation produces an overlapping mismatch-repair-deficiency
# phenotype that would confound MBD4-specific signature inference, so any
# tumour carrying both is dropped.
processed_tcga_data <- list()

for (i in names(tcga_data)) {
  df <- tcga_data[[i]]

  # MBD4-mutant tumour barcodes
  mbd4_mutant_cases   <- unique(df$Tumor_Sample_Barcode[df$Hugo_Symbol == "MBD4"])
  mbd4_filtered_data  <- df[df$Tumor_Sample_Barcode %in% mbd4_mutant_cases, ]

  # Exclude any of those that also carry MLH1 mutations
  mlh1_mutant_cases    <- unique(mbd4_filtered_data$Tumor_Sample_Barcode[mbd4_filtered_data$Hugo_Symbol == "MLH1"])
  double_filtered_data <- mbd4_filtered_data[!mbd4_filtered_data$Tumor_Sample_Barcode %in% mlh1_mutant_cases, ]

  # GDC `CONTEXT` is 11 nt centred on the variant; chars 5-7 give the trinuc
  double_filtered_data$trinuc_context <- substr(double_filtered_data$CONTEXT, 5, 7)
  processed_tcga_data[[i]] <- double_filtered_data
  print(i)
}

# -----------------------------------------------------------------------------
# 3. Build sigfit input matrices (one per cohort)
# -----------------------------------------------------------------------------
tri_nuc_matrices <- list()
for (i in names(processed_tcga_data)) {
  if (nrow(processed_tcga_data[[i]]) > 1) {
    df <- processed_tcga_data[[i]]
    matrix <- matrix(0, ncol = 4, nrow = nrow(df))
    colnames(matrix) <- c("Sample", "Ref", "Alt", "Trinuc")
    matrix[, "Sample"] <- df$Tumor_Sample_UUID
    matrix[, "Ref"]    <- df$Reference_Allele
    matrix[, "Alt"]    <- df$Tumor_Seq_Allele2
    matrix[, "Trinuc"] <- df$trinuc_context
    tri_nuc_matrices[[i]] <- matrix
  }
  print(i)
}

# =============================================================================
#  Mutational signature analysis
# =============================================================================
library(sigfit)
data("cosmic_signatures_v2")

# -----------------------------------------------------------------------------
# 4. Build 96-channel trinucleotide catalogues
# -----------------------------------------------------------------------------
tri_nuc_counts <- list()
tri_nuc_counts <- lapply(tri_nuc_matrices, build_catalogues)

# -----------------------------------------------------------------------------
# 5. Sweep nsignatures = 2:8 per cohort to choose optimum (GoF)
# -----------------------------------------------------------------------------
MCMC_signatures <- list()
for (i in names(tri_nuc_counts)) {
  MCMC_signatures[[i]] <- extract_signatures(counts      = tri_nuc_counts[[i]],
                                             nsignatures = 2:8,
                                             iter        = 5000,
                                             warmup      = 1000,
                                             seed        = 1756)
  print(i)
}

# Save GoF plots per cohort
for (i in names(MCMC_signatures)) {
  name <- paste0("MBD4_MT_", i)
  pdf(file = paste(name, "_Signature_Estimate.pdf", sep = "_"), height = 5, width = 5)
  plot_gof(MCMC_signatures[[i]])
  dev.off()
  print(i)
}

# -----------------------------------------------------------------------------
# 6. Refit at the cohort-specific optimum with longer MCMC chain
# -----------------------------------------------------------------------------
# Optima chosen from the GoF plots above.
sigs_extracted <- c("COAD" = 4, "UCEC" = 3)

optimal_signatures   <- list()
signature_exposures  <- list()
extracted_signatures <- list()

for (i in names(MCMC_signatures)) {
  nsig    <- as.integer(sigs_extracted[i])
  counts  <- tri_nuc_counts[[i]]

  sig_fit  <- extract_signatures(counts      = counts,
                                 nsignatures = nsig,
                                 iter        = 50000,
                                 warmup      = 1000,
                                 seed        = 1756)
  exposure <- retrieve_pars(sig_fit, par = "exposures",  hpd_prob = 0.95)
  ext_sigs <- retrieve_pars(sig_fit, par = "signatures")

  optimal_signatures[[i]]   <- sig_fit
  signature_exposures[[i]]  <- exposure
  extracted_signatures[[i]] <- ext_sigs
  print(i)
}

# -----------------------------------------------------------------------------
# 7. Cosine similarity of extracted signatures vs COSMIC v2
# -----------------------------------------------------------------------------
cosine_values <- list()
for (i in names(extracted_signatures)) {
  df <- extracted_signatures[[i]]$mean

  df1 <- as.matrix(df)
  df2 <- as.matrix(cosmic_signatures_v2)

  norm_df1 <- sqrt(rowSums(df1^2))
  norm_df2 <- sqrt(rowSums(df2^2))

  df1_unit <- (df1 / norm_df1)
  df2_unit <- (df2 / norm_df2)

  cos_sim <- df1_unit %*% t(df2_unit)
  cosine_values[[i]] <- cos_sim
  print(i)
}

# Prefix row names with the cohort label so the combined heatmap is readable
for (i in names(cosine_values)) {
  name <- i
  data <- cosine_values[[i]]
  rownames(data) <- paste0(i, "_", rownames(data))
  cosine_values[[i]] <- data
}


# -----------------------------------------------------------------------------
# 8. Combined cosine-similarity heatmap (COAD + UCEC)
# -----------------------------------------------------------------------------
library(pheatmap)

combined_cosine_values <- do.call(rbind, cosine_values)
pheatmap(combined_cosine_values,
         cluster_rows       = FALSE,
         cluster_cols       = FALSE,
         clustering_method  = "ward.D2",
         filename           = "MBD4_Mutant_Tumors_cosine_matrix.pdf",
         height             = 7,
         width              = 15)

save(combined_cosine_values, file = "sep_cosine_values.rda")

# -----------------------------------------------------------------------------
# 9. Match extracted signatures to COSMIC and refit for accurate exposures
# -----------------------------------------------------------------------------
matched_signatures <- lapply(extracted_signatures, function(x) {
  match_signatures(x, cosmic_signatures_v2, stat = "cosine")
})

# Hand-curated COSMIC matches (chosen from the cosine matrix above)
matched_signature_data <- list()
matched_signature_data[["COAD"]] <- cosmic_signatures_v2[c("Signature 14", "Signature 10", "Signature 9", "Signature 1"), ]
matched_signature_data[["UCEC"]] <- cosmic_signatures_v2[c("Signature 10", "Signature 14", "Signature 1"), ]

# Refit the chosen COSMIC signatures back to the original cohort catalogues
matched_signature_refit           <- list()
matched_signature_refit_exposures <- list()

for (i in names(tri_nuc_counts)) {
  data <- fit_signatures(counts     = tri_nuc_counts[[i]],
                         signatures = matched_signature_data[[i]],
                         iter       = 50000,
                         warmup     = 1000,
                         seed       = 1756)
  matched_signature_refit_exposures[[i]] <- retrieve_pars(data, par = "exposures", hpd_prob = 0.95)
  matched_signature_refit[[i]]           <- data
  print(i)
}

# -----------------------------------------------------------------------------
# 10. Plot extracted spectra and refit exposures per cohort
# -----------------------------------------------------------------------------
for (i in names(matched_signature_refit)) {
  data <- extracted_signatures[[i]]$mean
  pdf(file = paste("MBD4_", i, "optimal_signature_spectrum.pdf", sep = "_"),
      height = 5, width = 20)
  plot_spectrum(data)
  dev.off()
  print(i)
}

for (i in names(matched_signature_refit)) {
  data <- matched_signature_refit[[i]]
  pdf(file = paste(i, "refit_signature_exposure.pdf", sep = "_"),
      height = 20, width = 5)
  plot_exposures(data)
  dev.off()
  print(i)
}

save.image("MBD4_mutation_Signature_analysis.rda")