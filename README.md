# Metabolite-Driven Succination of MBD4 Links FH Loss to CpG Mutational Bias in HLRCC

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language: R](https://img.shields.io/badge/Language-R-276DC3.svg)](https://www.r-project.org/)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey.svg)](#citation)

Analysis code and figure-generation scripts accompanying the manuscript:

> **Metabolite-Driven Succination of MBD4 Links FH Loss to CpG Mutational Bias in HLRCC**
> *<!-- TODO: author list -->*
> *<!-- TODO: journal / preprint server, year, volume, pages, DOI -->*

---

## Abstract

<!--
TODO: Paste the manuscript abstract here.
-->

> *Abstract to be inserted on acceptance.*

---

## Overview

This repository contains the R scripts used to perform the computational analyses and generate the figures reported in the study. The work integrates:

- **Whole-exome mutational signature analysis** of HLRCC (Hereditary Leiomyomatosis and Renal Cell Cancer) tumours.
- **TCGA pan-cancer analysis** of *MBD4*-mutant tumours (COAD, UCEC) as a genetic phenocopy of fumarate-driven MBD4 loss-of-function.
- **TCGA UCEC *TET*-mutant analysis** as a comparator for 5-methylcytosine deamination–driven mutagenesis.
- **In vitro T:G mismatch repair assays** in HLRCC cell lines (NCCFH1, UOK262) using a methylated / unmethylated molecular beacon substrate.

Together, these analyses test the hypothesis that fumarate-mediated succination of MBD4 impairs T:G mismatch repair at CpG sites and contributes to the CpG-biased mutational landscape of HLRCC.

---

## Repository structure

```
.
├── HLRCC_Exome_Analysis.r          # SBS signature extraction from HLRCC exome MAFs
├── MBD4_Mut_Sig_Analysis.r         # TCGA COAD + UCEC MBD4-mutant signature analysis
├── tet_mutant_UCEC_Cases.r         # TCGA UCEC TET1/2/3-mutant signature analysis
├── TG_Mismatch Repair_Molecular_Beacon.r  # In vitro T:G mismatch repair beacon assay
├── LICENSE                         # MIT license
└── README.md
```

### Script descriptions

| Script | Purpose |
|---|---|
| [HLRCC_Exome_Analysis.r](HLRCC_Exome_Analysis.r) | Reads HLRCC tumour MAFs, builds trinucleotide context catalogues, performs *de novo* SBS extraction with `sigfit` (MCMC), matches extracted signatures to COSMIC v2, refits to estimate exposures, and plots a cosine-similarity heatmap. |
| [MBD4_Mut_Sig_Analysis.r](MBD4_Mut_Sig_Analysis.r) | Aggregates TCGA MAFs across COAD and UCEC, retains *MBD4*-mutant tumours while excluding *MLH1* co-mutants, performs SBS extraction (4 sigs COAD, 3 sigs UCEC), refits COSMIC matches, and generates normalised per-substitution trinucleotide heatmaps. |
| [tet_mutant_UCEC_Cases.r](tet_mutant_UCEC_Cases.r) | Performs analogous SBS extraction on TCGA UCEC *TET1/TET2/TET3*-mutant cases as a 5mC-deamination comparator. |
| [TG_Mismatch Repair_Molecular_Beacon.r](TG_Mismatch%20Repair_Molecular_Beacon.r) | Plots real-time fluorescence traces from molecular-beacon T:G mismatch-repair assays in NCCFH1 and UOK262 cells (± doxycycline induction; methylated vs unmethylated probe). |

---

## Figure mapping

<!-- TODO: fill in figure / panel numbers once final on proof. -->

| Manuscript figure | Script | Output file(s) |
|---|---|---|
| Fig. <!-- TODO -->, HLRCC SBS spectrum & exposures | `HLRCC_Exome_Analysis.r` | `HLRCC_enriched_signatures.pdf`, `HLRCC_signature_exposure.pdf`, `HLRCC_cosine_sim_heatmap.pdf` |
| Fig. <!-- TODO -->, MBD4-mutant TCGA signatures (COAD/UCEC) | `MBD4_Mut_Sig_Analysis.r` | `MBD4_*_optimal_signature_spectrum.pdf`, `*_refit_signature_exposure.pdf`, `MBD4_Mutant_Tumors_cosine_matrix.pdf` |
| Fig. <!-- TODO -->, MBD4-mutant trinucleotide heatmaps | `MBD4_Mut_Sig_Analysis.r` | `UCEC_MBD4_MT_Heatmap.pdf`, `UCEC_COAD_MT_Heatmap.pdf` |
| Fig. <!-- TODO -->, TET-mutant UCEC comparator | `tet_mutant_UCEC_Cases.r` | `tet_enriched_signatures.pdf`, `tet_exposures.pdf`, `UCEC_Tet_Heatmap.pdf` |
| Fig. <!-- TODO -->, T:G mismatch repair beacon assay | `TG_Mismatch Repair_Molecular_Beacon.r` | `NCCFH1_UOK262_Ctrl_Dox_Meth_Unmeth.pdf` |

---

## Requirements

- **R** ≥ 4.2
- **Stan** (installed automatically as a dependency of `sigfit`; requires a working C++ toolchain)
- R packages:
  - [`sigfit`](https://github.com/kgori/sigfit) (mutational signature extraction & fitting)
  - `pheatmap`
  - `ggplot2`
  - `ggpubr`

Install dependencies in R:

```r
install.packages(c("pheatmap", "ggplot2", "ggpubr", "devtools"))
devtools::install_github("kgori/sigfit", build_vignettes = TRUE)
```

> The mutational-signature scripts use `sigfit::extract_signatures()` with up to 50,000 MCMC iterations and will benefit from a multi-core machine. A typical TCGA-scale run takes 30 min–several hours depending on cohort size and `nsignatures`.

---

## Data availability

| Dataset | Source | Access |
|---|---|---|
| TCGA MAFs (COAD, UCEC) | NCI Genomic Data Commons | Public — [https://portal.gdc.cancer.gov](https://portal.gdc.cancer.gov) (ensemble-masked MAFs) |
| COSMIC SBS v2 reference signatures | bundled with `sigfit` | Public — `data("cosmic_signatures_v2")` |
| HLRCC tumour whole-exome MAFs | This study | Available from the corresponding author on reasonable request |
| Molecular-beacon fluorescence raw data (`NCCFH1_*`, `UOK262_*` CSVs) | This study | Available from the corresponding author on reasonable request |

Paths inside the scripts (e.g. `/blue/aikseng.ooi/...`, `setwd(...)`) reflect the original analysis environment and will need to be adjusted to match your local data location.

---

## Reproducing the analyses

1. Clone the repository.
   ```bash
   git clone https://github.com/<!-- TODO: org/user -->/Metabolite-Driven-Succination-of-MBD4-Links-FH-Loss-to-CpG-Mutational-Bias-in-HLRCC.git
   cd Metabolite-Driven-Succination-of-MBD4-Links-FH-Loss-to-CpG-Mutational-Bias-in-HLRCC
   ```
2. Obtain the input data (see *Data availability*) and update file paths near the top of each script.
3. Run each script independently from R / RStudio:
   ```r
   source("HLRCC_Exome_Analysis.r")
   source("MBD4_Mut_Sig_Analysis.r")
   source("tet_mutant_UCEC_Cases.r")
   source("TG_Mismatch Repair_Molecular_Beacon.r")
   ```

A fixed random seed (`seed = 1756`) is used in all `sigfit` calls so that signature extraction is deterministic.

---

## Citation

If you use code from this repository, please cite:

```bibtex
@article{<!-- TODO: key -->,
  title   = {Metabolite-Driven Succination of MBD4 Links FH Loss to CpG Mutational Bias in HLRCC},
  author  = {<!-- TODO: authors -->},
  journal = {<!-- TODO: journal -->},
  year    = {<!-- TODO: year -->},
  doi     = {<!-- TODO: DOI -->}
}
```

---

## Contact

Questions about the analyses or requests for data should be directed to the corresponding author (<!-- TODO: name + email -->).

For issues specific to the code in this repository, please open a [GitHub issue](../../issues).

---

## License

This code is released under the [MIT License](LICENSE).
