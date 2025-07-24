# ──────────────────────────────────────────────────────────────────────────────
# Project configuration: directories, file paths, and analysis parameters
# ──────────────────────────────────────────────────────────────────────────────

# Load core packages
library(tidyverse)  # Data import, wrangling, ggplot2
library(DESeq2)     # Differential expression analysis
library(ggrepel)    # Better label placement
library(ggsci)      # Scientific palettes for ggplot2

# ─── Directories ───────────────────────────────────────────────────────────────
base_dir   <- "/path/to/RNAseq_project"
counts_dir <- file.path(base_dir, "FeatureCounts")
meta_file  <- file.path(base_dir, "DEseq", "DESeq2", "samples_information.csv")
output_dir <- file.path(base_dir, "DEseq", "DESeq2", "DEwork", "Results")
plots_dir  <- file.path(output_dir, "Plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir,    recursive = TRUE, showWarnings = FALSE)

# ─── Analysis thresholds ───────────────────────────────────────────────────────
min_count        <- 10    # Minimum per‐gene count
min_samples      <- 5     # Number of samples meeting count threshold
pvalue_threshold <- 0.05  # Adjusted p‐value cutoff
log2fc_threshold <- 1     # Absolute log2 fold-change cutoff