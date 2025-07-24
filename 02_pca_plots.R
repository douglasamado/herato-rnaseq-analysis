# ──────────────────────────────────────────────────────────────────────────────
# Generates PCA plots from variance‐stabilized counts
# ──────────────────────────────────────────────────────────────────────────────

# Load configuration
source("config.R")

# Load DESeq2 dataset
dds <- readRDS(file.path(output_dir, "dds.rds"))

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA on all samples
pca_all <- plotPCA(vsd, intgroup = "Diet", returnData = TRUE)
pct_all <- round(100 * attr(pca_all, "percentVar"), 1)

p_all <- ggplot(pca_all, aes(PC1, PC2, color = Diet)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", pct_all[1], "%)"),
    y = paste0("PC2 (", pct_all[2], "%)")
  ) +
  scale_color_simpsons() +
  theme_bw(base_size = 16)

ggsave(
  filename = file.path(plots_dir, "PCA_all_samples.png"),
  plot    = p_all,
  width   = 6, height = 4, dpi = 300
)

# PCA on larvae only
meta <- as.data.frame(colData(dds))
larvae_samples <- rownames(meta)[meta$Stage == "Late.Larva" | meta$Stage == "Early.Larva"]
pca_larvae <- plotPCA(vsd[, larvae_samples], intgroup = "Diet", returnData = TRUE)
pct_larvae <- round(100 * attr(pca_larvae, "percentVar"), 1)

p_larvae <- ggplot(pca_larvae, aes(PC1, PC2, color = Diet)) +
  geom_point(size = 4, alpha = 0.8) +
  labs(
    x = paste0("PC1 (", pct_larvae[1], "%)"),
    y = paste0("PC2 (", pct_larvae[2], "%)")
  ) +
  scale_color_simpsons() +
  theme_bw(base_size = 16)

ggsave(
  filename = file.path(plots_dir, "PCA_larvae_only.png"),
  plot    = p_larvae,
  width   = 8, height = 5, dpi = 300
)