# ──────────────────────────────────────────────────────────────────────────────
# Creates volcano plots and saves filtered DEG tables
# ──────────────────────────────────────────────────────────────────────────────

# Load configuration
source("config.R")

# Load DESeq2 dataset
dds <- readRDS(file.path(output_dir, "dds.rds"))

# Volcano plotting helper
plot_volcano <- function(res, label, stage) {
  res_df <- as_tibble(res, rownames = "Gene") %>%
    filter(!is.na(pvalue)) %>%
    mutate(
      Regulation = case_when(
        pvalue < pvalue_threshold & log2FoldChange >= log2fc_threshold  ~ "Up",
        pvalue < pvalue_threshold & log2FoldChange <= -log2fc_threshold ~ "Down",
        TRUE                                                         ~ "NotSig"
      )
    )
  
  # Save DEG table
  write_csv(
    res_df,
    file.path(output_dir,
              paste0("DEGs_", label, "_", stage, "_filtered.csv"))
  )
  
  # Plot
  volcano <- ggplot(res_df, aes(log2FoldChange, -log10(pvalue))) +
    geom_point(aes(fill = Regulation, size = Regulation, alpha = Regulation),
               colour = "black", shape = 21) +
    scale_fill_manual(values = c(Up = "#ffad73", Down = "#26b3ff", NotSig = "grey")) +
    scale_size_manual(values = c(Up = 5, Down = 5, NotSig = 4)) +
    scale_alpha_manual(values = c(Up = 1, Down = 1, NotSig = 0.5)) +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed") +
    theme_minimal(base_size = 18) +
    labs(
      title = paste(label, stage),
      x     = "log2(Fold Change)",
      y     = "-log10(p-value)"
    )
  
  ggsave(
    filename = file.path(plots_dir,
                         paste0("volcano_", label, "_", stage, ".png")),
    plot   = volcano,
    width  = 12, height = 8, dpi = 300, bg = "white"
  )
}

# Define all contrasts
contrasts <- list(
  LateLarva = list(
    Aur_vs_Bif    = c("Group", "Auriculata.Late.Larva", "Biflora.Late.Larva"),
    BiflPaint_vs_Bif = c("Group", "Biflora.Painted.Late.Larva", "Biflora.Late.Larva"),
    BiflPaint_vs_Aur  = c("Group", "Biflora.Painted.Late.Larva", "Auriculata.Late.Larva")
  ),
  Pupa = list(
    Aur_vs_Bif    = c("Group", "Auriculata.Pupa", "Biflora.Pupa"),
    BiflPaint_vs_Bif = c("Group", "Biflora.Painted.Pupa", "Biflora.Pupa"),
    BiflPaint_vs_Aur  = c("Group", "Biflora.Painted.Pupa", "Auriculata.Pupa")
  ),
  BifloraStages = list(
    Late_vs_Early = c("Group", "Biflora.Late.Larva", "Biflora.Early.Larva"),
    Pupa_vs_Early = c("Group", "Biflora.Pupa", "Biflora.Early.Larva"),
    Late_vs_Pupa  = c("Group", "Biflora.Late.Larva", "Biflora.Pupa")
  )
)

# Loop through and plot
for (stage in names(contrasts)) {
  for (label in names(contrasts[[stage]])) {
    res <- results(dds, contrast = contrasts[[stage]][[label]])
    plot_volcano(res, label, stage)
  }
}