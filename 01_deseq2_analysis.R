# ──────────────────────────────────────────────────────────────────────────────
# Performs differential expression analysis and saves results
# ──────────────────────────────────────────────────────────────────────────────

# Load configuration
source("config.R")

# 1. READ RAW COUNTS
counts_raw <- read_tsv(
  file.path(counts_dir,
            "Herato.RNAseq.diets.ecp.featurecounts.anno-v3.txt"),
  skip = 1
)

# Strip annotation columns, set gene IDs as row names
count_matrix <- counts_raw %>%
  select(-c(Chr, Start, End, Strand, Length)) %>%
  column_to_rownames("Geneid") %>%
  as.matrix()

# 2. READ AND FORMAT METADATA
sample_info <- read_csv(meta_file) %>%
  mutate(
    Diet  = factor(Diet, levels = c("Biflora", "Auriculata", "Biflora.Painted")),
    Stage = factor(Stage, levels = c("Early.Larva", "Late.Larva", "Pupa")),
    Group = factor(paste(Diet, Stage, sep = "."))
  )

count_matrix <- count_matrix[, sample_info$ID]

# 3. BUILD DESEQ2 DATASET & FILTER LOW‐COUNT GENES
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = sample_info,
  design    = ~ Group
)

keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds  <- dds[keep, ]

# 4. RUN DIFFERENTIAL EXPRESSION
dds <- DESeq(dds)

# 5. SAVE DESEQ2 OBJECT & RAW RESULTS
saveRDS(dds, file = file.path(output_dir, "dds.rds"))

# Define contrasts for late-larva stage
contrasts_larva <- list(
  "Auriculata_vs_Biflora"      = c("Group", "Auriculata.Late.Larva", "Biflora.Late.Larva"),
  "BifloraPainted_vs_Biflora"  = c("Group", "Biflora.Painted.Late.Larva", "Biflora.Late.Larva"),
  "BifloraPainted_vs_Auriculata" = c("Group", "Biflora.Painted.Late.Larva", "Auriculata.Late.Larva")
)

# Save raw results tables
for (label in names(contrasts_larva)) {
  res <- results(dds, contrast = contrasts_larva[[label]])
  write_csv(
    as.data.frame(res) %>% rownames_to_column("Gene"),
    file.path(output_dir, paste0("DEGs_", label, "_LateLarva_raw.csv"))
  )
}