# =============================================================================
# LRP1 Expression Analysis — Litviňuková et al. 2020 Adult Human Heart Atlas
# Dataset: Litviňuková M et al. Cells of the adult human heart.
#          Nature 588: 466–472, 2020. doi: 10.1038/s41586-020-2797-4
#
# Produces:
#   - LRP1_UMAP_adult_heart.pdf                  : FeaturePlot coloured by LRP1 expression
#   - LRP1_pseudobulk_boxplot_adult_heart.pdf    : Pseudobulk DESeq2-normalised LRP1 counts
#
# Usage:
#   1. Set the two paths in the "User-defined paths" section below.
#   2. Run the script end-to-end in a fresh R session.
#
# Requirements: Seurat, Matrix, ggplot2, dplyr, DESeq2, tidyverse, rhdf5
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(tidyverse)
library(rhdf5)

# ── User-defined paths ────────────────────────────────────────────────────────
seurat_rds <- "path/to/litvinukova_global_seurat.rds"  # pre-converted Seurat object

# ── Gene of interest ──────────────────────────────────────────────────────────
GENE       <- "ENSG00000123384"
GENE_LABEL <- "LRP1"

# =============================================================================
# 1. Load data and rename metadata columns to match analysis conventions
# =============================================================================

obj <- readRDS(seurat_rds)

obj$author_cell_type <- obj$cell_type
obj$donor_id         <- obj$donor

table(obj$author_cell_type)
table(obj$donor_id)

# =============================================================================
# 2. Colour palette
# =============================================================================

cell_colours <- c(
  "Endothelial cell"            = "#E87D72",
  "Lymphatic Endothelial cell"  = "#F0C05A",
  "Fibroblast"                  = "#7FC97F",
  "Lymphoid"                    = "#AEC7E8",
  "Mast cell"                   = "#FCCDE5",
  "Neural cell"                 = "#CCEBC5",
  "Myeloid"                     = "#C39BD3",
  "Mural cell"                  = "#3A85A8",
  "Mesothelial cell"            = "#BEBEBE",
  "Atrial Cardiomyocyte"        = "#FFED6F",
  "Ventricular Cardiomyocyte"   = "#FDB462",
  "Adipocyte"                   = "#80B1D3"
)

# Fallback colours for any cell types not listed above
present_types <- unique(as.character(obj$author_cell_type))
missing_types <- setdiff(present_types, names(cell_colours))
if (length(missing_types) > 0) {
  extra_cols <- c("#F0C05A", "#56B4D9", "#BEBEBE", "#D95F02", "#7570B3")
  names(extra_cols) <- missing_types[seq_len(min(length(missing_types), length(extra_cols)))]
  cell_colours <- c(cell_colours, extra_cols)
}

# =============================================================================
# 3. Normalise
# =============================================================================

obj <- NormalizeData(obj)

# =============================================================================
# 4. UMAP — LRP1 expression
# =============================================================================

p_umap <- FeaturePlot(
  obj,
  features  = GENE,
  reduction = "umap",
  pt.size   = 0.04,
  cols      = c("lightgrey", "blue"),
  raster    = FALSE,
  alpha     = 0.7
) +
  labs(
    title = paste0(GENE_LABEL, " Expression"),
    x     = "UMAP 1",
    y     = "UMAP 2"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("LRP1_UMAP_adult_heart.png", p_umap, width = 7, height = 6, dpi = 300)
ggsave("LRP1_UMAP_adult_heart.pdf", p_umap, width = 7, height = 6)

# =============================================================================
# 5. Pseudobulk DESeq2 analysis
# =============================================================================

# ── 5.1 Aggregate counts per donor x cell type ───────────────────────────────
counts_mat <- GetAssayData(obj, layer = "counts")
meta       <- obj@meta.data
meta$cell  <- colnames(obj)

# Remove cells with missing donor or cell type annotation
meta       <- meta[!is.na(meta$donor_id) & !is.na(meta$author_cell_type), ]
counts_mat <- counts_mat[, rownames(meta)]

meta$group <- paste(meta$donor_id, meta$author_cell_type, sep = "__")
groups     <- unique(meta$group)

pb_counts <- sapply(groups, function(g) {
  cells <- meta$cell[meta$group == g]
  Matrix::rowSums(counts_mat[, cells, drop = FALSE])
})
pb_counts <- as.matrix(pb_counts)

# ── 5.2 Build sample metadata ─────────────────────────────────────────────────
pb_meta <- data.frame(
  sample    = groups,
  donor     = sub("__.*", "", groups),
  cell_type = sub(".*__", "", groups),
  stringsAsFactors = FALSE
)
rownames(pb_meta) <- groups

# ── 5.3 Filter groups with fewer than 10 cells ───────────────────────────────
cell_counts <- table(meta$donor_id, meta$author_cell_type)
keep_groups <- sapply(groups, function(g) {
  d  <- sub("__.*", "", g)
  ct <- sub(".*__", "", g)
  cell_counts[d, ct] >= 10
})

pb_counts <- pb_counts[, keep_groups, drop = FALSE]
pb_meta   <- pb_meta[keep_groups, , drop = FALSE]

# ── 5.4 Run DESeq2 ────────────────────────────────────────────────────────────
pb_counts        <- round(pb_counts)
mode(pb_counts)  <- "integer"

pb_meta$cell_type <- factor(pb_meta$cell_type)
if ("Mural cell" %in% levels(pb_meta$cell_type)) {
  pb_meta$cell_type <- relevel(pb_meta$cell_type, ref = "Mural cell")
}

dds <- DESeqDataSetFromMatrix(
  countData = pb_counts,
  colData   = pb_meta,
  design    = ~ donor + cell_type
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

resultsNames(dds)

# =============================================================================
# 6. Pseudobulk boxplot
# =============================================================================

norm_counts <- counts(dds, normalized = TRUE)

plot_df <- data.frame(
  sample     = colnames(norm_counts),
  norm_count = as.numeric(norm_counts[GENE, ]),
  cell_type  = pb_meta[colnames(norm_counts), "cell_type"]
)

p_box <- ggplot(plot_df, aes(x = cell_type, y = norm_count, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
  geom_jitter(width = 0.15, size = 2.2, alpha = 0.8) +
  scale_fill_manual(values = cell_colours, drop = FALSE) +
  labs(
    title = paste0(GENE_LABEL, " Pseudobulk Expression"),
    x     = NULL,
    y     = "Normalised counts"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave("LRP1_pseudobulk_boxplot_adult_heart.png", p_box, width = 9, height = 5, dpi = 300)
ggsave("LRP1_pseudobulk_boxplot_adult_heart.pdf", p_box, width = 9, height = 5)


# =============================================================================
# 7. Pairwise cell-type-vs-cell-type comparisons for LRP1
#    results(dds, contrast = c("cell_type", A, B)) gives log2(A / B).
#    Use RAW factor level names (with spaces, e.g. "Mural cell").
# =============================================================================

cell_levels <- levels(pb_meta$cell_type)
pairs       <- combn(cell_levels, 2, simplify = FALSE)

pairwise_results <- lapply(pairs, function(p) {
  res <- results(dds, contrast = c("cell_type", p[1], p[2]))
  row <- res[GENE, ]
  data.frame(
    cell_type_A = p[1],
    cell_type_B = p[2],
    log2FC      = row$log2FoldChange,   # log2(A / B); >0 = higher in A
    lfcSE       = row$lfcSE,
    pvalue      = row$pvalue,
    padj_gene   = row$padj              # BH across all genes (DESeq2 default)
  )
}) |> bind_rows()

# For a single pre-specified gene across many pairs, re-adjust across the
# SET of pairwise tests rather than genome-wide:
pairwise_results <- pairwise_results %>%
  mutate(padj_pairwise = p.adjust(pvalue, method = "BH")) %>%
  arrange(padj_pairwise)

print(pairwise_results)
write_csv(pairwise_results, "LRP1_pairwise_results_adult_heart.csv")