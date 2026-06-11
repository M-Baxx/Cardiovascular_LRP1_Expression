# =============================================================================
# LRP1 Expression Analysis — Zhao et al. 2025 Human Arterial scRNA-seq Atlas
# Dataset: Zhao Q et al. A cell and transcriptome atlas of human arterial
#          vasculature. Cell Genomics 5: 101034, 2025.
#          https://datasets.cellxgene.cziscience.com/c8f274de-bd5a-4604-b3ca-caf24fa8107a.h5ad
#
# Produces:
#   - LRP1_UMAP.pdf               : FeaturePlot coloured by LRP1 expression
#   - LRP1_UMAP_labelled.pdf      : DimPlot coloured by cell type with labels
#   - LRP1_boxplot.pdf            : Pseudobulk DESeq2-normalised LRP1 counts
#   - LRP1_dotplot_region.pdf     : Dot plot of LRP1 across regions x cell types
#   - LRP1_pseudobulk_results.csv : DESeq2 results table for LRP1
#
# Usage:
#   1. Set the two paths in the "User-defined paths" section below.
#   2. Run the script end-to-end in a fresh R session.
#
# Requirements: Seurat, reticulate, Matrix, ggplot2, dplyr, DESeq2,
#               tidyverse, rhdf5, ggrepel
#               Python virtualenv "r-anndata" with anndata installed
# =============================================================================

# ── Libraries ─────────────────────────────────────────────────────────────────
library(Seurat)
library(reticulate)
library(Matrix)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(tidyverse)
library(rhdf5)
library(ggrepel)

# ── User-defined paths ────────────────────────────────────────────────────────
h5ad_file  <- "path/to/zhao2025_arterial.h5ad"   # downloaded h5ad file
seurat_rds <- "path/to/zhao2025_seurat.rds"       # pre-converted Seurat object

# ── Gene of interest ──────────────────────────────────────────────────────────
GENE       <- "ENSG00000123384"
GENE_LABEL <- "LRP1"

# ── Cell type colour palette ──────────────────────────────────────────────────
cell_colours <- c(
  "Endothelial1"     = "#E87D72",
  "Endothelial2"     = "#F0C05A",
  "Fibroblast"       = "#7FC97F",
  "Lymphocyte"       = "#AEC7E8",
  "Macrophage"       = "#C39BD3",
  "Microvasculature" = "#56B4D9",
  "Smooth Muscle"    = "#3A85A8",
  "Uncharacterized"  = "#BEBEBE"
)

# =============================================================================
# 1. Download data (skip if already downloaded)
# =============================================================================

options(timeout = 600)  # 10 minutes

download.file(
  url      = "https://datasets.cellxgene.cziscience.com/c8f274de-bd5a-4604-b3ca-caf24fa8107a.h5ad",
  destfile = h5ad_file,
  mode     = "wb"
)

# =============================================================================
# 2. Inspect h5ad structure
# =============================================================================

h5closeAll()
h5ls(h5ad_file) |> print(max = 200)

# =============================================================================
# 3. Load Seurat object and transfer annotations from h5ad
# =============================================================================

obj <- readRDS(seurat_rds)

# Extract cell IDs and cell type labels
cell_ids      <- h5read(h5ad_file, "/obs/_index")
ct_categories <- h5read(h5ad_file, "/obs/author_cell_type/categories")
ct_codes      <- h5read(h5ad_file, "/obs/author_cell_type/codes")
author_cell_type <- ct_categories[ct_codes + 1]  # +1: codes are 0-indexed

# Transfer cell type labels
obj$author_cell_type <- author_cell_type[match(colnames(obj), cell_ids)]

# Extract and transfer UMAP coordinates
umap_coords <- h5read(h5ad_file, "/obsm/X_umap")
umap_coords <- t(umap_coords)  # transpose to cells x dims
rownames(umap_coords) <- cell_ids
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

umap_matched <- umap_coords[match(colnames(obj), cell_ids), ]
obj[["umap"]] <- CreateDimReducObject(
  embeddings = umap_matched,
  key        = "UMAP_",
  assay      = DefaultAssay(obj)
)

# Transfer donor IDs
donor_categories <- h5read(h5ad_file, "/obs/donor_id/categories")
donor_codes      <- h5read(h5ad_file, "/obs/donor_id/codes")
obj$donor_id     <- donor_categories[donor_codes + 1][match(colnames(obj), cell_ids)]

# Transfer tissue/region labels
tissue_categories <- h5read(h5ad_file, "/obs/tissue/categories")
tissue_codes      <- h5read(h5ad_file, "/obs/tissue/codes")
obj$tissue        <- tissue_categories[tissue_codes + 1][match(colnames(obj), cell_ids)]

# Verify transfers
cat("Shared cells:", length(intersect(colnames(obj), cell_ids)), "\n")
table(obj$author_cell_type)
table(obj$donor_id, obj$author_cell_type)
table(obj$tissue)

Idents(obj) <- "author_cell_type"

# =============================================================================
# 4. UMAP plots
# =============================================================================

# 4a. FeaturePlot — LRP1 expression
p_feature <- FeaturePlot(
  obj,
  features  = GENE,
  reduction = "umap",
  pt.size   = 0.04,
  cols      = c("lightgrey", "blue"),
  raster    = FALSE,
  alpha     = 0.7
) +
  labs(
    title = "LRP1 Expression in Adult Arterial Tissue Cells",
    x     = "UMAP 1",
    y     = "UMAP 2"
  ) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("LRP1_UMAP.png", p_feature, width = 7, height = 6, dpi = 300)
ggsave("LRP1_UMAP.pdf", p_feature, width = 7, height = 6)

# 4b. DimPlot — cell type labels
umap_df <- data.frame(
  UMAP_1    = Embeddings(obj, "umap")[, 1],
  UMAP_2    = Embeddings(obj, "umap")[, 2],
  cell_type = obj$author_cell_type
)

centroids <- umap_df %>%
  group_by(cell_type) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

p_labelled <- DimPlot(
  obj,
  reduction  = "umap",
  group.by   = "author_cell_type",
  pt.size    = 0.04,
  cols       = cell_colours,
  label      = FALSE,
  raster     = FALSE
) +
  geom_label_repel(
    data        = centroids,
    aes(x = UMAP_1, y = UMAP_2, label = cell_type),
    size        = 7,
    fontface    = "bold",
    box.padding = 0.5,
    fill        = "white",
    alpha       = 0.7,
    inherit.aes = FALSE
  ) +
  labs(
    title = "Cell Types in Adult Arterial Tissue",
    x     = "UMAP 1",
    y     = "UMAP 2"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    legend.position = "none"
  )

ggsave("LRP1_UMAP_labelled.png", p_labelled, width = 7, height = 6, dpi = 300)
ggsave("LRP1_UMAP_labelled.pdf", p_labelled, width = 7, height = 6)

# =============================================================================
# 5. Pseudobulk DESeq2 analysis
# =============================================================================

# ── 5.1 Aggregate counts per donor x cell type ───────────────────────────────
counts_mat <- GetAssayData(obj, layer = "counts")
meta       <- obj@meta.data
meta$cell  <- colnames(obj)
meta$group <- paste(meta$donor_id, meta$author_cell_type, sep = "__")

groups    <- unique(meta$group)
pb_counts <- sapply(groups, function(g) {
  cells <- meta$cell[meta$group == g]
  Matrix::rowSums(counts_mat[, cells, drop = FALSE])
})
pb_counts <- as.matrix(pb_counts)

# ── 5.2 Build sample metadata ─────────────────────────────────────────────────
pb_meta <- data.frame(
  sample    = groups,
  donor     = sub("__.*", "", groups),
  cell_type = sub(".*__", "", groups)
)
rownames(pb_meta) <- groups

pb_meta$cell_type <- factor(pb_meta$cell_type)
pb_meta$cell_type <- relevel(pb_meta$cell_type, ref = "Smooth Muscle")

# ── 5.3 Filter groups with fewer than 10 cells ───────────────────────────────
cell_counts <- table(meta$donor_id, meta$author_cell_type)
keep <- sapply(groups, function(g) {
  d  <- sub("__.*", "", g)
  ct <- sub(".*__", "", g)
  cell_counts[d, ct] >= 10
})

pb_counts <- pb_counts[, keep]
pb_meta   <- pb_meta[keep, ]

# ── 5.4 Run DESeq2 ────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = pb_counts,
  colData   = pb_meta,
  design    = ~ donor + cell_type
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

resultsNames(dds)

# ── 5.5 Extract LRP1 results for all cell type contrasts ─────────────────────
contrasts <- c(
  "Endothelial1"     = "cell_type_Endothelial1_vs_Smooth.Muscle",
  "Endothelial2"     = "cell_type_Endothelial2_vs_Smooth.Muscle",
  "Fibroblast"       = "cell_type_Fibroblast_vs_Smooth.Muscle",
  "Lymphocyte"       = "cell_type_Lymphocyte_vs_Smooth.Muscle",
  "Macrophage"       = "cell_type_Macrophage_vs_Smooth.Muscle",
  "Microvasculature" = "cell_type_Microvasculature_vs_Smooth.Muscle",
  "Uncharacterized"  = "cell_type_Uncharacterized_vs_Smooth.Muscle"
)

lrp1_results <- lapply(names(contrasts), function(ct) {
  res      <- results(dds, name = contrasts[ct])
  lrp1_row <- res[GENE, ]
  data.frame(
    cell_type = ct,
    baseMean  = lrp1_row$baseMean,
    log2FC    = lrp1_row$log2FoldChange,
    lfcSE     = lrp1_row$lfcSE,
    pvalue    = lrp1_row$pvalue,
    padj      = lrp1_row$padj
  )
}) |> bind_rows()

# Add Smooth Muscle as reference row (log2FC = 0 by definition)
lrp1_results <- bind_rows(
  data.frame(cell_type = "Smooth Muscle", baseMean = NA,
             log2FC = 0, lfcSE = 0, pvalue = NA, padj = NA),
  lrp1_results
)

print(lrp1_results)
write_csv(lrp1_results, "LRP1_pseudobulk_results.csv")

# =============================================================================
# 6. Pseudobulk boxplot
# =============================================================================

norm_counts <- counts(dds, normalized = TRUE)
lrp1_norm   <- norm_counts[GENE, ]

plot_df <- data.frame(
  sample     = names(lrp1_norm),
  norm_count = as.numeric(lrp1_norm),
  cell_type  = pb_meta[names(lrp1_norm), "cell_type"]
)

p_box <- ggplot(plot_df, aes(x = cell_type, y = norm_count, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 2.5, alpha = 0.8) +
  scale_fill_manual(values = cell_colours) +
  labs(
    title = "LRP1 Pseudobulk Expression Across Adult Arterial Cells",
    x     = NULL,
    y     = "Normalised Counts"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave("LRP1_boxplot.png", p_box, width = 9, height = 5, dpi = 300)
ggsave("LRP1_boxplot.pdf", p_box, width = 9, height = 5)

# =============================================================================
# 7. Dot plot — LRP1 expression across arterial regions and cell types
# =============================================================================

expr_vec <- as.numeric(GetAssayData(obj, layer = "data")[GENE, ])
meta_df  <- data.frame(
  expression = expr_vec,
  cell_type  = obj$author_cell_type,
  tissue     = obj$tissue
)

dot_df <- meta_df %>%
  group_by(tissue, cell_type) %>%
  summarise(
    avg_exp = mean(expression),
    pct_exp = mean(expression > 0) * 100,
    .groups = "drop"
  ) %>%
  mutate(tissue = tools::toTitleCase(as.character(tissue)))

p_dot_region <- ggplot(dot_df, aes(
  x      = tissue,
  y      = cell_type,
  size   = pct_exp,
  colour = avg_exp
)) +
  geom_point() +
  scale_colour_gradient(low = "lightgrey", high = "#b2182b",
                        name = "Avg\nExpression") +
  scale_size(range = c(0.5, 8), name = "% Expressing") +
  labs(
    title = "LRP1 Expression Across Arterial Regions and Cell Types",
    x     = NULL,
    y     = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3)
  )

ggsave("LRP1_dotplot_region.png", p_dot_region, width = 12, height = 6, dpi = 300)
ggsave("LRP1_dotplot_region.pdf", p_dot_region, width = 12, height = 6)

#Pairwise Comparison

# All cell-type levels present in the model
cell_levels <- levels(pb_meta$cell_type)

# Every unique pair (order: A vs B = log2(A / B))
pairs <- combn(cell_levels, 2, simplify = FALSE)

pairwise_results <- lapply(pairs, function(p) {
  A <- p[1]; B <- p[2]
  res <- results(dds, contrast = c("cell_type", A, B))
  row <- res[GENE, ]
  data.frame(
    cell_type_A = A,
    cell_type_B = B,
    log2FC      = row$log2FoldChange,   # log2(A / B)
    lfcSE       = row$lfcSE,
    pvalue      = row$pvalue,
    padj_gene   = row$padj              # BH across genes, per contrast
  )
}) |> bind_rows()

print(pairwise_results)

