#!/usr/bin/env Rscript

# Run from the repository root. Put the 10x folders in data/AD1 and
# data/Control2, or set DATA_DIR to the folder that contains both directories.
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(clustree)
})

input_dir <- Sys.getenv("DATA_DIR", "data")
output_dir <- file.path("results")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

sample_dirs <- c(AD1 = file.path(input_dir, "AD1"),
                 Control2 = file.path(input_dir, "Control2"))
missing_dirs <- sample_dirs[!dir.exists(sample_dirs)]
if (length(missing_dirs)) {
  stop("Missing input folder(s): ", paste(missing_dirs, collapse = ", "),
       ". See data/README.md.")
}

read_sample <- function(sample_name, folder) {
  counts <- Seurat::Read10X(data.dir = folder)
  if (is.list(counts)) {
    if (!"Gene Expression" %in% names(counts)) {
      stop("Expected a Gene Expression matrix in ", folder)
    }
    counts <- counts[["Gene Expression"]]
  }
  object <- Seurat::CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3,
    min.features = 100
  )
  object$sample <- sample_name
  object$condition <- if (sample_name == "AD1") "AD" else "Control"
  object$percent.mito <- Seurat::PercentageFeatureSet(
    object, pattern = "^MT-"
  )
  object$percent.ribo <- Seurat::PercentageFeatureSet(
    object, pattern = "^RP[SL]"
  )
  object$percent.globin <- Seurat::PercentageFeatureSet(
    object, pattern = "^HB[^(P)]"
  )
  object
}

objects <- lapply(names(sample_dirs), function(nm) read_sample(nm, sample_dirs[[nm]]))
names(objects) <- names(sample_dirs)

qc_before <- do.call(rbind, lapply(names(objects), function(nm) {
  data.frame(sample = nm, cells_before_filtering = ncol(objects[[nm]]))
}))

objects <- lapply(objects, function(object) {
  subset(object, subset = nFeature_RNA > 200 &
           nFeature_RNA < 5000 & percent.mito < 10)
})
qc_before$cells_after_filtering <- vapply(objects, ncol, integer(1))
qc_before$cells_removed <- qc_before$cells_before_filtering -
  qc_before$cells_after_filtering
qc_before$percent_removed <- round(
  100 * qc_before$cells_removed / qc_before$cells_before_filtering, 1
)
write.csv(qc_before, file.path(output_dir, "cell_filtering_summary.csv"),
          row.names = FALSE)

combined <- merge(
  x = objects[["AD1"]],
  y = objects[["Control2"]],
  add.cell.ids = c("AD1", "Control2"),
  project = "OlfactoryMucosa_AD"
)

combined <- Seurat::NormalizeData(
  combined, normalization.method = "LogNormalize", scale.factor = 10000
)
combined <- Seurat::FindVariableFeatures(
  combined, selection.method = "vst", nfeatures = 2000
)
combined <- Seurat::ScaleData(combined, features = Seurat::VariableFeatures(combined))
combined <- Seurat::RunPCA(combined, npcs = 40, verbose = FALSE)
combined <- Seurat::RunUMAP(combined, dims = 1:25, verbose = FALSE)
combined <- Seurat::FindNeighbors(combined, dims = 1:25, verbose = FALSE)
combined <- Seurat::FindClusters(combined, resolution = seq(0.1, 0.8, by = 0.1),
                                 verbose = FALSE)

p_pca <- Seurat::DimPlot(combined, reduction = "pca", group.by = "condition") +
  ggplot2::ggtitle("PCA by sample condition")
p_umap <- Seurat::DimPlot(combined, reduction = "umap", group.by = "condition") +
  ggplot2::ggtitle("UMAP by sample condition")
p_clusters <- Seurat::DimPlot(combined, reduction = "umap",
                              group.by = "RNA_snn_res.0.8", label = TRUE) +
  ggplot2::ggtitle("Clusters at resolution 0.8")
p_elbow <- Seurat::ElbowPlot(combined, ndims = 40)
resolution_columns <- grep("^RNA_snn_res\\.", colnames(combined@meta.data), value = TRUE)
p_tree <- clustree::clustree(
  combined@meta.data[, resolution_columns, drop = FALSE],
  prefix = "RNA_snn_res."
)

ggplot2::ggsave(file.path(figure_dir, "pca_by_condition.png"), p_pca,
                width = 7, height = 5, dpi = 300)
ggplot2::ggsave(file.path(figure_dir, "umap_by_condition.png"), p_umap,
                width = 7, height = 5, dpi = 300)
ggplot2::ggsave(file.path(figure_dir, "umap_clusters_res08.png"), p_clusters,
                width = 7, height = 5, dpi = 300)
ggplot2::ggsave(file.path(figure_dir, "pca_elbow.png"), p_elbow,
                width = 7, height = 5, dpi = 300)
ggplot2::ggsave(file.path(figure_dir, "clustree_resolutions.png"), p_tree,
                width = 9, height = 7, dpi = 300)

# The resolution 0.8 cluster markers are descriptive candidates for manual
# annotation. They are not disease biomarkers.
Seurat::Idents(combined) <- "RNA_snn_res.0.8"
cluster_markers <- Seurat::FindAllMarkers(
  combined, min.pct = 0.25, only.pos = TRUE, test.use = "wilcox"
)
cluster_markers <- subset(cluster_markers, !is.na(p_val_adj) & p_val_adj < 0.05)
write.csv(cluster_markers,
          file.path(output_dir, "cluster_markers_exploratory.csv"),
          row.names = FALSE)

# With only one sample per condition, cells are not independent biological
# replicates. Do not run condition-level FindMarkers or claim disease DE here.
summary_lines <- c(
  "Comparison labels: AD1 (one sample) and Control2 (one sample).",
  paste("Cells retained after QC:", sum(qc_before$cells_after_filtering)),
  paste("AD1 cells:", qc_before$cells_after_filtering[qc_before$sample == "AD1"]),
  paste("Control2 cells:", qc_before$cells_after_filtering[qc_before$sample == "Control2"]),
  paste("Clusters at resolution 0.8:", length(unique(combined$RNA_snn_res.0.8))),
  "Condition-level differential expression was not tested: there is one sample per condition.",
  "Cluster marker tests describe differences between cell clusters, not AD biomarkers."
)
writeLines(summary_lines, file.path(output_dir, "analysis_summary.txt"))

cat(paste(summary_lines, collapse = "\n"), "\n")
