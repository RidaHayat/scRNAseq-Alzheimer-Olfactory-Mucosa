# Single-cell RNA-seq analysis of Alzheimer’s disease olfactory mucosa

This project compares olfactory mucosa samples labeled `AD1` and `Control2`. The R script uses Seurat for quality control, normalization, dimensionality reduction, clustering, and gene-expression analysis. Figures from the analysis are included below.

## 1. Methods and parameters

| Step | Method | Parameters used in the script |
| --- | --- | --- |
| Data input | `Seurat::Read10X`, `CreateSeuratObject` | Input folders: `AD1` and `Control2`; `min.cells = 3`; `min.features = 100` |
| QC measurements | `PercentageFeatureSet` | Mitochondrial `^MT-`; ribosomal `^RP[SL]`; hemoglobin `^HB[^(P)]` |
| Cell filtering | `subset()` | `nFeature_RNA > 200`; `nFeature_RNA < 5000`; `percent.mito < 10` |
| Normalization | `NormalizeData` | `normalization.method = "LogNormalize"`; `scale.factor = 10000` |
| Variable genes and scaling | `FindVariableFeatures`, `ScaleData` | `selection.method = "vst"`; `nfeatures = 2000`; scales all genes |
| PCA and UMAP | `RunPCA`, `ElbowPlot`, `RunUMAP` | Elbow plot through 40 PCs; dimensions 1–25 used for UMAP and neighbor finding |
| Clustering | `FindNeighbors`, `FindClusters`, `clustree` | Resolutions 0.1–0.8 by 0.1; resolution 0.8 is shown in the final cluster UMAP |
| Cluster markers | `FindAllMarkers` | `min.pct = 0.25`; `only.pos = TRUE`; script filters to `p_val_adj < 0.05` |
| AD1 vs Control2 | `FindMarkers` | Script compares the two sample identities. A later call specifies `test.use = "wilcox"` and filters adjusted p-values below 0.05; confirm its grouping and identity settings before interpreting its output. |

## 2. Results and figures

### Cell numbers after filtering

The script records **2,003 retained AD1 cells** and **1,715 retained Control2 cells**, for **3,718 cells total**. It does not include the initial cell totals in its analysis output, so removal percentages are not reported here.

### Variable and abundant genes

The variable-feature plot labels the genes selected as the most variable in the merged analysis. The separate top-gene plots show the most abundant transcripts within each sample. The script comments that `MALAT1` ranks first and `MT2A` third in Control2. These are abundance rankings; the figures do not establish a statistically significant AD1-versus-Control2 change.

### Clusters and cell types

The script tests clustering resolutions from 0.1 to 0.8 and plots the resolution 0.8 result. The included UMAP and clustree show the cluster layout and how assignments change across resolutions. The script and figures included here do not provide cell-type labels or a SingleR annotation workflow, so specific cell identities are not assigned in this project summary.

## 3. Figures

### Quality control and gene expression

The feature scatter plots compare the number of detected genes and RNA counts per cell. The violin plots show QC distributions, including detected genes/counts and transcript fractions. The top-gene plots show the most abundant transcripts in each sample.

| AD1 | Control2 |
| --- | --- |
| ![AD1 feature scatter](figures/01_feature_scatter_ad1.jpeg) | ![Control2 feature scatter](figures/01_feature_scatter_control2.jpeg) |
| ![AD1 QC violins](figures/02_qc_violins_ad1.jpeg) | ![Control2 QC violins](figures/02_qc_violins_control2.jpeg) |
| ![AD1 mitochondrial, ribosomal, and globin fractions](figures/03_mito_ribo_ad1.jpeg) | ![Control2 mitochondrial, ribosomal, and globin fractions](figures/03_mito_ribo_control2.jpeg) |
| ![Most abundant genes in AD1](figures/04_top_genes_ad1.jpeg) | ![Most abundant genes in Control2](figures/04_top_genes_control2.jpeg) |

### Variable genes and dimensionality reduction

The variable-feature plot shows selected variable genes. PCA displays the principal-component structure, the feature plot colors the PCA by detected gene count, and the elbow plot shows the variance drop-off used to choose dimensions.

![Variable features](figures/05_variable_features.png)

![PCA plot](figures/06_pca_plot.jpeg)

![PCA feature plot](figures/07_feature_plot_pca.jpeg)

![Elbow plot](figures/08_elbow_plot.jpeg)

### UMAP and clustering

The condition UMAP colors cells by sample. The clustree plot compares cluster assignments at different resolutions. The final plot shows cluster labels at resolution 0.8.

![UMAP by sample](figures/09_umap_condition.jpeg)

![Clustree across resolutions](figures/10_clustree.jpeg)

![UMAP clusters at resolution 0.8](figures/11_umap_clusters_res08.jpeg)

## Files and data

- `scRNAseq_Alzheimers_Olfactory_Mucosa.R` — analysis script; it expects 10X input folders named `AD1` and `Control2` alongside the script.
- `figures/` — figures displayed in this README.

The count matrices and archives are not included.
