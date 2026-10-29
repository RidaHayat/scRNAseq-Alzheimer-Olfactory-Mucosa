#Project2

library(Seurat)
#creating the Sparse matrices
sparse_matrix_AD1 <- Seurat::Read10X(data.dir = "AD1" )
sparse_matrix_Control2 <- Seurat::Read10X(data.dir = "Control2" )

#creating the Seurat objects
seu_AD1 <- Seurat::CreateSeuratObject(counts = sparse_matrix_AD1,
                                  project = "SCRNASeq_AD1",
                                  min.cells = 3,
                                  min.features = 100)
seu_Control2 <- Seurat::CreateSeuratObject(counts = sparse_matrix_Control2,
                                      project = "SCRNASeq_Control2",
                                      min.cells = 3,
                                      min.features = 100)
#plotting to visualize the datasets
Seurat::FeatureScatter(seu_AD1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
Seurat::FeatureScatter(seu_Control2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#Visualizing QC per cell and gene

Seurat::VlnPlot(seu_AD1, features = c("nCount_RNA",
                                  "nFeature_RNA"))
Seurat::VlnPlot(seu_Control2, features = c("nCount_RNA",
                                  "nFeature_RNA"))

#Calculating the percentage of counts coming from transcript types
# mitochondrial genes for AD1
seu_AD1 <- Seurat::PercentageFeatureSet(seu_AD1, 
                                    pattern = "^MT-", 
                                    col.name = "percent.mito")
# ribosomal genes for AD1
seu_AD1 <- Seurat::PercentageFeatureSet(seu_AD1, 
                                    pattern = "^RP[SL]",
                                    col.name = "percent.ribo")
# hemoglobin genes (but not HBP) for AD1
seu_AD1 <- Seurat::PercentageFeatureSet(seu_AD1,
                                    pattern = "^HB[^(P)]",
                                    col.name = "percent.globin")

# mitochondrial genes for Control2
seu_Control2 <- Seurat::PercentageFeatureSet(seu_Control2, 
                                        pattern = "^MT-", 
                                        col.name = "percent.mito")
# ribosomal genes for Control2
seu_Control2 <- Seurat::PercentageFeatureSet(seu_Control2, 
                                        pattern = "^RP[SL]",
                                        col.name = "percent.ribo")
# hemoglobin genes (but not HBP) for Control2
seu_Control2 <- Seurat::PercentageFeatureSet(seu_Control2,
                                        pattern = "^HB[^(P)]",
                                        col.name = "percent.globin")
#plotting the distributions

Seurat::VlnPlot(seu_AD1, features = c("percent.mito",
                                  "percent.ribo",
                                  "percent.globin"))
Seurat::VlnPlot(seu_Control2, features = c("percent.mito",
                                  "percent.ribo",
                                  "percent.globin"))

#using a boxplot to see which genes are the most highly expressed
library(ggplot2)
library(Matrix)

most_expressed_boxplot <- function(object, ngenes = 20){
  
  # matrix of raw counts
  cts <- Seurat::GetAssayData(object, assay = "RNA", slot = "counts")
  
  # get percentage/cell
  cts <- t(cts)/colSums(cts)*100
  medians <- apply(cts, 2, median)
  
  # get top n genes
  most_expressed <- order(medians, decreasing = T)[ngenes:1]
  most_exp_matrix <- as.matrix((cts[,most_expressed]))
  
  # prepare for plotting
  most_exp_df <- stack(as.data.frame(most_exp_matrix))
  colnames(most_exp_df) <- c("perc_total", "gene")
  
  # boxplot with ggplot2
  boxplot <- ggplot(most_exp_df, aes(x=gene, y=perc_total)) +
    geom_boxplot() +
    coord_flip()
  return(boxplot)
}
most_expressed_boxplot(seu_AD1, 20)
most_expressed_boxplot(seu_Control2, 20)
#MALAT1 and MT2A switched: in the Control dataset MALAT1 is the most expressed while MT2A is the third

#Cell filtering
seu_AD1 <- subset(seu_AD1, subset = nFeature_RNA > 200 & 
                nFeature_RNA < 5000 &
                percent.mito < 10)

seu_Control2 <- subset(seu_Control2, subset = nFeature_RNA > 200 & 
                    nFeature_RNA < 5000 &
                    percent.mito < 10)
#Plotting to see the results of the filtering
Seurat::VlnPlot(seu_AD1, features = c("percent.mito",
                                      "percent.ribo",
                                      "percent.globin"))
Seurat::VlnPlot(seu_Control2, features = c("percent.mito",
                                           "percent.ribo",
                                           "percent.globin"))
# seu_AD1
#An object of class Seurat 
#17331 features across 2003 samples within 1 assay 
#Active assay: RNA (17331 features, 0 variable features)

#Seu_Control2
#An object of class Seurat 
#17521 features across 1715 samples within 1 assay 
#Active assay: RNA (17521 features, 0 variable features)

##Normalization_Single_Cell
seu_merge<-merge(seu_AD1,y=seu_Control2,add.cell.ids = c("AD1", "Control2"))

#Normalization
seu_merge <- Seurat::NormalizeData(seu_merge,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)


#Calculating the subset of features with high cell-to-cell variation
seu_merge <- Seurat::FindVariableFeatures(seu_merge,
                                    selection.method = "vst",
                                    nfeatures = 2000)
top10_merge <- head(Seurat::VariableFeatures(seu_merge), 10)
top10_merge

#plotting in scatter plot
vf_plot_merge <- Seurat::VariableFeaturePlot(seu_merge)
Seurat::LabelPoints(plot = vf_plot_merge,
                    points = top10_merge, repel = TRUE)



#scaling

seu_merge <- Seurat::ScaleData(seu_merge,
                         features = rownames(seu_merge))


##Dimensionality reduction
#Merging the two datasets
seu_merge <- Seurat::RunPCA(seu_merge)
#Running the pca analysis
Seurat::DimPlot(seu_merge, reduction = "pca")
Seurat::FeaturePlot(seu_merge, reduction = "pca", features = "nFeature_RNA")
#heatmaps
Seurat::DimHeatmap(seu_merge, dims = 1:12, cells = 500, balanced = TRUE)
#elbow plot
Seurat::ElbowPlot(seu_merge, ndims = 40)
#plotting the umap plot
seu_merge <- Seurat::RunUMAP(seu_merge, dims = 1:25)
Seurat::DimPlot(seu_merge, reduction = "umap")
seu_merge <- Seurat::FindNeighbors(seu_merge, dims = 1:25)
#Clustering

seu_merge <- Seurat::FindClusters(seu_merge, resolution = seq(0.1, 0.8, by=0.1))
library(clustree)
clustree::clustree(seu_merge@meta.data[,grep("RNA_snn_res", colnames(seu_merge@meta.data))],
                   prefix = "RNA_snn_res.")

Seurat::DimPlot(seu_merge, group.by = "RNA_snn_res.0.8")
#Finding the markers for all clusters
de_genes <- Seurat::FindAllMarkers(seu_merge,  min.pct = 0.25,
                                   only.pos = TRUE)

de_genes <- subset(de_genes, de_genes$p_val_adj<0.05)

write.csv(de_genes, "de_genes_FindAllMarkers.csv", row.names = F, quote = F)

# Select the top 10 biomarkers for each cluster
library(dbplyr)
library(dplyr)

top_specific_markers10 <- de_genes %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)
top_specific_markers10
#top_specific_markers <- de_genes %>%
#group_by(cluster) %>%
 # top_n(20, avg_log2FC)
#top_specific_markers

#top_specific_markers3 <- de_genes %>%
 # group_by(cluster) %>%
  #top_n(3, avg_log2FC)
#top_specific_markers3


#Finding the biomarkers differentiating by control group and test group

de_genes_cluster <- Seurat::FindMarkers(seu_merge,group.by = "orig.ident",ident.1 ="SCRNASeq_AD1" ,ident.2 ="SCRNASeq_Control2" )

# Select the top 10 biomarkers

top_specific_markers_group10 <- de_genes_cluster %>%
  top_n(10, avg_log2FC)
top_specific_markers_group10


seu_merge <- Seurat::SetIdent(seu_merge, value = "SingleR_annot")
deg <- Seurat::FindMarkers(seu_merge,
                                   ident.1 = "SCRNASeq_AD1",
                                   ident.2 = "SCRNASeq_Control2",
                                   group.by = seu_merge$orig.ident,
                                   test.use = "wilcox")
deg <- subset(deg, deg$p_val_adj<0.05)


top_specific_markers_group10 <- de_genes_cluster %>%
  top_n(10, avg_log2FC)
top_specific_markers_group10

# TODO: specify valid cluster identities before running a cluster-level FindMarkers comparison.
#umap->clustering->find cell population in clusters
#find genes expressed more in that specific cluster

#merge(seu1,seu2)<-before normalization
#goals: celltype annotation (markers, diff expr btw classes, check if genes are biomarkers)
#id of biomarkers <-label based on condition
#findmarkers(by="condition",ad,control)<-find diff expr btw clusters

saveRDS(seu_merge, "seu_merge.rds")

