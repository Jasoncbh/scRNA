rm(list=ls())

library(dplyr)
library(Seurat)
library(patchwork)

## =============1.Load the dataset
pbmc.data <- read.table(gzfile("2nd_AML1_ETO_out_gene_exon_tagged.dge.txt.gz"),sep = "",header = T)
pbmc.data
# Initialize the Seurat object with the raw (non-normalized data)
# min.cell：每个feature至少在多少个细胞中表达
# min.features：每个细胞中至少有多少个feature被检测到
a1 <- pbmc.data[1:5000,]
pbmc <- CreateSeuratObject(counts = a1, 
                           project = "pbmc", 
                           min.cells = 3,
                           min.features = 200)
pbmc
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#nFeature_RNA：每个细胞所检测到的基因数目
#nCount_RNA：每个细胞的UMI数目
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
######normalization标准化
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off
##########归一化##########
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
########################################
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

## =============8.确定使用PC个数
# each PC essentially representing a ‘metafeature’
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

#9.对细胞聚类
# 首先基于PCA空间构建一个基于欧氏距离的KNN图
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

pbmc <- RunUMAP(pbmc, dims = 1:10)

## =============10.将细胞在低维空间可视化UMAP/tSNE
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)

# 可视化
DimPlot(pbmc, reduction = "umap", label = T, label.size = 5)
DimPlot(pbmc, reduction = "tsne", label = T, label.size = 5)
saveRDS(pbmc, file = "data/pbmc_tutorial.rds")

## =============11.差异表达分析
# 在cluster2 vs else中差异表达
cluster1.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster1.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)


# 所有类的差异表达基因
# only.pos：只保留上调差异表达的基因
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(pbmc.markers,file = "data/pbmc.markers.csv")
head(pbmc.markers)


# 选择一些基因进行可视化，作者这里根据自己的知识背景选择的相关基因
VlnPlot(pbmc, features = c("BCL6", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

new.cluster.ids <- c("Naive CD4 T", 
                     "CD14+ Mono", 
                     "Memory CD4 T",
                     "B", 
                     "CD8 T",
                     "FCGR3A+ Mono",
                     "NK", 
                     "DC", 
                     "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 1.2) + NoLegend()
# 保存结果
pbmc@meta.data$cell_anno <- Idents(pbmc)
write.csv(pbmc@meta.data,file = "data/metadata.csv")
saveRDS(pbmc, file = "data/pbmc3k_final.rds")
