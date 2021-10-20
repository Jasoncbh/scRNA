rm(list = ls())
library(monocle3)
library(tidyverse)
library(imputeTS)
library(Seurat)
library(RColorBrewer)
library(Matrix)
library(cowplot)
library(future)
require(scales)
library(ggthemes)
library(sctransform)
library(patchwork)


values <- c(brewer.pal(9,"Set1"),
            brewer.pal(8,"Dark2"))


tem_raw <- list.files(path = "../sc_RNA/",pattern="*_gene_exon_tagged.dge.txt.gz")
tem_raw
temp <- tem_raw
temp
name <- character()
filter_umi <- 450
#########批量导入数据#######
for(i in 1:length(temp)){
  name[i] <- unlist(strsplit(temp[i],"_out_gene_exon_tagged.dge.txt.gz"))[1]
  message(paste(name[i], "is loading"))
  
  tmpvalue<-read.table(paste0("../sc_RNA/", temp[i]), sep = "\t", quote = "", row.names = 1, header = T)
  message(paste(name[i], "is loaded, now is adding name"))
  
  colnames(tmpvalue) <- paste0(name[i], "-", colnames(tmpvalue))
  message(paste0(name[i], "'s name added, now filtering ", filter_umi))
  
  tmpvalue <- tmpvalue[,colSums(tmpvalue) >= filter_umi]
  message(paste(name[i], "cells above", filter_umi, "filtered"))
  
  assign(name[i], tmpvalue)
  rm(tmpvalue)
}
message("data loading done, and strat merge counts file")
############################

# tmpvalue<-read.table(paste0("../sc_RNA/", tem_raw[1]), sep = "\t", quote = "", row.names = 1, header = T)
# name[1] <- unlist(strsplit(tem_raw[1],"_out_gene_exon_tagged.dge.txt.gz"))[1]
# colnames(tmpvalue) <- paste0(name[1], "-", colnames(tmpvalue))
# #######过滤#########
# tmpvalue <- tmpvalue[, colSums(tmpvalue) >= filter_umi]
# message(paste(name[1], "cells above", filter_umi, "filtered"))
# assign(name[1], tmpvalue)

############################################################################
summary <- NULL
summary <- cbind(summary, as.matrix(summary(colSums(get(name[1])))))
colnames(summary) <- paste0(name, "'s UMI")
summary

############################################################################..............data-merge
# dge <- get(name[1])
# for(p in 2:length(temp)) {
#   dge_temp <- get(name[p])
#   dge_merge <- dplyr::full_join(dge %>% mutate(name = rownames(dge)),
#                                 dge_temp %>% mutate(name = rownames(dge_temp)))
#   rownames(dge_merge) <- dge_merge$name
#   dge <-  dplyr::select(dge_merge, -(name)) %>% na_replace(0)
#   rm(dge_temp)
#   rm(dge_merge)
# }
# message("data merge done")
# saveRDS(dge, file = "merge.rds")

############################################################################...........单个数据导入
# dge <- get(name[1])
# dge_temp <- get(name[1])
# dge_merge <- dplyr::full_join(dge %>% mutate(name = rownames(dge)),
#                               dge_temp %>% mutate(name = rownames(dge_temp)))
# rownames(dge_merge) <- dge_merge$name
# dge <-  dplyr::select(dge_merge, -(name)) %>% na.replace(0)
# rm(dge_temp)
# rm(dge_merge)
# message("data merge done")
# saveRDS(dge, file = "AE.rds")
# rm(dge)
#################################################################
# dge <- readRDS("merge.rds")
########过滤条件########
filter_gene = 350
filter_cell = 5
########################

metadata_AE <- data.frame(
  matrix(unlist(strsplit(colnames(`2nd_AML1_ETO`),"-")), ncol = 2, byrow = T),
  row.names = colnames(`2nd_AML1_ETO`)
)
colnames(metadata_AE) <- c("tech","UMI")
metadata_AE$celltype <- "undefined"
message("Metadata_AE done")

metadata_Mig <- data.frame(
  matrix(unlist(strsplit(colnames(`2nd_MigR1`),"-")), ncol = 2, byrow = T),
  row.names = colnames(`2nd_MigR1`)
)
colnames(metadata_Mig) <- c("tech","UMI")
metadata_Mig$celltype <- "undefined"
message("Metadata_Mig done")

##########构建seurat对象###############
pbmc_AE <- CreateSeuratObject(counts = `2nd_AML1_ETO`, 
                           project = "AE", 
                           min.features = filter_gene, 
                           min.cells = filter_cell, 
                           meta.data = metadata_AE)

pbmc_Mig <- CreateSeuratObject(counts = `2nd_MigR1`, 
                           project = "Mig", 
                           min.features = filter_gene, 
                           min.cells = filter_cell, 
                           meta.data = metadata_Mig)

class(object.list)
object.list <- as.list(c(pbmc_AE,pbmc_Mig))
for (i in 1:length(object.list)){
  names(object.list)[i] <- paste("data_", i, sep = "")
}

# object.list <- c("pbmc_AE","pbmc_Mig")
 


##########过滤线粒体基因###############
pbmc_AE[["percent.mt"]] <- PercentageFeatureSet(pbmc_AE, pattern = "^MT-")
#nFeature_RNA：每个细胞所检测到的基因数目
#nCount_RNA：每个细胞的UMI数目
VlnPlot(pbmc_AE, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc_AE, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "tech")
plot2 <- FeatureScatter(pbmc_AE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "tech")
plot1 + plot2

####################################
hist(pbmc$percent.mt, breaks = 100)
hist(pbmc$nCount_RNA, breaks = 400)
hist(pbmc$nFeature_RNA, breaks = 400)

####################################
#############过滤###################
pbmc <- subset(pbmc, subset = nFeature_RNA > 350 & nFeature_RNA < 2000 & percent.mt < 10)

####################################---------normalization标准化
####################################
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
dev.off
dev.new()
plot1 + plot2
##############......................正常scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#####################################----------应用 sctransform 归一化
# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = c("percent.mt", "tech"))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#####################################------------确定使用PC个数
# each PC essentially representing a ‘metafeature’
# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc, ndims = 25)

#####################################------------对细胞聚类
# 首先基于PCA空间构建一个基于欧氏距离的KNN图
pbmc <- FindNeighbors(pbmc, dims = 1:21)
pbmc <- FindClusters(pbmc, resolution = 0.3)

pbmc <- RunUMAP(pbmc,
                dims = 1:21,
                reduction = "pca")

DimPlot(object = pbmc, pt.size = 1, reduction = "umap", group.by = "tech",label = T, label.size = 3)
DimPlot(object = pbmc, pt.size = 1, reduction = "umap",label = T, label.size = 5)

pbmc <- RunTSNE(pbmc,
                reduction = "pca",
                cells = NULL,
                dims = 1:21)

DimPlot(object = pbmc, pt.size = 1, reduction = "tsne", cols = "Set1")

DimPlot(object = pbmc, pt.size = 1, reduction = "tsne", label = T,  cols = "Paired")
DimPlot(object = pbmc, pt.size = 1, reduction = "tsne", split.by = "tech", label = T, repel = T, label.size = 5)
saveRDS(pbmc, file = "pbmc_test.rds")
VlnPlot(pbmc, features = c("NTS", "GZMB"))

cluster1.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster1.markers, n = 5)
VlnPlot(pbmc, features = c("TMSB10", "SAMSN1 "))


# only.pos：只保留上调差异表达的基因
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
head(pbmc.markers)
tail(pbmc.markers)
write.csv(pbmc.markers,file = "AE_sctransform.markers.csv")

# 选择一些基因进行可视化
VlnPlot(pbmc, features = c("SCP2", "SH3BGRL3"))
FeaturePlot(pbmc, features = c("LGALS1", "HPGDS", "TM4SF1", "CDK4", "SMC4", "LCP1", "LYZ", "CST3"))

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(pbmc, features = top5$gene) + NoLegend()

head(Idents(pbmc))

pbmc@meta.data$cell_anno <- Idents(pbmc)
write.csv(pbmc@meta.data,file = "metadata_sctransform.csv")
saveRDS(pbmc, file = "AE-Mig-sctransform_final.rds")
################################################################################
###########################....scCATCH........##################################

library(scCATCH)
?scCATCH
clu_ann2 <- scCATCH(object = pbmc.markers,species = "Human",tissue = c('Peripheral blood'))
new.cluster.ids2 <- clu_ann2$cell_type
names(new.cluster.ids2) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids2)

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 2)


################################################################################
###########################....SingleR........##################################

library(SingleR)
# library(remotes)
# library(celldex)
# refdata <- HumanPrimaryCellAtlasData()
# save(refdata,file = "refdata.RData")
# rm(refdata)
load("refdata.RData")
refdata
head(colnames(refdata))
head(rownames(refdata))

# 查看共有多少种细胞类型
unique(refdata@colData@listData[["label.main"]])

testdata <- GetAssayData(pbmc, slot="data")
dim(testdata)
testdata[1:30,1:4]
clusters <- pbmc@meta.data$seurat_clusters
table(clusters)

cellpred <- SingleR(test = testdata,  
                    ref = refdata, 
                    labels = refdata$label.main,
                    method = "cluster", 
                    clusters = clusters,
                    assay.type.test = "logcounts", 
                    assay.type.ref = "logcounts")

str(cellpred,max.level = 3)
metadata <- cellpred@metadata
head(metadata)
celltype = data.frame(ClusterID = rownames(cellpred), 
                      celltype = cellpred$labels, 
                      stringsAsFactors = F)
celltype
############################################################################
############################################################################

stat <- data.frame(cellident = pbmc@meta.data$tech, clusterident = pbmc@meta.data$seurat_clusters) %>%
  group_by(clusterident, cellident) %>%
  summarise(count=n()) %>%
  group_by(clusterident)

stat %>%  
  summarise(n()) %>%
  dplyr::filter(`n()` == 1)

table <- data.frame(table(pbmc@meta.data$tech))
info <- data.frame(clusterident =as.factor(0:(length(stat$clusterident)/length(levels)-1)))

for (i in 1:length(table$Var1)) {
  percentage <- left_join(stat %>% dplyr::filter(cellident == table$Var1[i]) %>% summarise(count / table$Freq[i] * 100),
                          stat %>% dplyr::filter(cellident == table$Var1[i]) %>% summarise(count), 
                          by = "clusterident")
  colnames(percentage) <- c("clusterident", paste0(table$Var1[i], "_ratio"), paste0(table$Var1[i], "_counts"))
  info <- left_join(info, percentage, by = "clusterident")
}

result <- info %>% mutate(ratio = info[,2]/info[,3])

summary(pbmc@meta.data$nFeature_RNA)
summary(pbmc@meta.data$nCount_RNA)




##########################............................Monocle3
rm(list = ls())
pbmc <- readRDS("AE-Mig-sctransform_final.rds")
library(Seurat)
library(monocle3)
expression_matrix <- as(as.matrix(pbmc@assays$RNA@counts),'sparseMatrix')
cell_metadata <- pbmc@meta.data 
gene_annotation <- data.frame(gene_short_name = row.names(pbmc@assays$RNA@counts),
                              row.names = row.names(pbmc@assays$RNA@counts))
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
##########预处理
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)

plot_pc_variance_explained(cds)

##########数据降维
##################

###########...tSNE 
cds <- reduce_dimension(cds, reduction_method="tSNE")
plot_cells(cds, reduction_method="tSNE", color_cells_by="cell_anno",label_cell_groups = T,cell_size = 1)

##########cluster cell

cds = cluster_cells(cds,
                    reduction_method = "tSNE",
                    k = 150)
plot_cells(cds, reduction_method = "tSNE", cell_size = 1)
plot_cells(cds, reduction_method = "tSNE", group_cells_by = "cluster")

##########...UMAP
cds <- reduce_dimension(cds,
                        cores = 1,
                        reduction_method = "UMAP",
                        preprocess_method = "PCA")
plot_cells(cds, 
           reduction_method = "UMAP", 
           color_cells_by = "cell_anno", 
           label_cell_groups = FALSE,
           cell_size = 1)
plot_cells(cds,
           reduction_method = "UMAP",
           genes = c("TM4SF1"),
           cell_size = 1)
plot_cells(cds, genes = c("TM4SF1", "CRHBP", "GNG11", "RUNX1T1"),cell_size = 1)
#####再次以UMAP进行聚类
cds <- cluster_cells(cds,reduction_method = "UMAP",k = 20)
plot_cells(cds, color_cells_by = "partition")


###########..........trajectory graph
#####################################
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cell_anno",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <-  order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)











