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
library(harmony)
values <- c(brewer.pal(9,"Set1"),
            brewer.pal(8,"Dark2"))
abc <- values[1:14]
tem_raw <- list.files(path = "../sc_RNA/Single-cell RNA-seq Rawdata/",pattern="*_gene_exon_tagged.dge.txt.gz")
tem_raw
temp <- tem_raw
temp
name <- character()
filter_umi <- 450
#########批量导入数据#######
for(i in 1:length(temp)){
  name[i] <- unlist(strsplit(temp[i],"_out_gene_exon_tagged.dge.txt.gz"))[1]
  message(paste(name[i], "is loading"))
  
  tmpvalue<-read.table(paste0("../sc_RNA/Single-cell RNA-seq Rawdata/", temp[i]), sep = "\t", quote = "", row.names = 1, header = T)
  message(paste(name[i], "is loaded, now is adding name"))
  
  colnames(tmpvalue) <- paste0(name[i], "-", colnames(tmpvalue))
  message(paste0(name[i], "'s name added, now filtering ", filter_umi))
  
  tmpvalue <- tmpvalue[,colSums(tmpvalue) >= filter_umi]
  message(paste(name[i], "cells above", filter_umi, "filtered"))
  
  assign(name[i], tmpvalue)
  rm(tmpvalue)
}
message("data loading done, and strat merge counts file")


###########################################...metadata
metadata_AE_1 <- data.frame(
  matrix(unlist(strsplit(colnames(`1st_AML1_ETO`),"-")), ncol = 2, byrow = T),
  row.names = colnames(`1st_AML1_ETO`)
)
colnames(metadata_AE_1) <- c("tech","UMI")
metadata_AE_1$celltype <- "undefined"
metadata_AE_1$type <- "AE"
metadata_AE_1$batch <- "1st"
metadata_AE_1$orig.ident <- "AE_1"
message("metadata_AE_1 done")
###
metadata_Mig_1 <- data.frame(
  matrix(unlist(strsplit(colnames(`1st_MigR1`),"-")), ncol = 2, byrow = T),
  row.names = colnames(`1st_MigR1`)
)
colnames(metadata_Mig_1) <- c("tech","UMI")
metadata_Mig_1$celltype <- "undefined"
metadata_Mig_1$type <- "MigR1"
metadata_Mig_1$orig.ident <- "MigR1_1"
metadata_Mig_1$batch <- "1st"
message("Metadata_Mig_1 done")
###
metadata_AE_2 <- data.frame(
  matrix(unlist(strsplit(colnames(`2nd_AML1_ETO`),"-")), ncol = 2, byrow = T),
  row.names = colnames(`2nd_AML1_ETO`)
)
colnames(metadata_AE_2) <- c("tech","UMI")
metadata_AE_2$celltype <- "undefined"
metadata_AE_2$type <- "AE"
metadata_AE_2$orig.ident <- "AE_2"
metadata_AE_2$batch <- "2nd"

message("Metadata_AE_2 done")
###
metadata_Mig_2 <- data.frame(
  matrix(unlist(strsplit(colnames(`2nd_MigR1`),"-")), ncol = 2, byrow = T),
  row.names = colnames(`2nd_MigR1`)
)
colnames(metadata_Mig_2) <- c("tech","UMI")
metadata_Mig_2$celltype <- "undefined"
metadata_Mig_2$type <- "MigR1"
metadata_Mig_2$orig.ident <- "MigR1_2"
metadata_Mig_2$batch <- "2nd"

message("Metadata_Mig_2 done")
###
metadata_AE_3 <- data.frame(
  matrix(unlist(strsplit(colnames(`3rd_AML1_ETO`),"-")), ncol = 2, byrow = T),
  row.names = colnames(`3rd_AML1_ETO`)
)
colnames(metadata_AE_3) <- c("tech","UMI")
metadata_AE_3$celltype <- "undefined"
metadata_AE_3$type <- "AE"
metadata_AE_3$orig.ident <- "AE_3"
metadata_AE_3$batch <- "3rd"

message("Metadata_AE_3 done")
###
metadata_Mig_3 <- data.frame(
  matrix(unlist(strsplit(colnames(`3rd_MigR1`),"-")), ncol = 2, byrow = T),
  row.names = colnames(`3rd_MigR1`)
)
colnames(metadata_Mig_3) <- c("tech","UMI")
metadata_Mig_3$celltype <- "undefined"
metadata_Mig_3$type <- "MigR1"
metadata_Mig_3$orig.ident <- "MigR1_3"
metadata_Mig_3$batch <- "3rd"

message("Metadata_Mig_3 done")
###
metadata_AE_4 <- data.frame(
  matrix(unlist(strsplit(colnames(`4th_AML1_ETO`),"-")), ncol = 2, byrow = T),
  row.names = colnames(`4th_AML1_ETO`)
)
colnames(metadata_AE_4) <- c("tech","UMI")
metadata_AE_4$celltype <- "undefined"
metadata_AE_4$type <- "AE"
metadata_AE_4$orig.ident <- "AE_4"
metadata_AE_4$batch <- "4th"

message("Metadata_AE_4 done")
###
metadata_Mig_4 <- data.frame(
  matrix(unlist(strsplit(colnames(`4th_MigR1`),"-")), ncol = 2, byrow = T),
  row.names = colnames(`4th_MigR1`)
)
colnames(metadata_Mig_4) <- c("tech","UMI")
metadata_Mig_4$celltype <- "undefined"
metadata_Mig_4$type <- "MigR1"
metadata_Mig_4$orig.ident <- "MigR1_4"
metadata_Mig_4$batch <- "4th"

message("Metadata_Mig_4 done")
#######################################
#######################################
filter_gene = 350
filter_cell = 5
##########构建seurat对象###############
pbmc_AE_1 <- CreateSeuratObject(counts = `1st_AML1_ETO`, 
                              project = "AE", 
                              min.features = filter_gene, 
                              min.cells = filter_cell, 
                              meta.data = metadata_AE_1)

pbmc_Mig_1 <- CreateSeuratObject(counts = `1st_MigR1`, 
                               project = "Mig", 
                               min.features = filter_gene, 
                               min.cells = filter_cell, 
                               meta.data = metadata_Mig_1)
#############################################################
pbmc_AE_2 <- CreateSeuratObject(counts = `2nd_AML1_ETO`, 
                                project = "AE", 
                                min.features = filter_gene, 
                                min.cells = filter_cell, 
                                meta.data = metadata_AE_2)

pbmc_Mig_2 <- CreateSeuratObject(counts = `2nd_MigR1`, 
                                 project = "Mig", 
                                 min.features = filter_gene, 
                                 min.cells = filter_cell, 
                                 meta.data = metadata_Mig_2)
#############################################################
pbmc_AE_3 <- CreateSeuratObject(counts = `3rd_AML1_ETO`, 
                                project = "AE", 
                                min.features = filter_gene, 
                                min.cells = filter_cell, 
                                meta.data = metadata_AE_3)

pbmc_Mig_3 <- CreateSeuratObject(counts = `3rd_MigR1`, 
                                 project = "Mig", 
                                 min.features = filter_gene, 
                                 min.cells = filter_cell, 
                                 meta.data = metadata_Mig_3)
#############################################################
pbmc_AE_4 <- CreateSeuratObject(counts = `4th_AML1_ETO`, 
                                project = "AE", 
                                min.features = filter_gene, 
                                min.cells = filter_cell, 
                                meta.data = metadata_AE_4)

pbmc_Mig_4 <- CreateSeuratObject(counts = `4th_MigR1`, 
                                 project = "Mig", 
                                 min.features = filter_gene, 
                                 min.cells = filter_cell, 
                                 meta.data = metadata_Mig_4)
##############################################################
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc_AE_1[["percent.mt"]] <- PercentageFeatureSet(pbmc_AE_1, pattern = "^MT-")
pbmc_AE_1 <- subset(pbmc_AE_1, subset = nFeature_RNA > 350 & nFeature_RNA < 2000 & percent.mt < 7)
########################
pbmc_AE_2[["percent.mt"]] <- PercentageFeatureSet(pbmc_AE_2, pattern = "^MT-")
pbmc_AE_2 <- subset(pbmc_AE_2, subset = nFeature_RNA > 350 & nFeature_RNA < 2000 & percent.mt < 6)
########################
pbmc_AE_3[["percent.mt"]] <- PercentageFeatureSet(pbmc_AE_3, pattern = "^MT-")
pbmc_AE_3 <- subset(pbmc_AE_3, subset = nFeature_RNA > 350 & nFeature_RNA < 2000 & percent.mt < 10)
########################
pbmc_AE_4[["percent.mt"]] <- PercentageFeatureSet(pbmc_AE_4, pattern = "^MT-")
pbmc_AE_4 <- subset(pbmc_AE_4, subset = nFeature_RNA > 350 & percent.mt < 7)
########################
pbmc_Mig_1[["percent.mt"]] <- PercentageFeatureSet(pbmc_Mig_1, pattern = "^MT-")
pbmc_Mig_1 <- subset(pbmc_Mig_1, subset = nFeature_RNA > 350 & nFeature_RNA < 2000 & percent.mt < 7)
########################
pbmc_Mig_2[["percent.mt"]] <- PercentageFeatureSet(pbmc_Mig_2, pattern = "^MT-")
pbmc_Mig_2 <- subset(pbmc_Mig_2, subset = nFeature_RNA > 350 & nFeature_RNA < 2000 & percent.mt < 6)
########################
pbmc_Mig_3[["percent.mt"]] <- PercentageFeatureSet(pbmc_Mig_3, pattern = "^MT-")
pbmc_Mig_3 <- subset(pbmc_Mig_3, subset = nFeature_RNA > 350 & nFeature_RNA < 2000 & percent.mt < 10)
########################
pbmc_Mig_4[["percent.mt"]] <- PercentageFeatureSet(pbmc_Mig_4, pattern = "^MT-")
pbmc_Mig_4 <- subset(pbmc_Mig_4, subset = nFeature_RNA > 350 & percent.mt < 7)

#########################################################################################################
object.list <- as.list(c(pbmc_AE_1,pbmc_Mig_1,pbmc_AE_2,pbmc_Mig_2,pbmc_AE_3,pbmc_Mig_3,pbmc_AE_4,pbmc_Mig_4))

sample_name <- c("AE1","Mig1","AE2","Mig2","AE3","Mig3","AE4","Mig4")
names(object.list) <- sample_name
# saveRDS(object.list,file = "object_list.rds")
object.list <- readRDS("object_list.rds")
##PCA降维
scRNA_harmony <- merge(object.list[[1]], y=c(object.list[[2]], object.list[[3]], object.list[[4]], object.list[[5]], 
                                             object.list[[6]], object.list[[7]], object.list[[8]]))
table(scRNA_harmony$batch)
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)

# scRNA_harmony <- JackStraw(scRNA_harmony, num.replicate = 100)
##整合
# system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident", kmeans_init_nstart=1, kmeans_init_iter_max=30)

#降维聚类
scRNA_harmony <- RunTSNE(scRNA_harmony,reduction = "harmony",dims = 1:50,tsne.method = "Rtsne",reduction.name = "tsne")
# scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:50) %>% FindClusters(resolution = 0.5)

##作图
#group_by_cluster
 DimPlot(scRNA_harmony, reduction = "tsne", label=T) 
#group_by_sample
DimPlot(scRNA_harmony, reduction = "tsne", group.by='orig.ident',cols = values) 
DimPlot(scRNA_harmony, reduction = "tsne", group.by='type') 
DimPlot(scRNA_harmony, reduction = "tsne", split.by ='orig.ident',label = T,ncol = 4,cols = values) 
DimPlot(scRNA_harmony, reduction = "tsne", split.by ='type',label = T,pt.size = 1,cols = values,ncol = 1) 
par(mfrow=c(2,4))
FeaturePlot(scRNA_harmony, features = c("TM4SF1"),split.by = "type",reduction = "tsne",label = T,ncol = 1,pt.size = 1)
cluster3.markers <- FindMarkers(scRNA_harmony, ident.1 = 3, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
marker_gene<- FindConservedMarkers(scRNA_harmony, ident.1 = 3, grouping.var = "type", verbose = FALSE)
write.csv(cluster3.markers,file = "marker.csv")
#combinate
plotc <- plot1+plot2
ggsave("scRNA_harmony_batch.png", plot = plotc, width = 10, height = 5)
saveRDS(scRNA_harmony, 'scRNA_harmony_batch.rds')
# plot3 = DimPlot(scRNA_harmony, reduction = "umap", group.by=) 
#combinate
scRNA_harmony_type <- RunHarmony(scRNA_harmony, group.by.vars = "type")

scRNA_harmony_type <- RunUMAP(scRNA_harmony_type, reduction = "harmony", dims = 1:30)
scRNA_harmony_type <- FindNeighbors(scRNA_harmony_type, reduction = "harmony", dims = 1:30) %>% FindClusters()
plot3 = DimPlot(scRNA_harmony, reduction = "umap", group.by="type")
plot =DimPlot(scRNA_harmony, reduction = "umap", label=T)
markergene_harmony <- FindConservedMarkers(scRNA_harmony_type, ident.1 = 3, grouping.var = "type", verbose = FALSE)
head(markergene_harmony)
###################################################################
#........................................................................去除第一次实验
object.list_none1 <- as.list(c(pbmc_AE_2,pbmc_Mig_2,pbmc_AE_3,pbmc_Mig_3,pbmc_AE_4,pbmc_Mig_4))
scRNA_harmony_none1 <- merge(object.list_none1[[1]], y=c(object.list_none1[[2]], object.list_none1[[3]], object.list_none1[[4]], object.list_none1[[5]], 
                                                         object.list_none1[[6]]))
table(scRNA_harmony_none1$batch)
scRNA_harmony_none1 <- NormalizeData(scRNA_harmony_none1) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
##整合
# system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
scRNA_harmony_none1 <- RunHarmony(scRNA_harmony_none1, group.by.vars = "type")
#降维聚类
scRNA_harmony_none1 <- RunUMAP(scRNA_harmony_none1, reduction = "harmony", dims = 1:30)
scRNA_harmony_none1 <- FindNeighbors(scRNA_harmony_none1, reduction = "harmony", dims = 1:30) %>% FindClusters()

##作图
#group_by_cluster
plot4 = DimPlot(scRNA_harmony_none1, reduction = "umap", label=T) 
#group_by_sample
plot5 = DimPlot(scRNA_harmony_none1, reduction = "umap", group.by='type') 
#combinate
plotc_batch <- plot4+plot5

markergene_harmony <- FindConservedMarkers(batch_1.combined, ident.1 = 3, grouping.var = "orig.ident", verbose = FALSE)
head(markers_1)
###################################################################
#######seurat#####
for (i in 1:length(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]])
  object.list[[i]] <- FindVariableFeatures(object.list[[i]])
}

scRNA.anchors <- FindIntegrationAnchors(object.list = object.list)
scRNA_seurat <- IntegrateData(anchorset = scRNA.anchors)
scRNA_seurat <- ScaleData(scRNA_seurat) %>% RunPCA(verbose=FALSE)
scRNA_seurat <- RunUMAP(scRNA_seurat, dims = 1:30)
scRNA_seurat <- FindNeighbors(scRNA_seurat, dims = 1:30) %>% FindClusters()
#group_by_cluster
plot4 = DimPlot(scRNA_seurat, reduction = "tSNE", label=T) 
#group_by_sample
plot5 = DimPlot(scRNA_seurat, reduction = "umap", group.by='orig.ident') 
#combinate

plotcs <- plot4+plot5

marker_gene <- FindConservedMarkers(batch_1.combined, ident.1 = 3, grouping.var = "orig.ident", verbose = FALSE)
head(markers_1)
#######################################################################
#..........................................................................................Harmony批次间差异比较
######################################################################
object.list_1 <- as.list(c(pbmc_AE_1,pbmc_Mig_1))
object.list_2 <- as.list(c(pbmc_AE_2,pbmc_Mig_2))
object.list_3 <- as.list(c(pbmc_AE_3,pbmc_Mig_3))
object.list_4 <- as.list(c(pbmc_AE_4,pbmc_Mig_4))

#######################################################################
scRNA_harmony_batch1 <- merge(object.list[[1]], object.list[[2]])
scRNA_harmony_batch2 <- merge(object.list[[3]], object.list[[4]])
scRNA_harmony_batch3 <- merge(object.list[[5]], object.list[[6]])
scRNA_harmony_batch4 <- merge(object.list[[7]], object.list[[8]])
#######################################################################
table(scRNA_harmony_batch1$orig.ident)
scRNA_harmony_batch1 <- NormalizeData(scRNA_harmony_batch1) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
##整合
scRNA_harmony_batch1 <- RunHarmony(scRNA_harmony_batch1, group.by.vars = "orig.ident")
#降维聚类
scRNA_harmony_batch1 <- RunUMAP(scRNA_harmony_batch1, reduction = "harmony", dims = 1:30)
# scRNA_harmony_batch1 <- RunTSNE(scRNA_harmony_batch1, reduction = "harmony", dims = 1:30)
scRNA_harmony_batch1<- RunTSNE(scRNA_harmony_batch1,reduction = "harmony",dims = 1:30,seed.use = 1,tsne.method = "Rtsne",reduction.name = "tsne")
scRNA_harmony_batch1 <- FindNeighbors(scRNA_harmony_batch1, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 1)
#group_by_cluster
plot_batch1_c= DimPlot(scRNA_harmony_batch1, reduction = "tsne", label=T) 
   DimPlot(object = scRNA_harmony_batch1, pt.size = 1, reduction = "tsne", group.by = "orig.ident", cols =abc)
   DimPlot(object = scRNA_harmony_batch1, pt.size = 1, reduction = "tsne", cols = "Set1")
   FeaturePlot(scRNA_harmony_batch1, features = c("TM4SF1"),  reduction = "tsne",min.cutoff = "q9")
   DimPlot(scRNA_harmony_batch1, reduction = "tsne",pt.size = 1,cols = abc, split.by = "orig.ident",label=T)
   #group_by_sample
plot_batch1_s = DimPlot(scRNA_harmony_batch1, reduction = "tsne", group.by='orig.ident')

##..........................................................................................

scRNA_harmony_batch2 <- NormalizeData(scRNA_harmony_batch2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_harmony_batch2 <- RunHarmony(scRNA_harmony_batch2, group.by.vars = "orig.ident")
# scRNA_harmony_batch2 <- RunUMAP(scRNA_harmony_batch2, reduction = "harmony", dims = 1:30)
scRNA_harmony_batch2<- RunTSNE(scRNA_harmony_batch2,reduction = "harmony",dims = 1:30,seed.use = 1,tsne.method = "Rtsne",reduction.name = "tsne")
scRNA_harmony_batch2 <- FindNeighbors(scRNA_harmony_batch2, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 1)
plot_batch2_c= DimPlot(scRNA_harmony_batch2, pt.size = 1, reduction = "tsne", group.by = "orig.ident", cols =abc, label=T) 
DimPlot(object = scRNA_harmony_batch2, pt.size = 1, reduction = "tsne", group.by = "orig.ident", cols =abc)
DimPlot(object = scRNA_harmony_batch2, pt.size = 1, reduction = "tsne", cols =abc)
DimPlot(scRNA_harmony_batch2, reduction = "tsne",pt.size = 1,cols = abc, split.by = "orig.ident",label=T)
# plot_batch2_s = DimPlot(scRNA_harmony_batch2, reduction = "umap", group.by='orig.ident') 
FeaturePlot(scRNA_harmony_batch2, features = c("TM4SF1"),  reduction = "tsne",min.cutoff = "q9")
###########################################################################################


scRNA_harmony_batch3 <- NormalizeData(scRNA_harmony_batch3) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_harmony_batch3 <- RunHarmony(scRNA_harmony_batch3, group.by.vars = "orig.ident")
scRNA_harmony_batch3 <- RunUMAP(scRNA_harmony_batch3, reduction = "harmony", dims = 1:30)
scRNA_harmony_batch3<- RunTSNE(scRNA_harmony_batch3,reduction = "harmony",dims = 1:30,seed.use = 1,tsne.method = "Rtsne",reduction.name = "tsne")

scRNA_harmony_batch3 <- FindNeighbors(scRNA_harmony_batch3, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 1)
plot_batch3_c= DimPlot(scRNA_harmony_batch3, reduction = "tsne", label=T) 
plot_batch3_s = DimPlot(scRNA_harmony_batch3, pt.size = 1, reduction = "tsne", group.by = "orig.ident", cols =abc) 
FeaturePlot(scRNA_harmony_batch3, features = c("TM4SF1"),  reduction = "tsne")
DimPlot(scRNA_harmony_batch3, reduction = "tsne",pt.size = 1,cols = abc, split.by = "orig.ident",label=T)


###########################################################################################                              


scRNA_harmony_batch4 <- NormalizeData(scRNA_harmony_batch4) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_harmony_batch4 <- RunHarmony(scRNA_harmony_batch4, group.by.vars = "orig.ident")
# scRNA_harmony_batch4 <- RunUMAP(scRNA_harmony_batch4, reduction = "harmony", dims = 1:30)
scRNA_harmony_batch4<- RunTSNE(scRNA_harmony_batch4,reduction = "harmony",dims = 1:30,seed.use = 1,tsne.method = "Rtsne",reduction.name = "tsne")

scRNA_harmony_batch4 <- FindNeighbors(scRNA_harmony_batch4, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 1)
plot_batch4_c= DimPlot(scRNA_harmony_batch4, reduction = "tsne", label=T,cols =abc) 
plot_batch4_s = DimPlot(scRNA_harmony_batch4,reduction = "tsne", group.by='orig.ident',cols =abc)                            
FeaturePlot(scRNA_harmony_batch4, features = c("TM4SF1"),  reduction = "tsne")
DimPlot(scRNA_harmony_batch4, reduction = "tsne",pt.size = 1,cols = abc, split.by = "orig.ident",label=T)

############################################################################################
############################################################################################
#seurat 
object.list_1 <- c(object.list[[1]],object.list[[2]])
object.list_1 <- lapply(X = object.list_1, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = object.list_1)
batch_1.anchors <- FindIntegrationAnchors(object.list = object.list_1, anchor.features = features)
# this command creates an 'integrated' data assay
batch_1.combined <- IntegrateData(anchorset = batch_1.anchors)
DefaultAssay(batch_1.combined) <- "integrated"
batch_1.combined <- ScaleData(batch_1.combined, verbose = FALSE)
batch_1.combined <- RunPCA(batch_1.combined, npcs = 30, verbose = FALSE)
batch_1.combined<- RunTSNE(batch_1.combined,reduction = "pca",dims = 1:30,cells = NULL,features = NULL,tsne.method = "Rtsne",reduction.name = "tsne")
# batch_1.combined <- RunUMAP(batch_1.combined, reduction = "pca", dims = 1:30)
batch_1.combined <- FindNeighbors(batch_1.combined, reduction = "pca", dims = 1:30)
batch_1.combined <- FindClusters(batch_1.combined, resolution = 0.5)
p1 <- DimPlot(batch_1.combined, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(batch_1.combined, reduction = "tsne", label = TRUE, repel = TRUE)
p1 + p2
DimPlot(batch_1.combined,reduction = "tsne", group.by='orig.ident',cols =abc)
DimPlot(batch_1.combined,reduction = "tsne", split.by ='orig.ident',cols =abc,label=T)
FeaturePlot(batch_1.combined, features = c("TM4SF1"),  reduction = "tsne",label=T)


#..........................................................

object.list_2 <- c(object.list[[3]],object.list[[4]])

object.list_2 <- lapply(X = object.list_2, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = object.list_2)
batch_2.anchors <- FindIntegrationAnchors(object.list = object.list_2, anchor.features = features)
# this command creates an 'integrated' data assay
batch_2.combined <- IntegrateData(anchorset = batch_2.anchors)
DefaultAssay(batch_2.combined) <- "integrated"
batch_2.combined <- ScaleData(batch_2.combined, verbose = FALSE)
batch_2.combined <- RunPCA(batch_2.combined, npcs = 50, verbose = FALSE)
# batch_2.combined <- RunUMAP(batch_2.combined, reduction = "pca", dims = 1:30)
batch_2.combined<- RunTSNE(batch_2.combined,reduction = "pca",dims = 1:30,cells = NULL,features = NULL,tsne.method = "Rtsne",reduction.name = "tsne")

batch_2.combined <- FindNeighbors(batch_2.combined, reduction = "pca", dims = 1:30)
batch_2.combined <- FindClusters(batch_2.combined, resolution = 1)
p3 <- DimPlot(batch_2.combined, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(batch_2.combined, reduction = "tsne", label = TRUE, repel = TRUE)
p3 + p4
DimPlot(batch_2.combined,pt.size = 1,reduction = "tsne", split.by ='orig.ident',cols =abc,label=T)
FeaturePlot(batch_2.combined, features = c("TM4SF1"),  reduction = "tsne")


object.list_2 <- lapply(X = object.list_2, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = object.list_2, nfeatures = 3000)
object.list_2 <- PrepSCTIntegration(object.list = object.list_2, anchor.features = features)
object.list_2 <- lapply(X = object.list_2, FUN = RunPCA, features = features)

batch_2.anchors <- FindIntegrationAnchors(object.list = object.list_2, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
batch_2.combined <- IntegrateData(anchorset = batch_2.anchors, normalization.method = "SCT", dims = 1:30)
batch_2.combined <- RunPCA(batch_2.combined, verbose = FALSE)
batch_2.combined <- RunTSNE(batch_2.combined, reduction = "pca", dims = 1:30)
batch_2.combined <- FindNeighbors(batch_2.combined, dims = 1:30, verbose = FALSE)
batch_2.combined <- FindClusters(batch_2.combined, resolution = 0.3,verbose = FALSE)
DimPlot(batch_2.combined, label = TRUE) + NoLegend()
 DimPlot(batch_2.combined,pt.size = 1, reduction = "FItSNE",label = T)
 DimPlot(batch_2.combined, reduction = "FItSNE",pt.size = 1, split.by = "orig.ident", label = TRUE)
 FeaturePlot(batch_2.combined, features = c("TM4SF1"),  reduction = "FItSNE",label = T)
 
 batch_2.combined <- RunTSNE(object = batch_2.combined, 
                 reduction.name = "FItSNE", 
                 reduction.key = "FItSNE_",
                 tsne.method = "FIt-SNE", 
                 fast_tsne_path = "../sc_RNA/Flt-tSNE/FItSNE.exe")

#..........................................................
object.list_3 <- lapply(X = object.list_3, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = object.list_3)
batch_3.anchors <- FindIntegrationAnchors(object.list = object.list_3, anchor.features = features)
# this command creates an 'integrated' data assay
batch_3.combined <- IntegrateData(anchorset = batch_3.anchors)
DefaultAssay(batch_3.combined) <- "integrated"
batch_3.combined <- ScaleData(batch_3.combined, verbose = FALSE)
batch_3.combined <- RunPCA(batch_3.combined, npcs = 30, verbose = FALSE)
batch_3.combined <- RunUMAP(batch_3.combined, reduction = "pca", dims = 1:30)
batch_3.combined <- FindNeighbors(batch_3.combined, reduction = "pca", dims = 1:30)
batch_3.combined <- FindClusters(batch_3.combined, resolution = 0.5)
p5 <- DimPlot(batch_3.combined, reduction = "umap", group.by = "orig.ident")
p6 <- DimPlot(batch_3.combined, reduction = "umap", label = TRUE, repel = TRUE)
p5 + p6

object.list_3 <- c(object.list[[5]],object.list[[6]])
object.list_3 <- lapply(X = object.list_3, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = object.list_3, nfeatures = 3000)
object.list_3 <- PrepSCTIntegration(object.list = object.list_3, anchor.features = features)
object.list_3 <- lapply(X = object.list_3, FUN = RunPCA, features = features)

batch_3.anchors <- FindIntegrationAnchors(object.list = object.list_3, normalization.method = "SCT",
                                          anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
batch_3.combined <- IntegrateData(anchorset = batch_3.anchors, normalization.method = "SCT", dims = 1:30)
batch_3.combined <- RunPCA(batch_3.combined, verbose = FALSE)
batch_3.combined <- RunTSNE(batch_3.combined, reduction = "pca", dims = 1:30)
batch_3.combined <- FindNeighbors(batch_3.combined, dims = 1:30, verbose = FALSE)
batch_3.combined <- FindClusters(batch_3.combined, resolution = 0.5,verbose = FALSE)
DimPlot(batch_3.combined, label = TRUE) + NoLegend()
DimPlot(batch_3.combined,pt.size = 1, reduction = "tsne",label = T)
DimPlot(batch_3.combined, reduction = "tsne", split.by = "orig.ident", label = TRUE)
FeaturePlot(batch_3.combined, features = c("TM4SF1"),  reduction = "tsne",label = T)

















#..........................................................
object.list_4 <- lapply(X = object.list_4, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = object.list_4)
batch_4.anchors <- FindIntegrationAnchors(object.list = object.list_4, anchor.features = features)
# this command creates an 'integrated' data assay
batch_4.combined <- IntegrateData(anchorset = batch_4.anchors)
DefaultAssay(batch_4.combined) <- "integrated"
batch_4.combined <- ScaleData(batch_4.combined, verbose = FALSE)
batch_4.combined <- RunPCA(batch_4.combined, npcs = 30, verbose = FALSE)
batch_4.combined <- RunUMAP(batch_4.combined, reduction = "pca", dims = 1:30)
batch_4.combined <- FindNeighbors(batch_4.combined, reduction = "pca", dims = 1:30)
batch_4.combined <- FindClusters(batch_4.combined, resolution = 0.5)
p7 <- DimPlot(batch_4.combined, reduction = "umap", group.by = "orig.ident")
p8 <- DimPlot(batch_4.combined, reduction = "umap", label = TRUE, repel = TRUE)
p7 + p8
merge_1<- DimPlot(batch_1.combined, reduction = "umap", split.by = "orig.ident")
merge_2<- DimPlot(batch_2.combined, reduction = "umap", split.by = "orig.ident")
merge_3<- DimPlot(batch_3.combined, reduction = "umap", split.by = "orig.ident")
merge_4<- DimPlot(batch_4.combined, reduction = "umap", split.by = "orig.ident")

markers_1 <- FindConservedMarkers(batch_1.combined, ident.1 = 3, grouping.var = "orig.ident", verbose = FALSE)
head(markers_1)







object.list_4 <- c(object.list[[7]],object.list[[8]])
object.list_4 <- lapply(X = object.list_4, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = object.list_4, nfeatures = 3000)
object.list_4 <- PrepSCTIntegration(object.list = object.list_4, anchor.features = features)
object.list_4 <- lapply(X = object.list_4, FUN = RunPCA, features = features)

batch_4.anchors <- FindIntegrationAnchors(object.list = object.list_4, normalization.method = "SCT",
                                          anchor.features = features, dims = 1:30, reduction = "rpca")
batch_4.combined <- IntegrateData(anchorset = batch_4.anchors, normalization.method = "SCT", dims = 1:30)
batch_4.combined <- RunPCA(batch_4.combined, verbose = FALSE)
batch_4.combined <- RunTSNE(batch_4.combined, reduction = "pca", dims = 1:30)
batch_4.combined <- FindNeighbors(batch_4.combined, dims = 1:30, verbose = FALSE)
batch_4.combined <- FindClusters(batch_4.combined, resolution = 0.5,verbose = FALSE)
DimPlot(batch_4.combined, label = TRUE) + NoLegend()
DimPlot(batch_4.combined,pt.size = 1, reduction = "tsne")
DimPlot(batch_4.combined, reduction = "tsne", split.by = "orig.ident", pt.size = 1,label = TRUE,cols = values)
FeaturePlot(batch_4.combined)

