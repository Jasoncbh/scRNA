library(ggsci)
rm(pbmc_1)
pbmc_1 <- merge(pbmc_AE_1,pbmc_Mig_1)
table(pbmc_1$orig.ident)

# store mitochondrial percentage in object meta data
pbmc_1<- PercentageFeatureSet(pbmc_1, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
pbmc_1 <- SCTransform(pbmc_1, vars.to.regress =c("percent.mt","orig.ident"), verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc_1 <- RunPCA(pbmc_1, verbose = FALSE)
ElbowPlot(pbmc_1, ndims = 30)
dims = 16
pbmc_1 <- RunTSNE(object = pbmc_1, 
                  dims = 1:dims, 
                  reduction.name = "FItSNE", 
                  reduction.key = "FItSNE_",
                  tsne.method = "FIt-SNE", 
                  fast_tsne_path = "./Flt-tSNE/FItSNE.exe")
# pbmc_1 <- RunTSNE(pbmc_1, dims = 1:dims, verbose = FALSE)
pbmc_1 <- FindNeighbors(pbmc_1, dims = 1:dims, verbose = FALSE)
pbmc_1 <- FindClusters(pbmc_1, verbose = FALSE,resolution = 0.5)
a <- DimPlot(pbmc_1, label = F,reduction = "FItSNE",pt.size =1.5,
             cols =values,
             label.size = 8)
a+theme(axis.line = element_line(size = 1.15,color = "black"), 
        axis.ticks.y = element_line(size = 1.2,),
        axis.ticks.length=unit(0.5,"lines"),axis.text.y = element_text(size = 15),
        axis.ticks.x =element_line(size = 1.2,),
        legend.position = "right")+
  ggtitle("流式分选后第0天")+  
  theme(plot.title = element_text(size=20,hjust = 0),legend.text=element_text(size=25))
b<- FeaturePlot(pbmc_1,features = "TM4SF1",pt.size = 2,label = T,ncol = 1,cols = c("grey", "red"),label.size = 8)
b+theme(axis.line = element_line(size = 1.15,color = "black"), 
        axis.ticks.y = element_line(size = 1.2,),
        axis.ticks.length=unit(0.5,"lines"),axis.text.y = element_text(size = 15),
        axis.ticks.x =element_line(size = 1.2,),
        legend.position = "none", )+
  ggtitle("流式分选后第0天")+  
  theme(plot.title = element_text(size=20,hjust = 0))
c <- DimPlot(pbmc_1,label = F,reduction = "FItSNE",pt.size = 1.25,split.by = "orig.ident",ncol = 1,
             cols =values,
             label.size = 8)
c+theme(axis.line = element_line(size = 1.15,color = "black"), 
       axis.ticks.y = element_line(size = 1.2,),
       axis.ticks.length=unit(0.5,"lines"),axis.text.y = element_text(size = 15),
       axis.ticks.x =element_line(size = 1.2,),
       legend.position = "right" )+
  ggtitle("流式分选后第0天")+  
  theme(plot.title = element_text(size=20,hjust = 0),legend.text=element_text(size=25))
pbmc1.markers <- FindAllMarkers(pbmc_1, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(pbmc1.markers,file = "pbmc1.marker.csv")
pbmc_1@active.ident <- pbmc_1$seurat_clusters
pbmc_1 <- RenameIdents(pbmc_1, 
                     `0` = "Neutrophil progenitor", 
                     `1` = "MEP", 
                     `2` = "Mast cell", 
                     `3` = "MEP", 
                     `4` = "Monocyte", 
                     `5` = "Mono/DC progenitors", 
                     `6` = "Neutrophil progenitor", 
                     `7` = "HSC-like", 
                     `8` = "Erythroid", 
                     `9` = "Natural killer T cell", 
                     `10` = "Basophil progenitor")


# as.data.frame(prop.table(table(Idents(pbmc_1),pbmc_1$type)))
as.data.frame(table(Idents(pbmc_1), pbmc_1$orig.ident))
cell1<-as.data.frame(table(Idents(pbmc_1), pbmc_1$type))
colnames(cell1)<-c("clu","type","value")

cell1$clu = factor(cell1$clu, levels=c('Mono/DC progenitors','Monocyte','Erythroid','MEP','Basophil progenitor','Mast cell','Natural killer T cell','Neutrophil progenitor','HSC-like')) ## 设置柱条的顺序
ggplot(data = cell1)+
  geom_bar(mapping = aes(x=clu,y=value,fill=type),position = "fill",stat = "identity")+
  theme(
    plot.margin = unit(c(1,1,1,1), "cm"),##########调整边界大小#############
    panel.background = element_blank(),
    plot.title = element_text(size = 20, face = "bold",hjust = 0.5, margin = margin(b = 15)),
    axis.line = element_line(size = 1.15,color = "black"),
    axis.title= element_blank(),
    axis.ticks.y = element_line(size = 1.2,),
    axis.ticks.length=unit(0.5,"lines"),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 18,colour = col,angle = 90,vjust = 0.5,hjust = 1),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key=element_rect(color=NA,fill=NA))+
  ggtitle("流式分选后第0天")+  
  scale_fill_npg()
col <- c("black","black","black","black","black","black","black","black","red")





class(id$value)
class(diamonds$price)
abc <- diamonds[,2:7]
ggplot(id,aes(x=clu,fill=type))+
  geom_bar(position = "fill")

class(abc$clu)
cell1$clu <- as.factor(cell1$clu)
cell1$type <- as.factor(cell1$type)
class(cell1$clu)
str(abc)
#####################################################################################################################################

pbmc_2 <- merge(pbmc_AE_2,pbmc_Mig_2)
pbmc_2 <- PercentageFeatureSet(pbmc_2, pattern = "^MT-", col.name = "percent.mt")

pbmc_2 <- SCTransform(pbmc_2,vars.to.regress = c("percent.mt","orig.ident"), verbose = FALSE)

# These are now standard steps in the Seurat workflow for visualization and clustering

pbmc_2 <- RunPCA(pbmc_2, verbose = FALSE)
ElbowPlot(pbmc_2, ndims = 30)
dims = 20
# pbmc_2 <- RunTSNE(pbmc_2, dims = 1:dims, verbose = FALSE)
pbmc_2 <- RunTSNE(object = pbmc_2, 
                  dims = 1:dims, 
                  reduction.name = "FItSNE", 
                  reduction.key = "FItSNE_",
                  tsne.method = "FIt-SNE", 
                  fast_tsne_path = "./Flt-tSNE/FItSNE.exe")
pbmc_2 <- FindNeighbors(pbmc_2, dims = 1:20, verbose = FALSE)
pbmc_2 <- FindClusters(pbmc_2, verbose = FALSE,resolution = 0.3)
a <- DimPlot(pbmc_2, label = F,reduction = "FItSNE",pt.size = 1.5,
             cols =c("#377EB8","#4DAF4A", "#984EA3" ,"#F781BF" ,"#E41A1C","#FF7F00" ,"#FFFF33" ,"#A65628", "#999999","#1B9E77" ,"#D95F02"  ),
             label.size = 6.5) 
a+theme(axis.line = element_line(size = 1.15,color = "black"), 
        axis.ticks.y = element_line(size = 1.2,),
        axis.ticks.length=unit(0.5,"lines"),axis.text.y = element_text(size = 15),
        axis.ticks.x =element_line(size = 1.2,),
        legend.position = "right", )+
  ggtitle("流式分选后第1天")+  
  theme(plot.title = element_text(size=20,hjust = 0),legend.text=element_text(size=25))
b <- FeaturePlot(pbmc_2,features = "TM4SF1",pt.size = 2,cols = c("grey", "red"),label = T,label.size = 7)
b+theme(axis.line = element_line(size = 1.15,color = "black"), 
        axis.ticks.y = element_line(size = 1.2,),
        axis.ticks.length=unit(0.5,"lines"),axis.text.y = element_text(size = 15),
        axis.ticks.x =element_line(size = 1.2,),
        legend.position = "none", )+
  ggtitle("流式分选后第1天")+  
  theme(plot.title = element_text(size=20,hjust = 0))
c <- DimPlot(pbmc_2,label = F,reduction = "FItSNE",pt.size = 1.25,split.by = "orig.ident",ncol = 1,
             label.size = 7.5,
             cols=c("#377EB8","#4DAF4A", "#984EA3" ,"#F781BF" ,"#E41A1C","#FF7F00" ,"#FFFF33" ,"#A65628", "#999999","#1B9E77" ,"#D95F02"  ))
c+theme(axis.line = element_line(size = 1.15,color = "black"), 
       axis.ticks.y = element_line(size = 1.2,),
       axis.ticks.length=unit(0.5,"lines"),axis.text.y = element_text(size = 15),
       axis.ticks.x =element_line(size = 1.2,),
       legend.position = "right" )+
  ggtitle("流式分选后第1天")+  
  theme(plot.title = element_text(size=20,hjust = 0),legend.text=element_text(size=25))
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc2.markers <- FindAllMarkers(pbmc_2, only.pos = TRUE)
pbmc2.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(pbmc2.markers,file = "pbmc_2_2marker.csv")
pbmc_2@active.ident <- pbmc_2$seurat_clusters
pbmc_2 <- RenameIdents(pbmc_2, 
                       `0` = "Natural killer T cell", 
                       `1` = "Mast cell", 
                       `2` = "Mk/Ery progenitor", 
                       `3` = "MEP", 
                       `4` = "HSC-like", 
                       `5` = "Progenitor", 
                       `6` = "Granulocyte-monocyte progenitor", 
                       `7` = "Neutrophil progenitor", 
                       `8` = "Megakaryocyte progenitor cell", 
                       `9` = "pDC/NK", 
                       `10` = "Megakaryocyte progenitor cell")

cell.prop2<-as.data.frame(table(Idents(pbmc_2), pbmc_2$type))
colnames(cell.prop2)<-c("clu","type","value")
# cell.prop2$clu = factor(cell.prop2$clu, levels=c('Megakaryocyte progenitor cell','Neutrophil progenitor','Erythroid','MEP','Basophil progenitor','Mast cell','Natural killer T cell','Neutrophil progenitor','HSC-like')) ## 设置柱条的顺序

ggplot(data = cell.prop2)+
  geom_bar(mapping = aes(x=clu,y=value,fill=type),position = "fill",stat = "identity")+
  theme(
    plot.margin = unit(c(1,1,1,1), "cm"),##########调整边界大小#############
    panel.background = element_blank(),
    plot.title = element_text(size = 20, face = "bold",hjust = 0.5, margin = margin(b = 15)),
    axis.line = element_line(size = 1.15,color = "black"),
    axis.title= element_blank(),
    axis.ticks.y = element_line(size = 1.2,),
    axis.ticks.length=unit(0.5,"lines"),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 18,colour = col,angle = 90,vjust = 0.5,hjust = 1),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key=element_rect(color=NA,fill=NA))+
  ggtitle("流式分选后第1天")+  
  scale_fill_npg()











saveRDS(pbmc_3,file = "batch3.rds")
# rm(pbmc_3)
pbmc_3 <- merge(pbmc_AE_3,pbmc_Mig_3)
pbmc_3 <- PercentageFeatureSet(pbmc_3, pattern = "^MT-", col.name = "percent.mt")

pbmc_3 <- SCTransform(pbmc_3,vars.to.regress = c("percent.mt","orig.ident"), verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc_3 <- RunPCA(pbmc_3, verbose = FALSE)
ElbowPlot(pbmc_3, ndims = 25)
dims =15
pbmc_3 <- RunTSNE(object = pbmc_3, 
                dims = 1:dims, 
                reduction.name = "FItSNE", 
                reduction.key = "FItSNE_",
                tsne.method = "FIt-SNE", 
                fast_tsne_path = "./Flt-tSNE/FItSNE.exe")
# pbmc_3 <- RunTSNE(pbmc_3, dims = 1:dims, verbose = FALSE)
pbmc_3 <- FindNeighbors(pbmc_3, dims = 1:15, verbose = FALSE)
pbmc_3 <- FindClusters(pbmc_3, verbose = FALSE,resolution = 0.5)
a <- DimPlot(pbmc_3, label = TRUE,reduction = "FItSNE",pt.size = 1.25,
             label.size = 6,
             cols = c("#377EB8","#4DAF4A", "#984EA3" ,"#F781BF" ,"#FF7F00" ,"#FFFF33" ,"#E41A1C","#A65628", "#999999","#1B9E77" ,"#D95F02" ,"#7570B3" ))
a+theme(axis.line = element_line(size = 1.15,color = "black"), 
        axis.ticks.y = element_line(size = 1.2,),
        axis.ticks.length=unit(0.5,"lines"),axis.text.y = element_text(size = 15),
        axis.ticks.x =element_line(size = 1.2,),
        legend.position = "right", )+
  ggtitle("流式分选后第3天")+  
  theme(plot.title = element_text(size=20,hjust = 0),legend.text=element_text(size=15))


b <- FeaturePlot(pbmc_3,features = "TM4SF1",pt.size = 2,label = T,label.size = 7,cols = c("grey", "red"),max.cutoff = 2,ncol = 1)
b+theme(axis.line = element_line(size = 1.15,color = "black"), 
        axis.ticks.y = element_line(size = 1.2,),
        axis.ticks.length=unit(0.5,"lines"),axis.text.y = element_text(size = 15),
        axis.ticks.x =element_line(size = 1.2,),
        legend.position = "none", )+
  ggtitle("流式分选后第3天")+  
  theme(plot.title = element_text(size=20,hjust = 0))
c <- DimPlot(pbmc_3,label = F,reduction = "FItSNE",pt.size = 1.25,split.by = "orig.ident",ncol = 1,
        label.size = 6,
             cols = c("#377EB8","#4DAF4A", "#984EA3" ,"#F781BF" ,"#FF7F00" ,"#FFFF33" ,"#E41A1C","#A65628", "#999999","#1B9E77" ,"#D95F02" ,"#7570B3" ))
c+theme(axis.line = element_line(size = 1.15,color = "black"), 
        axis.ticks.y = element_line(size = 1.2,),
        axis.ticks.length=unit(0.5,"lines"),axis.text.y = element_text(size = 15),
        axis.ticks.x =element_line(size = 1.2,),
        legend.position = "right" )+
  ggtitle("流式分选后第3天")+  
  theme(plot.title = element_text(size=20,hjust = 0),legend.text=element_text(size=20),axis.title = element_text(size = 20))

FindMarkers(pbmc_3, ident.1 = 0,min.pct = 0.25)
pbmc3.markers <- FindAllMarkers(pbmc_3, only.pos = TRUE)
pbmc3.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(pbmc3.markers,file = "pbmc_3marker.csv")
pbmc_3@active.ident <- pbmc_3$seurat_clusters
pbmc_3 <- RenameIdents(pbmc_3, 
                       `0` = "Granulocyte-monocyte progenitor", 
                       `1` = "Natural killer T (NKT) cell", 
                       `2` = "Cells in G2/M phase", 
                       `3` = "Monocyte", 
                       `4` = "Mast cell", 
                       `5` = "Progenitor", 
                       `6` = "HSC-like", 
                       `7` = "Eosinophil progenitor", 
                       `8` = "Granulocyte-monocyte progenitor", 
                       `9` = "Erythroid", 
                       `10` = "Neutrophil progenitor", 
                       `11` = "Natural killer T (NKT) cell")



cell.prop3<-as.data.frame(prop.table(table(Idents(pbmc_3), pbmc_3$type)))
colnames(cell.prop3)<-c("cluster","type","proportion")
ggplot(cell.prop3,aes(type,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  scale_fill_manual(values = c("#377EB8","#4DAF4A", "#984EA3" ,"#F781BF" ,"#FF7F00" ,"#FFFF33" ,"#E41A1C","#A65628", "#999999","#1B9E77" ,"#D95F02" ,"#7570B3" ))+
  theme(axis.ticks.length=unit(0.5,'cm'))+
  guides(fill=guide_legend(title=NULL))









library(ggplot2)
library(Seurat)
dims=20
rm(pbmc_4)
pbmc_4 <- merge(pbmc_AE_4,pbmc_Mig_4)
# pbmc_4[["percent.mt"]] <- PercentageFeatureSet(pbmc_4, pattern = "^MT-")
# 
# pbmc_4 <- subset(pbmc_4, subset = nFeature_RNA > 350 & percent.mt < 7)
pbmc_4 <- PercentageFeatureSet(pbmc_4, pattern = "^MT-", col.name = "percent.mt")
table(pbmc_4$type)
pbmc_4 <- SCTransform(pbmc_4,vars.to.regress = c("percent.mt","tech"))
pbmc_4 <- RunPCA(pbmc_4,features =VariableFeatures(object = pbmc_4)[1:500], npcs = 50)
ElbowPlot(pbmc_4, ndims = 25)
pbmc_4 <- RunTSNE(object = pbmc_4, 
                dims = 1:dims, 
                reduction.name = "FItSNE", 
                reduction.key = "FItSNE_",
                tsne.method = "FIt-SNE", 
                fast_tsne_path = "./Flt-tSNE/FItSNE.exe")
pbmc_4 <- FindNeighbors(pbmc_4, dims = 1:dims, verbose = FALSE)
pbmc_4 <- FindClusters(pbmc_4,resolution = 0.6)
a <- DimPlot(pbmc_4, label = F,reduction = "FItSNE",pt.size = 2,label.size = 6.5,
             cols =c("#377EB8","#4DAF4A", "#984EA3" ,"#F781BF" ,"#FF7F00" ,"#FFFF33" ,"#E41A1C","#A65628", "#999999","#1B9E77" ,"#D95F02"  )) 
a+theme(axis.line = element_line(size = 1.15,color = "black"), 
        axis.ticks.y = element_line(size = 1.2,),
        axis.ticks.length=unit(0.5,"lines"),axis.text.y = element_text(size = 15),
        axis.ticks.x =element_line(size = 1.2,),
        legend.position = "right", )+
  ggtitle("流式分选后第7天")+  
  theme(plot.title = element_text(size=20,hjust = 0),legend.text=element_text(size=25))

b <- FeaturePlot(pbmc_4,features = "TM4SF1",pt.size =2,label = T,cols = c("grey", "red"),label.size = 7)
b+theme(axis.line = element_line(size = 1.2,color = "black"), 
       axis.ticks.y = element_line(size = 1.2,),
       axis.ticks.length=unit(0.5,"lines"),axis.text.y = element_text(size = 15),
       axis.ticks.x =element_line(size = 1.2,),
       legend.position = "none", )+
  ggtitle("流式分选后第7天")+  
  theme(plot.title = element_text(size=20,hjust = 0))
c <- DimPlot(pbmc_4,label = F,reduction = "FItSNE",pt.size = 1.25,split.by = "orig.ident",ncol = 1,
        cols =c("#377EB8","#4DAF4A", "#984EA3" ,"#F781BF" ,"#FF7F00" ,"#FFFF33" ,"#E41A1C","#A65628", "#999999","#1B9E77" ,"#D95F02"  ))

c+theme(axis.line = element_line(size = 1.15,color = "black"), 
        axis.ticks.y = element_line(size = 1.2,),
        axis.ticks.length=unit(0.5,"lines"),axis.text.y = element_text(size = 15),
        axis.ticks.x =element_line(size = 1.2,),
        legend.position = "right" )+
  ggtitle("流式分选后第7天")+  
  theme(plot.title = element_text(size=20,hjust = 0),legend.text=element_text(size=25))

FindMarkers(pbmc_4, ident.1 = 3)
pbmc4.markers <- FindAllMarkers(pbmc_4, only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.1)
pbmc4.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(pbmc4.markers,file = "pbmc_4marker_4.csv")
pbmc_4@active.ident <- pbmc_4$seurat_clusters
pbmc_4 <- RenameIdents(pbmc_4, 
                       `0` = "CMP (undetermined)", 
                       `1` = "Mast cell progenitor", 
                       `2` = "Natural killer T cell", 
                       `3` = "Mk/Ery progenitor", 
                       `4` = "Progenitor cell", 
                       `5` = "Neutrophil progenitor", 
                       `6` = "HSC-like", 
                       `7` = "Progenitor cell", 
                       `8` = "Natural killer T cell", 
                       `9` = "Monocyte", 
                       `10` = "Eosinophil progenitor",
                       `11` = "Monocyte",
                       `12` = "Neutrophil progenitor",
                       `13` = "Microglial cell",
                       `14` = "Cells in G2/M phase",
                       `15` = "Natural killer T cell")

cell4$clu = factor(cell4$clu, levels=c('Neutrophil progenitor','Monocyte',"Mk/Ery progenitor",'Mast cell',"Cells in G2/M phase",'Natural killer T cell',"Microglial cell","Eosinophil progenitor", "Progenitor cell",'HSC-like')) ## 设置柱条的顺序

cell4<-as.data.frame(table(Idents(pbmc_4), pbmc_4$type))
colnames(cell4)<-c("clu","type","value")
# cell.prop2$clu = factor(cell.prop2$clu, levels=c('Megakaryocyte progenitor cell','Neutrophil progenitor','Erythroid','MEP','Basophil progenitor','Mast cell','Natural killer T cell','Neutrophil progenitor','HSC-like')) ## 设置柱条的顺序
col <- c("black","black","black","black","black","black","black","black","black","red")
ggplot(data = cell4)+
  geom_bar(mapping = aes(x=clu,y=value,fill=type),position = "fill",stat = "identity")+
  theme(legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key=element_rect(color=NA,fill=NA),
    plot.margin = unit(c(1,1,1,1), "cm"),##########调整边界大小#############
    panel.background = element_blank(),
    plot.title = element_text(size = 20, face = "bold",hjust = 0.5, margin = margin(b = 15)),
    axis.line = element_line(size = 1.15,color = "black"),
    axis.title= element_blank(),
    axis.ticks.y = element_line(size = 1.2,),
    axis.ticks.length=unit(0.5,"lines"),
    axis.text.x = element_text(size = 18,color = col,angle = 90,vjust = 0.5,hjust = 1),
    axis.ticks.x = element_blank())+
  scale_fill_npg()+
ggtitle("流式分选后第7天")
  markers.to.plot <- c( "TPSAB1","TPSB2", "CPA3", "CCL2",  "GZMB",
                       "H2AFZ", "HMGB2","RRM2", "MALAT1","SOX4", "RAB33A", "CTSG", "CXCL8", "CMA1", "TM4SF1", "CRHBP", 
                       "GNG11","SOD2","CHI3L1", "LGMN", "ISG15", "IFI6", "MX1","C1QB","C1QA","CDC20","CCNB1","PLK1")
  DotPlot(pbmc_4, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
    RotatedAxis()+  theme(
      plot.margin = unit(c(1,1,1,1), "cm"),##########调整边界大小#############
      panel.background = element_blank(),
      plot.title = element_text(size = 20, face = "bold",hjust = 0.5, margin = margin(b = 15)),
      axis.line = element_line(size = 1.15,color = "black"),
      axis.title= element_blank(),
      axis.ticks.y = element_line(size = 1.2,),
      axis.ticks.length=unit(0.5,"lines"),
      axis.text.y = element_text(size = 15),
      axis.text.x = element_text(size = 18,colour = col,angle = 90,vjust = 0.5,hjust = 1),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 12),
      legend.key=element_rect(color=NA,fill=NA))+
    ggtitle("流式分选后第7天")


saveRDS(pbmc_AE_4,file = "ae4.rds")
# rm(pbmc_4)
saveRDS(pbmc_AE_3,file = "AE3.rds")
saveRDS()
##############################################################################
#..................................seurat_sctransform_爆内存
pbmc_merge <- merge(pbmc_1,y=c(pbmc_2,pbmc_3,pbmc_4),add.cell.ids = c("1st", "2nd","3rd","4th"),project = "AE_MIG_merge")
pbmc_merge.list <- SplitObject(pbmc_merge, split.by = "batch")
table(pbmc_merge$batch)
saveRDS(pbmc_merge.list,file = "pbmc_merge.list.rds")
pbmc_merge.list <- readRDS("object_list.rds")
pbmc_merge.list <- lapply(X = pbmc_merge.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = pbmc_merge.list, nfeatures = 3000)
pbmc_merge.list <- PrepSCTIntegration(object.list = pbmc_merge.list, anchor.features = features)
pbmc_merge.list <- lapply(X = pbmc_merge.list, FUN = RunPCA, features = features)

immune.anchors <- FindIntegrationAnchors(object.list = pbmc_merge.list, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)



###############################################################################

