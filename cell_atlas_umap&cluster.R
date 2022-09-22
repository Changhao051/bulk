library(harmony)
library(Seurat)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
setwd("~/Desktop/matrix/")
FLYT<-read.table('FLYT-filter-matrix.csv',sep = ",", check.names = F, comment.char = "#", header = TRUE, row.names = 1)
SXT<-read.table('SXT-filter-matrix.csv',sep = ",", check.names = F, comment.char = "#", header = TRUE, row.names = 1)
HLQ<-read.table('HLQ-filter-matrix.csv',sep = ",", check.names = F, comment.char = "#", header = TRUE, row.names = 1)
HLH<-read.table('HLH-filter-matrix.csv',sep = ",", check.names = F, comment.char = "#", header = TRUE, row.names = 1)
DZT<-read.table('DZT-filter-matrix.csv',sep = ",", check.names = F, comment.char = "#", header = TRUE, row.names = 1)
SMT<-read.table('SMT-filter-matrix.csv',sep = ",", check.names = F, comment.char = "#", header = TRUE, row.names = 1)
head(SXT)
data<-list(FLYT,SXT,HLQ,HLH,DZT,SMT)
name<-c("FLYT","SXT","HLQ","HLH","DZT","SMT")
scrna.list<-list()
for (i in name) {
  scrna.list[[i]] <- data[[i]]
  scrna.list[[i]]<-CreateSeuratObject(scrna.list[[i]],project = i,min.cells = 3,       # 1% of total cells
                                      min.features = 200)#过滤细胞
}
exp <- merge(scrna.list[[1]],y=c(scrna.list[[2]],scrna.list[[3]],scrna.list[[4]],scrna.list[[5]],scrna.list[[6]]))

scrna <- exp
##==harmony
scrna_harmony <- NormalizeData(scrna) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
system.time({scrna_harmony <- RunHarmony(scrna_harmony, group.by.vars = "orig.ident")})

scrna_harmony <- RunUMAP(scrna_harmony, reduction = "harmony", dims = 1:25)
scrna_harmony <- FindNeighbors(scrna_harmony, reduction = "harmony", dims = 1:25)%>% FindClusters(resolution=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0))
##clustree
library(clustree)
clustree(scrna_harmony@meta.data, prefix = "RNA_snn_res.")
#umap
DimPlot(scrna_harmony, reduction = "umap",group.by = 'RNA_snn_res.1.4',label = T)

##find marker
Idents(scrna_harmony)<-"orig.ident"
marker <- FindAllMarkers(scrna_harmony, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write(marker,file = "marker.csv")


# marker dotplot
seurat.object_copy <-scrna_harmony
Idents(seurat.object_copy) <- "RNA_snn_res.1.4"
table(seurat.object_copy$RNA_snn_res.1.4)
my_levels <- c("7","9","13","15","17","18","20","21","23","32","39",
               "1","6","10","22","28","31","35","37","40",
               "24",
               "3","4","5","8","19",
               "2","11","16","27","33","34","38",
               "0","12","29","30","41",
               "14","25","26","36")
Idents(seurat.object_copy) <- factor(Idents(seurat.object_copy), levels= my_levels)
marker<-c("Gene_ID")
col4<-colorRampPalette(brewer.pal(6,"Blues"))(20)

DotPlot(seurat.object_copy,features = marker2,dot.scale = 9)+coord_flip()+theme_bw()+
    theme(panel.grid = element_blank(),axis.title.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
    scale_color_gradientn(values = seq(0,1,0.2),colours = col4)+
    labs(x= NULL,y= NULL)+guides(size = guide_legend(order = 3))


# heatmap
  subobj <- subset(seurat.object_copy, downsample = 100)
  Idents(subobj) <- factor(Idents(subobj), levels= my_levels)
  library(Seurat)
  top_n = Br.markers %>% 
    filter(p_val_adj  <= 0.05) %>%
    filter(pct.1  >= 0.5) %>%
    group_by(cluster) %>% 
    top_n(n = x, wt = avg_log2FC)
  
  DoHeatmap(subobj,features = Ne_marekr$gene,hjust = 0.5,size = 4,slot = "scale.data"
  )+scale_fill_gradientn(colors = c("white", "gray", "blue"))
  
# subcluster analysis
Idents(scrna_harmony)<-"cell_type"
cell_subset<-scrna_harmony[,scrna_harmony$cell_type%in%c("cellsubset")]
cell_subset <- NormalizeData(cell_subset, normalization.method = "LogNormalize", scale.factor = 1e4) 
cell_subset <- FindVariableFeatures(cell_subset, selection.method = 'vst', nfeatures = 2000)
cell_subset <- ScaleData(cell_subset)
cell_subset <- RunPCA(cell_subset, features = VariableFeatures(object = cell_subset)) 
cell_subset <- FindNeighbors(cell_subset, dims = 1:30)
cell_subset <- FindClusters(cell_subset, resolution = 0.7 )
cell_subset <- RunUMAP(cell_subset, dims = 1:30)
DimPlot(cell_subset, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()

##Histograms of cell ration

library(ggplot2)
Idents(scrna_harmony)<-"cell_type"
table(Idents(scrna_harmony))
prop.table(table(Idents(scrna_harmony)))
table(Idents(scrna_harmony),scrna_harmony$orig.ident)
cellratio<- prop.table(table(Idents(scrna_harmony),scrna_harmony$orig.ident),margin = 2)
cellratio<- as.data.frame(cellratio)
mycolor<-c("#F5A118","#B2D9AD","#3D83C5","#2EB066","#674699","#114B84", "red")
cellratio$Var2 <- factor(cellratio$Var2, c("life cycle number"))
ggplot(cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.5,size = 0.5)+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = mycolor)+coord_polar(theta = 'y', direction = 1)+theme(axis.text = element_blank()) 

#donut chart of cell ratio
library(lessR)

celltype_D<-FetchData(sample,"cell_type")
PieChart(cell_type,data = celltype_D,hole =0.9,fill = mycolor,labels_cex = 0.001,main = "sample1",values = 'prop',values_size = 0.8)
celltype_ES<-FetchData(sample,"cell_type")
PieChart(cell_type,data = celltype_ES,hole =0.9,fill = mycolor,labels_cex = 0.001,main = "sample2",main_cex = 1,values = 'prop',values_size = 0.8)
celltype_AS<-FetchData(sample,"cell_type")
PieChart(cell_type,data = celltype_AS,hole =0.9,fill = mycolor,labels_cex = 0.001,main = "sample3",main_cex = 1,values = 'prop',values_size = 0.8)



















