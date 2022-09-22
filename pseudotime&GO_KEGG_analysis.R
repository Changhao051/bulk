#monocle & GO analysis

library(dplyr)
library(Seurat)
library(patchwork) 
library(monocle)
setwd("~/Desktop/~")
data<-scrna_harmony[,scrna_harmony$cell_type%in%c("cell_type")]
HL<-scrna_harmony[,scrna_harmony$orig.ident%in%c("SXT","HLQ","HLH","DZT")]
data<-HL[,HL$cell_type%in%c("Epidermal/muscle")]
Idents(data)<-"RNA_snn_res.1.4"
expr_matrix <- as(as.matrix(data@assays$RNA@counts), 'sparseMatrix')
p_data <-data@meta.data 
f_data <- data.frame(gene_short_name = row.names(data),row.names = row.names(data))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())  
f_data <- data.frame(gene_short_name = row.names(data),row.names = row.names(data))
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1) 
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) 
deg.cluster <- FindAllMarkers(data)
express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds) 
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds) # 这里是可视化
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 4)

p1=plot_cell_trajectory(cds, color_by = "RNA_snn_res.1.4",cell_size = 0.1)
p2=plot_cell_trajectory(cds, color_by = "orig.ident",cell_size = 0.1)
p3=plot_cell_trajectory(cds, color_by = "State",cell_size = 0.1)
p4=plot_cell_trajectory(cds, color_by = "cell_type",cell_size = 0.1)
p1+p2+p3+p4
p5 = plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                  color_by = "RNA_snn_res.1.4",
                                  cell_size = 0.3)
p6 = plot_complex_cell_trajectory(cds, x = 1, y = 2,
                                  color_by = "orig.ident",
                                  cell_size = 0.3)

p5+p6
BEAM_res <- BEAM(cds, branch_point = 3, cores = 8)
BEAM_res1 <- BEAM_res[order(BEAM_res$qval),]
BEAM_deg<- cds[row.names(subset(BEAM_res,
                                qval < 1e-4)),]
BEAM_res1 <- BEAM_res1[,c("gene_short_name", "pval", "qval")]
head(BEAM_res1)
tmp1 = plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
                                                        qval < 1e-4)),],
                                   branch_point = 2,
                                   num_clusters = 4,
                                   cores = 1,
                                   use_gene_short_name = T,
                                   hmcols = colorRampPalette(rev(brewer.pal(9,"PRGn")))(62),
                                   show_rownames = T,
                                   return_heatmap = T
)
tmp1$ph_res




BEAM_res<-top_n(BEAM_res,n = 100)


##GO and KEGG_analysis

BEAM_genes <- top_n(BEAM_res1,n = 100, dplyr::desc(qval)) %>% pull(gene_short_name) %>% as.character()
write.csv(BEAM_genes,file = "HL_data_branch3_deg_top100.csv")

tmp1 = plot_genes_branched_heatmap(cds[BEAM_genes,],  branch_point = 3, 
                                   num_clusters = 3,
                                   cores = 1,
                                   hmcols = colorRampPalette(rev(brewer.pal(9,"PRGn")))(62),
                                   show_rownames = T,return_heatmap = T
)



tmp1$ph_res

gene_group=tmp1$annotation_row
gene_group$gene=rownames(gene_group)
X= gene_group
library(clusterProfiler)

info <- read.csv("~/Desktop/GO_dataset.csv")
info = info[,-1] 
cluster<-split(X,X$Cluster)
TERM2GENE <- data.frame(info[,2], info[,1])
TERM2NAME <- data.frame(info[,2], info[,3])
TERM2NAME <- unique(TERM2NAME[order(TERM2NAME[,1]),])
reslut<-lapply(cluster, function(x){
  enrich <- enricher(x$gene,
                     TERM2GENE=TERM2GENE,
                     TERM2NAME=TERM2NAME,
                     pvalueCutoff=10,
                     qvalueCutoff=10,
                     pAdjustMethod = "fdr"
                     
  )
  return(enrich)
  # arplarplot(enrich)%>%ggsave(filename = paste(unique(x$cluster),'barplot.eps'))
  # dotplot(enrich)%>%ggsave(filename = paste(unique(x$cluster),'dotplot.eps'))
  
})

barplot(reslut$`32`,title = "XXX",font.size = 8,showCategory=20)


p0 <- barplot(reslut$`1`,title = "cluster1",font.size = 8,showCategory=20)
p1 <- barplot(reslut$`2`,title = "cluster2",font.size = 8)
p2 <- barplot(reslut$`3`,title = "cluster3",font.size = 8)
#p3 <- barplot(reslut$DZT,title = "DZT",font.size = 8)
p0+p1+p2
