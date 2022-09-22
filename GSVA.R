
Idents(scrna_harmony) <- "orig.ident"
data<-scrna_harmony[,scrna_harmony$orig.ident%in%c("X1","X2")]
#data<-scrna_harmony
Idents(data) <- "orig.ident"
Idents(data) <- "RNA_snn_res.1.4"
exp<- AverageExpression(data,assays = "RNA",slot = "data")$RNA
go <- read.csv('~/Desktop/GO_dataset.csv')
ref<-data.frame(des=go$pathway_description,genes=go$gene)
ref<-split(ref$genes,f = ref$des)
res <- GSVA::gsva(expr = exp,ref)
view(res)
var <- apply(res, 1,var) %>% as.data.frame()

pheatmap::pheatmap(res[order(var$.,decreasing = T)[1:100],],cluster_cols = F,cluster_rows = T)
res<-data.frame(res)
tma <- c(order(res$X0,decreasing = T)[1:10],
         order(res$X1,decreasing = T)[1:10],
         order(res$X2,decreasing = T)[1:10],
         order(res$X3,decreasing = T)[1:10],
         order(res$X4,decreasing = T)[1:10],
         order(res$X1,decreasing = T)[1:10],
         order(res$X6,decreasing = T)[1:10],
         order(res$X7,decreasing = T)[1:10],
         order(res$X8,decreasing = T)[1:10],
         order(res$X9,decreasing = T)[1:10],
         order(res$X10,decreasing = T)[1:10],
         order(res$X11,decreasing = T)[1:10],
         order(res$X12,decreasing = T)[1:10],
         order(res$X13,decreasing = T)[1:10],
         order(res$X14,decreasing = T)[1:10],
         order(res$X11,decreasing = T)[1:10],
         order(res$X16,decreasing = T)[1:10],
         order(res$X17,decreasing = T)[1:10],
         order(res$X18,decreasing = T)[1:10],
         order(res$X19,decreasing = T)[1:10],
         order(res$X20,decreasing = T)[1:10],
         order(res$X21,decreasing = T)[1:10],
         order(res$X22,decreasing = T)[1:10],
         order(res$X23,decreasing = T)[1:10],
         order(res$X24,decreasing = T)[1:10],
         order(res$X21,decreasing = T)[1:10],
         order(res$X26,decreasing = T)[1:10],
         order(res$X27,decreasing = T)[1:10],
         order(res$X28,decreasing = T)[1:10],
         order(res$X29,decreasing = T)[1:10],
         order(res$X30,decreasing = T)[1:10],
         order(res$X31,decreasing = T)[1:10],
         order(res$X32,decreasing = T)[1:10],
         order(res$X33,decreasing = T)[1:10],
         order(res$X34,decreasing = T)[1:10],
         order(res$X31,decreasing = T)[1:10],
         order(res$X36,decreasing = T)[1:10],
         order(res$X37,decreasing = T)[1:10],
         order(res$X38,decreasing = T)[1:10],
         order(res$X39,decreasing = T)[1:10],
         order(res$X40,decreasing = T)[1:10],
         order(res$X41,decreasing = T)[1:10])

tma <- c(order(res$FLYT,decreasing = T)[1:20],
         order(res$SXT,decreasing = T)[1:20])

table(res[tma,])

write.csv(res[tma,], file = "XX.csv")
pheatmap::pheatmap(res[tma,],cluster_cols = F,cluster_rows = F,angle_col = 45,color =colorRampPalette(c( "white","#FFF8DC","#483D8B"))(70))
cv<- read.csv(file = "XX.csv",row.names = 1,header =T,sep = ",")
cv=t(cv)
write.csv(cv, file = "横裂生殖GOtop10.csv")
pheatmap::pheatmap(cv,cluster_cols = F,cluster_rows = F,angle_col =90,color =c("#D28370","#7496BA",'#674699'),fontsize_col = 10)




