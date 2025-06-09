sce = Seurat_tmp[,Seurat_tmp@meta.data$cellclass %in% c("PODO","PEC","GEC","MC","LEUK")]
sce <- SCTransform(sce, vars.to.regress = c("percent.mt","percent.rpl", "percent.rps","nCount_RNA"), verbose = TRUE)
sce <- RunPCA(sce, verbose = TRUE)
ElbowPlot(sce, ndims = 50)
sce <- FindNeighbors(sce, dims = 1:20, verbose = TRUE)

Seurat_tmp  = sce
network  <- as.matrix(Seurat_tmp@graphs$SCT_snn)
snn=as(network, "dgCMatrix")
snn_df=as.data.frame(summary(snn))

colnames(snn_df) <- c("from", "to", "weight")
g <- graph_from_data_frame(snn_df, directed = FALSE)
set.seed(0)

leiden<-leiden.community(g,resolution=0.6)
membership<-leiden[["membership"]]
index<-match(c(1:ncol(Seurat_tmp)),as.numeric(names(membership)))
leidengroup<-as.numeric(membership[index])

leidengroup


Seurat_tmp$leidengroup <- leidengroup
Seurat_tmp <- RunUMAP(Seurat_tmp, dims = 1:20, verbose = F)
p4 = DimPlot(Seurat_tmp, reduction = "umap", group.by = "cellclass",label = T,label.size = 5)+
  theme_bw()& 
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text = element_blank(), 
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5,size=14) 
  )&
  theme(aspect.ratio = 1)
p4
ggsave("p4.png", plot = p4, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)
ggsave("p4_leiden.pdf", plot = p4, path = "C:/R/DKD/pictures_pdf/", width = 8, height = 6, dpi = 600)

saveRDS(Seurat_tmp, "C:/R/DKD/Seurat_tmp_qiu.rds")
Seurat_tmp <- readRDS("C:/R/DKD/Seurat_tmp_qiu.rds")


vst_data <- read.csv("C:/R/DKD/exprset1.csv", row.names = 1)


Seurat_tmp@assays$SCT$data


shared_genes<-intersect(rownames(vst_data),rownames(Seurat_tmp))
sc_exprs <- as.matrix(Seurat_tmp@assays$SCT$data)
correlation_matrix<-cor(vst_data[shared_genes,],sc_exprs[shared_genes,])
dim(correlation_matrix)



phenotype <- ifelse(grepl("DKD", rownames(correlation_matrix)), 1, 0)
data = list(x = correlation_matrix, y = phenotype)
fit = SGL(data, leidengroup, type = "logit",alpha =0.5)
lam<-fit[['lambdas']]

cvfit<-cvSGL(data,leidengroup,type='logit',nfold = 5,alpha =0.5,lambdas = lam)
error<-cvfit$lldiff
m<-min(error)
h<-which(error==m)



a<-fit[["beta"]]
b<-a[,h]
LP_SGL_pos<-colnames(Seurat_tmp)[which(b>0)]
LP_SGL_neg<-colnames(Seurat_tmp)[which(b<0)]
Background<-colnames(Seurat_tmp)[which(b==0)]


Seurat_tmp$cell_type <- "Background"
Seurat_tmp$cell_type[LP_SGL_pos] <- "LP_SGL_pos"
Seurat_tmp$cell_type[LP_SGL_neg] <- "LP_SGL_neg"



p5 = DimPlot(Seurat_tmp, reduction = "umap", group.by = "cell_type", label = T,label.size = 5,
             cols = c("Background" = "gray", "LP_SGL_pos" = "red", "LP_SGL_neg" = "blue")) + 
  theme(legend.title = element_blank())  +
  theme_bw()&
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5,size=14) 
  )&
  theme(aspect.ratio = 1)
p5
ggsave("p5.png", plot = p5, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)
ggsave("p5.pdf", plot = p5, path = "C:/R/DKD/pictures_pdf/", width = 8, height = 6, dpi = 600)

saveRDS(Seurat_tmp, "C:/R/DKD/Seurat_leidengroup.rds")