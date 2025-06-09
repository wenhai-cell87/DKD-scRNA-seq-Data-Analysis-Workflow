library(ggVennDiagram)
library(Seurat) 
library(ggplot2) 
library(harmony)  
library(dplyr) 
library(stringr)  
library(tibble)  
library(here)  
library(future) 
library(hdf5r)
library(leidenAlg)
library(DESeq2)
library(SGL)
library(igraph)
library(ggvenn)
library(randomForest)
library(car)
library(EnhancedVolcano)
library(rms)
library(pROC)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(WGCNA)
set.seed(1234)  
aggr_input_dir <- "C:/R/DKD/sc_RNA-data/"


h5_files <- list.files(aggr_input_dir, pattern = "*DN", full.names = TRUE)
seurat_list <- lapply(h5_files, function(file) {
  
  data <- Read10X_h5(file)
  sample_group <-  if (grepl("DN", file)) {
    sub(".*(DN[0-9]+).*", "\\1", basename(file))  
  } else {
    "Unknown"
  }
  seurat_object <- CreateSeuratObject(counts = data, project = basename(file))
  seurat_object$group <- sample_group  
  return(seurat_object)
})


group_ids <- sapply(seurat_list, function(x) unique(x$group))

combined <- merge(seurat_list[[1]], y =seurat_list[-1], add.cell.ids = group_ids)

head(combined@meta.data)
rnaAggr = combined


rnaAggr <- AddMetaData(object = rnaAggr, metadata = data.frame(orig.ident = combined@meta.data$group, row.names = rownames(rnaAggr@meta.data)))
head(rnaAggr@meta.data)
rnaAggr <- PercentageFeatureSet(rnaAggr, pattern = "^MT-", col.name = "percent.mt")
rnaAggr <- PercentageFeatureSet(rnaAggr, pattern = "^RPL", col.name = "percent.rpl")
rnaAggr <- PercentageFeatureSet(rnaAggr, pattern = "^RPS", col.name = "percent.rps")
Idents(rnaAggr) = rnaAggr$group
pdf(file.path("C:\\R\\DKD\\pictures", "step1_qc.pdf"))
VlnPlot(object = rnaAggr, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0)
VlnPlot(object = rnaAggr, features = c("percent.mt", "percent.rps", "percent.rpl"), ncol = 2, pt.size = 0, y.max = 1)
VlnPlot(object = rnaAggr, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "orig.ident", pt.size = 0)
VlnPlot(object = rnaAggr, features = c("percent.mt", "percent.rps", "percent.rpl"), ncol = 2, group.by = "orig.ident", pt.size = 0, y.max = 1)
dev.off()

rnaAggr <- subset(rnaAggr, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA < 16000 & percent.mt < 0.5 & percent.rps < 0.3 & percent.rpl < 0.3)


rnaAggr <- SCTransform(rnaAggr, vars.to.regress = c("percent.mt","percent.rpl", "percent.rps","nCount_RNA"), verbose = TRUE)
rnaAggr <- RunPCA(rnaAggr, verbose = TRUE)
ElbowPlot(rnaAggr, ndims = 50)
#rnaAggr <- RunHarmony(rnaAggr, "orig.ident", plot_convergence = TRUE, assay.use="SCT")
rnaAggr <- FindNeighbors(rnaAggr, dims = 1:30, verbose = TRUE)


Seurat_tmp  = rnaAggr
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



Seurat_tmp$leidengroup <- leidengroup
Seurat_tmp <- RunUMAP(Seurat_tmp, dims = 1:30, verbose = F)
p1 = DimPlot(Seurat_tmp, reduction = "umap", group.by = "leidengroup",label = T)+
  theme_bw()& 
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),axis.text = element_blank(),
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5,size=14) 
  )&
  theme(aspect.ratio = 1)
p1
ggsave("p1.png", plot = p1, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)


genes = c(
  
  "PTPRC", # LEUK 白细胞
  "NPHS2", # PODO 足细胞
  "AQP2", #PC 主细胞
  "SLC12A3",#DCT 远曲小管
  "SLC26A7", #ICA A 型插入细胞
  "CFH","CLDN1",# PEC  肾小球壁上皮细胞
  "SLC26A4", #ICB B 型插入细胞
  "PDGFRB","ITGA8", # 系膜细胞 (MC) 
  "PECAM1","KDR", # GEC 肾小球内皮细胞
  "CUBN", "LRP2", # PCT 近端曲小管
  "SLC8A1", #CNT 连接小管
  "SLC12A1"#LOH/TAL 粗升支
  
)



FEA = FeaturePlot(Seurat_tmp,features = genes ,order = T, #raster = T,
                  #keep.scale = NULL,
                  reduction = "umap")&
  #scale_x_continuous("")&scale_y_continuous("")&
  theme_bw()& 
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),axis.text = element_blank(),
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5,size=14) 
  )&
  theme(aspect.ratio = 1)
FEA
ggsave("FEA.png", plot = FEA, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)


Idents(Seurat_tmp) = Seurat_tmp$leidengroup
options(repr.plot.width=16, repr.plot.height=10)
dot1 = DotPlot(Seurat_tmp, features = genes)+ 
  theme_bw()&
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    plot.title = element_text(hjust = 0.5,size=14) 
  )&
  theme(aspect.ratio = 1) +RotatedAxis()
dot1
ggsave("dot1.png", plot = dot1, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)



Seurat_tmp$cellclass <- case_when(
  Seurat_tmp$leidengroup %in% c("1") ~ "CNT",
  Seurat_tmp$leidengroup %in% c("2") ~ "DCT",
  Seurat_tmp$leidengroup %in% c("3") ~ "LOH/TAL",
  Seurat_tmp$leidengroup %in% c("4") ~ "PC",
  Seurat_tmp$leidengroup %in% c("5") ~ "PCT",
  Seurat_tmp$leidengroup %in% c("6") ~ "PCT",
  Seurat_tmp$leidengroup %in% c("7") ~ "PCT",
  Seurat_tmp$leidengroup %in% c("8") ~ "LOH/TAL",
  Seurat_tmp$leidengroup %in% c("9") ~ "LOH/TAL",
  Seurat_tmp$leidengroup %in% c("10") ~ "GEC",
  Seurat_tmp$leidengroup %in% c("11") ~ "ICA",
  Seurat_tmp$leidengroup %in% c("12") ~ "PCT",
  Seurat_tmp$leidengroup %in% c("13") ~ "PODO",
  Seurat_tmp$leidengroup %in% c("14") ~ "LOH/TAL",
  Seurat_tmp$leidengroup %in% c("15") ~ "PEC",
  Seurat_tmp$leidengroup %in% c("16") ~ "ICB",
  Seurat_tmp$leidengroup %in% c("17") ~ "DCT",
  Seurat_tmp$leidengroup %in% c("18") ~ "MC",
  Seurat_tmp$leidengroup %in% c("19") ~ "PCT",
  Seurat_tmp$leidengroup %in% c("20") ~ "LEUK"
)


p3 = DimPlot(Seurat_tmp, reduction = "umap", group.by = "cellclass",label = T,label.size = 5)+
  theme_bw()& 
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text = element_blank(), 
    legend.position = "none", 
    
    plot.title = element_text(hjust = 0.5,size=14) 
  )&
  theme(aspect.ratio = 1)
p3
ggsave("p3.png", plot = p3, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)
ggsave("p3.pdf", plot = p3, path = "C:/R/DKD/pictures_pdf/", width = 8, height = 6, dpi = 600)

Idents(Seurat_tmp) = Seurat_tmp$cellclass
dot2 = DotPlot(Seurat_tmp, features = genes)+ 
  theme_bw()& 
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    plot.title = element_text(hjust = 0.5,size=18) 
  )&
  theme(aspect.ratio = 1) +RotatedAxis()
dot2
ggsave("dot2.png", plot = dot2, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)
ggsave("dot2.pdf", plot = dot2, path = "C:/R/DKD/pictures_pdf/", width = 8, height = 6, dpi = 600)



H1 = DoHeatmap(Seurat_tmp, features =genes)+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))+
  theme(plot.title = element_text(size = 20, hjust = 0.5))

H1
ggsave("H1.png", plot = H1, path = "C:/R/DKD/pictures/", width = 12, height = 6, dpi = 600)
ggsave("H1.pdf", plot = H1, path = "C:/R/DKD/pictures_pdf/", width = 12, height = 6, dpi = 600)



saveRDS(Seurat_tmp, "C:/R/DKD/Seurat_tmp.rds")
Seurat_tmp <- readRDS("C:/R/DKD/Seurat_tmp.rds")

dot2 = DotPlot(Seurat_tmp, features = c("NRG1","SLC39A14"))+ 
  theme_bw()&
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    plot.title = element_text(hjust = 0.5,size=14)
  )&
  theme(aspect.ratio = 1) +RotatedAxis()
dot2

FEA1 = FeaturePlot(Seurat_tmp,features = genes ,order = T, #raster = T,
                   #keep.scale = NULL,
                   reduction = "umap")&
  #scale_x_continuous("")&scale_y_continuous("")&
  theme_bw()& 
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),axis.text = element_blank(), 
    legend.position = "none", #去掉图例
    plot.title = element_text(hjust = 0.5,size=14) 
  )&
  theme(aspect.ratio = 1)
FEA1

ggsave("FEA1.png", plot = FEA1, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)