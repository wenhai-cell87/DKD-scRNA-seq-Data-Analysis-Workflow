library(ggVennDiagram)
library(Seurat)  # 用于创建和处理Seurat对象
library(ggplot2)  # 用于绘图
library(harmony)  # 用于批次效应校正
#library(DoubletFinder)  # 用于双重细胞检测
library(dplyr)  # 数据操作
library(stringr)  # 字符串操作
library(tibble)  # 数据框处理
library(here)  # 用于定义相对路径
library(future)  # 用于并行计算
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
library(org.Hs.eg.db)  # 人类基因注释数据库
library(WGCNA)
packageVersion("Seurat")

#BiocManager"Seurat"#BiocManager::install("ggVenn")
#remotes::install_github("gaospecial/ggVennDiagram")
set.seed(1234)  # 设置随机种子确保结果可复现

aggr_input_dir <- "C:/R/DKD/sc_RNA-data/"


h5_files <- list.files(aggr_input_dir, pattern = "*DN", full.names = TRUE)
seurat_list <- lapply(h5_files, function(file) {
  # 读取数据
  data <- Read10X_h5(file)
  sample_group <-  if (grepl("DN", file)) {
    sub(".*(DN[0-9]+).*", "\\1", basename(file))  # 提取DN后面的数字
  } else {
    "Unknown"
  }
  seurat_object <- CreateSeuratObject(counts = data, project = basename(file))
  seurat_object$group <- sample_group  # 将组别（如 Control1, DN2）添加到 Seurat 对象的元数据
  return(seurat_object)
})


# 提取每个样本的分组信息作为cell ID前缀
group_ids <- sapply(seurat_list, function(x) unique(x$group))

# 合并
combined <- merge(seurat_list[[1]], y =seurat_list[-1], add.cell.ids = group_ids)

head(combined@meta.data)
rnaAggr = combined

# 将样本信息添加到元数据中
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

# 数据过滤
rnaAggr <- subset(rnaAggr, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA < 16000 & percent.mt < 0.5 & percent.rps < 0.3 & percent.rpl < 0.3)

# visualize doublets before proceeding with snRNA preprocessing
rnaAggr <- SCTransform(rnaAggr, vars.to.regress = c("percent.mt","percent.rpl", "percent.rps","nCount_RNA"), verbose = TRUE)
rnaAggr <- RunPCA(rnaAggr, verbose = TRUE)
ElbowPlot(rnaAggr, ndims = 50) # to determine number of dimensions for clustering
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


# 将 leidengroup 添加到 Seurat 对象的元数据中，绘制UMAP图，使用 leidengroup 进行分组
Seurat_tmp$leidengroup <- leidengroup
Seurat_tmp <- RunUMAP(Seurat_tmp, dims = 1:30, verbose = F)
p1 = DimPlot(Seurat_tmp, reduction = "umap", group.by = "leidengroup",label = T)+
  theme_bw()& #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
    axis.ticks = element_blank(),axis.text = element_blank(), #去掉坐标轴刻度和数字
    legend.position = "none", #去掉图例
    plot.title = element_text(hjust = 0.5,size=14) #改变标题位置和字体大小
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
  theme_bw()& #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
    axis.ticks = element_blank(),axis.text = element_blank(), #去掉坐标轴刻度和数字
    legend.position = "none", #去掉图例
    plot.title = element_text(hjust = 0.5,size=14) #改变标题位置和字体大小
  )&
  theme(aspect.ratio = 1)
FEA
ggsave("FEA.png", plot = FEA, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)


Idents(Seurat_tmp) = Seurat_tmp$leidengroup
options(repr.plot.width=16, repr.plot.height=10)
dot1 = DotPlot(Seurat_tmp, features = genes)+ 
theme_bw()& #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
    #axis.ticks = element_blank(),axis.text = element_blank(), #去掉坐标轴刻度和数字
    #legend.position = "none", #去掉图例
    plot.title = element_text(hjust = 0.5,size=14) #改变标题位置和字体大小
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
  theme_bw()& #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
    axis.ticks = element_blank(),axis.text = element_blank(), #去掉坐标轴刻度和数字
    legend.position = "none", #去掉图例
    
    plot.title = element_text(hjust = 0.5,size=14) #改变标题位置和字体大小
  )&
  theme(aspect.ratio = 1)
p3
ggsave("p3.png", plot = p3, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)
ggsave("p3.pdf", plot = p3, path = "C:/R/DKD/pictures_pdf/", width = 8, height = 6, dpi = 600)

Idents(Seurat_tmp) = Seurat_tmp$cellclass
dot2 = DotPlot(Seurat_tmp, features = genes)+ 
  theme_bw()& #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
    #axis.ticks = element_blank(),axis.text = element_blank(), #去掉坐标轴刻度和数字
    #legend.position = "none", #去掉图例
    plot.title = element_text(hjust = 0.5,size=18) #改变标题位置和字体大小
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
  theme_bw()& #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
    #axis.ticks = element_blank(),axis.text = element_blank(), #去掉坐标轴刻度和数字
    #legend.position = "none", #去掉图例
    plot.title = element_text(hjust = 0.5,size=14) #改变标题位置和字体大小
  )&
  theme(aspect.ratio = 1) +RotatedAxis()
dot2

FEA1 = FeaturePlot(Seurat_tmp,features = genes ,order = T, #raster = T,
                  #keep.scale = NULL,
                  reduction = "umap")&
  #scale_x_continuous("")&scale_y_continuous("")&
  theme_bw()& #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
    axis.ticks = element_blank(),axis.text = element_blank(), #去掉坐标轴刻度和数字
    legend.position = "none", #去掉图例
    plot.title = element_text(hjust = 0.5,size=14) #改变标题位置和字体大小
  )&
  theme(aspect.ratio = 1)
FEA1

ggsave("FEA1.png", plot = FEA1, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)




sce = Seurat_tmp[,Seurat_tmp@meta.data$cellclass %in% c("PODO","PEC","GEC","MC","LEUK")]
sce <- SCTransform(sce, vars.to.regress = c("percent.mt","percent.rpl", "percent.rps","nCount_RNA"), verbose = TRUE)
sce <- RunPCA(sce, verbose = TRUE)
ElbowPlot(sce, ndims = 50) # to determine number of dimensions for clustering
#rnaAggr <- RunHarmony(rnaAggr, "orig.ident", plot_convergence = TRUE, assay.use="SCT")
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

# 将 leidengroup 添加到 Seurat 对象的元数据中，绘制UMAP图，使用 leidengroup 进行分组
Seurat_tmp$leidengroup <- leidengroup
Seurat_tmp <- RunUMAP(Seurat_tmp, dims = 1:20, verbose = F)
p4 = DimPlot(Seurat_tmp, reduction = "umap", group.by = "cellclass",label = T,label.size = 5)+
  theme_bw()& #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
    axis.ticks = element_blank(),axis.text = element_blank(), #去掉坐标轴刻度和数字
    legend.position = "none", #去掉图例
    plot.title = element_text(hjust = 0.5,size=14) #改变标题位置和字体大小
  )&
  theme(aspect.ratio = 1)
p4
ggsave("p4.png", plot = p4, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)
ggsave("p4_leiden.pdf", plot = p4, path = "C:/R/DKD/pictures_pdf/", width = 8, height = 6, dpi = 600)

saveRDS(Seurat_tmp, "C:/R/DKD/Seurat_tmp_qiu.rds")
Seurat_tmp <- readRDS("C:/R/DKD/Seurat_tmp_qiu.rds")





vst_data <- read.csv("C:/R/DKD/exprset1.csv", row.names = 1)


Seurat_tmp@assays$SCT$data
R.version
# 找共享基因并计算相关性矩阵
shared_genes<-intersect(rownames(vst_data),rownames(Seurat_tmp))
sc_exprs <- as.matrix(Seurat_tmp@assays$SCT$data)
correlation_matrix<-cor(vst_data[shared_genes,],sc_exprs[shared_genes,])
dim(correlation_matrix)


# 创建SGL模型所需的数据列表、训练SGL模型
phenotype <- ifelse(grepl("DKD", rownames(correlation_matrix)), 1, 0)
data = list(x = correlation_matrix, y = phenotype)
fit = SGL(data, leidengroup, type = "logit",alpha =0.5)
lam<-fit[['lambdas']]
# 使用交叉验证选择最优正则化参数
cvfit<-cvSGL(data,leidengroup,type='logit',nfold = 5,alpha =0.5,lambdas = lam)
error<-cvfit$lldiff
m<-min(error)
h<-which(error==m)


# 提取模型的系数并分类
a<-fit[["beta"]]
b<-a[,h]
LP_SGL_pos<-colnames(Seurat_tmp)[which(b>0)]
LP_SGL_neg<-colnames(Seurat_tmp)[which(b<0)]
Background<-colnames(Seurat_tmp)[which(b==0)]

# 将LP_SGL_pos和LP_SGL_neg放入Seurat对象
Seurat_tmp$cell_type <- "Background"
Seurat_tmp$cell_type[LP_SGL_pos] <- "LP_SGL_pos"
Seurat_tmp$cell_type[LP_SGL_neg] <- "LP_SGL_neg"



p5 = DimPlot(Seurat_tmp, reduction = "umap", group.by = "cell_type", label = T,label.size = 5,
             cols = c("Background" = "gray", "LP_SGL_pos" = "red", "LP_SGL_neg" = "blue")) + 
  theme(legend.title = element_blank())  +
  theme_bw()& #改变ggplot2的主题
  theme( #进一步修改主题
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), #去掉背景线
    axis.ticks = element_blank(),axis.text = element_blank(), #去掉坐标轴刻度和数字
    legend.position = "none", #去掉图例
    plot.title = element_text(hjust = 0.5,size=14) #改变标题位置和字体大小
  )&
  theme(aspect.ratio = 1)
p5
ggsave("p5.png", plot = p5, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)
ggsave("p5.pdf", plot = p5, path = "C:/R/DKD/pictures_pdf/", width = 8, height = 6, dpi = 600)

saveRDS(Seurat_tmp, "C:/R/DKD/Seurat_leidengroup.rds")



# 计算差异基因并可视化
Seurat_tmp <- PrepSCTFindMarkers(Seurat_tmp)
DEG_LP_SGL <- FindMarkers(Seurat_tmp, 
                          ident.1 = "LP_SGL_pos", 
                          #ident.2 = "LP_SGL_neg", 
                          group.by = "cell_type",
                          logfc.threshold = 1,  # 设置log fold change的阈值
                          min.pct = 0.1         # 设定在多少比例的细胞中表达该基因
)   
dim(DEG_LP_SGL)
significant_genes <- subset(DEG_LP_SGL, p_val_adj < 0.05)
dim(significant_genes)
options(repr.plot.width=10, repr.plot.height=8)
vol1 <- EnhancedVolcano(significant_genes,
                        lab= rownames(significant_genes),
                        #selectLab = genes,
                        x='avg_log2FC',
                        y='p_val_adj',
                        
                        pointSize=2,
                        #title = "LP_SGL_pos VS LP_SGL_neg",
                        subtitle = NULL,
                        labSize = 5,
                        pCutoff = 0.001,      # p值截断值：水平线
                        FCcutoff = 1,         # FC截断值：垂直线
                        #boxedLabels = T,
                        parseLabels = F,labFace = "bold",drawConnectors = T,widthConnectors = 1.0,lengthConnectors = unit(0.01, "npc"),
                        legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", expression(p - value ~ and
                                                                                              ~ log[2] ~ FC)),
                        gridlines.major = F,
                        gridlines.minor = F,
)+theme(aspect.ratio = 1)
vol1
ggsave("vol1.png", plot = vol1, path = "C:/R/DKD/pictures/", width = 8, height = 6, dpi = 600)
ggsave("vol1.pdf", plot = vol1, path = "C:/R/DKD/pictures_pdf/", width = 8, height = 6, dpi = 600)


# KEGG
diff_entrez<-bitr(
  common,
  fromType='SYMBOL',
  toType='ENTREZID',
  OrgDb='org.Hs.eg.db'
)
head(diff_entrez)


KEGG_enrich <- clusterProfiler::enrichKEGG(gene = diff_entrez$ENTREZID,
                                           organism = "hsa", #物种Homo sapiens 
                                           pvalueCutoff = 0.05,#pvalue阈值
                                           qvalueCutoff = 0.05,#qvalue阈值
                                           pAdjustMethod = "BH",#p值矫正方法，one of "holm", 
                                           #"hochberg", "hommel", 
                                           #"bonferroni", "BH", 
                                           #"BY", "fdr", "none"
                                           minGSSize = 10,#富集分析中考虑的最小基因集合大小
                                           maxGSSize = 500)#富集中考虑的最大基因集合大小
#将RNTREZ转换为Symbo
KEGG_enrich<-setReadable(KEGG_enrich,
                         OrgDb = org.Hs.eg.db,
                         keyType = 'ENTREZID')
#提取KEGG富集结果表格
KEGG_result<-KEGG_enrich@result
KEGG_result

library(enrichplot)#R包加载
KEGG = dotplot(KEGG_enrich,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory = 10)#展示top10通路
KEGG
ggsave("KEGG.png", plot = KEGG, path = "C:/R/DKD/pictures/", width = 8, height = 10, dpi = 600)
ggsave("KEGG.pdf", plot = KEGG, path = "C:/R/DKD/pictures_pdf/", width = 8, height = 10, dpi = 600)

# GO

go_enrich<-clusterProfiler::enrichGO(gene = diff_entrez$ENTREZID,
                                     ont = 'all',#可选'BP','CC','MF' or 'all'
                                     keyType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db,
                                     pAdjustMethod = "BH",#p值矫正方法
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)
#将RNTREZ转换为Symbol
go_enrich<-DOSE::setReadable(go_enrich,
                             OrgDb = org.Hs.eg.db,
                             keyType = 'ENTREZID')

#去除冗余的GO term
go_geo<- simplify(go_enrich, cutoff=0.7, by="p.adjust",
                  select_fun=min)
#提取goG富集结果表格 
go_result<-go_geo@result
go_result


GO = dotplot(go_geo,
        x = "GeneRatio",
        color = "p.adjust",
        showCategory=8,
        split='ONTOLOGY',
        label_format = Inf)+#不换行
  #分面
  facet_grid(ONTOLOGY~.,
             space = 'free_y',#面板大小根据y轴自行调整
             scale='free_y'#子图坐标轴根据y轴自行调整
  )
GO
ggsave("GO.png", plot = GO, path = "C:/R/DKD/pictures/", width = 8, height = 10, dpi = 600)
ggsave("GO.pdf", plot = GO, path = "C:/R/DKD/pictures_pdf/", width = 8, height = 10, dpi = 600)

jihe = (p3+dot2)/(p4+p5)
ggsave("jihe.png", plot = jihe, path = "C:/R/DKD/pictures/", width = 12, height = 10, dpi = 600)


jihe1 = KEGG+GO
ggsave("jihe1.png", plot = jihe1, path = "C:/R/DKD/pictures/", width = 16, height = 10, dpi = 600)

head(vst_data)


# 分泌蛋白
# 读取.tsv文件
data <- read.table("C:/R/DKD/protein_class_SPOCTOPUS.tsv", header = TRUE, sep = "\t")

protein <- intersect(rownames(significant_genes), data$Gene)

x <- list(differential_gene=rownames(significant_genes), secreted_proteins=data$Gene)

venn = ggVennDiagram(x, label_percent_digit = 3  ,set_size = 3)+ scale_fill_gradient(low="grey90",high = "red")
venn

pUp <- ggvenn(x,c("differential_gene", "secreted_proteins"), stroke_linetype = 1, stroke_size = 1.1,
              stroke_color="#00332c" ,digits = 0,
              set_name_color = "black",
              show_percentage = F,
              set_name_size = 6,
              fill_color = c("#fef6e4","#8bd3dd"),fill_alpha=0.8,
              text_color="#00332c",text_size=5) +
  theme_minimal() +
  theme(panel.grid = element_blank(),text=element_blank(),
        plot.title = element_text(color="black", size=16, face="bold",hjust=0,vjust=-15))
pUp





ggsave("pUp.png", plot = pUp, path = "C:/R/DKD/pictures/", width = 6, height = 6, dpi = 600)

ggsave("pUp.pdf", plot = pUp, path = "C:/R/DKD/pictures_pdf/", width = 6, height = 6, dpi = 600)





common <- intersect(rownames(significant_genes), data$Gene)

common1 <- intersect(rownames(vst_data), common)

common2 <- intersect(rownames(test_data), common)


common3 <- intersect(common1, common2)



# train数据集  在bulk转录组数据中提取出差异基因的表达矩阵
exp_brca_final <- vst_data[rownames(vst_data) %in% common3, ]
dim(exp_brca_final)
# 去除重复的行
exp_brca_final <- exp_brca_final[!duplicated(rownames(exp_brca_final)), ]
dim(exp_brca_final)

# 数据整理
exp_brca_final <- as.data.frame(t(exp_brca_final))
exp_brca_final$sample <- rownames(exp_brca_final)
exp_brca_final$condition <- ifelse(grepl("DKD", exp_brca_final$sample), "DKD", "Con")
rownames(exp_brca_final) <- exp_brca_final$sample
lasso_data <- exp_brca_final[ , -226]
dim(lasso_data)
#write.csv(lasso_data, "C:/vscode/pictures_DKD/lasso_data.csv", row.names = TRUE)


# test数据集 
# 读入bulk转录组数据并进行数据整理 (使用 DESeq2包对数据进行方差稳定变换 (VST) 归一化)
bulk_dataset = read.xlsx("C:/R/DKD/Raw_counts_glom.xlsx")
# 创建样本信息
colData <- data.frame(row.names = colnames(bulk_dataset)[-1])
# 创建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = bulk_dataset[, -1],
                              colData = colData, 
                              design = ~ 1)

vsd <- vst(dds, blind = FALSE)
test_data <- assay(vsd)
rownames(test_data) <- bulk_dataset$Gene.symbol
dim(test_data)
head(test_data)
#test_data <- read.csv("C:/R/DKD/exprset1.csv", row.names = 1)


# 在bulk转录组数据中提取出差异基因的表达矩阵
test_exp_brca_final <- test_data[rownames(test_data) %in% common3, ]
test_exp_brca_final <- test_exp_brca_final[!duplicated(rownames(test_exp_brca_final)), ]
dim(test_exp_brca_final)


# 数据整理
test_exp_brca_final <- as.data.frame(t(test_exp_brca_final))
test_exp_brca_final$sample <- rownames(test_exp_brca_final)
test_exp_brca_final$condition <- ifelse(grepl("DN", test_exp_brca_final$sample), "DKD", "Con")
rownames(test_exp_brca_final) <- test_exp_brca_final$sample
test_lasso_data <- test_exp_brca_final[ , -226]
dim(test_lasso_data)
#write.csv(lasso_data, "C:/vscode/pictures_DKD/lasso_data.csv", row.names = TRUE)


# 设置自变量和因变量
set.seed(1234)
x <- as.matrix(lasso_data[ , c(1:225)])
y <- as.matrix(lasso_data$condition)


# 设置自变量和因变量
x.test <- as.matrix(test_lasso_data[ , c(1:225)])
y.test <- as.matrix(test_lasso_data$condition)









                                                                                                                   

#使用glmnet()建模
library(glmnet)
alpha1_fit <- glmnet(x, y, alpha = 1, family = "binomial", nlambda = 100)
# family参数用于指定模型所假设的概率分布族，这将影响到模型的形式和假设。不同的问题类型和数据类型需要选择不同的分布族
# "binomial"：用于二元分类问题，假设因变量是二元的，表示患病与正常或1与0的概率，这是逻辑回归（logistic distribution）的常见用途。
# nlambda: 这个参数表示要拟合的lambda（正则化参数）的数量。Lambda是一个控制正则化强度的参数，
#它控制着Lasso回归中特征的稀疏性程度。nlambda 指定了在不同lambda值上进行模型拟合的次数，
#以便找到合适的正则化强度。一般使用默认值100
pdf("C:/R/DKD/pictures_pdf/lasso.pdf", width = 6, height = 6)
#png("C:/R/DKD/pictures/lasso.png", width = 6*600, height = 6*600, res = 600)
plot(alpha1_fit, xvar = "lambda", label = TRUE)
dev.off()


# 交叉验证
alpha1.fit.cv <- cv.glmnet(x, y, type.measure = "deviance", alpha = 1, family = "binomial")
#type.measure 参数用于指定在进行交叉验证时用于度量模型性能的指标。
#ype.measure 被设置为 "deviance"，这是用于分类问题的一种常见指标，特别是对于逻辑回归等模型。
#它表示对数似然的负二倍，可用于评估分类模型的拟合。


# 绘图
pdf("C:/R/DKD/pictures_pdf/lasso1.pdf", width = 6, height = 6)
#png("C:/R/DKD/pictures/lasso1.png", width = 6*600, height = 6*600, res = 600)
plot(alpha1.fit.cv, xvar = "lambda", label = TRUE)
dev.off()
# 左边虚线为λ min，意思是偏差最小时的λ，代表在该lambda取值下，模型拟合效果最高。变量数是29，相比λ-se，保留下来的变量更多
print(alpha1.fit.cv)
# 提取特征
feature_all <- as.data.frame(as.matrix(coef(alpha1.fit.cv, s = alpha1.fit.cv$lambda.1se)))
colnames(feature_all) <- "coff"
feature_opt <- feature_all %>% filter(abs(coff) > 0)
rownames(feature_opt)

feature<- feature_opt[-1, , drop = FALSE]
feature

x_df = as.data.frame(x)
y <- as.factor(lasso_data[,226])  # 将分类标签转为因子类型
x_df$score <- predict(alpha1.fit.cv, newx = as.matrix(x_df), s = "lambda.1se", type = "response")

ddist <- datadist(x_df)
options(datadist="ddist")


roc1<-roc(y,as.numeric(x_df$MALT1))
auc_value <- auc(roc1)
ci_auc <- ci.auc(roc1)  # 计算AUC的置信区间
# 创建图像并保存
pdf("C:/R/DKD/Supplementary figs1/MALT1.pdf", width = 6, height = 6)
#png("C:/R/DKD/pictures/ADAM15.png", width = 6*600, height = 6*600, res = 600)
plot(roc1, col = "red", lwd = 2)
# 添加AUC值到图例
legend_text <- paste("AUC =", sprintf("%.4f", auc_value), 
                       "\n95% CI:", 
                       sprintf("%.4f", ci_auc[1]), "-", sprintf("%.4f", ci_auc[3]))
legend("bottomright", legend = legend_text, col = "red", lwd = 2, bty = "n")
# 添加基因名字（例如，ADAM15）到图上
text(0.5, 0.8, labels = "MALT1", cex = 1.5, col = "black")
dev.off()






# 开始进行模型检验
#首先是训练集
x_df = as.data.frame(x)
predCV.train <- predict(alpha1.fit.cv, newx = as.matrix(x_df),
                                                 s = "lambda.1se",
                                                 type = "response")
actuals.train <- ifelse(y == "DKD", 1, 0)

#其次是测试集
x.test_df = as.data.frame(x.test)
predCV.test  <- predict(alpha1.fit.cv, newx = as.matrix(x.test),
                                                 s = "lambda.1se",
                                                 type = "response")
actuals.test <- ifelse(y.test == "DKD", 1, 0)


predCV.train <- as.vector(predCV.train)
actuals.train <- as.vector(actuals.train)
## plot ROC 
# 保存为PNG文件
png("C:/R/DKD/pictures/ROC_curve.png",width = 6*600, height = 6*600, res = 600)
pdf("C:/R/DKD/pictures_pdf/ROC_curve.pdf",width = 6, height =6)
# 生成图像
x1 <- plot.roc(actuals.train, predCV.train,
               smooth = F,
               lwd = 2,
               ylim = c(0, 1),
               xlim = c(1, 0),
               legacy.axes = T,
               main = "",
               col = "red")

actuals.test <- as.vector(actuals.test)
predCV.test <- as.vector(predCV.test)
x2 <- plot.roc(actuals.test, predCV.test,
               smooth = F,
               add = T,
               lwd = 2,
               ylim = c(0, 1),
               xlim = c(1, 0),
               legacy.axes = T,
               main = "",
               col = "seagreen3")

# 添加AUC值到图例
legend.name <- c(paste("Train: AUC", sprintf("%.4f", x1[["auc"]])),
                 paste("Test: AUC", sprintf("%.4f", x2[["auc"]])))
legend("bottomright", 
       legend = legend.name,
       lwd = 2,
       col = c("red", "seagreen3"),
       bty = "n")

# 关闭设备以保存图像
dev.off()






