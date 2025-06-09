
Seurat_tmp <- PrepSCTFindMarkers(Seurat_tmp)
DEG_LP_SGL <- FindMarkers(Seurat_tmp, 
                          ident.1 = "LP_SGL_pos", 
                          
                          group.by = "cell_type",
                          logfc.threshold = 1,  
                          min.pct = 0.1        
)   
dim(DEG_LP_SGL)
significant_genes <- subset(DEG_LP_SGL, p_val_adj < 0.05)
dim(significant_genes)
options(repr.plot.width=10, repr.plot.height=8)
vol1 <- EnhancedVolcano(significant_genes,
                        lab= rownames(significant_genes),
                        
                        x='avg_log2FC',
                        y='p_val_adj',
                        
                        pointSize=2,
                        
                        subtitle = NULL,
                        labSize = 5,
                        pCutoff = 0.001,     
                        FCcutoff = 1,         
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
                                           organism = "hsa", 
                                           pvalueCutoff = 0.05,
                                           qvalueCutoff = 0.05,
                                           pAdjustMethod = "BH",
                                           
                                           minGSSize = 10,#富集分析中考虑的最小基因集合大小
                                           maxGSSize = 500)#富集中考虑的最大基因集合大小

KEGG_enrich<-setReadable(KEGG_enrich,
                         OrgDb = org.Hs.eg.db,
                         keyType = 'ENTREZID')

KEGG_result<-KEGG_enrich@result
KEGG_result

library(enrichplot)
KEGG = dotplot(KEGG_enrich,
               x = "GeneRatio",
               color = "p.adjust",
               showCategory = 10)
KEGG
ggsave("KEGG.png", plot = KEGG, path = "C:/R/DKD/pictures/", width = 8, height = 10, dpi = 600)
ggsave("KEGG.pdf", plot = KEGG, path = "C:/R/DKD/pictures_pdf/", width = 8, height = 10, dpi = 600)

# GO

go_enrich<-clusterProfiler::enrichGO(gene = diff_entrez$ENTREZID,
                                     ont = 'all',
                                     keyType = "ENTREZID",
                                     OrgDb = org.Hs.eg.db,
                                     pAdjustMethod = "BH",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff = 0.05)

go_enrich<-DOSE::setReadable(go_enrich,
                             OrgDb = org.Hs.eg.db,
                             keyType = 'ENTREZID')


go_geo<- simplify(go_enrich, cutoff=0.7, by="p.adjust",
                  select_fun=min)

go_result<-go_geo@result
go_result


GO = dotplot(go_geo,
             x = "GeneRatio",
             color = "p.adjust",
             showCategory=8,
             split='ONTOLOGY',
             label_format = Inf)+
  
  facet_grid(ONTOLOGY~.,
             space = 'free_y',
             scale='free_y'
  )
GO
ggsave("GO.png", plot = GO, path = "C:/R/DKD/pictures/", width = 8, height = 10, dpi = 600)
ggsave("GO.pdf", plot = GO, path = "C:/R/DKD/pictures_pdf/", width = 8, height = 10, dpi = 600)

jihe = (p3+dot2)/(p4+p5)
ggsave("jihe.png", plot = jihe, path = "C:/R/DKD/pictures/", width = 12, height = 10, dpi = 600)


jihe1 = KEGG+GO
ggsave("jihe1.png", plot = jihe1, path = "C:/R/DKD/pictures/", width = 16, height = 10, dpi = 600)

head(vst_data)



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



exp_brca_final <- vst_data[rownames(vst_data) %in% common3, ]
dim(exp_brca_final)

exp_brca_final <- exp_brca_final[!duplicated(rownames(exp_brca_final)), ]
dim(exp_brca_final)


exp_brca_final <- as.data.frame(t(exp_brca_final))
exp_brca_final$sample <- rownames(exp_brca_final)
exp_brca_final$condition <- ifelse(grepl("DKD", exp_brca_final$sample), "DKD", "Con")
rownames(exp_brca_final) <- exp_brca_final$sample
lasso_data <- exp_brca_final[ , -226]
dim(lasso_data)




bulk_dataset = read.xlsx("C:/R/DKD/Raw_counts_glom.xlsx")
colData <- data.frame(row.names = colnames(bulk_dataset)[-1])
dds <- DESeqDataSetFromMatrix(countData = bulk_dataset[, -1],
                              colData = colData, 
                              design = ~ 1)

vsd <- vst(dds, blind = FALSE)
test_data <- assay(vsd)
rownames(test_data) <- bulk_dataset$Gene.symbol
dim(test_data)
head(test_data)


test_exp_brca_final <- test_data[rownames(test_data) %in% common3, ]
test_exp_brca_final <- test_exp_brca_final[!duplicated(rownames(test_exp_brca_final)), ]
dim(test_exp_brca_final)


test_exp_brca_final <- as.data.frame(t(test_exp_brca_final))
test_exp_brca_final$sample <- rownames(test_exp_brca_final)
test_exp_brca_final$condition <- ifelse(grepl("DN", test_exp_brca_final$sample), "DKD", "Con")
rownames(test_exp_brca_final) <- test_exp_brca_final$sample
test_lasso_data <- test_exp_brca_final[ , -226]
dim(test_lasso_data)


set.seed(1234)
x <- as.matrix(lasso_data[ , c(1:225)])
y <- as.matrix(lasso_data$condition)


x.test <- as.matrix(test_lasso_data[ , c(1:225)])
y.test <- as.matrix(test_lasso_data$condition)



library(glmnet)
alpha1_fit <- glmnet(x, y, alpha = 1, family = "binomial", nlambda = 100)
pdf("C:/R/DKD/pictures_pdf/lasso.pdf", width = 6, height = 6)
plot(alpha1_fit, xvar = "lambda", label = TRUE)
dev.off()


alpha1.fit.cv <- cv.glmnet(x, y, type.measure = "deviance", alpha = 1, family = "binomial")


pdf("C:/R/DKD/pictures_pdf/lasso1.pdf", width = 6, height = 6)
plot(alpha1.fit.cv, xvar = "lambda", label = TRUE)
dev.off()

print(alpha1.fit.cv)

feature_all <- as.data.frame(as.matrix(coef(alpha1.fit.cv, s = alpha1.fit.cv$lambda.1se)))
colnames(feature_all) <- "coff"
feature_opt <- feature_all %>% filter(abs(coff) > 0)
rownames(feature_opt)

feature<- feature_opt[-1, , drop = FALSE]
feature

x_df = as.data.frame(x)
y <- as.factor(lasso_data[,226]) 
x_df$score <- predict(alpha1.fit.cv, newx = as.matrix(x_df), s = "lambda.1se", type = "response")

ddist <- datadist(x_df)
options(datadist="ddist")


roc1<-roc(y,as.numeric(x_df$MALT1))
auc_value <- auc(roc1)
ci_auc <- ci.auc(roc1)  

pdf("C:/R/DKD/Supplementary figs1/MALT1.pdf", width = 6, height = 6)

plot(roc1, col = "red", lwd = 2)

legend_text <- paste("AUC =", sprintf("%.4f", auc_value), 
                     "\n95% CI:", 
                     sprintf("%.4f", ci_auc[1]), "-", sprintf("%.4f", ci_auc[3]))
legend("bottomright", legend = legend_text, col = "red", lwd = 2, bty = "n")
text(0.5, 0.8, labels = "MALT1", cex = 1.5, col = "black")
dev.off()


x_df = as.data.frame(x)
predCV.train <- predict(alpha1.fit.cv, newx = as.matrix(x_df),
                        s = "lambda.1se",
                        type = "response")
actuals.train <- ifelse(y == "DKD", 1, 0)


x.test_df = as.data.frame(x.test)
predCV.test  <- predict(alpha1.fit.cv, newx = as.matrix(x.test),
                        s = "lambda.1se",
                        type = "response")
actuals.test <- ifelse(y.test == "DKD", 1, 0)


predCV.train <- as.vector(predCV.train)
actuals.train <- as.vector(actuals.train)

png("C:/R/DKD/pictures/ROC_curve.png",width = 6*600, height = 6*600, res = 600)
pdf("C:/R/DKD/pictures_pdf/ROC_curve.pdf",width = 6, height =6)

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


legend.name <- c(paste("Train: AUC", sprintf("%.4f", x1[["auc"]])),
                 paste("Test: AUC", sprintf("%.4f", x2[["auc"]])))
legend("bottomright", 
       legend = legend.name,
       lwd = 2,
       col = c("red", "seagreen3"),
       bty = "n")

dev.off()
