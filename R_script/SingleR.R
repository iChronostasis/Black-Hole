#####   Chronostasis   #####
# Cell_typing
# SingleR细胞类型注释
  Seurat.data <- readRDS("Markers/Markers_0.5.rds")
# SCTransform()函数
  Seurat.data <- SCTransform(Seurat.data, vars.to.regress = "percent.mt")
# Load the SingleRdata  
  Seurat.data_for_SingleR <- GetAssayData(Seurat.data, slot="data")
  clusters=Seurat.data@meta.data$seurat_clusters
# 根据需求选择上述网站的数据集作为SingleR的参考数据集
  library(SingleR)
  library(scRNAseq)
  library(celldex)
#### Human being ####
#####################
## Using multiple references
  library(celldex)
  hpca <- celldex::HumanPrimaryCellAtlasData(ensembl=TRUE)# if you have the genename,choose FALSE
  bpe <- celldex::BlueprintEncodeData(ensembl=TRUE)

## annotation
  pred.human <- SingleR(test = Seurat.data_for_SingleR, assay.type.test=1,
                      ref = list(BPE=bpe, HPCA=hpca), 
                      labels = list(bpe$label.main, hpca$label.main))
  
## Check the final label from the combined assignment.
  table(pred.human$labels) 
  head(pred.human$orig.results$BPE$labels)
  head(pred.human$orig.results$HPCA$labels)

## visualization
  Seurat.data$hpca_type <- pred.human$orig.results$HPCA$labels
  Seurat.data$bpe_type <- pred.human$orig.results$BPE$labels
  
## visualization
  library(ggsci)
  p1 <- DimPlot(Seurat.data, group.by='hpca_type',pt.size = 3)+scale_color_ucscgb()
  p2 <- DimPlot(Seurat.data, group.by='bpe_type',pt.size = 3)+scale_color_ucscgb()
  library(cowplot)
  p<-plot_grid(p1, p2, nrow=1, ncol=2)
  ggsave("singler.pdf", plot = p, width = 15, height = 10) 
  
  
####### Mouse #######
#####################
  library(celldex)
  #immgen <- ImmGenData()
  immgen <- celldex::ImmGenData()
  mouseRNA <- celldex::MouseRNAseqData()
  monacoImm <- celldex::MonacoImmuneData()

## annotation
  pred.mouse <- SingleR(test = Seurat.data_for_SingleR, assay.type.test=1,
                        ref = list(Immgen=immgen, MouseRNA=mouseRNA,MonacoImmune=monacoImm), 
                        labels = list(immgen$label.main, mouseRNA$label.main,monacoImm$label.main))
  
## Check the final label from the combined assignment.
  table(pred.mouse$labels) 
  
## Check the 'winning' reference for each cell.
  table(pred.mouse$reference)   
  head(pred.mouse$orig.results$Immgen$labels)  
  head(pred.mouse$orig.results$MouseRNA$labels)
  head(pred.mouse$orig.results$MonacoImmune$labels)
    
## visualization
  Seurat.data$immgen <- pred.mouse$orig.results$Immgen$labels
  Seurat.data$MouseRNA <- pred.mouse$orig.results$MouseRNA$labels
  Seurat.data$MonacoImmun <- pred.mouse$orig.results$MonacoImmune$labels
  
  library(ggsci)
  p1 <- DimPlot(Seurat.data, group.by='immgen',pt.size = 3)+scale_color_ucscgb()
  p2 <- DimPlot(Seurat.data, group.by='MouseRNA',pt.size = 3)+scale_color_ucscgb()
  p3 <- DimPlot(Seurat.data, group.by='MonacoImmun',pt.size = 3)+scale_color_ucscgb()
  library(cowplot)
  p<-plot_grid(p1, p2,p3,nrow=1, ncol=3)
  ggsave("singler.pdf", plot = p, width = 15, height = 10)   
  
  
  
  
  
  
