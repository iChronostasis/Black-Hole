########################   Chronostasis   ############################
######################################################################
# Slingshot: Trajectory Inference for Single-Cell Data
## Load R packages
  library(slingshot) 
  library(Rtsne)
  library(mclust)
# Load the data(you wanna analyze)
## e.g. After selected the part of cell type for further analysis
  ident_df <- data.frame(cell=names(Idents(Seurat.data.combined)),celltype = Seurat.data.combined@active.ident)
  cl<-ident_df[ident_df$celltype=="SynTII"|ident_df$celltype=="SynTII Precursor"|ident_df$celltype=="LaTP",]
  part_celltype <- subset(Seurat.data.combined, cells=as.vector(ident_df[ident_df$celltype=="SynTII"|ident_df$celltype=="SynTII Precursor"|ident_df$celltype=="LaTP",1]))
#  prepare the input   
  sdata <- as.sparse(part_celltype@assays$RNA@data)
  sdata <- CreateSeuratObject(counts = sdata)
  sce <- as.SingleCellExperiment(sdata)

# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
  geneFilter <- apply(assays(sce)$counts,1,function(x){
    sum(x >= 3) >= 10
  })
  sce <- sce[geneFilter, ]

# Normalization    
  FQnorm <- function(counts){
    rk <- apply(counts,2,rank,ties.method='min')
    counts.sort <- apply(counts,2,sort)
    refdist <- apply(counts.sort,1,median)
    norm <- apply(rk,2,function(r){ refdist[r] })
    rownames(norm) <- rownames(counts)
    return(norm)
  }
  assays(sce)$norm <- FQnorm(assays(sce)$counts)

# Dimensionality Reduction
  pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
  rd1 <- pca$x[,1:2]
  
  plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

  rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
  colnames(rd2) <- c('UMAP1', 'UMAP2')
  
  plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
  
  reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)
  
# Clustering Cells
  library(mclust, quietly = TRUE)
  
  cl1 <- Mclust(rd1)$classification
  colData(sce)$GMM <- cl1
  
  library(RColorBrewer)
  plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
  
  cl2 <- kmeans(rd1, centers = 4)$cluster
  colData(sce)$kmeans <- cl2
  
  plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)
  
# Using Slingshot
  sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')
  
  summary(sce$slingPseudotime_1)
  
  library(grDevices)
  colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
  plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
  
  plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
  lines(SlingshotDataSet(sce), lwd=2, col='black')
  
  plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
  lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
  
#####################################################
# Downstream Analysis
## Identifying temporally dynamic genes
  library(tradeSeq)
  
# fit negative binomial GAM
  sce <- fitGAM(sce)
  
# test for dynamic expression
  ATres <- associationTest(sce)

  topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
  pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
  heatdata <- as.matrix(mycds@assayData[["exprs"]])
  heatmap(log1p(heatdata), Colv = NA,
          ColSideColors = brewer.pal(9,"Set1")[heatclus])
  
### Figures(According to the PMID:33141023)  
  library(pheatmap)
  DimPlot(part_celltype, reduction = "umap", repel = TRUE)
  pheatmap(as.matrix(log1p(heatdata)),annotation_col = all_label,cluster_rows = F,show_colnames = F,cluster_cols = F)
  
  
######################################################################  
########################### End ######################################  
######################################################################   
  
  
  
  
