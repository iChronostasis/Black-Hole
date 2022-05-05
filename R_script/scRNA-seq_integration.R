########################   Chronostasis   ############################
######################################################################
### 6个数据进行合并
  setwd("/Users/hecate/研一/Term_2/Project/out")
# 创建每个单独的Seurat对象对于每个样本
  seurat_obj <- list()
  for (file in c("outs-N-uninfected", "outs-1day","outs-3day","outs-5day","outs-7day","outs-12day")){
    seurat_data <- Read10X(data.dir = paste0(file,"/filtered_gene_bc_matrices/mm10flu_v1"))
    seurat_obj[[file]] <- CreateSeuratObject(counts = seurat_data, 
                                             min.features = 100, 
                                             project = file)
    #assign(file, seurat_obj)
  }

# 合并数据集  
  Seurat.data.combined <- merge(seurat_obj[[1]], y=c(seurat_obj[[2]], seurat_obj[[3]], 
                                                     seurat_obj[[4]], seurat_obj[[5]], seurat_obj[[6]]))

# split the dataset into a list of two seurat objects (stim and CTRL)
  ifnb.list <- SplitObject(Seurat.data.combined, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
  ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

# select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = ifnb.list)

  immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

# this command creates an 'integrated' data assay
  Seurat.data.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
  DefaultAssay(Seurat.data.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
  Seurat.data.combined <- ScaleData(Seurat.data.combined, verbose = FALSE)
  Seurat.data.combined <- RunPCA(Seurat.data.combined, npcs = 30, verbose = FALSE)
  Seurat.data.combined <- RunUMAP(Seurat.data.combined, reduction = "pca", dims = 1:15)
  Seurat.data.combined <- FindNeighbors(Seurat.data.combined, reduction = "pca", dims = 1:15)
  Seurat.data.combined <- FindClusters(Seurat.data.combined, resolution = 0.6) 
  
# Visualization
  dir.create("Visualization")
  p1 <- DimPlot(Seurat.data.combined, reduction = "umap", group.by = "orig.ident")
  p2 <- DimPlot(Seurat.data.combined, reduction = "umap", split.by = "orig.ident")
  #p2 <- DimPlot(Seurat.data.combined, reduction = "umap", label = TRUE, repel = TRUE)
  p<-p1 + p2
  ggsave("Visualization/cluster.pdf", plot = p, width = 9, height = 8) 
  
# Finding differentially expressed features (cluster biomarkers)
  DefaultAssay(Seurat.data.combined) <- "integrated"
  Seurat.data.combined.markers <- FindAllMarkers(Seurat.data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  Seurat.data.combined.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  
# 
  
  
  
  
  