########################   Chronostasis   ##########################################
####################################################################################
  rm(list=ls())
## Setup the Seurat Object
  library(dplyr)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(data.table)
# Load the dataset(some types of dataset)
  setwd("~/Project/result")
### 1. input:filtered_gene_bc_matrices/~/barcodes.tsv，genes.tsv，matrix.mtx  
  Seurat.data <- Read10X(data.dir = "~/aggregation/filtered_gene_bc_matrices_mex/mm10flu_v1/")
# Initialize the Seurat object with the raw (non-normalized data).
  Seurat.data <- CreateSeuratObject(counts = Seurat.data, project = "mm10flu_v1", min.cells = 3, min.features = 200)

### 2. when the file is csv or txt(when you need)
  #data<-fread("GSE152248_AllStages_AllNuclei_datamatrix.txt")[,-1]
  #gene<-fread("GSE152248_AllStages_AllNuclei_datamatrix.txt")[,1]
  #row.names(data)<-gene$V1
  #colnames(data)[1]<-c("Gene")
  data<-read.csv("GSE127136_project_IgA_nephropathy_counts.csv",header=T,row.names=1)
  Seurat.data <- CreateSeuratObject(counts = data)
  #meta.data<-fread("../GSE152248_AllStages_AllNuclei_clusters.txt")
  #row.names(meta.data) <- meta.data$cell_id
  #Seurat.data <- AddMetaData(object = Seurat.data, metadata = meta.data)

### 3. when the file is .h5(when you need)
#### e.g.
  #download.file("https://cf.10xgenomics.com/samples/cell-exp/4.0.0/Parent_NGSC3_DI_PBMC/Parent_NGSC3_DI_PBMC_filtered_feature_bc_matrix.h5",
              #destfile = "3p_pbmc10k_filt.h5")
  Seurat.data <- Read10X_h5("data/update/3p_pbmc10k_filt.h5",use.names = T)
  Seurat.data <- CreateSeuratObject(Seurat.data,project = " ")

#################################################################################

## Standard pre-processing workflow
  dir.create("QC")
###### QC and selecting cells for further analysis
### 计算质控指标
# 计算细胞中线粒体基因比例
  Seurat.data[["percent.mt"]] <- PercentageFeatureSet(Seurat.data, pattern = "^MT-")
# 计算细胞中核糖体基因比例
  Seurat.data[["percent.rb"]] <- PercentageFeatureSet(Seurat.data, pattern = "^RP[LS]")
# 计算红细胞比例(For need)
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB.genes <- CaseMatch(HB.genes, rownames(Seurat.data))
  Seurat.data[["percent.HB"]]<-PercentageFeatureSet(Seurat.data, features=HB.genes) 
  
### 查看质控指标
# 设置绘图元素
  theme.set2 = theme(axis.title.x=element_blank())
  plot.featrures = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.rb", "percent.HB")
  group = "orig.ident"
  # 质控前小提琴图
  plots = list()
  for(i in seq_along(plot.featrures)){
    plots[[i]] = VlnPlot(Seurat.data, group.by=group, pt.size = 0,
                         features = plot.featrures[i]) + theme.set2 + NoLegend()}
  violin <- wrap_plots(plots = plots, nrow=2)    
  ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 9, height = 8) 
  
# 设置质控指标
  plots = list()
  for(i in seq_along(plot.featrures)){
    plots[[i]] = VlnPlot(Seurat.data, group.by=group, pt.size = 0,
                         features = plot.featrures[i]) + theme.set2 + NoLegend()}
  quantile(Seurat.data$nFeature_RNA, seq(0.01, 0.1, 0.01))
  quantile(Seurat.data$nFeature_RNA, seq(0.9, 1, 0.01))
  plots[[1]] + geom_hline(yintercept = 500) + geom_hline(yintercept = 4500)
  quantile(Seurat.data$nCount_RNA, seq(0.9, 1, 0.01))
  plots[[2]] + geom_hline(yintercept = 22000)
  #quantile(Seurat.data$percent.mt, seq(0.9, 1, 0.01))
  #plots[[3]] + geom_hline(yintercept = 20)
  #quantile(Seurat.data$percent.HB, seq(0.9, 1, 0.01))
  #plots[[5]] + geom_hline(yintercept = 1)
  
### 质控
# 设置质控标准
  minGene=200
  maxGene=2500
  #maxUMI=22000
  pctMT=5
  #pctHB=1
# 数据质控并绘制小提琴图
  Seurat.data <- subset(Seurat.data, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & percent.mt < pctMT )
  plots = list()
  for(i in seq_along(plot.featrures)){
    plots[[i]] = VlnPlot(Seurat.data, group.by=group, pt.size = 0.01,
                         features = plot.featrures[i]) + theme.set2 + NoLegend()}
  violin <- wrap_plots(plots = plots, nrow=2)    
  ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 9, height = 8) 
  
# Normalizing the data
  Seurat.data <- NormalizeData(Seurat.data, normalization.method = "LogNormalize", scale.factor = 10000)
  
# Identification of highly variable features
  Seurat.data <- FindVariableFeatures(Seurat.data, selection.method = "vst", nfeatures = 2000)
  
# Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(Seurat.data), 10)
  
# plot variable features with and without labels
  plot1 <- VariableFeaturePlot(Seurat.data)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  
# Scaling the data
  all.genes <- rownames(Seurat.data)
  Seurat.data <- ScaleData(Seurat.data, features = all.genes)

# SCTransform()函数 ====== NormalizeData,ScaleData,FindVariableFeatures
  Seurat.data <- SCTransform(Seurat.data, vars.to.regress = "percent.mt")  ###c("percent.mt","percent.rb")
################################################################################################################################
#############################################Optional###########################################################################
# Cell cycle scoring (in order to mitigate the effects of cell cycle heterogeneity in the data)(when you need ==>> optional) ###
  Seurat.data <- CellCycleScoring(Seurat.data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)             ### 
                                                                                                                             ###
# view cell cycle scores and phase assignments                                                                               ###
  head(Seurat.data[[]])                                                                                                      ###
                                                                                                                             ###
# Visualize the distribution of cell cycle markers across                                                                    ###
  RidgePlot(Seurat.data, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)                                           ###
                                                                                                                             ###
# Re-Scale the data                                                                                                          ###
  Seurat.data <- ScaleData(Seurat.data, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Seurat.data))       ###
                                                                                                                             ###
# Now, a PCA on the variable genes no longer returns components associated with cell cycle(compare to the next line)         ###
  #Seurat.data <- RunPCA(Seurat.data, features = VariableFeatures(Seurat.data), nfeatures.print = 10)                        ###
                                                                                                                             ###
# When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase                                  ###
  Seurat.data <- RunPCA(Seurat.data, features = c(s.genes, g2m.genes))                                                       ###
  DimPlot(Seurat.data)                                                                                                       ###
################################################################################################################################
################################################################################################################################


# Perform linear dimensional reduction
  Seurat.data <- RunPCA(Seurat.data, features = VariableFeatures(object = Seurat.data))
# Examine and visualize PCA results a few different ways
  print(Seurat.data[["pca"]], dims = 1:5, nfeatures = 5)
  VizDimLoadings(Seurat.data, dims = 1:2, reduction = "pca")
  DimPlot(Seurat.data, reduction = "pca")
### DimHeatmap()允许我们探索数据中的最初的异质性，然后决定下游使用多少个PC进行分析。
  DimHeatmap(Seurat.data, dims = 1, cells = 500, balanced = TRUE)
  DimHeatmap(Seurat.data, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the ‘dimensionality’ of the dataset
  Seurat.data <- JackStraw(Seurat.data, num.replicate = 100)
  Seurat.data <- ScoreJackStraw(Seurat.data, dims = 1:20)
  JackStrawPlot(Seurat.data, dims = 1:15)
  ElbowPlot(Seurat.data)

# Cluster the cells
  dir.create("Cluster")
  
  Seurat.data <- FindNeighbors(Seurat.data, dims = 1:10)
  Seurat.data <- FindClusters(Seurat.data, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
  head(Idents(Seurat.data), 5)
  
# Run non-linear dimensional reduction (UMAP/tSNE)  
  Seurat.data <- RunUMAP(Seurat.data, dims = 1:10)
  #Seurat.data <- RunTSNE(Seurat.data, npcs = 30, verbose = FALSE)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
  DimPlot(Seurat.data, reduction = "umap")
  p <- DimPlot(Seurat.data, reduction = "umap")
  ggsave("Cluster/Cluster_0.5.pdf", plot = p, width = 9, height = 8) 
  
# Save the ClusterData
## resolution = 0.5
  saveRDS(Seurat.data, file = "Cluster/clustering_0.5.rds")
  Seurat.data <- readRDS("Cluster/clustering_0.5.rds")
  
# Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
  #cluster2.markers <- FindMarkers(Seurat.data, ident.1 = 2, min.pct = 0.25)
  #head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
  #cluster5.markers <- FindMarkers(Seurat.data, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
  #head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
  Seurat.data.markers <- FindAllMarkers(Seurat.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  Seurat.data.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
  #VlnPlot(Seurat.data, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
  #VlnPlot(Seurat.data, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
  #FeaturePlot(Seurat.data, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
# DoHeatmap()为给定的细胞和特征生成一个表达热图。在这种情况下，我们为每个簇绘制前 20 个标记（或所有标记，如果小于 20）。
  Seurat.data.markers %>% group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC) -> top10
  DoHeatmap(Seurat.data, features = top10$gene) + NoLegend()

# Save the MarkersData
  dir.create("Markers")
  saveRDS(Seurat.data, file = "Markers/Markers_0.5.rds")
  write.csv(Seurat.data.markers,file = "Markers/Markers_0.5")
  
# Assigning cell type identity to clusters'
## Or use (SingleR.R) or CellTyping.R .etc  
## Seurat : Use canonical markers to easily match the unbiased clustering to known cell type
### Featureplot：according to the representative marker genes in order to assign cell type identity to clusters'
  FeaturePlot(Seurat.data, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
### when you make sure , assign cell type identity to clusters'
  dir.create("Cell_typing")
  new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                       "NK", "DC", "Platelet")
  names(new.cluster.ids) <- levels(Seurat.data)
  Seurat.data <- RenameIdents(Seurat.data, new.cluster.ids)
  p <- DimPlot(Seurat.data, reduction = "umap", label = TRUE, label.size = 8) + NoLegend()
  ggsave("Cell_typing/cell_typing_0.5.pdf", plot = p, width = 9, height = 8) 
# save the data  
  saveRDS(Seurat.data, file = "Seurat.data_final.rds")
  
####################################################################################
########### After Assigning cell type identity to clusters #########################
############### Select the cluster for the further analysis ########################  
################### e.g. monocle slingshot .etc ####################################  
# 读取包含细胞类型的数据
  Seurat.data <- readRDS("Monocle/Seurat.data_cell_typing.rds")
# Select the cluster for the further analysis  
# e.g. NK or cluster id
  ident_df <- data.frame(cell=names(Idents(Seurat.data)), cluster=Seurat.data$cluster)
  ## NK
  Seurat.data.subcluster <- subset(Seurat.data, cells=as.vector(ident_df[ident_df$cluster=="1",1]))
  ## cluster id
  Seurat.data.subcluster <- subset(Seurat.data, cells=as.vector(ident_df[ident_df$cluster=="2",1]))
  saveRDS(Seurat.data.subcluster, file = "Seurat.data.subcluster.rds")
  
# re-clustering(if you need)
  Seurat.data.subcluster <- RunPCA(Seurat.data.subcluster, features = VariableFeatures(object =Seurat.data.subcluster),npcs = 10)
  #Seurat.data.subcluster <- RunTSNE(Seurat.data.subcluster,dims = 1:15)
  Seurat.data.subcluster <- RunUMAP(Seurat.data.subcluster, reduction = "pca", dims = 1:10)
  Seurat.data.subcluster<- FindNeighbors(Seurat.data.subcluster,reduction = "pca", dims = 1:10)
  Seurat.data.subcluster<- FindClusters(Seurat.data.subcluster, resolution = 0.5) 
  
  p <- DimPlot(Seurat.data.subcluster, reduction = "umap")  
  ggsave("images/cluster2_reclustering.pdf", p, width = 8, eight = 6)
  
# Verify the cell type
# mouse
## macrophages
  FeaturePlot(Seurat.data,features = c("Csf1r","Cd14","Cd163"), reduction = "umap")
### M1
  FeaturePlot(Seurat.data,features = c("Cd80","Fcgr1","Cd86","Fcgr2b","Fcgr3","Cd40"), reduction = "umap")
### M2
  FeaturePlot(Seurat.data,features = c("Cd163","Mrc1","Ccl3"), reduction = "umap")

# human
## macrophages
  FeaturePlot(Seurat.data,features = c(), reduction = "umap")
### M1
  FeaturePlot(Seurat.data,features = c(), reduction = "umap")
### M2
  FeaturePlot(Seurat.data,features = c(), reduction = "umap")  

############### Select the celltype for the further analysis #######################
################### e.g. monocle slingshot .etc ####################################  
# Select the celltype for the further analysis
  ident_df <- data.frame(cell=names(Idents(Seurat.data)),celltype = Seurat.data@active.ident)
  part_celltype <- subset(Seurat.data, cells=as.vector(ident_df[ident_df$celltype=="NK"|ident_df$celltype=="Macrophages",1]))
####################################################################################


####################################################################################
############## The code of other Figures refer to Visualization.R ##################
####################################################################################
##################### Now You can do the monocle (Monocle.R) #######################
#################################################################################### 

##########################################################################################################################
