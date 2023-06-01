########################   Chronostasis   ############################
######################################################################
# Method: 1.Seurat; 2.Harmnoy; 3.LIGER
## 1. Seurat
# Enter commands in R (or R studio, if installed)
  #install.packages('Seurat')
# scRNA-seq integration can be used both to correct for technical differences between datasets 
# (i.e. batch effect correction), and to perform comparative scRNA-seq analysis of across experimental conditions
  rm(list = ls())
# Load R packages
  library(Seurat)
  library(tidyverse)
  library(patchwork)
# For multi- datasets/samples
# Load the outputs of cellranger count rather than cellranger aggr, and combine into one Seurat object
  in_data_dir <- "/Users/sth/"
  setwd(in_data_dir)
  dir.create("figures")
  dir.create("out_data")
  fig_dir <- "~/figures/"
  out_data_dir <- "~/out_data/"
# get the sample name  
  samples <- dir(paste0(in_data_dir,"/data"))
# get the sample group
  sampleGroup <- read.table("sampleGroup.txt",header=T.sep="\t")
# create seurat object:
  seurat_list <- lapply(samples, function(sample){
    cur_data <- Read10X(paste0(in_data_dir,"/data/",sample,'/outs/filtered_feature_bc_matrix/'))
    cur_seurat <- CreateSeuratObject(
      counts = cur_data,
      min.cells = 3,
      min.features = 200,
      project = sampleGroup$Group
    )
    cur_seurat$SampleID <- sample
    return(cur_seurat)
  })
# assign the sample name & sample group
  for (i in 1:length(seurat_list)){
    seurat_list[[i]]$group <- sampleGroup[i,2]
    print(i)
  }
  names(seurat_list) <- samples   
#############################################
############################################# 
########1.QC after merge multi-samples#######
############################################# 
############################################# 
# Preprocess(QC after merge multi-samples)
  dir.create(paste0(fig_dir,"/QC/","QC_after_merge"))
## merge seurat object
  seurat_obj <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])
## Calculate mito gene (ribosome gene or else) abundance percentage  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^MT-")
  
### If you need other quality control indicators,use the following code 
### Ribosome gene
  #seurat_obj[["percent.rb"]] <- PercentageFeatureSet(object = seurat_obj, pattern = "^RP[LS]")  
  ### Erythrocyte gene
  #HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  #HB.genes <- CaseMatch(HB.genes, rownames(seurat_list[[i]])) 
  #seurat_obj[["percent.HB"]] <- PercentageFeatureSet(object = seurat_obj, features = HB.genes)

## Visualize QC metrics as a violin plot   
  pdf(paste0(fig_dir,"/QC/QC_after_merge/violin_before_QC.pdf"), width=10, height=10)
  # png(paste0(fig_dir, "pngs/qc_violin_plot.png"), width=10, height=10, res=250, units='in')
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2, pt.size=0, )
  dev.off()
  
## Filter the outliers(according to the violinplot,or use the standard QC metrics)  
### use the standard QC metrics to preprocess   
  seurat_obj <- subset(seurat_obj, nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

## Visualize the QC result as a violin plot  
  pdf(paste0(fig_dir,"/QC/QC_after_merge/violin_after_QC.pdf"), width=10, height=10)
  # png(paste0(fig_dir, "pngs/qc_violin_plot.png"), width=10, height=10, res=250, units='in')
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 2, pt.size=0, )
  dev.off()
#############################################  
# Normalizing the data
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
# Identification of highly variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
# Run the standard workflow for visualization and clustering
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30)
  #seurat_obj <- RunTSNE(seurat_obj, npcs = 30, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5) 
# visualize the cluster before integrate
  dir.create(paste0(fig_dir,"Visualization"))
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident")
  ggsave(paste0(fig_dir,"Visualization/QC_after_merge_Before_integrate_cluster.pdf"), plot = p, width = 9, height = 8) 
#############################################

#############################################
############################################# 
######2.QC individually for each sample######
############################################# 
############################################# 
# Preprocess(QC individually for each sample)
  dir.create(paste0(fig_dir,"QC"))
  dir.create(paste0(fig_dir,"QC/","QC_for_each_sample"))
## Extract mitochondria gene list(some datasets have none mito genes,because of the cellranger preprocess)
  #rownames(seurat_list[[' ']]@assays$RNA@counts)[grep(pattern = "^MT",x = rownames(seurat_list[[' ']]@assays$RNA@counts))]
  ALL_mtgene <- lapply(samples, function(sample){
    grep(pattern = "^MT",x = rownames(seurat_list[[sample]]@assays$RNA@counts),value = TRUE)
  })
  names(ALL_mtgene) <- samples

## Calculate mito gene (ribosome gene or else) abundance percentage
  #seurat_list[[' ']][["percent.mt"]] <- PercentageFeatureSet(object = seurat_list[[' ']], pattern = "^MT-")
  for (i in samples){
    seurat_list[[i]][["percent.mt"]] <- PercentageFeatureSet(object = seurat_list[[i]], pattern = "^MT-")
    print(i)
  }
  
### If you need other quality control indicators,add the following to the for{}  
### Ribosome gene
  #seurat_list[[i]][["percent.rb"]] <- PercentageFeatureSet(object = seurat_list[[i]], pattern = "^RP[LS]")  
### Erythrocyte gene
  #HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  #HB.genes <- CaseMatch(HB.genes, rownames(seurat_list[[i]])) 
  #seurat_list[[i]][["percent.HB"]] <- PercentageFeatureSet(object = seurat_list[[i]], features = HB.genes)
  
## Visualize QC metrics as a violin plot
  plot.featrures = c("nFeature_RNA", "nCount_RNA","percent.mt") ###c("percent.rb", "percent.HB")
  for (i in samples){
    violin <- VlnPlot(object = seurat_list[[i]], features = plot.featrures, pt.size=0.1)
    ggsave(paste0(fig_dir,"/QC/","QC_for_each_sample/",'vlnplot_before_QC_',i,'.pdf') , plot = violin , width = 10, height = 8)
    print(i)
  }
  
## Visulization of feature-feature, Pearson correlation is on the top
  for (i in samples){
    plot1 <- FeatureScatter(object = seurat_list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(object = seurat_list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    plot <- plot1 + plot2
    ggsave(paste0(fig_dir,"/QC/","QC_for_each_sample/",'featurescatter_before_QC_',i,'.pdf') , plot = plot , width = 10, height = 8)
    print(i)
  }
  
## Filter the outliers(according to the violinplot,or use the standard QC metrics)  
### use the standard QC metrics to preprocess 
  for (i in samples){
    seurat_list[[i]] <- subset(x = seurat_list[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    print(i)
  }
### use customize QC metrics to preprocess   
  #seurat_list[['']] <- subset(x = seurat_list[['']], subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 5)

## Visualize the QC result as a violin plot  
  for (i in samples){
    violin <- VlnPlot(object = seurat_list[[i]], features = plot.featrures, pt.size=0.1)
    ggsave(paste0(fig_dir,"/QC/","QC_for_each_sample/",'vlnplot_after_QC_',i,'.pdf') , plot = violin , width = 10, height = 8)
    print(i)
  }
# merge seurat object
  seurat_obj <- merge(x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)])
#############################################  
# Run the standard workflow for visualization and clustering
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:30)
  #seurat_obj <- RunTSNE(seurat_obj, npcs = 30, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5) 
# visualize the cluster before integrate
  dir.create(paste0(fig_dir,"Visualization"))
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident")
  ggsave(paste0(fig_dir,"Visualization/QC_for_each_sample_Before_integrate_cluster.pdf"), plot = p, width = 9, height = 8) 


#############################################
############################################# 
#############################################
############################################# 
# clean up
  rm(seurat_list)
  gc()
# add metadata(when you have the file or you could make one if you need)
  #sample_metadata <- read.csv(file = "~/metaData.csv")
  #rownames(sample_metadata) <- sample_metadata$Sample.ID
  #cell_metadata <- sample_metadata[seurat_obj$SampleID,]
  #for(meta in names(cell_metadata)){
    #seurat_obj[[meta]] <- cell_metadata[[meta]]
  #}
#############################################################    
######################### Integrate #########################
#############################################################    
# split the dataset into a list of two seurat objects (EXPERIMENT and CTRL)
  seurat_obj.list <- SplitObject(seurat_obj, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
  seurat_obj.list <- lapply(X = seurat_obj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

# select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = seurat_obj.list)

  seurat_obj.anchors <- FindIntegrationAnchors(object.list = seurat_obj.list, anchor.features = features)

# this command creates an 'integrated' data assay
  Seurat.data.combined <- IntegrateData(anchorset = seurat_obj.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
  DefaultAssay(Seurat.data.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
  Seurat.data.combined <- ScaleData(Seurat.data.combined, verbose = FALSE)
  Seurat.data.combined <- RunPCA(Seurat.data.combined, npcs = 30, verbose = FALSE)
  Seurat.data.combined <- RunUMAP(Seurat.data.combined, reduction = "pca", dims = 1:30)
  #Seurat.data.combined <- RunTSNE(Seurat.data.combined, npcs = 30, verbose = FALSE)
  Seurat.data.combined <- FindNeighbors(Seurat.data.combined, reduction = "pca", dims = 1:30)
  Seurat.data.combined <- FindClusters(Seurat.data.combined, resolution = 0.6) 
  
# Visualization
  p1 <- DimPlot(Seurat.data.combined, reduction = "umap", group.by = "orig.ident")
  p2 <- DimPlot(Seurat.data.combined, reduction = "umap", split.by = "orig.ident")
  #p2 <- DimPlot(Seurat.data.combined, reduction = "umap", label = TRUE, repel = TRUE)
  ggsave(paste0(fig_dir,"Visualization/After_Integrate_cluster.pdf"), plot = p1, width = 9, height = 8) 
  ggsave(paste0(fig_dir,"Visualization/After_Integrate_cluster_splitby_sample.pdf"), plot = p2, width = 20, height = 18)
  
# Finding differentially expressed features (cluster biomarkers)
  DefaultAssay(Seurat.data.combined) <- "integrated"
  Seurat.data.combined.markers <- FindAllMarkers(Seurat.data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  Seurat.data.combined.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
## save the markers  
  saveRDS(Seurat.data.combined.markers,paste0(out_data_dir,"Seurat.data.combined.markers.rds"))
# save the Seurat object
  saveRDS(Seurat.data.combined,paste0(out_data_dir,"Seurat.data.combined.rds"))

####################################################################################
########### After Assigning cell type identity to clusters #########################
############### Select the cluster for the further analysis ########################  
################### e.g. monocle slingshot .etc ####################################  
# 读取包含细胞类型的数据
  Seurat.data.combined <- readRDS("Monocle/Seurat.data.combined_cell_typing.rds")
# Select the cluster for the further analysis  
# e.g. NK or cluster id
  ident_df <- data.frame(cell=names(Idents(Seurat.data.combined)), cluster=Seurat.data.combined$cluster)
## NK
  Seurat.data.combined.subcluster <- subset(Seurat.data.combined, cells=as.vector(ident_df[ident_df$cluster=="1",1]))
## cluster id
  Seurat.data.combined.subcluster <- subset(Seurat.data.combined, cells=as.vector(ident_df[ident_df$cluster=="2",1]))
  saveRDS(Seurat.data.combined.subcluster, file = paste0(out_data_dir,"Seurat.data.combined.subcluster.rds"))
  
# re-clustering(if you need)
  dir.create(paste0(fig_dir,"Re-clustering"))
  Seurat.data.combined.subcluster <- RunPCA(Seurat.data.combined.subcluster, features = VariableFeatures(object =Seurat.data.combined.subcluster),npcs = 10)
  #Seurat.data.combined.subcluster <- RunTSNE(Seurat.data.combined.subcluster,dims = 1:15)
  Seurat.data.combined.subcluster <- RunUMAP(Seurat.data.combined.subcluster, reduction = "pca", dims = 1:10)
  Seurat.data.combined.subcluster<- FindNeighbors(Seurat.data.combined.subcluster,reduction = "pca", dims = 1:10)
  Seurat.data.combined.subcluster<- FindClusters(Seurat.data.combined.subcluster, resolution = 0.5) 
  
  p <- DimPlot(Seurat.data.combined.subcluster, reduction = "umap")  
  ggsave(paste0(fig_dir,"Re-clustering/cluster2_reclustering.pdf"), p, width = 8, eight = 6)
  
# Verify the cell type
  DefaultAssay(Seurat.data.combined) <- "integrated"
# mouse
## macrophages
  FeaturePlot(Seurat.data.combined,features = c("Csf1r","Cd14","Cd163"), reduction = "umap")
### M1
  FeaturePlot(Seurat.data.combined,features = c("Cd80","Fcgr1","Cd86","Fcgr2b","Fcgr3","Cd40"), reduction = "umap")
### M2
  FeaturePlot(Seurat.data.combined,features = c("Cd163","Mrc1","Ccl3"), reduction = "umap")
  
# human
## macrophages
  FeaturePlot(Seurat.data.combined,features = c(), reduction = "umap")
### M1
  FeaturePlot(Seurat.data.combined,features = c(), reduction = "umap")
### M2
  FeaturePlot(Seurat.data.combined,features = c(), reduction = "umap")  
  
############### Select the celltype for the further analysis #######################
################### e.g. monocle slingshot .etc ####################################  
# Select the celltype for the further analysis
  ident_df <- data.frame(cell=names(Idents(Seurat.data.combined)),celltype = Seurat.data.combined@active.ident)
  part_celltype <- subset(Seurat.data.combined, cells=as.vector(ident_df[ident_df$celltype=="NK"|ident_df$celltype=="Macrophages",1]))
  
####################################################################################
####################################################################################
############## The code of other Figures refer to Visualization.R ##################
####################################################################################
##################### Now You can do the monocle (Monocle.R) #######################
#################################################################################### 
  
# Identify differential expressed genes across conditions 
## Seurat: e.g.(label : stim == stimulated)
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())
###  start to do comparative analyses and look at the differences induced by stimulation
### One way to look broadly at these changes is to plot the average expression of both the stimulated and control cells and look for genes that are visual outliers on a scatter plot.
### select the CD4 Naive T cells    
  t.cells <- subset(Seurat.data.combined, idents = "CD4 Naive T")
  Idents(t.cells) <- "stim"
  avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
  avg.t.cells$gene <- rownames(avg.t.cells)
### select the CD14 Mono cells    
  cd14.mono <- subset(Seurat.data.combined, idents = "CD14 Mono")
  Idents(cd14.mono) <- "stim"
  avg.cd14.mono <- as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))
  avg.cd14.mono$gene <- rownames(avg.cd14.mono)

### Here, we take the average expression of both the stimulated and control naive T cells and CD14 monocyte populations and generate the scatter plots, 
### highlighting genes that exhibit dramatic responses to interferon stimulation.  
  genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
  p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
  p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
  p1 + p2
  
### identified common cell types across condition ,and mark them
  Seurat.data.combined$celltype.stim <- paste(Idents(Seurat.data.combined), Seurat.data.combined$stim, sep = "_")
  Seurat.data.combined$celltype <- Idents(Seurat.data.combined)
  Idents(Seurat.data.combined) <- "celltype.stim"
### use FindMarkers() to find the genes that are different between stimulated and control B cells 
  b.interferon.response <- FindMarkers(Seurat.data.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
  head(b.interferon.response, n = 15)

### visualize the result    
  FeaturePlot(immune.combined, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3,
              cols = c("grey", "red"))
  
  plots <- VlnPlot(immune.combined, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "celltype",
                   pt.size = 0, combine = FALSE)
  wrap_plots(plots = plots, ncol = 1)
  
##########################################################################################################################  
#############################################
############################################# 
##############Other Methods##################
#############################################
############################################# 
## 2. Harmnoy (From the Seurat V3 of Harmony)[https://portals.broadinstitute.org/harmony/SeuratV3.html]
### Install the packages
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("harmony")
### Load the packages
  library(harmony)
  dir.create(paste0(fig_dir,"Harmony"))
### Load the data(same as before) [After QC]  
  scRNA_harmony<- seurat_obj
### Before running Harmony, make a Seurat object and following the standard pipeline through PCA.
  scRNA_harmony <- scRNA_harmony %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = scRNA_harmony@var.genes, npcs = 20, verbose = FALSE)
### Or use the SCTransform
  #sce.all <- SCTransform(seu_obj)%>% RunPCA(verbose=FALSE)
##### Before run harmony  
  options(repr.plot.height = 5, repr.plot.width = 12)
  p1 <- DimPlot(object = scRNA_harmony, reduction = "pca", pt.size = .1, group.by = "orig.ident")
  p2 <- VlnPlot(object = scRNA_harmony, features = "PC_1", group.by = "orig.ident", pt.size = .1)
  p3 <- plot_grid(p1,p2)
  ggsave(paste0(fig_dir,"Harmony/Before_Harmony.pdf"), p3, width = 9, eight = 8)
### integrate the dataset by harmony 
#### set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round
  scRNA_harmony <- scRNA_harmony %>% RunHarmony("orig.ident", plot_convergence = T)
### Check the generated embeddings
  harmony_embeddings <- Embeddings(scRNA_harmony, 'harmony')
  harmony_embeddings[1:5, 1:5]
### After run harmony  
  options(repr.plot.height = 5, repr.plot.width = 12)
  p1 <- DimPlot(object = scRNA_harmony, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
  p2 <- VlnPlot(object = scRNA_harmony, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
  p3 <- plot_grid(p1,p2)
  ggsave(paste0(fig_dir,"Harmony/After_Harmony.pdf"), p3, width = 9, eight = 8)

### Downstream analysis    
#### To use the corrected Harmony embeddings rather than PCs, set reduction = 'harmony'  
  scRNA_harmony <- scRNA_harmony %>% 
    RunUMAP(reduction = "harmony", dims = 1:20) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
    FindClusters(resolution = 0.5) %>% 
    identity()
### Visualization(similar to the Seurat)
### When save the result,Remember to change the directory name
  options(repr.plot.height = 4, repr.plot.width = 10)
  DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident", pt.size = .1, split.by = 'orig.ident')
  
  options(repr.plot.height = 4, repr.plot.width = 6)
  DimPlot(scRNA_harmony, reduction = "umap", label = TRUE, pt.size = .1)
  
#### plot the distribution among clusters [https://github.com/hemberg-lab/scRNA.seq.course/blob/master/course_files/utils/custom_seurat_functions.R]
  plot_integrated_clusters = function (srat) { 
    ## take an integrated Seurat object, plot distributions over orig.ident
    library(Seurat)
    library(patchwork)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    
    count_table <- table(srat@meta.data$seurat_clusters, srat@meta.data$orig.ident)
    count_mtx   <- as.data.frame.matrix(count_table)
    count_mtx$cluster <- rownames(count_mtx)
    melt_mtx    <- melt(count_mtx)
    melt_mtx$cluster <- as.factor(melt_mtx$cluster)
    
    cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)
    
    sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
    cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
    melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)
    
    colnames(melt_mtx)[2] <- "dataset"
    
    
    p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") + 
      theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("")
    p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) + 
      geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() + 
      scale_fill_brewer(palette = "Set2") +
      ylab("Fraction of cells in each dataset") + xlab("Cluster number") + theme(legend.position="top")
    
    p2 + p1 + plot_layout(widths = c(3,1))
    
  }
  
  plot_integrated_clusters(scRNA_harmony) 
  
  
#############################################
############################################# 
############# Other Methods #################
#############################################
############################################# 
## 3. LIGER [https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/liger.html] 
### Install the packages
  install.packages('rliger')
  remotes::install_github('satijalab/seurat-wrappers')
  library(rliger)
  library(SeuratWrappers)
### Load the data(same as before) [After QC]  
  scRNA_liger <- seurat_obj
#### normalize expression data to account for differences in sequencing depth and efficiency between cells
#### identify variably expressed genes and scale the data so that each gene has the same variance
  scRNA_liger <- NormalizeData(scRNA_liger)
  scRNA_liger <- FindVariableFeatures(scRNA_liger)
  scRNA_liger <- ScaleData(scRNA_liger, split.by = "orig.ident", do.center = F)
### Joint Matrix Factorization  
  scRNA_liger <- RunOptimizeALS(scRNA_liger, k = 30, lambda = 5, split.by = "orig.ident") ## this one takes a while
### Quantile Normalization and Joint Clustering  
  scRNA_liger <- RunQuantileNorm(scRNA_liger, split.by = "stim")
### You can optionally perform Louvain clustering (`FindNeighbors` and `FindClusters`) after
### `RunQuantileNorm` according to your needs
  scRNA_liger <- FindNeighbors(scRNA_liger, reduction = "iNMF", dims = 1:20)
  scRNA_liger <- FindClusters(scRNA_liger, resolution = 0.55)
### Dimensional reduction and plotting
  scRNA_liger <- RunUMAP(scRNA_liger, dims = 1:ncol(scRNA_liger[["iNMF"]]), reduction = "iNMF")
  scRNA_liger    <- SetIdent(scRNA_liger,value = "orig.ident")
  DimPlot(scRNA_liger, group.by = c("orig.ident", "seurat_annotations"), ncol = 3)
  
  
##########################################################################################################################    
  
  
  
  
