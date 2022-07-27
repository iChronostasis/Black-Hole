########################   Chronostasis   #######################################
#################################################################################
# Some Figures(1,2,3,4,5一般用marker gene进行画图，可以确定细胞分型正确与否；6，7在细胞注释之后画图)

### Multi-plots
  library(cowplot)
  pdf("markers1234.pdf", width = 18, height =  12)
  plot_grid(p1,p2,p3,p4,align = "v", nrow = 3,labels =c("cluster1","cluster2", "cluster3","cluster4"))
  dev.off()

### 1.Featureplot：可视化关注的marker gene所富集的簇，可以通过每一簇区别于其他簇的marker gene验证后续的细胞注释的结果
  FeaturePlot(Seurat.data, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
#### more color(rainbow)  
  FeaturePlot(Seurat.data,c("MS4A1","LYZ","NKG7","PPBP","LTF","HBA1","FCER1A","IL7R","FCGR3B")) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
  
### 2. dotplot 气泡图
  Idents(Seurat.data) <- factor(Idents(Seurat.data), levels = c("Monocytes","Stellate cells","Myofibroblasts", "T cells", "Plasma cells"))
  markers.to.plot <- c("SPP1","FABP5","SELENOP","COL1A1","COL3A1","COL1A2","SPARC","RPL23P2","LINC02154","CD2","TRAC","MZB1","IGHG1","JCHAIN","XBP1")
  marker.to.plot  <- unique(markers.to.plot)
  DotPlot(Seurat.data, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
    RotatedAxis()
  
### 3. 每个簇的基因表达值的热图展示 Heatmap for gene expression level per cluster
  cluster.averages.Seurat.data <- AverageExpression(Seurat.data, return.seurat = TRUE)
  head(cluster.averages.Seurat.data[["RNA"]][, 1:5])
  
  pdf("pbm_newheatmap.pdf", width = 4, height = 5)
  DoHeatmap(cluster.averages.Seurat.data, features = unlist(TopFeatures(Seurat.data[["pca"]], balanced = TRUE)), label=F, size = 3)
  dev.off()

### 4. VlnPlot 小提琴图
  features.plot.markers = c("GPX3","MIOX","CD3E","CD86")
  pdf("violin_featrues.pdf", width = 12, height = 5)
  VlnPlot(object =Seurat.data,features = features.plot.markers,pt.size=0, ncol = 4)
  dev.off()

### 5. RidgePlot 山脊图
  features.plot.markers = c("GPX3","MIOX","CD3E","CD86")
  pdf("Ridge rb_featrues.pdf", width = 8, height =8 )
  RidgePlot(Seurat.data, features = features.plot.markers, ncol= 2)
  dev.off()

### 6. 各类细胞组成成分图
  cell.prop <- as.data.frame(table(Idents(Seurat.data)))
  colnames(cell.prop)<-c("cluster","proportion")
  library(ggplot2)
  ggplot(cell.prop,aes(label,proportion,fill = cluster))+
    geom_bar(stat="identity",position="fill")+
    scale_fill_brewer(palette = "Set1")+
    guides(fill=guide_legend(ncol=2))+
    coord_flip()+
    ggtitle("")+
    theme_classic() +
    theme(axis.text = element_text(size = 15),axis.title= element_text(size = 20))
  
### 7. 冲击图（细胞组成成分图）
  library(RColorBrewer)
  library(ggalluvial) 
  
  cell.prop<-as.data.frame(prop.table(table(Seurat.data@meta.data$celltype,Seurat.data@meta.data$orig.ident)))
  colnames(cell.prop)<-c("celltype","orig.ident","proportion")
  
  colourCount = length(unique(cell.prop$celltype))
  colourCount
  #9
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  
  p<-ggplot(cell.prop,aes(orig.ident,proportion,fill=celltype,alluvium=celltype,stratum=celltype))+
    geom_bar(stat="identity",position="fill",width=1/2)+
    geom_stratum(width=1/2,size=0,colour="grey")+geom_alluvium(aes(fill=celltype),width=1/2,curve_type="linear")+
    scale_fill_manual(values = getPalette(colourCount))+ 
    ggtitle(paste0("cell proportion of ",name,""))+
    xlab("")+
    ylab("")+
    theme_bw()+
    guides(fill=guide_legend(title="celltype"))+
    theme( axis.title=element_text(size=11,face="bold"),
           axis.ticks.length=unit(0.2,'cm'), #刻度长短
           axis.ticks.y=element_blank(),
           axis.text.y=element_text(colour="black",size=10,face="bold"),
           axis.line.y=element_blank(),
           axis.text.x=element_text(colour="black",size=10,face="bold"), 
           axis.title.x=element_text(size = 12,face="bold"),
           panel.border = element_blank(),
           axis.line = element_line(colour = "black",size=0.8), 
           legend.text=element_text(colour="black",size=12),  
           legend.title=element_text(colour="black", size=12,face="bold"),
           #legend.position = c(0.9,0.85),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           plot.title = element_text(hjust = 0.5))+
    coord_flip()+
    scale_y_continuous(limits =c(0, 1) ,expand = c(0,0))+
    scale_x_discrete(expand=c(0,0.3)) #控制x轴到图形的距离
  #+geom_text(aes(label = Frequency), size = 5, hjust = 0.5, vjust = 1, position = "stack")  
  
  ggsave(paste0("",name,"_group_celltype_proportion2.png"), plot = p, width = 9, height = 3)


### 7.1 细胞成分图（堆积）
  ggplot(cell.prop,aes(location,prop_each_loc,fill=cluster,alluvium=cluster))+
  geom_bar(stat="identity",width=0.5)+
  geom_alluvium()+
  labs(x="Menstrual_Phase",y="celltype", title = "cell_prop_in_diff_Menstrual_Phase", color = "") +
  scale_fill_d3("category20")+theme_bw() +
  theme(panel.grid.major=element_line(colour=NA),#图片背景#
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 14,color="black"),#X轴标度字体#
    axis.text.y = element_text(size = 14,color="black"),#Y轴标度字体#
    axis.title.x = element_text(size = 14),#x轴名字字体#
    axis.title.y = element_text(size = 14),#y轴名字字体#
    legend.text = element_text(size = 14))+#图例字体的大小#
  theme(axis.ticks.length=unit(0.5,'cm'))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill=guide_legend(override.aes = list(size=4)))


### 8. Find particular type of interest and highlighting ----(from nbisweden.github.io)
  findCells <- function(obj, column, values, name=NULL) {
    stopifnot(is(obj, "Seurat"))
    stopifnot(is.character(column))
    stopifnot(column %in% names(obj@meta.data))
    col <- obj@meta.data[[column]]
    stopifnot(is.character(col) || is.factor(col))
    values <- unique(values)
    stopifnot(is.character(values) || is.factor(values))
    if (length(values)>1 && is.null(name))
      stop("findCells: specify a name to be used for the selection")
    if(is.null(name))
      name <- values
    stopifnot(is.character(name))
    rem <- setdiff(c(values), col)
    if(length(rem)>0)stop("findCells: requested value(s) never occurs in this column: ", rem)
    l <- list(colnames(obj)[ col %in% values ])
    names(l) <- name
    l
  }                                       #findCells
  
##########################################################
## e.g.
## Highlighting: compare T cells in hpca and bpe groups
  p_t1 <- DimPlot(Seurat.data,group.by = "bpe_type",reduction = "umap",
                  cells.highlight = findCells(Seurat.data,"bpe_type",c("CD8+ T-cells","CD4+ T-cells"), name = "T cells"))+
    theme(plot.margin = unit(c(3,0,1,0),'lines'))
  p_t2 <- DimPlot(Seurat.data,group.by = "hpca_type",reduction = "umap",
                  cells.highlight = findCells(Seurat.data,"hpca_type","T_cells"))+
    theme(plot.margin = unit(c(3,0,1,0),'lines'))
  pdf("findcells_Tcells.pdf", width = 14, height = 6)
  plot_grid(p_t2,p_t1, nrow = 1,labels = c("hpca_T cells","bpe_T cells"),label_size= 18, hjust = -1)
  dev.off()
  ## find intersected NK cells
  Seurat.data@meta.data[which(Seurat.data@meta.data$hpca_type == "NK_cell" & Seurat.data@meta.data$bpe_type == "NK cells"),]
  

#### 9. plot the distribution among clusters (integrate the datasets)
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

  plot_integrated_clusters(Seurat.data)
  plot_integrated_clusters(scRNA_harmony) 
#################################################################################

