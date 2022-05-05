########################   Chronostasis   ############################
######################################################################
# Load the R packages
  library(monocle)
  library(Seurat)
  library(ggplot2)
  library(ggsci)
# 1. create the directory
  dir.create("Monocle")  
# 2. create the CellDataSet(including the cell type or part celltype)
  Seurat.data <- readRDS("")
  
  data <- as(as.matrix(Seurat.data@assays$RNA@data), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data = Seurat.data@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
  cds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())

# 估计尺寸因子和分散
  cds <- estimateSizeFactors(cds)
  cds<- estimateDispersions(cds)

##选择排序基因
### 选择seurat确定的高变基因
  order.genes <- VariableFeatures(Seurat.data.partcelltype)
  saveRDS(order.genes, "Monocle/order.genes.seurat.rds")
# 设置排序基因
  mycds <- setOrderingFilter(cds, order.genes)
  #mycds <- setOrderingFilter(mycds, order.genes.adj)  #如果调整了排序基因，用这行代码。
  p <- plot_ordering_genes(mycds)
  ggsave("Monocle/order.genes.seurat.pdf", p, width = 8, height = 6)

### monocle选择高变基因
  disp_table <- dispersionTable(mycds)
  order.genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit) %>% 
    pull(gene_id) %>% as.character()
  mycds <- setOrderingFilter(mycds, order.genes)
  p <- plot_ordering_genes(mycds)
  ggsave("Monocle/order.genes.monocle.pdf", p, width = 8, height = 6)

##降维
  mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')

##排序
  mycds <- orderCells(mycds)
  mycds <- orderCells(mycds, root_state = 7)  #调整排序时用
##结果可视化
#State轨迹分布图
  plot1 <- plot_cell_trajectory(mycds, color_by = "State")
  ggsave("Monocle/Trajectory_State.pdf", plot = plot1, width = 10, height = 6.5)

#Celltype轨迹分布图
  plot2 <- plot_cell_trajectory(mycds, color_by = "cell_type")
  ggsave("Monocle/Trajectory_Celltype.pdf", plot = plot2, width = 10, height = 6.5)

#Pseudotime轨迹图
  plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
  ggsave("Monocle/Trajectory_Pseudotime.pdf", plot = plot3, width = 10, height = 6.5)

#合并作图
  plotc <- plot1|plot2|plot3
  ggsave("Monocle/Trajectory_Combination.pdf", plot = plotc, width = 10, height = 3.5)

#分面图
  p <- plot_cell_trajectory(mycds, color_by = "cell_type") + facet_wrap(~cell_type, nrow = 3)
  ggsave("Monocle/Trajectory_Facet_cell_type.pdf", plot = p, width = 10, height = 10)

  p <- plot_cell_trajectory(mycds, color_by = "State") +facet_wrap(~State, nrow = 1)
  ggsave("Monocle/Trajectory_Facet_State.pdf", plot = p, width = 10, height = 10)
# Save the data
  saveRDS(mycds,"Monocle/mycds_NK.rds")  