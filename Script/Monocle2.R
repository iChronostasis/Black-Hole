########################   Chronostasis   ############################
######################################################################
# Install the packages
  source("http://bioconductor.org/biocLite.R")
  biocLite()  
  biocLite("monocle")
# Load the R packages
  library(monocle)
  library(Seurat)
  library(ggplot2)
  library(ggsci)
# 1. create the directory
  dir.create("Monocle")  
# 2. create the CellDataSet(including the cell type or part celltype)
  Seurat.data <- readRDS("")
  
  expr_matrix <- as.matrix(Seurat.data@assays$RNA@data)
  sample_sheet <- data.frame(Seurat.data@meta.data)
  #sample_sheet <- data.frame(row.names = colnames(Seurat.data),
                             #id = colnames(Seurat.data),
                             #Time = Seurat.data@meta.data$Time,
                             #cell_type = Seurat.data$cell_type,
                             #stringsAsFactors = FALSE)
  gene_annotation <- data.frame(row.names = rownames(Seurat.data),
                                id = rownames(Seurat.data),
                                gene_short_name = rownames(Seurat.data),
                                stringsAsFactors = FALSE)
  
  pd <- new("AnnotatedDataFrame", data = sample_sheet)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
# Construct monocle cds  
  monocle <- newCellDataSet(expr_matrix,  
                            phenoData = pd, 
                            featureData = fd,
                            #lowerDetectionLimit=0.1,
                            expressionFamily = negbinomial() )
  #monocle$Time <- factor(monocle$Time,levels = c("0day","1day","3day","5day","7day","12day"))
# 估计尺寸因子和分散
  monocle <- estimateSizeFactors(monocle)
  monocle <- estimateDispersions(monocle)

## The distribution of mRNA totals across the cells
  pData(monocle)$Total_mRNAs <- Matrix::colSums(exprs(monocle))
  monocle <- monocle[,pData(monocle)$Total_mRNAs < 1e6]
  upper_bound <- 10^(mean(log10(pData(monocle)$Total_mRNAs)) + 2*sd(log10(pData(monocle)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(pData(monocle)$Total_mRNAs)) - 2*sd(log10(pData(monocle)$Total_mRNAs)))
  
  pdf('Monocle/Total_mRNAs.pdf')
  qplot(Total_mRNAs, data = pData(monocle),color = Seurat.data$Time, geom = "density" ) + 
    geom_vline(xintercept = lower_bound) +
    geom_vline(xintercept = upper_bound)
  dev.off()
  
# Filtering low-quality cells
  monocle <- detectGenes(monocle, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(monocle),
                                      num_cells_expressed >= 10))
  length(expressed_genes)
  
  
########################################################  
### monocle选择高变基因 [无marker gene时]
######################################################## 
  disp_table <- dispersionTable(monocle)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.01 & dispersion_empirical >= 1 * dispersion_fit) %>% 
    pull(gene_id) %>% as.character()
  monocle <- setOrderingFilter(monocle, unsup_clustering_genes)
  pdf("Monocle/select_genes_monocle.pdf")
  plot_ordering_genes(monocle)
  dev.off()
  
######################################################## 
### 选择seurat确定的高变基因
######################################################## 
  #unsup_clustering_genes <- VariableFeatures(Seurat.data.partcelltype)
### Ordering based on genes that differ between clusters 
  #monocle <- detectGenes(monocle, min_expr = 0.1)
  #fData(monocle)$use_for_ordering <-
    #fData(monocle)$num_cells_expressed > 0.05 * ncol(monocle)
  #plot_pc_variance_explained(monocle, return_all = F)
  
########################################################  
# choose genes that define a cell's progress[dpfeature]
######################################################## 
## One effective way to isolate a set of ordering genes is to simply compare the cells collected 
## at the beginning of the process to those at the end and find the differentially expressed genes
  diff_test_res <- differentialGeneTest(monocle[expressed_genes,],
                                        fullModelFormulaStr = "~Time")
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  monocle <- setOrderingFilter(monocle, ordering_genes)
  pdf("Monocle/select_genes_dpfeature.pdf")
  plot_ordering_genes(monocle)
  dev.off()
########################################################  
## 降维 
  monocle <- reduceDimension(monocle,max_components = 2,method = 'DDRTree')
## 排序
  monocle <- orderCells(monocle, reverse = TRUE)
  
  monocle <- orderCells(monocle, root_state = 7)  #调整排序时用
  
### The function below is handy for identifying the State which contains most of the cells from time zero  
  monocle_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
      T0_counts <- table(pData(cds)$State, pData(cds)$Time)[,"0day"]
      return(as.numeric(names(T0_counts)[which
                                         (T0_counts == max(T0_counts))]))
    } else {
      return (1)
    }
  }
  monocle <- orderCells(monocle, root_state = monocle_state(monocle)) 
## 结果可视化
# State轨迹分布图
  plot1 <- plot_cell_trajectory(monocle, color_by = "State")
  ggsave("Monocle/Trajectory_State.pdf", plot = plot1, width = 10, height = 6.5)
  
# Celltype轨迹分布图
  plot2 <- plot_cell_trajectory(monocle, color_by = "cell_type")
  ggsave("Monocle/Trajectory_Celltype.pdf", plot = plot2, width = 10, height = 6.5)
  
# Pseudotime轨迹图
  plot3 <- plot_cell_trajectory(monocle, color_by = "Pseudotime")
  ggsave("Monocle/Trajectory_Pseudotime.pdf", plot = plot3, width = 10, height = 6.5)

# Time轨迹图    
  plot4 <- plot_cell_trajectory(monocle, color_by = "Time")
  ggsave("Monocle/Trajectory_Time.pdf", plot = plot4, width = 10, height = 6.5)
  
# 合并作图
  plotc <- plot1|plot2|plot3|plot4
  ggsave("Monocle/Trajectory_Combination.pdf", plot = plotc, width = 10, height = 3.5)
  
# 分面图
  p <- plot_cell_trajectory(monocle, color_by = "cell_type") + facet_wrap(~Time, nrow = 3)
  ggsave("Monocle/Trajectory_Facet_cell_type_Time.pdf", plot = p, width = 10, height = 10)
  
# Save the data
  saveRDS(monocle,"Monocle/monocle.rds")    
  
  
## if you don't have a timeseries, you might need to set the root based on where certain marker genes are expressed
## using your biological knowledge of the system
  blast_genes <- row.names(subset(fData(monocle),
                                  gene_short_name %in% c("CCNB2", "MYOD1", "MYOG")))
  plot_genes_jitter(monocle[blast_genes,],
                    grouping = "State",
                    min_expr = 0.1)
  
  expressed_genes <-  row.names(subset(fData(monocle),
                                            num_cells_expressed >= 10))
  filtered <- monocle[expressed_genes,]
  my_genes <- row.names(subset(fData(filtered),
                               gene_short_name %in% c("Cd27","Cd244a","Kit","Il2rb","Il7r","Klrk1","Klrb1c","Ncr1","Klra","Itga2","Klrc1","Itgav","Itgax","Itgam","Spn","Klrg1")))
  cds_subset <- filtered[my_genes,]
  plot_genes_in_pseudotime(cds_subset, color_by = "Time")


### Finding Genes that Distinguish Cell Type or State http://cole-trapnell-lab.github.io/monocle-release/images/vignette/plot_diff_res_multi-1.png
  to_be_tested <- row.names(subset(fData(monocle),
                                   gene_short_name %in% c("UBC", "NCAM1", "ANPEP")))
  cds_subset <- monocle[to_be_tested,]
  diff_test_res <- differentialGeneTest(cds_subset,
                                        fullModelFormulaStr = "~cellType")
  diff_test_res[,c("gene_short_name", "pval", "qval")]
  plot_genes_jitter(cds_subset,
                    grouping = "cellType",
                    color_by = "cellType",
                    nrow= 1,
                    ncol = NULL,
                    plot_trend = TRUE)
  
  
### Visualize the expression changes of each cluster under pseudo-time (ridge plot)
  library(ggridges)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggpubr)
  df <- pData(cds) ## pData(cds)取出的是cds对象中cds@phenoData@data的内容
  p4<-ggplot(df, aes(Pseudotime, y = cluster2, fill=cluster2)) +
    geom_density_ridges(alpha=0.6,bins=20) +
    geom_vline(xintercept = c(0,5,10),linetype=2)+
    theme_ridges() + 
    theme(legend.position="none",panel.spacing = unit(0.1, "lines"),strip.text.x = element_text(size = 8))+
    theme(
      panel.grid = element_blank()
    )+
    scale_fill_d3("category20")+
    xlab("Pseudotime")
  ggsave("cell_trajectory_T_density.pdf", p4, width = 10, height = 10)

### Finding Genes that Change as a Function of Pseudotime  
  to_be_tested <- row.names(subset(fData(monocle),
                                   gene_short_name %in% c("MYH3", "MEF2C", "CCNB2", "TNNT1")))
  cds_subset <- monocle[to_be_tested,]
  diff_test_res <- differentialGeneTest(cds_subset,
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  diff_test_res[,c("gene_short_name", "pval", "qval")]
  plot_genes_in_pseudotime(cds_subset, color_by = "Time")
  
  
## Analyzing Branches in Single-Cell Trajectories
  BEAM_res <- BEAM(monocle, branch_point = 1, cores = 1)
  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  
  plot_genes_branched_heatmap(monocle[row.names(subset(BEAM_res,
                                                    qval < 1e-4)),],
                              branch_point = 1,
                              num_clusters = 4,
                              cores = 1,
                              use_gene_short_name = T,
                              show_rownames = T)
  
  monocle_genes <- row.names(subset(fData(monocle),
                                 gene_short_name %in% c("Ccnd2", "Sftpb", "Pdpn")))
  plot_genes_branched_pseudotime(lung[lung_genes,],
                                 branch_point = 1,
                                 color_by = "Time",
                                 ncol = 1)
######################################################################  
  
