celltype <- read.table('/Users/hecate/研一/Term 2/2022.03.21/PanglaoDB_markers_27_Mar_2020.tsv',
                       sep = '\t',header = T,stringsAsFactors = F)

celltype <- read.table('/home/chengxin/scRNA_seq/PanglaoDB_markers_27_Mar_2020.tsv',
                       sep = '\t',header = T,stringsAsFactors = F)


iHyperGD <- function(backGround,
                     geneNum.pathway,
                     geneNum.check,
                     geneNum.located){
  1-phyper(geneNum.located-1,
           geneNum.pathway,
           backGround - geneNum.pathway,
           geneNum.check,
           lower.tail = T)
}

hyperType <- function(genes,type,backGround){
  geneNum.pathway <- length(which(celltype$cell.type == type))
  geneNum.check   <- length(genes)
  geneNum.located <- sum(genes %in% celltype$official.gene.symbol[which(celltype$cell.type == type)])
  
  iHyperGD(backGround,geneNum.pathway,geneNum.check,geneNum.located)
}

cellTypeEst_hyper <- function(cluster,m,backGround = 2000,type = NA){
  genes <- m$gene[which(m$cluster == cluster)]
  
  p.val <- c()
  if(is.na(type)){
    for (Type in unique(celltype$cell.type)) {
      p.val <- c(p.val,hyperType(genes,Type,backGround))
    }
  }
  names(p.val) <- unique(celltype$cell.type)
  sort(p.val,decreasing = F)
}
# 一次一个分簇
cellTypeEst_hyper(2,m=pbmc.markers) 
for (i in (0:26)){
  print(cellTypeEst_hyper(i,m=pbmc.markers)[1])
}

