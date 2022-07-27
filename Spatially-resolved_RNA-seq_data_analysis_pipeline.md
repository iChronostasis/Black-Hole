# Spatially-resolved RNA-seq data analysis pipeline

<!-- TOC -->

- [Spatially-resolved RNA-seq data analysis pipeline](#spatially-resolved-rna-seq-data-analysis-pipeline)
  - [上游分析](#上游分析)
    - [Space Ranger使用步骤](#space-ranger使用步骤)
  - [下游分析](#下游分析)
    - [10x Visium](#10x-visium)
    - [Slide-seq](#slide-seq)
  - [参考](#参考)

<!-- /TOC -->

## 上游分析

1. 取样，实验设计，测序

   ![img](https://www.genenergy.cn/upload/images/2020/7/1717927829.jpg)

2. Space ranger软件处理测序文件

首先，space ranger将前面产生的fastq测序数据比对到参考基因组上，进行基因表达的定量，即**UMI计数**。接着，依靠图像处理算法确定组织位置，即**spots定量**。这样就可以得到一个spots和基因的表达矩阵。最后，基于这个表达矩阵，进行下一步**降维，聚类及差异分析**。

主要分为四步：

**第一步基因组比对**，采用star将reads比对到参考基因组上，根据比对位置对reads进行分类：分为外显子区，内含子区及基因间隔区。

**第二步MAPQ调整**：对于有些reads同时比对到外显子区及其它区域，则优先认为比对到外显子；若该reads在外显子区域的MAPQ为255，则更可信

**第三步转录组比对**：进一步将比对到外显子区的reads与己知转录本进行比对，如果比对上并且链相同，则认为比对上该转录本。只有unique比对到转录本的reads才会作为UMI计数

**第四步UMI计数**：主要有两个校正过程(1)对UMI碱基进行校正：主要校正UMI中的测序错误(2)对基因进行校正：保留多数reads支持的基因注释

**Space Ranger 如何识别组织区域？**

主要解决两个关键问题：一是基准点对齐，二是确定组织位置；

图中红色点就是基准点，这些基准点组成一个基准框，基准框的角和边都是独特的，左上角是沙漏型，右上角是实心六边形；右下角是空心六边形，左下角是三角形。

首先提取“看起来”像基准点的spots，将这些候选基准点与已知的基准点模式对齐，识别出基准线独特的角和边，坐标转换，将spatial barcode与组织图像联系起来。

接着，计算和比较组织切片放置的多个估计值，这些估计值用于训练分类器，将捕获区域内的每个像素标记为组织还是背景。这样就可以识别出组织区域。

![img](https://www.genenergy.cn/upload/images/2020/7/1717834798.jpg)

也就是说每个捕获区域最多捕获4992个spots。拿两个spots举例说明一下：行坐标是4，偶数行，第一个spot的列坐标0，2，4，这个spot的坐标为4，4这个spot行坐标是5，奇数行，第一个spot的列坐标从1开始，1，3，5，7，坐标信息就是5，7。这个坐标信息就代表spatial barcode在切片中的空间位置。

经过基因组比对，识别组织区域，可以得到spots和基因的一个表达矩阵。基于表达矩阵，进行下一步降维，聚类及差异分析。降维过程中有一个非常重要的参数就是主成分个数。这个参数对后续spots聚类影响非常大的。这里space ranger 默认为10。

**space ranger主要有四个命令**：

1.mkfastq实现的功能将测序产生的bcl格式数据转换成fastq格式

2.mkgtf/mkref用于构建参考基因组索引文件

3.有了原始fastq测序数据及参考基因组索引文件，可以用count命令对单样本进行spot和基因进行定量，定量结果储存在features，barcodes和matrix三个文件里。

4.mat2csv进行格式转换，将spot和基因进行定量转换为一个行为基因列为barcode的表达矩阵，csv格式。

另外，还有testrun:测试软件是否安装成功；upload：上传日志文件；sitecheck：查看配置，输出系统信息

### Space Ranger使用步骤

1. 下载Space Ranger - 1.0.0

[下载链接]: https://support.10xgenomics.com/spatial-gene-expression/software/overview/welcome

```bash
curl -o spaceranger-1.0.0.tar.gz "http://cf.10xgenomics.com/releases/spatial-exp/spaceranger-1.0.0.tar.gz?Expires=1575402715&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9zcGF0aWFsLWV4cC9zcGFjZXJhbmdlci0xLjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE1NzU0MDI3MTV9fX1dfQ__&Signature=fACB1rzbHv1rwUicNqL8SheRe6FkFOKxow5cTXcZPfOPBOTBEplElFMnOi4Xv4A2X3kydX45B-JnIaRj7I6a2doGEMTyqv84BnM5LxHAVBtWrXJyQqXbKKtgl9Dxe4BDnM9rPKhs6o2UbmWWAHX8Xu4J3~vgP3yXbhovuyl6OqCxu5p82oxTeOfN0bONqZdZ33svlAXJhatUTdpse2YCSRJZzov69NSHF6gE5DXl6iu5RWU7AgnjFgCuEFkQMwyn-FoYi2~i0s2fOFK0RCVI07~YKNDsjz3eXgOoHjWGPtWw5DAbPpTB2~32xkGzYeIYeZjH6m5JEgNGuvfWEyj~Aw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
```

2. 解压reference

   ```bash
   tar -xzvf refdata-cellranger-GRCh38-3.0.0.tar.gz
   ```

3. 运行spaceranger

   命令：

   （1）**mkfastq命令**实现功能是将测序仪产生的bcl格式转换为fastq格式，生成三类fastq压缩文件R1,R2,I1;R1read1端测序数据，R2reads2短测序数据，I储存sample index 信息；

   ```bash
   spaceranger mkfastq --id=Library1 --run=/path/to/data --csv=sample.csv --use-bases-mask =Y28,I8,Y151
   ```

   ![img](https://www.genenergy.cn/upload/images/2020/7/1717927829.jpg)

   （2）**mkgtf/ref** ：构建参考基因组索引文件；人和小鼠参考基因组可从官网下载，其它物种从Ensemble，UCSC等数据库下载fasta基因组文件和gtf注释文件。mkgtf实现的功能只过滤GTF文件，只保留感兴趣的基因注释 (gtf第三列特征类型必须有exon)；

   ```bash
   cellranger mkgtf hg38.gtf hg38.filtered.gtf -attribute= gene_ biotype:protein_ coding
   ```

   mkref---输入FASTA 和过滤后的GTF文件构建索引；如果基因组文件比较大，mkref这个过程可能需要几个小时时间。

   ```bash
   cellranger mkref --genome-refdata- ·cellranger-GRCh38-3.0.0 --fasta=hg38.fa l --genes= =hg38.filtered, gtf
   ```

   (3)**count**:实现的功能是spots和基因定量，输入是mkfastq得到的fastq文件及组织切片的图像，输出结果是一个outs目录，里面储存spots gene的定量结果及初步的降维聚类结果。

```bash
spaceranger count --id= =case --transcriptome=refdata-cellranger-GRCh38-3.0.0
--fastqs=./case --r2-length=91 --r1-length=28 --image Case-V19N13-090-C1.jpg --slide 
V19N13-090 --area C1 --loupe-alignment case_V19N13-090-C1.json
```

![img](https://www.genenergy.cn/upload/images/2020/7/1717106986.jpg)

这张芯片最上面有个serial  number ，这是每个visium波片上印刷的唯一标识符，四个捕获区域，从上倒下一次标记为A1,B1,C1,D1；如果space ranger无法自动识别出组织区域，需要manual alignment导出json文件，人工识别组织区域。这里特别提示以下：space ranger count命令对fastq测序数据命名要求非常严格。

（4）**mat2csv实现的功能进行格式转换**

```bash
spaceranger mat2csv \
filtered_feature_bc_matrix.h5|filtered_feature_bc_matrix/ \
gene_cell_exprs_table_raw.csv

输出：行为基因列为spots的表达矩阵，csv文件
```

目前，spaceranger 还未发布多个空转样本合并分析命令；但可以用cellranger aggr合并spotsXgene表达矩阵。特别注意：barcodeID后面数字与cellranger aggr输入aggr.csv文件顺序是一致的

TAAGAGATCACCTTAT-1：细胞ID加后缀-1表示Sample1

GGCACTAGACTACAA-2：细胞ID加后缀-2表示Sample2

* 结果解读：

![img](https://www.genenergy.cn/upload/images/2020/7/17171742614.jpg)

space ranger count命令输出结果保存在一个outs目录下，文件非常多。

1. 一个所有结果汇总的html页面，即web_summary.html。

web_summary.html的结果分成了summary和analysis两部分，summary主要是一些描述信息，比如识别出来的组织区域spots数，检测到的基因数目，测序质量，reads比对情况等信息。 如下图：

![img](https://www.genenergy.cn/upload/images/2020/7/17171838317.jpg)

第一部分是**spots＆基因数目的评估结果**，样本2总共检测到2698个被组织覆盖的spots，每个spot平均测序reads为115K左右，每个spot检测到的基因中位数是5861；

第二部分**指标反应的是碱基测序质量的信息**，一般情况下，基本上不会发生测序质量问题

第三部分描述是**识别出来的组织区域及spots详细信息**，图中红色框就是基准框，捕获区域内蓝色区域即是自动识别出来的组织区域。

![img](https://www.genenergy.cn/upload/images/2020/7/17171852583.jpg)

第四部分是**reads比对情况**，比对率一般不会低于70%，如果比对率过低，考虑参考基因组物种不对或者样品污染； 

第五部分是**样本信息，ID，试剂型号，slide编号，参考基因组及软件版本号等**

2. spatial目录储存着spatial barcode 空间位置信息；

![img](https://www.genenergy.cn/upload/images/2020/7/17171951896.jpg)

spatial 目录里面有个文件叫tissue_positions_lists.csv，里面存着每个spatial barcode ID的空间位置信息，总共六列。

**第一列：spatial barcode ID；第二列：是否覆盖组织区域；第三列：行坐标；第四列：列坐标；第五列：每个spot中心列像素坐标；第六列：每个spot中心行像素坐标；**

有了这些位置信息，我们就可以将spatial barcodes分析结果与它在组织中的空间位置联系起来。

3. filtered_freature_bc_matrix储存着spotXgene定量结果；

   ![img](https://www.genenergy.cn/upload/images/2020/7/17171940599.jpg)

filtered_feature_bc_matrix目录包括barcodes.tsv.gz，features.tsv.gz，matrix.mtx.gz三个文件。

barcodes.tsv.gz里面储存spatial barcode ID

features.tsv.gz文件储存检测到的基因id和symbol

matrix.mtx.gz是关于基因，Spot表达矩阵。第一行三个数字分别表示这个样本检测到了基因，spots和UMI总数。从第二行开始，第一列表示基因序号，二列表示spot序号，第三列表示检测到UMI数目。举例说一下，33509 1 53表示在features文件里第33509那个基因在第一个spot检测到的UMI数目是53。第一个spot对应的就是features文件第一个spatial barcode ID。

4. analysis目录储存着数据降维，聚类及差异分析结果；

5. cloupe.cloupe用于可视化，可以用loupe browser打开，查看基因的表达，辅助鉴定细胞类型，创建和修改子群，差异分析等操作。

## 下游分析

* 使用 Seurat 对空间数据集进行分析、可视化和集成

### 10x Visium

```R
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
```

Dataset： Visium v1 化学生成的矢状小鼠脑切片数据集，有两个连续的前部和两个（匹配的）连续的后部

```R
# 方法1
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

# 方法2
brain<-Load10X_Spatial(
  data.dir = 'E:\\BioInfo\\DATA\\Visium\\test',  #该目录包含10X提供的matrix.mtx，genes.tsv（或features.tsv）和barcodes.tsv文件。
  assay = "Spatial",
  slice = "slice1",  #组织切片存储图像的名称
  filter.matrix = TRUE, #仅保留已确定位于组织上方的spot
  to.upper = FALSE, #将所有功能名称转换为大写。 例如，当分析需要在人类和小鼠基因名称之间进行比较时，该功能将非常有用。
)
# 查看Seurat对象信息
brain

# 方法3
expr <- "151673_filtered_feature_bc_matrix.h5"
expr.mydata <- Seurat::Read10X_h5(filename = expr)
mydata <- Seurat::CreateSeuratObject(counts = expr.mydata, project = '151673', assay = 'Spatial')
mydata$slice <- 1
mydata$region <- '151673'
imgpath = "/public/home/fengyang/Spatial/brain/cortex/151673/data/"
img <- Seurat::Read10X_Image("/public/home/fengyang/Spatial/brain/cortex/151673/data")
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = mydata)]
mydata[['image']] <- img

# 该方法主要是读取以下三种文件
# h5文件：151673_filtered_feature_bc_matrix.h5；
# image文件：tissue_lowres_image.png；
# location文件：tissue_positions_list.txt
```

数据存储：

来自 10x 的 visium 数据包含以下数据类型：

* 基因表达矩阵点

* 组织切片的图像（在数据采集期间从 H&amp;E 染色中获得）

* 将原始高分辨率图像与此处用于可视化的较低分辨率图像相关联的缩放因子。

  在 Seurat 对象中，逐点基因表达矩阵类似于典型的“RNA”，<code style="box-sizing: border-box; font-family: Menlo, Monaco, Consolas, &quot;Courier New&quot;, monospace; font-size: 13.5px; padding: 2px 4px; color: rgb(51, 51, 51); background-color: rgb(248, 248, 248); border-radius: 4px;">Assay</code><font style="box-sizing: border-box; vertical-align: inherit;"><font style="box-sizing: border-box; vertical-align: inherit;">包含点水平，而不是单细胞水平数据。</font><font style="box-sizing: border-box; vertical-align: inherit;">图像本身存储在</font></font><code style="box-sizing: border-box; font-family: Menlo, Monaco, Consolas, &quot;Courier New&quot;, monospace; font-size: 13.5px; padding: 2px 4px; color: rgb(51, 51, 51); background-color: rgb(248, 248, 248); border-radius: 4px;">images</code><font style="box-sizing: border-box; vertical-align: inherit;"><font style="box-sizing: border-box; vertical-align: inherit;">Seurat 对象</font><font style="box-sizing: border-box; vertical-align: inherit;">的新</font><font style="box-sizing: border-box; vertical-align: inherit;">槽中。</font><font style="box-sizing: border-box; vertical-align: inherit;">该</font></font><code style="box-sizing: border-box; font-family: Menlo, Monaco, Consolas, &quot;Courier New&quot;, monospace; font-size: 13.5px; padding: 2px 4px; color: rgb(51, 51, 51); background-color: rgb(248, 248, 248); border-radius: 4px;">images</code><font style="box-sizing: border-box; vertical-align: inherit;"><font style="box-sizing: border-box; vertical-align: inherit;">槽还存储将斑点与其在组织图像上的物理位置相关联所需的信息。</font></font>



1. 数据预处理

首先需要对数据进行归一化，以解决数据点之间测序深度的差异。我们注意到，对于空间数据集，分子计数/点的差异可能很大，特别是如果整个组织的细胞密度存在差异。我们在这里看到了大量的异质性，这需要有效的标准化。

```R
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
```



![qc-1](images/qc-1.png)

这些图表明，不同点的分子计数差异不仅是技术性的，而且还取决于组织解剖结构。例如，神经元耗尽的组织区域（如皮质白质）可重现地表现出较低的分子计数。因此，`LogNormalize()`强制每个数据点在标准化后具有相同的基础“大小”的标准方法（例如函数）可能会出现问题。

作为替代方案，我们建议使用<mark> sctransform</mark>（Hafemeister 和 Satija，Genome Biology 2019），它构建了基因表达的正则化负二项式模型，以便在保留生物差异的同时考虑技术伪影。

```R
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
```

![norm.test3-1](images/norm.test3-1.png)

计算每个特征（基因）与 UMI 数量（`nCount_Spatial`此处为变量）的相关性。然后我们根据基因的平均表达将基因分组，并生成这些相关性的箱线图。您可以看到对数标准化未能充分标准化前三组中的基因，这表明技术因素继续影响高表达基因的标准化表达估计。相比之下，sctransform 归一化大大减轻了这种影响。

2. 基因表达可视化

在 Seurat 中，我们具有探索空间数据固有的视觉特性并与之交互的功能。`SpatialFeaturePlot()`Seurat 中的功能扩展`FeaturePlot()`，并且可以在组织组织学的顶部叠加分子数据。

例如，在小鼠大脑的这个数据集中，基因 Hpca 是一个强大的海马标记，而 Ttr 是一个脉络丛的标记。

```R
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
```

![](images/featureplot-1.png)

Seurat 中的默认参数强调分子数据的可视化。可以通过更改以下参数来调整斑点的大小（及其透明度）以改善组织学图像的可视化：

- `pt.size.factor`- 这将缩放点的大小。默认值为 1.6
- `alpha`- 最小和最大透明度。默认值为 c(1, 1)。
- 尝试设置为`alpha`c(0.1, 1)，降低表达较低的点的透明度

```R
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2
```

![](images/fpe1-1.png)

3. 降维、聚类和可视化

继续对 RNA 表达数据进行降维和聚类，与 scRNA-seq 分析相同的工作流程。

```R
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
```

接下来可以在 UMAP 空间（用`DimPlot()`）或用 覆盖在图像上可视化聚类结果`SpatialDimPlot()`。

```R
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2
```

![dim.plots-1](images/dim.plots-1.png)

可以使用该`cells.highlight`参数在`SpatialDimPlot()`. 这对于区分单个集群的空间定位非常有用：

```R
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3, 5, 8)), facet.highlight = TRUE, ncol = 3)
```

![facetdim-1](images/facetdim-1.png)

交互式绘图：

```R
SpatialDimPlot(brain, interactive = TRUE)
SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)
LinkedDimPlot(brain)
```



4. 空间变量特征基因的识别

Seurat 提供两种工作流程来识别与组织内空间位置相关的分子特征。

第一个是<mark>基于组织内预先注释的解剖区域执行差异表达</mark>，这可以<mark>从无监督聚类或先验知识中确定</mark>。这种策略在这种情况下会起作用，因为上面的集群表现出明显的空间限制。

```R
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
```

![](images/de-1.png)

另一种方法`FindSpatiallyVariables()`是<mark>搜索在没有预注释的情况下表现出空间模式的特征基因</mark>。默认方法（`method = 'markvariogram`），是由激发[Trendsceek](https://www.nature.com/articles/nmeth.4634)，该模型空间转录数据作为标记点处理并计算“变差函数”，它识别基因，其表达水平依赖于它们的空间位置。更具体地说，此过程计算 gamma(r) 值，该值测量相距特定“r”距离的两个点之间的相关性。默认情况下，我们在这些分析中使用 '5' 的 r 值，并且仅计算可变基因的这些值（其中变异的计算独立于空间位置）以节省时间。

```R
brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
    selection.method = "markvariogram")

# 可视化这个度量识别的前 6 个特征基因的表达
top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))
```

![](images/spatial.vf.plot-1.png)

5. 子集解剖区域

与单细胞对象一样，可以将对象子集化以关注数据的子集。在这里，我们大约对额叶皮层进行了子集化。此过程还有助于将这些数据与下一节中的皮质 scRNA-seq 数据集集成。首先，我们采用集群的子集，然后根据确切位置进一步细分。子集化后，我们可以在完整图像或裁剪图像上可视化皮层细胞。

```R
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2
```



![](images/subset1.plot-1.png)



6. 与单细胞数据集成

在~50um 处，visium 检测的斑点将包含多个细胞的表达谱。对于越来越多的可以使用 scRNA-seq 数据的系统，用户可能有兴趣对每个空间体素进行“解卷积”以预测细胞类型的基本组成。我们始终使用积分方法（与反卷积方法相反）发现了卓越的性能，这可能是因为表征空间和单细胞数据集的噪声模型大不相同，并且积分方法专门设计为对这些差异具有鲁棒性。

我们首先加载数据（可[在此处](https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1)下载），预处理 scRNA-seq 参考，然后执行标签转移。该过程为每个点输出每个 scRNA-seq 派生类的概率分类。我们将这些预测添加为 Seurat 对象中的新分析。

```R
allen_reference <- readRDS("../data/allen_cortex.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30)

# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "subclass", label = TRUE)
```

![](images/sc.data3-1.png)

```R
anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
    weight.reduction = cortex[["pca"]], dims = 1:30)
# 得到每个class的每个点的预测分数
cortex[["predictions"]] <- predictions.assay

# 可以区分这些神经元亚型的不同连续层
DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
```

基于这些预测分数，我们还可以预测位置受空间限制的*细胞类型*。我们使用基于标记点过程的相同方法来定义空间可变特征，但使用细胞类型预测分数作为“标记”而不是基因表达。

```R
cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "markvariogram", features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)
```

![](images/sc.data8-1.png)



最后，我们表明我们的综合程序能够恢复神经元和非神经元子集的已知空间定位模式，包括层流兴奋性、第 1 层星形胶质细胞和皮质灰质。

```R
SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT","L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))
```

![](images/sc.data9-1.png)

7. 在 Seurat 中处理多个切片

```R
# 读取数据
brain2 <- LoadData("stxBrain", type = "posterior1")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)
# merge多个切片
brain.merge <- merge(brain, brain2)
# 联合降维&聚类分析
DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)
# 可视化
DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
```

![](images/joint.viz-1.png)

```R
### SpatialFeaturePlot()默认将所有切片绘制为列，将分组/特征绘制为行
SpatialDimPlot(brain.merge)
SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))
```



![](images/joint.viz2-1.png)



![](images/joint.viz3-1.png)

### Slide-seq

数据来源：小鼠海马体的[Slide-seq v2](https://www.biorxiv.org/content/10.1101/2020.03.12.989806v1)生成的数据集。

```R
InstallData("ssHippo")
slide.seq <- LoadData("ssHippo")
```

1. 数据预处理

基因表达数据对珠子的初始预处理步骤类似于其他空间 Seurat 分析和典型的 scRNA-seq 实验。在这里，我们注意到许多珠子包含特别低的 UMI 计数，但选择保留所有检测到的珠子用于下游分析。

```R
plot1 <- VlnPlot(slide.seq, features = "nCount_Spatial", pt.size = 0, log = TRUE) + NoLegend()
slide.seq$log_nCount_Spatial <- log(slide.seq$nCount_Spatial)
plot2 <- SpatialFeaturePlot(slide.seq, features = "log_nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
```

![](images/qc.ss-1.png)



```R
# 使用sctransform规范化数据并执行标准的 scRNA-seq 降维和聚类工作流程
slide.seq <- SCTransform(slide.seq, assay = "Spatial", ncells = 3000, verbose = FALSE)
slide.seq <- RunPCA(slide.seq)
slide.seq <- RunUMAP(slide.seq, dims = 1:30)
slide.seq <- FindNeighbors(slide.seq, dims = 1:30)
slide.seq <- FindClusters(slide.seq, resolution = 0.3, verbose = FALSE)
#  在 UMAP 空间（用DimPlot()）或珠坐标空间用可视化聚类结果SpatialDimPlot()
plot1 <- DimPlot(slide.seq, reduction = "umap", label = TRUE)
plot2 <- SpatialDimPlot(slide.seq, stroke = 0)
plot1 + plot2  在 UMAP 空间（用DimPlot()）或珠坐标空间用可视化聚类结果SpatialDimPlot()
```

![](images/dim.plots.ss-1.png)

```R
SpatialDimPlot(slide.seq, cells.highlight = CellsByIdentities(object = slide.seq, idents = c(1,6, 13)), facet.highlight = TRUE)
```

![](images/dim.plots.ss-2.png)

2. 与 scRNA-seq 参考集成

为了便于对 Slide-seq 数据集进行细胞类型注释，我们正在利用现有的小鼠单细胞 RNA-seq 海马数据集，该数据集由[Saunders*、Macosko* 等人制作。2018 年](https://doi.org/10.1016/j.cell.2018.07.028)。数据可[在此处](https://www.dropbox.com/s/cs6pii5my4p3ke3/mouse_hippocampus_reference.rds?dl=0)作为已处理的 Seurat 对象下载，原始计数矩阵可在[DropViz 网站](http://dropviz.org/)上[找到](http://dropviz.org/)。

```R
ref <- readRDS("../data/mouse_hippocampus_reference.rds")
```

该论文的原始注释在 Seurat 对象的单元格元数据中提供。这些注释以多种“分辨率”提供，从大类 ( `ref$class`) 到细胞类型内的子簇 ( `ref$subcluster`)。出于本小插图的目的，我们将对细胞类型注释 (`ref$celltype`)进行修改，我们认为这取得了良好的平衡。

我们将首先运行 Seurat 标签转移方法来预测每个珠子的主要细胞类型。

```R
anchors <- FindTransferAnchors(reference = ref, query = slide.seq, normalization.method = "SCT",npcs = 50)
predictions.assay <- TransferData(anchorset = anchors, refdata = ref$celltype, prediction.assay = TRUE,
    weight.reduction = slide.seq[["pca"]], dims = 1:50)
slide.seq[["predictions"]] <- predictions.assay
# 可视化预测分数
DefaultAssay(slide.seq) <- "predictions"
SpatialFeaturePlot(slide.seq, features = c("Dentate Principal cells", "CA3 Principal cells", "Entorhinal cortex",
    "Endothelial tip", "Ependymal", "Oligodendrocyte"), alpha = c(0.1, 1))
```

![](images/transfer.viz.ss-1.png)

```R
slide.seq$predicted.id <- GetTransferPredictions(slide.seq)
Idents(slide.seq) <- "predicted.id"
SpatialDimPlot(slide.seq, cells.highlight = CellsByIdentities(object = slide.seq, idents = c("CA3 Principal cells",
    "Dentate Principal cells", "Endothelial tip")), facet.highlight = TRUE)
```

![](images/max.idents.ss-1.png)

3. 空间变量特征的识别

正如 Visium 中所述，可以通过两种一般方式识别空间可变特征：预先注释的解剖区域之间的差异表达测试或测量特征对空间位置依赖性的统计数据。

在这里，我们`FindSpatiallyVariableFeatures()`通过设置`method = 'moransi'`. 

Moran's I 计算整体空间自相关并给出一个统计量（类似于相关系数），用于衡量特征对空间位置的依赖性。这使我们能够根据特征的表达在空间上的变化程度对特征进行排名。为了便于对该统计数据进行快速估计，我们实施了一种基本的分箱策略，该策略将在 Slide-seq puck 上绘制一个矩形网格，并对每个分箱内的特征和位置进行平均。x 和 y 方向的 bin 数量分别由`x.cuts`和`y.cuts`参数控制。此外，虽然不是必需的，但安装可选`Rfast2`包（`install.packages('Rfast2')`)，将通过更有效的实现显着减少运行时间。

```R
DefaultAssay(slide.seq) <- "SCT"
slide.seq <- FindSpatiallyVariableFeatures(slide.seq, assay = "SCT", slot = "scale.data", features = VariableFeatures(slide.seq)[1:1000],
    selection.method = "moransi", x.cuts = 100, y.cuts = 100)
# 可视化由 Moran's I 识别的前 6 个特征的表达
SpatialFeaturePlot(slide.seq, features = head(SpatiallyVariableFeatures(slide.seq, selection.method = "moransi"),
    6), ncol = 3, alpha = c(0.1, 1), max.cutoff = "q95")
```

![](images/spatial.vf.plot.ss-1.png)

## 参考

1. [spacer ranger](https://www.genenergy.cn/news/156/2066.html)
2. [Seurat](https://satijalab.org/seurat/articles/spatial_vignette.html)
