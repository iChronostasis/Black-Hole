## ε(*′･∀･｀)зﾞ
# About Black Hole
<mark>仅供本实验室</mark><br>
2022.04.16 上传关于单细胞数据分析的流程教程及部分代码脚本。<br>
2022.05.05 目前储存并更新关于单细胞数据分析的流程教程及平时常使用的代码脚本。<br>
2022.05.07 更新并上传了部分代码脚本，同时自我练习规范编写代码的习惯，可能以往的代码脚本的一些编写习惯会在实际使用过程中进行修改。<br>
2022.05.09 scRNA-seq_integration.R 更新了Harmony和LIGER的合并方法，但是LIGER的使用率相比Seurat和Harmnoy好像有点低...另外，Visualization.R 更新了遇到了值得存储的图。<br>
2022.05.10 更新Monocle2.R代码脚本，后续会在实际应用中继续更新和添加。更新了空间转录组数据分析教程的内容。<br>
2022.05.12 上传Cellranger和scrublet的使用代码脚本，后续如果搭建流程会使用到。<br>
2022.05.17 上传了与scanpy相关的代码脚本。<br>
### To be continued...

# User Guide

1. **The_workflow_of_single-cell_data_analysis.md** 包含了大部分的single cell RNA-seq数据处理分析的教程内容，常用及部分不常用的R包或者python包。
2. **xmind_The_workflow_of_single-cell_data_analysis.pdf** 是关于single cell RNA-seq数据处理分析的思维导图，建议与上一个文件一齐食用。
3. **The_workflow_of_single-cell_data_analysis.html** 以及 **The_workflow_of_single-cell_data_analysis.pdf** 以上两个可以支持本地有该存储库的images文件夹的情况下离线在本地查看本教程内容。
4. **images** 该文件夹储存了**The_workflow_of_single-cell_data_analysis.md**和**Spatially-resolved_RNA-seq_data_analysis_pipeline.md**教程中的相关素材。
5. **Spatially-resolved_RNA-seq_data_analysis_pipeline.md** 介绍了关于空间转录组scRNA-seq的分析教程。**xmind_Spatially-resolved_RNA-seq_data_analysis_pipeline.pdf** 是关于空间转录组scRNA-seq的分析的思维导图
6. **Cellranger** 该文件夹储存了Cellranger 针对10x Genomics单细胞数据分析的脚本。
    * Cellranger/Cellranger.sh shell命令主体
    * Cellranger/task_cellranger.sh 运行的命令脚本
7. **Script** 该文件夹储存了相关的代码脚本。建议可与教程内容相结合食用。
    * Script/Standard_processing_workflow.R 单个单细胞数据集基本分析流程
    * Script/scRNA-seq_integration.R 当合并多个数据集时的Seurat,Harmony,LIGER基本分析流程(推荐Seurat)
    * Script/Visualization.R 多个结果图的代码整理
    * Script/SingleR.R 使用R包SingleR包进行细胞类型注释
    * Script/cellTyping.R 使用PanglaoDB数据库进行细胞类型注释
    * Script/Monocle2.R 使用Monocle2包实现伪时间分析+ 轨迹分析
    * Script/Slingshot.R 伪时间分析 + 轨迹分析(之前遇到过，顺带上传)
    * Script/scrublet.ipynb scrublet包去除doublets的代码脚本
    * Script/scanpy_Preprocessing_and_clustering.py scanpy处理单细胞数据脚本
    * Script/scapny_Trajectory_inference.py scanpy伪时间分析+轨迹分析
    * Script/scanpy_integrating_data_using_ingest.py scanpy使用ingest合并数据集
    * Script/scanpy_Visualization.py scanpy 可视化相关的代码
    
# 预祝食用愉快 (｡･∀･)ﾉﾞ
