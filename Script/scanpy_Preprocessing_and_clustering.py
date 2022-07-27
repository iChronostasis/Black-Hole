########################   Chronostasis   ##########################################
####################################################################################
# scanpy tutorial
# Load the packages
    import numpy as np
    import pandas as pd
    import scanpy as sc
# Set the parameters
    sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3) 日志输出的详细程度，其中详细程度的含义如下：0='error', 1='warning', 2='info', 3='hint', 4=more details, 5=even more details, etc .
    sc.settings.verbosity = 2  # reduce the verbosity
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor='white')    
# The file to save 
    results_file = 'write/~.h5ad'  # the file that will store the analysis results
# Load the data
    adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)      
# Preprocessing
    sc.pl.highest_expr_genes(adata, n_top=20, )
## Basic filtering【similar to the parmeters of seurat::Createseuratobject 】
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
## With `pp.calculate_qc_metrics`, we can compute many metrics very efficiently
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
## A violin plot of some of the computed quality measures
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)
## Remove cells that have too many mitochondrial genes expressed or too many total counts
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')    
## Actually do the filtering by slicing the `AnnData` object   
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :] 
## Total-count normalize (library-size correct) the data matrix $\mathbf{X}$ to 10,000 reads per cell, so that counts become comparable among cells.    
    sc.pp.normalize_total(adata, target_sum=1e4)
## Logarithmize the data
    sc.pp.log1p(adata)
## Identify highly-variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
### extracting highly variable genes
    sc.pl.highly_variable_genes(adata)    
## Set the `.raw` attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. 
## This simply freezes the state of the AnnData object.    
    adata.raw = adata
## Actually do the filtering 
    adata = adata[:, adata.var.highly_variable]
## Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed.
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
## Scale each gene to unit variance. Clip values exceeding standard deviation 10
    sc.pp.scale(adata, max_value=10)
    
# PCA
## Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data
    sc.tl.pca(adata, svd_solver='arpack')
## make a scatter plot in the PCA coordinate (useless)
    sc.pl.pca(adata, color='CST3')
## This gives us information about how many PCs we should consider in order to compute the neighborhood relations of cells
    sc.pl.pca_variance_ratio(adata, log=True)
## save the result
    adata.write(results_file)
    
# Computing the neighborhood graph
## compute the neighborhood graph of cells using the PCA representation of the data matrix
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
## Embedding the neighborhood graph
    sc.tl.umap(adata)
## Plot the umap (similar to featureplot)    
    sc.pl.umap(adata, color=['CST3', 'NKG7', 'PPBP'])          
# Clustering the neighborhood graph [Cluster]
## we recommend the Leiden graph-clustering method (community detection based on optimizing modularity)     
    sc.tl.leiden(adata)
    sc.pl.umap(adata, color=['leiden', 'CST3', 'NKG7'])
## Save the result
    adata.write(results_file)    
## Finding marker genes
###  by default, the `.raw` attribute of AnnData is used in case it has been initialized before
### t-test
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)     
### sc.settings.verbosity = 2  # reduce the verbosity
###  wilcoxon
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
## Save the result
    adata.write(results_file)    
### As an alternative, let us rank genes using logistic regression    
    sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
    c.pl.rank_genes_groups(adata, n_genes=25, sharey=False)

# Cell_typing for clusters
## define a list of marker genes for later reference
    marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',  
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP']        
## Reload the object that has been save with the Wilcoxon Rank-Sum test result
    adata = sc.read(results_file)
## Show the 10 top ranked genes per cluster 0, 1, ..., 7 in a dataframe
    pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)
## Get a table with the scores and groups
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}).head(5)        
## Compare to a single cluster
    sc.tl.rank_genes_groups(adata, 'leiden', groups=['0'], reference='1', method='wilcoxon')
    sc.pl.rank_genes_groups(adata, groups=['0'], n_genes=20)     
### If we want a more detailed view for a certain group, use sc.pl.rank_genes_groups_violin
    sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)    
## Reload the object with the computed differential expression (i.e. DE via a comparison with the rest of the groups): 
    adata = sc.read(results_file)
    sc.pl.rank_genes_groups_violin(adata, groups='0', n_genes=8)
## compare a certain gene across groups, use the following (similar to the violinplot in every clusters)    
    sc.pl.violin(adata, ['CST3', 'NKG7', 'PPBP'], groupby='leiden')
## Actually mark the cell types
    new_cluster_names = [
    'CD4 T', 'CD14 Monocytes',
    'B', 'CD8 T', 
    'NK', 'FCGR3A Monocytes',
    'Dendritic', 'Megakaryocytes']
    adata.rename_categories('leiden', new_cluster_names)
## Assign the cell type to the clusters        
    sc.pl.umap(adata, color='leiden', legend_loc='on data', title='', frameon=False, save='.pdf')
## Dotplot
    sc.pl.dotplot(adata, marker_genes, groupby='leiden');
## violin plot
    sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', rotation=90);     
### save the file    
    adata.write(results_file, compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading
    adata.raw.to_adata().write('./write/~_withoutX.h5ad')
## If you want to export to "csv", you have the following options:
    # Export single fields of the annotation of observations
    # adata.obs[['n_counts', 'louvain_groups']].to_csv(
    #     './write/pbmc3k_corrected_louvain_groups.csv')

    # Export single columns of the multidimensional annotation
    # adata.obsm.to_df()[['X_pca1', 'X_pca2']].to_csv(
    #     './write/pbmc3k_corrected_X_pca.csv')

    # Or export everything except the data using `.write_csvs`.
    # Set `skip_data=False` if you also want to export the data.
    # adata.write_csvs(results_file[:-5], )    
    