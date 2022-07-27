########################   Chronostasis   ##########################################
####################################################################################
# Core plotting functions
## Scatter plots for embeddings
    import scanpy as sc
    import pandas as pd
    from matplotlib.pyplot import rc_context
    
    sc.set_figure_params(dpi=100, color_map = 'viridis_r')
    sc.settings.verbosity = 1
    sc.logging.print_header()
### Load the datasets
    adata = sc.datasets.pbmc68k_reduced()
### Visualization of gene expression and other variables
#### rc_context is used for the figure size, in this case 4x4
    with rc_context({'figure.figsize': (4, 4)}):
        sc.pl.umap(adata, color='CD79A')
    with rc_context({'figure.figsize': (3, 3)}):
        sc.pl.umap(adata, color=['CD79A', 'MS4A1', 'IGJ', 'CD3D', 'FCER1A', 'FCGR3A', 'n_counts', 'bulk_labels'], s=50, frameon=False, ncols=4, vmax='p99') 
#### The functions for scatterplots have many options that allow fine tuning of the images. For example, we can look at the clustering as follows:
##### compute clusters using the leiden method and store the results with the name `clusters`
    sc.tl.leiden(adata, key_added='clusters', resolution=0.5)    
    with rc_context({'figure.figsize': (5, 5)}):
    sc.pl.umap(adata, color='clusters', add_outline=True, legend_loc='on data', 
               legend_fontsize=12, legend_fontoutline=2,frameon=False,  
               title='clustering of cells', palette='Set1') 
### Identification of clusters based on known marker genes
    marker_genes_dict = {
    'B-cell': ['CD79A', 'MS4A1'],     
    'Dendritic': ['FCER1A', 'CST3'],    
    'Monocytes': ['FCGR3A'],
    'NK': ['GNLY', 'NKG7'],
    'Other': ['IGLL1'],
    'Plasma': ['IGJ'],
    'T-cell': ['CD3D'],
    }
    
#### dotplot
##### the color represents the mean expression within each of the categories (in this case in each cluster) and the dot size indicates the fraction of cells in the categories expressing a gene       
    sc.pl.dotplot(adata, marker_genes_dict, 'clusters', dendrogram=True)
##### This information can be used to manually annotate the cells as follows
##### create a dictionary to map cluster to annotation label
    cluster2annotation = {
        '0': 'Monocytes',
        '1': 'Dendritic',
        '2': 'T-cell',
        '3': 'NK',
        '4': 'B-cell',
        '5': 'Dendritic',
        '6': 'Plasma',
        '7': 'Other',
        '8': 'Dendritic',
    }
##### add a new `.obs` column called `cell type` by mapping clusters to annotation using pandas `map` function
    adata.obs['cell type'] = adata.obs['clusters'].map(cluster2annotation).astype('category')    
    sc.pl.dotplot(adata, marker_genes_dict, 'cell type', dendrogram=True)        
    sc.pl.umap(adata, color='cell type', legend_loc='on data',
           frameon=False, legend_fontsize=10, legend_fontoutline=2)
    
#### violin plot
    with rc_context({'figure.figsize': (4.5, 3)}):
        sc.pl.violin(adata, ['CD79A', 'MS4A1'], groupby='clusters' )
    with rc_context({'figure.figsize': (4.5, 3)}):
        sc.pl.violin(adata, ['n_genes', 'percent_mito'], groupby='clusters', stripplot=False, inner='box')  # use stripplot=False to remove the internal dots, inner='box' adds a boxplot inside violins
    
#### stacked-violin plot
    ax = sc.pl.stacked_violin(adata, marker_genes_dict, groupby='clusters', swap_axes=False, dendrogram=True)
    
#### matrixplot
    sc.pl.matrixplot(adata, marker_genes_dict, 'clusters', dendrogram=True, cmap='Blues', standard_scale='var', colorbar_title='column scaled\nexpression')    
##### scale and store results in layer
    adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X
    sc.pl.matrixplot(adata, marker_genes_dict, 'clusters', dendrogram=True, 
                 colorbar_title='mean z-score', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r')
   
##### Combining plots in subplots
    import matplotlib.pyplot as plt

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20,4), gridspec_kw={'wspace':0.9})

    ax1_dict = sc.pl.dotplot(adata, marker_genes_dict, groupby='bulk_labels', ax=ax1, show=False)
    ax2_dict = sc.pl.stacked_violin(adata, marker_genes_dict, groupby='bulk_labels', ax=ax2, show=False)
    ax3_dict = sc.pl.matrixplot(adata, marker_genes_dict, groupby='bulk_labels', ax=ax3, show=False, cmap='viridis')        
    
#### Heatmaps
    ax = sc.pl.heatmap(adata, marker_genes_dict, groupby='clusters', cmap='viridis', dendrogram=True)    
##### The heatmap can also be plotted on scaled data
    ax = sc.pl.heatmap(adata, marker_genes_dict, groupby='clusters', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r', dendrogram=True, swap_axes=True, figsize=(11,4))    

#### Tracksplot [the same information as the heatmap]   
    ax = sc.pl.tracksplot(adata, marker_genes_dict, groupby='clusters', dendrogram=True)
    
### Visualization of marker genes
    sc.tl.rank_genes_groups(adata, groupby='clusters', method='wilcoxon')
#### Visualize marker genes using dotplot
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=4)
##### In order to get a better representation we can plot log fold changes* instead of gene expression    
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=4, values_to_plot='logfoldchanges', min_logfoldchange=3, vmax=7, vmin=-7, cmap='bwr')
##### Focusing on particular groups
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=30, values_to_plot='logfoldchanges', min_logfoldchange=4, vmax=7, vmin=-7, cmap='bwr', groups=['1', '5'])
##### Visualize marker genes using matrixplot 
    sc.pl.rank_genes_groups_matrixplot(adata, n_genes=3, use_raw=False, vmin=-3, vmax=3, cmap='bwr', layer='scaled')
##### Visualize marker genes using stacked violin plots
    sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=3, cmap='viridis_r')
##### Visualize marker genes using heatmap
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=3, use_raw=False, swap_axes=True, vmin=-3, vmax=3, cmap='bwr', layer='scaled', figsize=(10,7), show=False);
###### Showing 10 genes per category, turning the gene labels off and swapping the axe
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, use_raw=False, swap_axes=True, show_gene_labels=False,
                                vmin=-3, vmax=3, cmap='bwr')   
#### Visualize marker genes using tracksplot
    sc.pl.rank_genes_groups_tracksplot(adata, n_genes=3)   
    
### Comparison of marker genes using split violin plots
    with rc_context({'figure.figsize': (9, 1.5)}):
        sc.pl.rank_genes_groups_violin(adata, n_genes=20, jitter=False)    
### Dendrogram options
##### compute hierarchical clustering using PCs (several distance metrics and linkage methods are available).
    sc.tl.dendrogram(adata, 'bulk_labels')
    ax = sc.pl.dendrogram(adata, 'bulk_labels')

### Plot correlation
    ax = sc.pl.correlation_matrix(adata, 'bulk_labels', figsize=(5,3.5))






    