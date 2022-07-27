########################   Chronostasis   ##########################################
####################################################################################
# The following tutorial describes a simple PCA-based method for integrating data we call [ingest]
# The ingest assumes an annotated reference dataset that captures the biological variability of interest
    import scanpy as sc
    import pandas as pd
    import seaborn as sns
## Set the parameters
    sc.settings.verbosity = 1             # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')
## Load the datasets
### an annotated reference dataset `adata_ref` and a dataset for which you want to query labels and embeddings `adata`
    adata = sc.datasets.pbmc68k_reduced()
    adata_ref_file = 'write/pbmc3k.h5ad' ### before
    adata_ref = sc.read(adata_ref_file)
## To use `sc.tl.ingest`, the datasets need to be defined on the same variables
    var_names = adata_ref.var_names.intersection(adata.var_names)
    adata_ref = adata_ref[:, var_names]
    adata = adata[:, var_names]    
## The model and graph (here PCA, neighbors, UMAP) trained on the reference data will explain the biological variation observed within it
    sc.pp.pca(adata_ref)
    sc.pp.neighbors(adata_ref)
    sc.tl.umap(adata_ref)  
    #sc.tl.louvain(adata_ref)
    sc.pl.umap(adata_ref, color='louvain')
##  Mapping the datasets using ingest       
    sc.tl.ingest(adata, adata_ref, obs='louvain')
    adata.uns['louvain_colors'] = adata_ref.uns['louvain_colors']  # fix colors
    sc.pl.umap(adata, color=['louvain', 'bulk_labels'], wspace=0.5)
## By comparing the 'bulk_labels' annotation with 'louvain', we see that the data has been reasonably mapped
    adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])
    adata_concat.obs.louvain = adata_concat.obs.louvain.astype('category')
    adata_concat.obs.louvain.cat.reorder_categories(adata_ref.obs.louvain.cat.categories, inplace=True)  # fix category ordering
    adata_concat.uns['louvain_colors'] = adata_ref.uns['louvain_colors']  # fix category colors     
    sc.pl.umap(adata_concat, color=['batch', 'louvain'])

## Using BBKNN
    sc.tl.pca(adata_concat)
    
    %%time
    sc.external.pp.bbknn(adata_concat, batch_key='batch')  # running bbknn 1.3.6
    
    sc.tl.umap(adata_concat)
    sc.pl.umap(adata_concat, color=['batch', 'louvain'])
    
# Pancreas [example]
## note that this collection of batches is already intersected on the genes
    adata_all = sc.read('data/pancreas.h5ad')     
    adata_all.shape   ### dim(adata_all)
## Inspect the cell types observed in these studies
    counts = adata_all.obs.celltype.value_counts()
    counts
### To simplify visualization, let's remove the 5 minority classes
    counts.index[-5:]   
    minority_classes = counts.index[-5:].tolist()        # get the minority classes
    adata_all = adata_all[                               # actually subset
        ~adata_all.obs.celltype.isin(minority_classes)]
    adata_all.obs.celltype.cat.reorder_categories(       # reorder according to abundance 
        counts.index[:-5].tolist(), inplace=True)  
## Seeing the batch effect
    sc.pp.pca(adata_all)
    sc.pp.neighbors(adata_all)
    sc.tl.umap(adata_all)
       
    sc.pl.umap(adata_all, color=['batch', 'celltype'], palette=sc.pl.palettes.vega_20_scanpy)

## It can be well-resolved using BBKNN
    %%time
    sc.external.pp.bbknn(adata_all, batch_key='batch')
### If one prefers to work more iteratively starting from one reference dataset, one can use ingest    
    sc.tl.umap(adata_all)
    sc.pl.umap(adata_all, color=['batch', 'celltype'])
## Mapping onto a reference batch using ingest
    adata_ref = adata_all[adata_all.obs.batch == '0']
    sc.pp.pca(adata_ref)
    sc.pp.neighbors(adata_ref)
    sc.tl.umap(adata_ref)
    sc.pl.umap(adata_ref, color='celltype')
    
## Iteratively map labels (such as 'celltype') and embeddings (such as 'X_pca' and 'X_umap') from the reference data onto the query batches
    adatas = [adata_all[adata_all.obs.batch == i].copy() for i in ['1', '2', '3']]    
    sc.settings.verbosity = 2  # a bit more logging
    for iadata, adata in enumerate(adatas):
        print(f'... integrating batch {iadata+1}')
        adata.obs['celltype_orig'] = adata.obs.celltype  # save the original cell type
        sc.tl.ingest(adata, adata_ref, obs='celltype')
## Each of the query batches now carries annotation that has been contextualized with `adata_ref`. By concatenating, we can view it together.
    adata_concat = adata_ref.concatenate(adatas)
    adata_concat.obs.celltype = adata_concat.obs.celltype.astype('category')
    adata_concat.obs.celltype.cat.reorder_categories(adata_ref.obs.celltype.cat.categories, inplace=True)  # fix category ordering
    adata_concat.uns['celltype_colors'] = adata_ref.uns['celltype_colors']  # fix category coloring
    sc.pl.umap(adata_concat, color=['batch', 'celltype'])
    
## Evaluating consistency
### subset the data to the query batches
    adata_query = adata_concat[adata_concat.obs.batch.isin(['1', '2', '3'])]
    
    sc.pl.umap(
    adata_query, color=['batch', 'celltype', 'celltype_orig'], wspace=0.4)
    
## Cell types conserved across batches
### focus on cell types that are conserved with the reference, to simplify reading of the confusion matrix
    obs_query = adata_query.obs
    conserved_categories = obs_query.celltype.cat.categories.intersection(obs_query.celltype_orig.cat.categories)  # intersected categories
    obs_query_conserved = obs_query.loc[obs_query.celltype.isin(conserved_categories) & obs_query.celltype_orig.isin(conserved_categories)]  # intersect categories
    obs_query_conserved.celltype.cat.remove_unused_categories(inplace=True)  # remove unused categoriyes
    obs_query_conserved.celltype_orig.cat.remove_unused_categories(inplace=True)  # remove unused categoriyes
    obs_query_conserved.celltype_orig.cat.reorder_categories(obs_query_conserved.celltype.cat.categories, inplace=True)  # fix category ordering
    pd.crosstab(obs_query_conserved.celltype, obs_query_conserved.celltype_orig)            
### All cell types
    pd.crosstab(adata_query.obs.celltype, adata_query.obs.celltype_orig)

## Visualizing distributions across batches    
### Density plot
    sc.tl.embedding_density(adata_concat, groupby='batch')
    sc.pl.embedding_density(adata_concat, groupby='batch')     
### Partial visualizaton of a subset of groups in embedding
    for batch in ['1', '2', '3']:
        sc.pl.umap(adata_concat, color='batch', groups=[batch])    
    
    
    
            
        
        
        
        
        
        
        
        