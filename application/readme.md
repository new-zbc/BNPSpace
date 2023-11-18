## MOB data

File *data/MOB_sce.RData* is a `SingleCellExperiment` object containing original count data and locations for each spot. 

Folder *script/* contains the source code to analyze MOB data via BNPSpace and other competitive methods and the plot function to generate figures in the paper. 

File *MCMCresults.RData* is the MCMC posterior sample of BNPSpace. 

File *clustering_result_sce.RData* is the `SingleCellExperiment` object with clustering results of MOB data. 

File *BIC_result.RData* is the pBIC result for MOB data. 

File *0-result-MOB.xlsx* is the enrichment analysis of discriminating genes.

## DLPFC data

File *data/151509_counts.RData* is a `SingleCellExperiment` object containing original count data and locations for each spot. 

Folder *script/* contains the source code to analyze MOB data via BNPSpace and other competitive methods and the plot function to generate figures in the paper. 

File *BNPSpace_results.RData* contains clustering results of BNPSpace and post-process results after hierarchical clustering.

File *different_K_ARI.RData* is clustering results of other methods given different numbers of clusters. 

File *result_other_method.RData* is clustering results of other competitive methods.

File *0-result-DLPFC.xlsx* is the enrichment analysis of discriminating genes.