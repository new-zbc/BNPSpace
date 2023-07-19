## MOB data

File *MOB_sce.RData* is a `SingleCellExperiment` object containing original count data and locations for each spot. 

File *MOB_analysis.R* is the source code to analyze MOB data via BNPSpace and other competitive methods. 

File "MCMCresults.RData" is the MCMC posterior samples of BNPSpace. 

File "result_other_method.RData" is the clustering results of MOB data for other competitive methods. 

File "result_sce.RData" is a `SingleCellExperiment` object containing original data, ground truth, and clustering results of each method. 

## DLPFC data

File *151509_counts.RData* is a `SingleCellExperiment` object containing original count data and locations for each spot. 

File *BNPSpace_results.RData* contains clustering results of BNPSpace and postprocess results after hierarchical clustering.

File *different_K_ARI.RData* is clustering results of other methods given different number of clusters. 

File *result_other_method.RData* is clustering results of other competitive methods.