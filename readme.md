## BNPSpace

### User manual

The data structure used for BNPSpace is `SingleCellExperiment` object containing assay "counts" and col data "row" and "col". "row" and "col" are the coordinates of each spot. 

```R
source("demo/dataSim.R")
> sce
class: SingleCellExperiment 
dim: 500 1600 
metadata(0):
assays(2): counts logcounts
rownames(500): gene_1 gene_2 ... gene_499 gene_500
rowData names(0):
colnames(1600): spot_1 spot_2 ... spot_1599 spot_1600
colData names(4): row col label sizeFactor
reducedDimNames(1): PCA
mainExpName: NULL
altExpNames(0):
```

The above simulation data contains two assays, original counts data and log-normalized data (logcounts). There are 500 genes and 1600 spots. PCA represents the dimension reduction features using principal component analysis. The col data names "row" and "col" are spatial coordinates, "label" is the cluster ground truth, and "sizeFactor" is the size factor of each spot.

Then, we can run BNPSpace.

```R
library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)
source("R/utils.R")
source("R/main.R")
### construct spots network
Adj = find_neighbors(sce, "ST", "image")
neighbors = find_neighbor_index(Adj, "ST")
result = NSCFS(sce, neighbors, n_clusters = NULL, f = 1, n_iters = 500, seed = 1)
ARI_value = randIndex(table(result$pred_label, colData(sce)$label))

### clustering result
result$pred_label
### feature selection result
result$gamma

```

