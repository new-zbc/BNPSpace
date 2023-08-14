library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)


source("R/utils.R")
source("R/main.R")

### construct spots network
Adj = find_neighbors(sce, "ST", "image")

neighbors = find_neighbor_index(Adj, "ST")
result = NSCFS(sce, neighbors, n_clusters = NULL, f = 1, n_iters = 100, seed = 1)
library(flexclust)
ARI_value = randIndex(table(result$pred_label, colData(sce)$label))


