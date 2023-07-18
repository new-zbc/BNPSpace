
library(SingleCellExperiment)
library(Matrix)
library(flexclust)
library(scater)
library(scran)
source("R/utils.R")
source("R/main.R")

set.seed(1)
load("application/MOBdata/MOB_sce.RData")
### quality control
sce = spatialFilter(sce, cutoff_max = 10)

### construct spots network
Adj = find_neighbors(sce, "ST", "lattice")
neighbors = find_neighbor_index(Adj, "ST")

seed = 20220222
####################################################################
######################### proposed method  #############################
###################################################################
#### run NSCFS
result1 = NSCFS(sce, neighbors, n_clusters = NULL, f = 1, n_iters = 500, seed = 2)

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = result1$pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))
ARI_value

save(result1, file = "application/MOBdata/MCMCresults.RData")


####################################################################
######################### SC.MEB 1  #################################
###################################################################
library(SC.MEB)
set.seed(seed)
K_set = 2:20
data_mat = reducedDim(sce)[,1:10]

selection = SC.MEB(data_mat, Adj_sp = Adj, parallel = FALSE, K_set = K_set)

res = selectK(selection, K_set = K_set, criterion = "BIC", c = 1)

est_K = res$best_K_BIC
pred_label = res$best_K_label

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))

SC.MEB1 = list(K = est_K, label = pred_label, ARI = ARI_value)



####################################################################
######################### SC.MEB 2  #################################
###################################################################
library(SC.MEB)
set.seed(seed)
K_set = 2:20
data_mat = reducedDim(sce)

selection = SC.MEB(data_mat, Adj_sp = Adj, parallel = FALSE, K_set = K_set)

res = selectK(selection, K_set = K_set, criterion = "BIC", c = 1)

est_K = res$best_K_BIC
pred_label = res$best_K_label

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))

SC.MEB2 = list(K = est_K, label = pred_label, ARI = ARI_value)


####################################################################
######################### DR.SC 1  #################################
###################################################################
library(DR.SC)
set.seed(seed)
K_set = 2:20

sce_drsc = spatialPreprocess(raw_sce)
sce_drsc = spatialFilter(sce_drsc, cutoff_feature = 0, cutoff_max = 0)

data_mat = t(assay(sce_drsc, "logcounts")[rowData(sce_drsc)$is.HVG, ])

reslist = DR.SC_fit(data_mat, q = 10, K= K_set, Adj_sp = Adj, coreNum = 1)

res = selectModel(reslist, criteria = 'MBIC', pen.const=1)

est_K = res$bestK
pred_label = res$cluster

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))

DR.SC1 = list(K = est_K, label = pred_label, ARI = ARI_value)


####################################################################
######################### DR.SC 2  #################################
###################################################################
library(DR.SC)
set.seed(seed)
K_set = 2:20

reslist = DR.SC_fit(data_mat, q = 15, K= K_set, Adj_sp = Adj, coreNum = 1)

res = selectModel(reslist, criteria = 'MBIC', pen.const=1)

est_K = res$bestK
pred_label = res$cluster

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))

DR.SC2 = list(K = est_K, label = pred_label, ARI = ARI_value)


####################################################################
######################### BayesSpace  #################################
###################################################################
library(BayesSpace)
set.seed(seed)

qPlot(sce, qs = seq(2, 20), burn.in = 100, nrep = 1000, use.dimred = "PCA", 
      platform = "ST", init.method = "kmeans", model = "normal", 
      precision = "equal", gamma = 2, alpha = 1, beta = 0.01)


set.seed(2)
sce_result = BayesSpace::spatialCluster(sce, q=4, use.dimred = "PCA", d = 15,
                                        platform = "ST", init.method = "kmeans", model = "normal", 
                                        precision = "equal", nrep = 5000, gamma = 2, 
                                        alpha = 1, beta = 0.01, save.chain = FALSE)

pred_label = colData(sce_result)$spatial.cluster

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))

BayesSpace = list(label = pred_label, ARI = ARI_value)



####################################################################
######################### k-means  #################################
###################################################################
set.seed(seed)
data_mat = reducedDim(sce)

K_set = 2:20
wss = rep(0, length(K_set))
for(i in 1:length(K_set)){
  wss[i] = sum(kmeans(data_mat, centers = K_set[i], nstart = 5)$withinss)
}

plot(K_set, wss, type = "b")

pred_label <- kmeans(data_mat, centers = 4, nstart = 3)$cluster

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))

k_means = list(label = pred_label, ARI = ARI_value)


####################################################################
######################### GMM method  #############################
###################################################################
set.seed(seed)
library(mclust)
data_mat = reducedDim(sce)

fit_int = Mclust(data_mat, G = 2:20)

est_K = fit_int$G

pred_label <- fit_int$classification

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))

GMM = list(K = est_K, label = pred_label, ARI = ARI_value)


#################################################################
############################# Louvain  ##########################
##################################################################
set.seed(seed)
library(Seurat)
x_object = as.Seurat(sce, counts = "counts", data = "logcounts", project =  "sce_to_seurat")

x_object = Seurat::FindNeighbors(x_object, reduction = "PCA", dim=1:15, k.param = 10)

x_object = Seurat::FindClusters(x_object)

pred_label <- x_object@meta.data$seurat_clusters

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))

est_K = length(unique(pred_label))

Louvain = list(K = est_K, label = pred_label, ARI = ARI_value)


res_summary = list(SC.MEB1 = SC.MEB1, SC.MEB2 = SC.MEB2, DR.SC1 = DR.SC1, DR.SC2 = DR.SC2,
                   BayesSpace = BayesSpace, kmeans = k_means, GMM = GMM, Louvain = Louvain)

save(res_summary, file = "real_data/MOB_data/result_other_method.RData")


