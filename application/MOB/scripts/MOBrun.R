library(SingleCellExperiment)
library(Matrix)
library(flexclust)
library(scater)
library(scran)
source("R/utils.R")
source("R/main.R")

set.seed(1)
load("application/MOB/data/MOB_sce.RData")
raw_sce = sce
### quality control
sce = spatialFilter(sce, cutoff_max = 10)

### construct spots network
Adj = find_neighbors(sce, "ST", "lattice")
neighbors = find_neighbor_index(Adj, "ST")

seed = 1

#load("real_data/MOB_data/clustering_result_sce.RData")

####################################################################
######################### SC.MEB 1  #################################
###################################################################
library(SC.MEB)
set.seed(seed)
K_set = 4
data_mat = reducedDim(sce)[,1:15]

timestart = Sys.time()
selection = SC.MEB(data_mat, Adj_sp = Adj, parallel = FALSE, K_set = K_set)
timeend = Sys.time()
timerun = timeend - timestart
print(timerun)

res = selectK(selection, K_set = K_set, criterion = "BIC", c = 1)

est_K = res$best_K_BIC
pred_label = res$best_K_label

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))

colData(sce)$SCMEB = pred_label


####################################################################
######################### DR.SC   #################################
###################################################################
library(DR.SC)
set.seed(seed)
K_set = 4

sce_drsc = spatialPreprocess(raw_sce)
sce_drsc = spatialFilter(sce_drsc, cutoff_feature = 0, cutoff_max = 0)


data_mat = t(assay(sce_drsc, "logcounts")[rowData(sce_drsc)$is.HVG, ])

timestart = Sys.time()
reslist = DR.SC_fit(data_mat, q = 10, K= K_set, Adj_sp = Adj, coreNum = 1)
res = selectModel(reslist, criteria = 'MBIC', pen.const=1)
timeend = Sys.time()
timerun = timeend - timestart
print(timerun)

est_K = res$bestK
pred_label = res$cluster

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))

colData(sce)$DRSC = pred_label

#DR.SC2 = list(K = est_K, label = pred_label, ARI = ARI_value)


####################################################################
######################### BayesSpace  #################################
###################################################################
library(BayesSpace)
set.seed(seed)

# qPlot(sce, qs = seq(2, 20), burn.in = 100, nrep = 1000, use.dimred = "PCA", 
#       platform = "ST", init.method = "kmeans", model = "normal", 
#       precision = "equal", gamma = 2, alpha = 1, beta = 0.01)

timestart = Sys.time()
sce_result = BayesSpace::spatialCluster(sce, q=4, use.dimred = "PCA", d = 15,
                                        platform = "ST", init.method = "kmeans", model = "normal", 
                                        precision = "equal", nrep = 5000, gamma = 2, 
                                        alpha = 1, beta = 0.01, save.chain = FALSE)
timeend = Sys.time()
timerun = timeend - timestart
print(timerun)

pred_label = colData(sce_result)$spatial.cluster

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))

colData(sce)$BayesSpace = pred_label
#BayesSpace = list(label = pred_label, ARI = ARI_value)


# ####################################################################
# ######################### k-means  #################################
# ###################################################################
# set.seed(seed)
# data_mat = reducedDim(sce)
# 
# # K_set = 2:20
# # wss = rep(0, length(K_set))
# # for(i in 1:length(K_set)){
# #   wss[i] = sum(kmeans(data_mat, centers = K_set[i], nstart = 5)$withinss)
# # }
# # 
# # plot(K_set, wss, type = "b")
# 
# pred_label <- kmeans(data_mat, centers = 4, nstart = 3)$cluster
# 
# ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
# pred_spot = pred_label[!is.na(colData(sce)$label)]
# ARI_value = randIndex(table(pred_spot, ground_truth))
# 
# colData(sce)$Kmeans = pred_label
# #k_means = list(label = pred_label, ARI = ARI_value)
# 
# 
# ####################################################################
# ######################### GMM method  #############################
# ###################################################################
# set.seed(seed)
# library(mclust)
# data_mat = reducedDim(sce)
# 
# fit_int = Mclust(data_mat, G = 4)
# 
# est_K = fit_int$G
# 
# pred_label <- fit_int$classification
# 
# ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
# pred_spot = pred_label[!is.na(colData(sce)$label)]
# ARI_value = randIndex(table(pred_spot, ground_truth))
# 
# colData(sce)$GMM = pred_label
# 
# #GMM = list(K = est_K, label = pred_label, ARI = ARI_value)



#################################################################
############################# Louvain  ##########################
##################################################################
set.seed(1)
library(Seurat)
x_object = as.Seurat(sce, counts = "counts", data = "logcounts", project =  "sce_to_seurat")

timestart = Sys.time()
x_object = Seurat::FindNeighbors(x_object, reduction = "PCA", dim=1:10, k.param = 20)

x_object = Seurat::FindClusters(x_object, resolution = 0.3)
timeend = Sys.time()
timerun = timeend - timestart
print(timerun)

pred_label <- x_object@meta.data$seurat_clusters

ground_truth = colData(sce)$label[!is.na(colData(sce)$label)]
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth))

est_K = length(unique(pred_label))

colData(sce)$Louvain = pred_label

#Louvain = list(K = est_K, label = pred_label, ARI = ARI_value)


source("R/main.R")
set.seed(1)
####################################################################
######################### proposed method  #############################
###################################################################
#### run NSCFS
timestart = Sys.time()
result1 = NSCFS(sce, neighbors, n_clusters = NULL, f = 1, K_init = 10, n_iters = 300, seed = 1)
timeend = Sys.time()
timerun = timeend - timestart
print(timerun)


ground_truth = colData(sce)$label
pred_label = result1$pred_label
pred_spot = pred_label[!is.na(colData(sce)$label)]
ARI_value = randIndex(table(pred_spot, ground_truth[!is.na(colData(sce)$label)]))
ARI_value


colData(sce)$BNPSpace = pred_label
rowData(sce)$gamma = result1$gamma

save(sce, file = "application/MOB/clustering_result_sce.RData")
