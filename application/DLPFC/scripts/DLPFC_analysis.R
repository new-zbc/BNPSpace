library(SingleCellExperiment)
library(SC.MEB)
library(Matrix)
library(flexclust)
library(scater)
library(scran)
source("R/utils.R")
source("R/main.R")

sampleID = 151509


data_folder = "application/DLPFC"

file_name = paste0(data_folder, "/data/", sampleID, "_counts.RData")

load(file_name)


### quality control
sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
### construct spots network
Adj = find_neighbors(sce2, "Visium", "lattice")
neighbors = find_neighbor_index(Adj, "Visium")

seed = 2022022

####################################################################
######################### proposed method  #############################
###################################################################
#### run NSCFS
timestart = Sys.time()
result1 = NSCFS(sce2, neighbors, n_clusters = NULL, f = 1.5, n_iters = 500, seed = 3)
timeend = Sys.time()
timerun = timeend - timestart
print(timerun)

ground_truth = colData(sce2)$label
pred_spot = result1$pred_label
ARI_value = randIndex(table(pred_spot, ground_truth))
print(ARI_value)

file_name = paste0("application/DLPFC/", sampleID, "/","MCMCresults.RData")
save(result1, file = file_name)




####################################################################
######################### SC.MEB  #################################
###################################################################
library(SC.MEB)
set.seed(seed)
K_set = 7
data_mat = reducedDim(sce2)

timestart = Sys.time()
selection = SC.MEB(data_mat, Adj_sp = Adj, parallel = FALSE, K_set = K_set)
timeend = Sys.time()
timerun = timeend - timestart
print(timerun)

res = selectK(selection, K_set = K_set, criterion = "BIC", c = 1)

est_K = res$best_K_BIC
pred_label = res$best_K_label

ground_truth = colData(sce2)$label
pred_spot = pred_label
ARI_value = randIndex(table(pred_spot, ground_truth))

SC.MEB = list(K = est_K, label = pred_label, ARI = ARI_value)




####################################################################
######################### DR.SC   #################################
###################################################################
library(DR.SC)
set.seed(seed)
K_set = 7

sce3 = spatialFilter(sce1, cutoff_sample = 100, cutoff_feature = 0, cutoff_max = 0)

data_mat = t(assay(sce3, "logcounts")[rowData(sce3)$is.HVG, ])

timestart = Sys.time()
reslist = DR.SC_fit(data_mat, q = 15, K= K_set, Adj_sp = Adj, coreNum = 1)
timeend = Sys.time()
timerun = timeend - timestart
print(timerun)

res = selectModel(reslist, criteria = 'MBIC', pen.const=1)

est_K = res$bestK
pred_label = res$cluster

ground_truth = colData(sce3)$label
pred_spot = pred_label
ARI_value = randIndex(table(pred_spot, ground_truth))

DR.SC = list(K = est_K, label = pred_label, ARI = ARI_value)


####################################################################
######################### BayesSpace  #################################
###################################################################
library(BayesSpace)
set.seed(seed)

# qPlot(sce2, qs = seq(2, 20), burn.in = 100, nrep = 1000, use.dimred = "PCA",
#       platform = "Visium", init.method = "kmeans", model = "normal",
#       precision = "equal", gamma = 3, alpha = 1, beta = 0.01)
timestart = Sys.time()
sce_result = BayesSpace::spatialCluster(sce2, q=7, use.dimred = "PCA",
                                        platform = "Visium", init.method = "kmeans", model = "normal",
                                        precision = "equal", nrep = 5000, gamma = 3,
                                        alpha = 1, beta = 0.01, save.chain = FALSE)
timeend = Sys.time()
timerun = timeend - timestart
print(timerun)

pred_label = colData(sce_result)$spatial.cluster

ground_truth = colData(sce2)$label
pred_spot = pred_label
ARI_value = randIndex(table(pred_spot, ground_truth))

BayesSpace = list(label = pred_label, ARI = ARI_value)



# ####################################################################
# ######################### k-means  #################################
# ###################################################################
# set.seed(seed)
# data_mat = reducedDim(sce2)
# pred_label <- kmeans(data_mat, centers = 7, nstart = 5)$cluster
# 
# ground_truth = colData(sce2)$label
# pred_spot = pred_label
# ARI_value = randIndex(table(pred_spot, ground_truth))
# 
# k_means = list(label = pred_label, ARI = ARI_value)
# 
# 
# ####################################################################
# ######################### GMM method  #############################
# ###################################################################
# set.seed(seed)
# library(mclust)
# data_mat = reducedDim(sce2)
# 
# fit_int = Mclust(data_mat, G = 2:20)
# 
# est_K = fit_int$G
# 
# pred_label <- fit_int$classification
# 
# ground_truth = colData(sce2)$label
# pred_spot = pred_label
# ARI_value = randIndex(table(pred_spot, ground_truth))
# 
# GMM = list(K = est_K, label = pred_label, ARI = ARI_value)


#################################################################
############################# Louvain  ##########################
##################################################################
set.seed(seed)
library(Seurat)
timestart = Sys.time()
x_object = as.Seurat(sce2, counts = "counts", data = "logcounts", project =  "sce_to_seurat")

x_object = Seurat::FindNeighbors(x_object, reduction = "PCA", dim=1:15)

x_object = Seurat::FindClusters(x_object)
timeend = Sys.time()
timerun = timeend - timestart
print(timerun)

pred_label <- x_object@meta.data$seurat_clusters

ground_truth = colData(sce2)$label
pred_spot = pred_label
ARI_value = randIndex(table(pred_spot, ground_truth))

est_K = length(unique(pred_label))

Louvain = list(K = est_K, label = pred_label, ARI = ARI_value)


res_summary = list(SC.MEB = SC.MEB,  DR.SC = DR.SC, 
                   BayesSpace = BayesSpace, Louvain = Louvain)

file_name = paste0("application/DLPFCdata/","result_other_method.RData")
save(res_summary, file = file_name)


