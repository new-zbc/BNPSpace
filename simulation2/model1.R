library(doParallel)
library(foreach)


data_folder = "simulation1/scenario1_1"
n.cluster = 10
n.data_set = 50

############################################
### proposed method
############################################

cl = makeCluster(n.cluster)
registerDoParallel(cl)
mydata1 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar%
  {
    library(Rcpp)
    library(SingleCellExperiment)
    library(RcppArmadillo)
    library(flexclust)

    filename = paste0(data_folder, "/data/", i, ".RData")
    load(file = filename)

    source("R/utils.R")

    ### construct spots network
    Adj = find_neighbors(sce, "ST", "image")
    
    neighbors = find_neighbor_index(Adj, "ST")
    result = NSCFS(sce, neighbors, n_clusters = NULL, f = 1, n_iters = 2000, seed = 1)
    library(flexclust)
    ARI_value = randIndex(table(result$pred_label, colData(sce)$label))

    file_name = paste0(data_folder, "/NSCFS/", i, ".RData")
    
    saving = list()
    saving$gamma = result$gamma
    saving$K_iter = result$MCMCList$K_iter
    saving$ARI = ARI_value
    
    save(saving, file = file_name)

    c(ARI_value, length(unique(result$pred_label)))
  }
stopCluster(cl)


file_name = paste0(data_folder, "/others/", "NSCFS", ".txt")
write.table(mydata1, file = file_name)


############################################
### SC.MEB(PCs = 15)
############################################

cl = makeCluster(n.cluster)
registerDoParallel(cl)
mydata3 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar%
  {
    filename = paste0(data_folder, "/data/", i, ".RData")
    load(file = filename)
    
    K_set = 2:10
    library(SC.MEB)
    library(SingleCellExperiment)
    Adj = find_neighbors2(sce, platform = "ST")
    
    data_mat = reducedDim(sce)
    
    selection = SC.MEB(data_mat, Adj_sp = Adj, parallel = FALSE)
    
    res = selectK(selection, K_set = K_set, criterion = "BIC", c = 1)
    
    est_K = res$best_K_BIC
    pred_label = res$best_K_label
    
    
    library(flexclust)
    ARI_value = randIndex(table(pred_label, colData(sce)$label))
    
    c(ARI_value, est_K)
  }
stopCluster(cl)


file_name = paste0(data_folder, "/others/", "SCMEB2", ".txt")
write.table(mydata3, file = file_name)



####################################################
######## DR.SC method (PCs = 15)
####################################################
cl = makeCluster(n.cluster)
registerDoParallel(cl)  
mydata5 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar% 
  {
    filename = paste0(data_folder, "/data/", i, ".RData")
    load(file = filename)
    
    library(DR.SC)
    library(SingleCellExperiment)
    library(SC.MEB)
    Adj = find_neighbors2(sce, platform = "ST")
    
    reslist = DR.SC_fit(t(assay(sce, "logcounts")), q = 15, K= 2:10, Adj_sp = Adj, coreNum = 1)
    
    res = selectModel(reslist, criteria = 'MBIC', pen.const=1)
    
    est_K = res$bestK
    pred_label = res$cluster
    
    library(flexclust)
    ARI_value = randIndex(table(pred_label, colData(sce)$label))
    
    c(ARI_value, est_K)
  }
stopCluster(cl)

file_name = paste0(data_folder, "/others/", "SC_RC2", ".txt")
write.table(mydata5, file = file_name)


############################################
### BayesSpace(K  = K)
############################################

cl = makeCluster(n.cluster)
registerDoParallel(cl)  
mydata6 <- foreach(i=1:n.data_set) %dopar% 
  {
    filename = paste0(data_folder, "/data/", i, ".RData")
    load(file = filename)
    library(flexclust)
    library(BayesSpace)
    K = length(unique(colData(sce)$label))
    sce_result = BayesSpace::spatialCluster(sce, q=K, use.dimred = "PCA", 
                                            platform = "ST", init.method = "kmeans", model = "normal", 
                                            precision = "equal", nrep = 5000, gamma = 2, 
                                            alpha = 1, beta = 0.01, save.chain = FALSE)
    
    ARI_value = randIndex(table(colData(sce_result)$spatial.cluster, colData(sce)$label))
    
    ARI_value
  }
stopCluster(cl)

mydata6 = unlist(mydata6)
file_name = paste0(data_folder, "/others/", "BayesSpace", ".txt")
write.table(mydata6, file = file_name)




######################################
### K-means
#####################################
cl = makeCluster(n.cluster)
registerDoParallel(cl)  
mydata7 <- foreach(i=1:n.data_set) %dopar% 
  {
    filename = paste0(data_folder, "/data/", i, ".RData")
    load(file = filename)
    
    library(SingleCellExperiment)
    data_mat = reducedDim(sce)
    K = length(unique(colData(sce)$label))
    x_gmm <- kmeans(data_mat, centers = K, nstart = 5)$cluster
    
    library(flexclust)
    ARI_value = randIndex(table(x_gmm, colData(sce)$label))
    
    ARI_value
  }
stopCluster(cl)

mydata7 = unlist(mydata7)

file_name = paste0(data_folder, "/others/", "kmeans", ".txt")
write.table(mydata7, file = file_name)


######################################
### GMM
#####################################
cl = makeCluster(n.cluster)
registerDoParallel(cl)  
mydata8 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar% 
  {
    filename = paste0(data_folder, "/data/", i, ".RData")
    load(file = filename)
    
    library(SingleCellExperiment)
    library(mclust)
    data_mat = reducedDim(sce)
    
    fit_int = Mclust(data_mat, G = 2:10)
    
    K_est = fit_int$G
    
    x_gmm <- fit_int$classification
    
    library(flexclust)
    ARI_value = randIndex(table(x_gmm, colData(sce)$label))
    
    c(ARI_value, K_est)
  }
stopCluster(cl)


file_name = paste0(data_folder, "/others/", "GMM", ".txt")
write.table(mydata8, file = file_name)



######################################
### Louvain
#####################################
cl = makeCluster(n.cluster)
registerDoParallel(cl)  
mydata9 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar% 
  {
    filename = paste0(data_folder, "/data/", i, ".RData")
    load(file = filename)
    
    library(SingleCellExperiment)
    library(Seurat)
    library(flexclust)
    x_object = as.Seurat(sce, counts = "counts", data = "logcounts", project =  "sce_to_seurat")
    
    x_object = Seurat::FindNeighbors(x_object, reduction = "PCA", dim=1:15)
    
    x_object = Seurat::FindClusters(x_object)
    ARI_value = randIndex(table(x_object@meta.data$seurat_clusters, colData(sce)$label))
    
    K_est = length(unique(x_object@meta.data$seurat_clusters))
    c(ARI_value, K_est)
  }
stopCluster(cl)

file_name = paste0(data_folder, "/others/", "Louvain", ".txt")
write.table(mydata9, file = file_name)



############################################
### proposed method f = 0
############################################

cl = makeCluster(n.cluster)
registerDoParallel(cl)
mydata10 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar%
  {
    library(Rcpp)
    library(SingleCellExperiment)
    library(RcppArmadillo)
    library(flexclust)
    
    filename = paste0(data_folder, "/data/", i, ".RData")
    load(file = filename)
    
    source("R/utils.R")
    
    ### construct spots network
    Adj = find_neighbors(sce, "ST", "image")
    neighbors = find_neighbor_index(Adj, "ST")
    result = NSCFS(sce, neighbors, n_clusters = NULL, f = 0, n_iters = 2000, seed = 1)
    library(flexclust)
    ARI_value = randIndex(table(result$pred_label, colData(sce)$label))
    
    file_name = paste0(data_folder, "/NSCFS_0/", i, ".RData")
    
    saving = list()
    saving$gamma = result$gamma
    saving$K_iter = result$MCMCList$K_iter
    saving$ARI = ARI_value
    save(saving, file = file_name)
    
    c(ARI_value, length(unique(result$pred_label)))
  }
stopCluster(cl)


file_name = paste0(data_folder, "/others/", "NSCFS_0", ".txt")
write.table(mydata10, file = file_name)


