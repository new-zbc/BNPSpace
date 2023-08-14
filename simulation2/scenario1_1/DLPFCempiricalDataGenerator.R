library(SingleCellExperiment)
library(Matrix)
library(flexclust)
library(scater)
library(scran)
source("R/utils.R")

sampleID = 151676

data_folder = paste0("application/DLPFCdata/", sampleID)

file_name = paste0(data_folder, "/", sampleID, "_counts.RData")

load(file_name)

sce1 = sce1[, !is.na(colData(sce1)$label)]

sce1 = spatialPreprocess(sce1)


data_folder = "simulation2/scenario1_1"
K = 7
prob = 0
p_gamma = 50 # the number of effective variables
p = 500 # total number of variables

for(data_index in 1:50){
  
  set.seed(data_index)
  
  ground_truth = as.character(colData(sce1)$label)
  s = colData(sce1)$sizeFactor
  N = length(ground_truth)
  
  
  # True parameters
  non_effect = rgamma(p-p_gamma, 2, rate = 1)
  effect = rgamma(p_gamma, 2, rate = 1)
  d1 = c(effect, non_effect)
  d2 = c(effect + 3, non_effect)
  d3 = c(effect + 6, non_effect)
  d4 = c(effect + 11, non_effect)
  
  index1 = sample(p_gamma, 25)
  index2 = base::setdiff(1:p_gamma, index1)
  
  effect5 = effect
  effect5[index1] = effect5[index1] + 3
  d5 = c(effect5, non_effect)
  
  effect6 = effect
  effect6[index2] = effect6[index2] + 3
  d6 = c(effect6, non_effect)
  
  index3 = sample(p_gamma, 25)
  effect7 = effect
  effect7[index3] = effect7[index3] + 6
  d7 = c(effect7, non_effect)
  
  d = cbind(d1, d2, d3, d4, d5, d6, d7)
  
  ##### generate count
  X = matrix(NA, N, p)
  for(i in 1:N){
    for(j in 1:p){
      ind = rbinom(1,1, prob)
      if( ind == 1){
        X[i, j] = 0
      }else{
        if(ground_truth[i] == "Layer1"){X[i, j] = rpois(1, lambda = s[i]*d[j, 1])}
        if(ground_truth[i] == "Layer2"){X[i, j] = rpois(1, lambda = s[i]*d[j, 2])}
        if(ground_truth[i] == "Layer3"){X[i, j] = rpois(1, lambda = s[i]*d[j, 3])}
        if(ground_truth[i] == "Layer4"){X[i, j] = rpois(1, lambda = s[i]*d[j, 4])}
        if(ground_truth[i] == "Layer5"){X[i, j] = rpois(1, lambda = s[i]*d[j, 5])}
        if(ground_truth[i] == "Layer6"){X[i, j] = rpois(1, lambda = s[i]*d[j, 6])}
        if(ground_truth[i] == "WM"){X[i, j] = rpois(1, lambda = s[i]*d[j, 7])}
      }
    }
  }
  
  colnames(X) = paste0("gene_", 1:dim(X)[2])
  rownames(X) = paste0("spot_", 1:dim(X)[1])
  
  print(mean(X == 0))
  
  sce <- SingleCellExperiment(list(counts = t(X)), 
                              colData = DataFrame(row = colData(sce1)$row, 
                                                  col = colData(sce1)$col, 
                                                  label = ground_truth))
  
  ## spatial preprocess
  ## conduct PCA and compute size factors
  sce = spatialPreprocess(sce)
  
  
  file_name = paste0(data_folder, "/", "data/", data_index, ".RData")
  save(sce, file = file_name)
  
}

