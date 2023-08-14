library(MASS)
library(SingleCellExperiment)
library(flexclust)
library(scater)
library(scran)

find_coordinate <- function(height, width){
  result = matrix(0, height*width, 2)
  for(i in 1:(height*width)){
    part1 = i %/% height
    part2 = i %% height
    
    if(part2 != 0){
      result[i, 1] = part1 + 1
      result[i, 2] = part2
    }
    else{
      result[i, 1] = part1
      result[i, 2] = part2 + height
    }
  }
  return(result)
}


data_folder = "simulation1/scenario3_2"
K = 7
p_gamma = 20 # the number of effective variables
p = 500 # total number of variables
prob = 0.2 # the drop out probability

height = 40
width = 40

load("simulation1/7_clusters_pattern.RData")

for(data_index in 1:50){
  
  set.seed(data_index)
  
  ground_truth = c(labels) + 1
  position = find_coordinate(height, width)
  
  ##### generate count
  N = height * width
  s = rep(0, N) 
  
  # True parameters
  non_effect = rgamma(p-p_gamma, 2, rate = 1)
  effect = rgamma(p_gamma, 2, rate = 1)
  d1 = c(effect, non_effect)
  d2 = c(effect + 3, non_effect)
  d3 = c(effect + 6, non_effect)
  d4 = c(effect + 10, non_effect)
  
  index1 = sample(p_gamma, 10)
  index2 = base::setdiff(1:p_gamma, index1)
  
  effect5 = effect
  effect5[index1] = effect5[index1] + 3
  d5 = c(effect5, non_effect)
  
  effect6 = effect
  effect6[index2] = effect6[index2] + 3
  d6 = c(effect6, non_effect)
  
  index3 = sample(p_gamma, 10)
  effect7 = effect
  effect7[index3] = effect7[index3] + 6
  d7 = c(effect7, non_effect)
  
  d = cbind(d1, d2, d3, d4, d5, d6, d7)
  
  X = matrix(NA, N, p)
  for(i in 1:N){
    s[i] = runif(1, 0.4, 1.6)
    for(j in 1:p){
      ind = rbinom(1,1, prob)
      if( ind == 1){
        X[i, j] = 0
      }else{
        if(ground_truth[i] == 1){X[i, j] = rpois(1, lambda = s[i]*d[j, 1])}
        if(ground_truth[i] == 2){X[i, j] = rpois(1, lambda = s[i]*d[j, 2])}
        if(ground_truth[i] == 3){X[i, j] = rpois(1, lambda = s[i]*d[j, 3])}
        if(ground_truth[i] == 4){X[i, j] = rpois(1, lambda = s[i]*d[j, 4])}
        if(ground_truth[i] == 5){X[i, j] = rpois(1, lambda = s[i]*d[j, 5])}
        if(ground_truth[i] == 6){X[i, j] = rpois(1, lambda = s[i]*d[j, 6])}
        if(ground_truth[i] == 7){X[i, j] = rpois(1, lambda = s[i]*d[j, 7])}
      }
    }
  }
  gamma_truth = c(rep(1, p_gamma), rep(0, p - p_gamma))
  
  library(SingleCellExperiment)
  size_factors = rowSums(X) / mean(rowSums(X))
  
  print(mean(X == 0))
  
  colnames(X) = paste0("gene_", 1:dim(X)[2])
  rownames(X) = paste0("spot_", 1:dim(X)[1])
  
  sce <- SingleCellExperiment(list(counts = t(X)), 
                              colData = DataFrame(row = position[, 1], 
                                                  col = position[, 2], 
                                                  label = ground_truth))
  
  sce <- logNormCounts(sce)
  sce <- scater::runPCA(sce, subset_row=1:500, ncomponents=15, 
                        exprs_values="logcounts")
  
  
  file_name = paste0(data_folder, "/", "data/", data_index, ".RData")
  save(sce, file = file_name)
  
}


