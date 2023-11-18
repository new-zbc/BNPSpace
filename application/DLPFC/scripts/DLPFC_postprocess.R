postProcess <- function(gene_exp, pred_label, ground_truth){
  # gene_exp should be K by p matrix
  library(flexclust)
  group_distance = dist(gene_exp)
  fit = hclust(group_distance, method = "average")
  group_megre = fit$merge
  
  step = dim(group_megre)[1]
  ARI = rep(0, step)
  
  N = length(ground_truth)
  label_step = matrix(NA, N, step)
  
  pred_label = - pred_label
  for(i in 1:step){
    pred_label[pred_label %in% group_megre[i, ]] = i
    ARI[i] = randIndex(table(pred_label, ground_truth))
    label_step[, i] = pred_label
  }
  result = cbind(K= step:1, group_megre, ARI)
  index = which.max(ARI)
  list(ARI = result, pred = label_step[, index], all = label_step)
}





DLPFC_ARI <- function(sampleID){
  
  data_folder = "application/DLPFCdata"
  
  file_name = paste0(data_folder, "/", sampleID, "_counts.RData")
  
  load(file_name)
  
  
  ### quality control
  sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
  ### construct spots network
  Adj = find_neighbors(sce2, "Visium", "lattice")
  neighbors = find_neighbor_index(Adj, "Visium")
  
  load(paste0(data_folder, "/", "MCMCresults.RData"))
  result1 = res_summary$NSCFS
  
  ground_truth = colData(sce2)$label
  pred_spot = result1$pred_label
  ARI_value = randIndex(table(pred_spot, ground_truth))
  temp = c(length(unique(pred_spot)), 0, 0, ARI_value)
  
  gamma = res_summary$NSCFS$gamma > 0.5
  
  clusters = unique(pred_spot)
  K = length(clusters)
  d_list = list()
  for(k in 1:K){
    d = rowMeans(res_summary$NSCFS$MCMCList$D_iter[gamma, k, 300:1000])
    d_list = append(d_list, list(d))
  }
  
  gene_exp = as.matrix(as.data.frame(d_list))
  colnames(gene_exp) = 1:K
  
  
  result = postProcess(t(gene_exp), pred_spot, ground_truth)
  
  result$ARI = rbind(temp, result$ARI)
  
  return(result)
  
}
