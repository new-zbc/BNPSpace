library(Rcpp)
library(RcppArmadillo)
library(SingleCellExperiment)

sourceCpp("R/model.cpp")
sourceCpp("R/fix_clusters_model.cpp")


#### MCMC procedure
NSCFS <- function(sce, Adj, n_clusters = NULL, d=1, K_init = 15, n_iters = 100, seed = seed){
  
  set.seed(seed)
  data_set = as.matrix(t(assay(sce, "counts")))
  N = dim(data_set)[1]
  p = dim(data_set)[2]
  # hyperparameters
  GAMMA = 1
  a = 1
  b = 1
  GAMMA = 1
  aw = 0.05
  bw = 1.95
  ar = 0.1
  br = 1.9
  
  if(is.null(n_clusters)){
    #### initialize parameters
    K_start = K_init
    group_t <- kmeans(reducedDim(sce), centers = K_start, nstart = 5)$cluster
    
    gamma_t = c(rep(1,5), rep(0, p-5))
    D_t = exp(matrix(rnorm(p*K_start), p, K_start))
    D_0_t = exp(rnorm(p))
    R_t = matrix(0, N, p)
    R_t[data_set == 0] = 1
    pi_t = 1 - mean(R_t)
    
    
    result = runMCMC(group_t-1, D_t, gamma_t, D_0_t, R_t, pi_t, f=d, a, b, GAMMA, 
                     aw, bw, ar, br, data_set, Adj, sizeFactors(sce), n_iters)
    
    max_class = max(result$K_iter)
    prob_pred = matrix(0, N, max_class)
    for(i in 1:N){
      for(j in 1:max_class){
        prob_pred[i, j] = mean((result$group_iter[i, ] + 1) == j)
      }
    }
    
    pred_label = as.vector(apply(prob_pred, 1,  function(x) which(x == max(x))[1]))
    
  }else{
    #### initialize parameters
    K_start = n_clusters
    group_t <- kmeans(reducedDim(sce), centers = K_start, nstart = 5)$cluster
    
    gamma_t = c(rep(1,5), rep(0, p-5))
    D_t = exp(matrix(rnorm(p*K_start), p, K_start))
    D_0_t = exp(rnorm(p))
    R_t = matrix(0, N, p)
    R_t[data_set == 0] = 1
    pi_t = 1 - mean(R_t)
    
    result = run_fix_clusters(group_t-1, D_t, gamma_t, D_0_t, R_t, pi_t, f=d, a, b, 
                              aw, bw, ar, br, data_set, Adj, sizeFactors(sce), n_iters)
    
    prob_pred = matrix(0, N, n_clusters)
    for(i in 1:N){
      for(j in 1:n_clusters){
        prob_pred[i, j] = mean((result$group_iter[i, ] + 1) == j)
      }
    }
    
    pred_label = as.vector(apply(prob_pred, 1,  function(x) which(x == max(x))[1]))
    
  }
  
  gamma_est = rowMeans(result$gamma_iter[, (floor(n_iters/2)+1):n_iters])
  
  #### refinement of predicted label
  max_neighbors = dim(Adj)[2]
  label_refine = pred_label
  for(i in 1:length(pred_label)){
    index = Adj[i, which(Adj[i, ] != 0)]
    if(length(index) >= max_neighbors-1){
      neighbor_label = unique(pred_label[index])
      if(length(neighbor_label) == 1){
        label_refine[i] = neighbor_label
      }
    }
  }
  
  
  output = list()
  output$pred_label = label_refine
  output$gamma = gamma_est
  output$K = length(unique(label_refine))
  output$MCMCList = result
  
  
  return(output)
}




select_f <- function(sce, Adj, n_clusters = NULL, f=1, K_init = 15, n_iters = 100, seed = seed){
  
  set.seed(seed)
  data_set = as.matrix(t(assay(sce, "counts")))
  N = dim(data_set)[1]
  p = dim(data_set)[2]
  # hyperparameters
  GAMMA = 1
  a = 1
  b = 1
  GAMMA = 1
  aw = 0.05
  bw = 1.95
  ar = 0.1
  br = 1.9
  
  #### initialize parameters
  K_start = K_init
  group_t <- kmeans(reducedDim(sce), centers = K_start, nstart = 5)$cluster
  
  gamma_t = c(rep(1,5), rep(0, p-5))
  D_t = exp(matrix(rnorm(p*K_start), p, K_start))
  D_0_t = exp(rnorm(p))
  R_t = matrix(0, N, p)
  R_t[data_set == 0] = 1
  pi_t = 1 - mean(R_t)
  
  
  result = runMCMC(group_t-1, D_t, gamma_t, D_0_t, R_t, pi_t, f, a, b, GAMMA, 
                   aw, bw, ar, br, data_set, Adj, sizeFactors(sce), n_iters)
  
  max_class = max(result$K_iter)
  prob_pred = matrix(0, N, max_class)
  for(i in 1:N){
    for(j in 1:max_class){
      prob_pred[i, j] = mean((result$group_iter[i, ] + 1) == j)
    }
  }
  
  pred_label = as.vector(apply(prob_pred, 1,  function(x) which(x == max(x))[1]))
  gamma_est = rowMeans(result$gamma_iter[, (floor(n_iters/2)+1):n_iters])
  
  Dev_mean = mean(result$log_prob_iter[(floor(n_iters/2)+1):n_iters])
  
  D0_est = rowMeans(result$D_0_iter, (floor(n_iters/2)+1):n_iters)
  D_est = matrix(0, nrow = dim(result$D_iter), ncol = 30)
  for(k in 1:30){
    D_est[, k] = rowMeans(result$D_iter[,k,(floor(n_iters/2)+1):n_iters])
  }
  
  R_est = matrix(as.numeric(result$R_est > 0.5), nrow = N)
  
  s = sizeFactors(sce)
  Dev_est = 0
  for(i in 1:N){
    for(j in 1:p){
      if(R_est[i,j] == 0){
        if(gamma_est[j] == 1){Dev_est = Dev_est + dpois(data_set[i,j], lambda = s[i]*D_est[j, pred_label[i]], log = T)}
        else{Dev_est = Dev_est + dpois(data_set[i,j], lambda = s[i]*D0_est[j], log = T)}
      }
    }
  }
  
  pd = Dev_est - Dev_mean
  
  mDIC = -2 * Dev_est + log(N) * pd
  
  
  output = list()
  output$pred_label = pred_label
  output$gamma = gamma_est
  output$K = length(unique(pred_label))
  output$mDIC = mDIC
  output$MCMCList = result
  
  
  return(output)
}
