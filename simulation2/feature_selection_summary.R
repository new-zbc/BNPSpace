#### Assessment function for variable selection #######
selection.perf <- function(gamma, gamma_truth, threshold = 0.5)
{
  library(ROCR)
  p = length(gamma)
  subset.selection = rep(0, p)
  subset.selection[gamma >= threshold] = 1
  if(sum(subset.selection) >= p)
  {
    if(sum(gamma) == p){
      return(c(sensitivity = 1, 
               Specificity = 0, 
               MCC = NA, 
               AUC = NA))
    }
    else{
      pred.auc = prediction(gamma, gamma_truth)
      perf <- performance(pred.auc,'auc')
      auc = perf@y.values[[1]]
      return(c(sensitivity = 1, 
               Specificity = 0, 
               MCC = NA, 
               AUC = auc))
    }
    
  }
  
  if(sum(subset.selection) == 0)
  {
    if(sum(gamma) == 0){
      return(c(sensitivity = 0, 
               Specificity = 1, 
               MCC = NA, 
               AUC = NA))
    }
    else{
      pred.auc = prediction(gamma, gamma_truth)
      perf <- performance(pred.auc,'auc')
      auc = perf@y.values[[1]]
      return(c(sensitivity = 1, 
               Specificity = 0, 
               MCC = NA, 
               AUC = auc))
    }
  }
  
  result = table(subset.selection, gamma_truth)
  sens = result[2,2]/sum(result[,2])
  spec = result[1,1]/sum(result[,1])
  mcc = (result[2,2]*result[1,1] - result[1,2]*result[2,1]) / sqrt(sum(result[1,]) * sum(result[2,]))
  mcc = mcc /  sqrt(sum(result[,1]) * sum(result[,2]))
  pred.auc = prediction(gamma, gamma_truth)
  perf <- performance(pred.auc,'auc')
  auc = perf@y.values[[1]]
  return(c(sensitivity = sens, 
           Specificity = spec, 
           MCC = mcc, 
           AUC = auc))
}


BayesianFDR <- function(PPI){
  ka = seq(0, 0.99, length.out = 100)
  result = rep(0, 100)
  for(i in 1:length(ka)){
    term1 = sum((1-PPI) * as.numeric(1 - PPI < ka[i]))
    term2 = sum(as.numeric(1 - PPI < ka[i]))
    if(term2 > 0){result[i] = term1/term2}
  }
  if(length(which(result > 0.05)) > 0){out = ka[min(which(result > 0.05))]}
  else{out = 0.99}
  return(out)
}


BayesianFDR(saving$gamma)




library(doParallel)
library(foreach)

data_folder = "simulation1/scenario1_1"
n.cluster = 5
n.data_set = 50

############################################
### proposed method
############################################

gamma_truth = rep(0, 500)
gamma_truth[1:20] = 1

mydata <- foreach(i=1:n.data_set, .combine = rbind) %do% 
  {
    
    library(SingleCellExperiment)
    
    filename1 = paste0(data_folder, "/data/", i, ".RData")
    load(file = filename1)
    
    filename2 = paste0(data_folder, "/NSCFS/", i, ".RData")
    load(file = filename2)
    
    selection_res = selection.perf(saving$gamma, gamma_truth = gamma_truth, threshold = 0.5)
    selection_res
  }


row.names(mydata) = 1:dim(mydata)[1]
file_name = paste0(data_folder, "/feature_selection_result.csv")
write.csv(mydata, file = file_name, quote = FALSE)

vec1 = colMeans(mydata, na.rm = T)
vec2 = apply(mydata, 2, sd, na.rm = T)

print(vec1)
print(vec2)



data_folder = "simulation1/scenario1_1"
n.cluster = 5
n.data_set = 50

############################################
### proposed method
############################################

gamma_truth = rep(0, 500)
gamma_truth[1:20] = 1

mydata <- foreach(i=1:n.data_set, .combine = rbind) %do% 
  {
    
    library(SingleCellExperiment)
    
    filename1 = paste0(data_folder, "/data/", i, ".RData")
    load(file = filename1)
    
    filename2 = paste0(data_folder, "/NSCFS/", i, ".RData")
    load(file = filename2)
    
    selection_res = selection.perf(saving$gamma, gamma_truth = gamma_truth, threshold = BayesianFDR(saving$gamma))
    selection_res
  }


row.names(mydata) = 1:dim(mydata)[1]
file_name = paste0(data_folder, "/feature_selection_result_2.csv")
write.csv(mydata, file = file_name, quote = FALSE)

vec1 = colMeans(mydata, na.rm = T)
vec2 = apply(mydata, 2, sd, na.rm = T)

print(vec1)
print(vec2)
