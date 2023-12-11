folders = c("simulation1/scenario1_1", "simulation1/scenario1_2", "simulation1/scenario1_3", 
            "simulation1/scenario2_1", "simulation1/scenario2_2", "simulation1/scenario2_3", 
            "simulation1/scenario3_1", "simulation1/scenario3_2", "simulation1/scenario3_3")

result_mean = matrix(0, length(folders), 4)
result_sd = matrix(0, length(folders), 4)

for(i in 1:length(folders)){
  
  folder = folders[i]
  
  file = paste0(folder, "/feature_selection_result_2.csv")
  data = read.csv(file = file, row.names = 1)
  
  
  result_mean[i, ] = colMeans(data, na.rm = T)
  result_sd[i, ] = apply(data, 2, sd, na.rm = T)
}

res = data.frame(pi = rep(c(0.1, 0.2, 0.3), 3), K = rep(c(3, 5, 7), each = 3))

res_mean = as.data.frame(cbind(res, result_mean))
res_sd = as.data.frame(cbind(res, result_sd))

colnames(res_mean) = c("pi", "K", "Sensitivity", "Specificity", "MCC", "AUC")

colnames(res_sd) = c("pi", "K", "Sensitivity", "Specificity", "MCC", "AUC")

write.csv(res_mean, file = "reproduce/figures_and_tables/tableS1mean.csv")

write.csv(res_sd, file = "reproduce/figures_and_tables/tableS1sd.csv")
