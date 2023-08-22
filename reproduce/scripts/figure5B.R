
plot_ARI = function(sampleID){
  library(ggplot2)
  library(ggrepel)
  source("R/utils.R")
  data_folder = "application/DLPFCdata"
  file_name = paste0(data_folder, "/", "different_K_ARI.RData")
  load(file_name)
  
  colnames(ARI_K_results) = c("K", "SCMEB(PCs=10)", "SCMEB(PCs=15)", "DR.SC(PCs=10)", "DR.SC(PCs=15)", "BayesSpace", "kmeans", "GMM")
  ARI_K_results = ARI_K_results[,-c(2, 4)]
  colnames(ARI_K_results) = c("K", "SCMEB", "DRSC", "BayesSpace", "Kmeans", "GMM")
  
  library(reshape)
  ARI_K_results_melt = melt(ARI_K_results, id.vars = "K", variable_name = "Method")
  colnames(ARI_K_results_melt) = c("K", "Method", "ARI")
  
  load(paste0(data_folder, "/", "BNPSpace_results.RData"))
  result_ARI = BNPSpace_res$post
  
  result_ARI = as.data.frame(result_ARI$ARI[, c(1, 4)])
  result_ARI$Method = rep("BNPSpace", dim(result_ARI)[1])
  result_ARI = result_ARI[-c(dim(result_ARI)[1], dim(result_ARI)[1]-1), ]
  colnames(result_ARI) = c("K", "ARI", "Method")
  
  result_ARI$ARI = round(result_ARI$ARI, 3)
  
  result_ARI = data.frame(K = result_ARI$K, Method = result_ARI$Method, ARI = result_ARI$ARI)
  
  ARI_K_results_melt = rbind(result_ARI, ARI_K_results_melt)
  
  ARI_K_results_melt = ARI_K_results_melt[!(ARI_K_results_melt$K %in% c(2, 16:20)),]
  
  index = which.max(result_ARI$ARI)
  
  K = length(unique(ARI_K_results_melt$K))
  df = data.frame(K = 3:15, y = rep(min(ARI_K_results_melt$ARI) - 0.02, K), xend = 3:15, yend = rep(0, K))
  for(k in unique(ARI_K_results_melt$K)){
    df$yend[which(df$K == k)] = max(ARI_K_results_melt$ARI[ARI_K_results_melt$K == k])
  }
  
  ARI_K_results_melt$Method = factor(ARI_K_results_melt$Method, levels = c("BNPSpace", "BayesSpace", 
                                                                           "DRSC", "SCMEB", 
                                                                           "Kmeans", "GMM"))
  
  p <- ggplot() + geom_point(data = ARI_K_results_melt, aes(x = K, y = ARI, color = Method, shape = Method), size = 5)
  #geom_point() + geom_line(linetype = "dashed")
  #p <- p + geom_segment(aes(x = K, y = y, xend = xend, yend = yend), data = df,  alpha = 0.4)
  
  p <- p + geom_line(data = ARI_K_results_melt[ARI_K_results_melt$Method == "BNPSpace",], aes(x = K, y = ARI, color = Method))
  
  p <- p + 
    theme_bw() + theme(axis.title = element_text(size = 15, face = "bold"),
                       axis.text = element_text(size = 15, face = "bold"),
                       legend.text = element_text(size = 15, face = "bold"),
                       legend.title = element_text(size = 18, face='bold'),
                       legend.justification = c(1,1),
                       legend.background = element_rect(colour = "black", fill = NA),
                       legend.position = c(1, 1)) + 
    scale_color_manual(values = c("red", "darkgreen", "hotpink", "darkgray", "orange", "purple"))
  p
}


pARI = plot_ARI(151509)
ggsave(pARI, filename = "reproduce/figures_and_tables/figure5B.png", width = 11, height = 8, bg = "white")
