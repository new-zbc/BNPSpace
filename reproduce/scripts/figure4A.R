source("R/utils.R")
color_assign <- function(pred, ground_truth, color_pal){
  R = length(unique(pred))
  C = length(unique(ground_truth))
  p_table = table(pred, ground_truth)
  size_order = order(rowSums(p_table))
  
  p_table = p_table / matrix(rep(rowSums(p_table), each = C), byrow = T,ncol = C)
  
  color_pred = rep(NA, R)
  layer_used = rep(FALSE, C)
  count = 0
  for(i in 1:R){
    #i = which(size_order == k)
    index = which.max(p_table[i, ])
    if(!layer_used[index] & p_table[i, index] > 0.3){
      layer_used[index] = TRUE
      color_pred[i] = color_pal[index]
    }
    else{
      color_pred[i] = color_pal[C+1+count]
      count = count+1
    }
  }
  return(color_pred)
}



plot_cluster <- function(sampleID){
  
  color_pal = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF", "#2CA02CFF", "#9467BDFF",
                "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
                "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
                "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
  #show_col(color_pal)
  data_folder = "application/DLPFC"
  file_name = paste0(data_folder, "/data/", sampleID, "_counts.RData")
  load(file_name)
  
  ### quality control
  sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
  ### construct spots network
  Adj = find_neighbors(sce2, "Visium", "lattice")
  neighbors = find_neighbor_index(Adj, "Visium")
  
  load(paste0(data_folder, "/", "result_other_method.RData"))
  
  library(flexclust)
  spaGCN_ARI = randIndex(colData(sce2)$spaGCN_7, colData(sce2)$label)
  
  
  library(BayesSpace)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  
  metadata(sce1)$BayesSpace.data <- list()
  metadata(sce1)$BayesSpace.data$platform <- "Visium"
  metadata(sce1)$BayesSpace.data$is.enhanced <- FALSE
  
  metadata(sce2)$BayesSpace.data <- list()
  metadata(sce2)$BayesSpace.data$platform <- "Visium"
  metadata(sce2)$BayesSpace.data$is.enhanced <- FALSE
  
  si <- 12; tsi <- 16
  
  colData(sce2)$label = as.character(colData(sce2)$label)
  colData(sce2)$label[colData(sce2)$label == "Layer1"] = "L1"
  colData(sce2)$label[colData(sce2)$label == "Layer2"] = "L2"
  colData(sce2)$label[colData(sce2)$label == "Layer3"] = "L3"
  colData(sce2)$label[colData(sce2)$label == "Layer4"] = "L4"
  colData(sce2)$label[colData(sce2)$label == "Layer5"] = "L5"
  colData(sce2)$label[colData(sce2)$label == "Layer6"] = "L6"
  colData(sce2)$label = factor(colData(sce2)$label, levels = c("L1", "L2", "L3", "L4", "L5", "L6", "WM"))
  
  p1 <- clusterPlot(sce2, label=colData(sce2)$label, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,1]) +
    labs(title="Mannual annotation", fill = "Domains") + scale_fill_manual(values = color_pal[1:length(unique(colData(sce2)$label))]) +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  
  colData(sce2)$BayesSpace = as.factor(res_summary$BayesSpace$label)
  
  p4 <- clusterPlot(sce2, label= colData(sce2)$BayesSpace, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("BayesSpace: ARI=", round(res_summary$BayesSpace$ARI, 3)), fill = "Domains") + 
    scale_fill_manual(values = color_assign(colData(sce2)$BayesSpace, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  
  
  spaGCN_7 = as.numeric(colData(sce2)$spaGCN_7) + 1
  spaGCN_7[spaGCN_7 == 1] = 0
  spaGCN_7[spaGCN_7 == 7] = 1
  spaGCN_7[spaGCN_7 == 0] = 7
  
  spaGCN_7[spaGCN_7 == 3] = 0
  spaGCN_7[spaGCN_7 == 4] = 3
  spaGCN_7[spaGCN_7 == 0] = 4
  
  spaGCN_7 = as.factor(spaGCN_7)
  p6 <- clusterPlot(sce2, label= spaGCN_7, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("SpaGCN: ARI=", round(spaGCN_ARI, 3)), fill = "Domains") + 
    scale_fill_manual(values = color_assign(spaGCN_7, colData(sce2)$label, color_pal)) +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) 
 
  
  load(paste0(data_folder, "/", "BNPSpace_results.RData"))
  result_ARI = BNPSpace_res$post
  result_label = result_ARI$pred
  clusters = as.numeric(names(table(result_label)))
  pred_label = rep(NA, length(result_label))
  for(i in 1:length(result_label)){
    for(k in 1:length(clusters)){
      if(result_label[i] == clusters[k]){pred_label[i] = k}
    }
  }
  
  colData(sce2)$BNPSpace = as.factor(pred_label)
  
  p8 <- clusterPlot(sce2, label= colData(sce2)$BNPSpace, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("BNPSpace: ARI=", round(max(result_ARI$ARI[, 4]), 3)), fill = "Domains") + 
    scale_fill_manual(values = color_assign(colData(sce2)$BNPSpace, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
  
  library(cowplot)
  p <- plot_grid(p1, p8,  p4, p6, byrow = T, nrow = 1, ncol = 4)
  p
  
}

cluster = plot_cluster(151509)

ggsave(cluster, filename = "reproduce/figures_and_tables/figure5A.png", width = 20, height = 5, bg = "white")


