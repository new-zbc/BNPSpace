library(ggsci)
mypal = pal_d3(palette = "category20")(20)
mypal
library(scales)
show_col(mypal)
new_pal = mypal
new_pal[3] = mypal[4]
new_pal[4] = mypal[3]
show_col(new_pal)
color_pal = new_pal

color_assign_1 <- function(pred, ground_truth, color_pal){
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


color_assign_2 <- function(pred, ground_truth, color_pal){
  R = length(unique(pred))
  C = length(unique(ground_truth))
  p_table = table(pred, ground_truth)
  p_table = p_table / matrix(rep(colSums(p_table), each = R), byrow = F,ncol = C)
  
  color_pred = rep(NA, R)
  #layer_used = rep(FALSE, C)
  count = 0
  for(i in 1:C){
    index = which.max(p_table[, i])
    if(is.na(color_pred[index]) & p_table[i, index] > 0.5){
      #layer_used[index] = TRUE
      color_pred[index] = color_pal[i]
    }
  }
  
  color_pred[which(is.na(color_pred))] = color_pal[C + (1:sum(is.na(color_pred)))]
  return(color_pred)
}

color_assign = color_assign_1


#sampleID = 151509

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
  #Adj = find_neighbors(sce2, "Visium", "lattice")
  #neighbors = find_neighbor_index(Adj, "Visium")
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
  
  si <- 8; tsi <- 10
  
  cluster_truth = unique(colData(sce2)$label)
  
  # if(length(cluster_truth) == 7){
  #   colData(sce2)$label = factor(colData(sce2)$label, levels = c("Layer1", "Layer2", "Layer3", "Layer4", "Layer5", "Layer6", "WM"))
  # }
  # else{
  #   colData(sce2)$label = factor(colData(sce2)$label, levels = c("Layer3", "Layer4", "Layer5", "Layer6", "WM"))
  # }
  
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
    labs(title="Ground Truth", fill = "Domains") + scale_fill_manual(values = color_pal[1:length(unique(colData(sce2)$label))]) +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  
  SC_MEB = list(res_summary$SC.MEB1, res_summary$SC.MEB2)
  index = which.max(c(res_summary$SC.MEB1$ARI, res_summary$SC.MEB2$ARI))
  
  colData(sce2)$SCMEB = as.factor(SC_MEB[[index]]$label)
  
  
  p2 <- clusterPlot(sce2, label= colData(sce2)$SCMEB, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("SCMEB: ARI=", round(SC_MEB[[index]]$ARI, 3))) + 
    scale_fill_manual(values = color_assign(colData(sce2)$SCMEB, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  
  DR_SC = list(res_summary$DR.SC1, res_summary$DR.SC2)
  index = which.max(c(res_summary$DR.SC1$ARI, res_summary$DR.SC2$ARI))
  
  colData(sce2)$DRSC = as.factor(DR_SC[[index]]$label)
  
  p3 <- clusterPlot(sce2, label= colData(sce2)$DRSC, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("DRSC: ARI=", round(DR_SC[[index]]$ARI, 3))) + 
    scale_fill_manual(values = color_assign(colData(sce2)$DRSC, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  colData(sce2)$BayesSpace = as.factor(res_summary$BayesSpace$label)
  
  p4 <- clusterPlot(sce2, label= colData(sce2)$BayesSpace, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("BayesSpace: ARI=", round(res_summary$BayesSpace$ARI, 3))) + 
    scale_fill_manual(values = color_assign(colData(sce2)$BayesSpace, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  
  colData(sce2)$kmeans = as.factor(res_summary$kmeans$label)
  p5 <- clusterPlot(sce2, label= colData(sce2)$kmeans, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("Kmeans: ARI=", round(res_summary$kmeans$ARI, 3))) + 
    scale_fill_manual(values = color_assign(colData(sce2)$kmeans, colData(sce2)$label, color_pal)) +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
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
    labs(title=paste0("spaGCN: ARI=", round(spaGCN_ARI, 3))) + 
    scale_fill_manual(values = color_assign(spaGCN_7, colData(sce2)$label, color_pal)) +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  res_summary$Louvain$label = as.numeric(res_summary$Louvain$label) 
  colData(sce2)$Louvain = as.factor(res_summary$Louvain$label)
  
  p7 <- clusterPlot(sce2, label= colData(sce2)$Louvain, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("Louvain: ARI=", round(res_summary$Louvain$ARI, 3))) + 
    scale_fill_manual(values = color_assign(colData(sce2)$Louvain, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  result_ARI = DLPFC_ARI(sampleID)
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
    labs(title=paste0("BNPSpace: ARI=", round(max(result_ARI$ARI[, 4]), 3))) + 
    scale_fill_manual(values = color_assign(colData(sce2)$BNPSpace, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  library(cowplot)
  p <- plot_grid(p1, p8,  p4, p6, byrow = T, nrow = 1, ncol = 4, rel_heights = c(0.5, 1,1,1))
  p
  
}


merge_list <- function(sampleID){
  
  color_pal = c("#E377C2FF", "#9467BDFF", "#1F77B4FF", "#D62728FF", "#FF7F0EFF", "#2CA02CFF",
                "#7F7F7FFF", "#BCBD22FF", "#17BECFFF", "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF",
                "#FF9896FF", "#C5B0D5FF", "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF" ,
                "#9EDAE5FF")
  
  data_folder = paste0("real_data/DLPFC_data/", sampleID)
  file_name = paste0(data_folder, "/", sampleID, "_counts.RData")
  load(file_name)
  
  ### quality control
  sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
  ### construct spots network
  Adj = find_neighbors(sce2, "Visium", "lattice")
  neighbors = find_neighbor_index(Adj, "Visium")
  load(paste0(data_folder, "/", "reults_second_try.RData"))
  pred_spot = res_summary$NSCFS$pred_label
  
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
  
  si <- 8; tsi <- 10
  
  result_ARI = DLPFC_ARI(sampleID)
  
  #####################################################
  result_label = result_ARI$pred
  clusters = as.numeric(names(table(result_label)))
  pred_label = rep(NA, length(result_label))
  for(i in 1:length(result_label)){
    for(k in 1:length(clusters)){
      if(result_label[i] == clusters[k]){pred_label[i] = k}
    }
  }
  
  label_rename = pred_spot
  sum_table = table(pred_label, pred_spot)
  for(k in 1:length(unique(pred_label))){
    index = which(sum_table[k, ] != 0)
    if(length(index) > 1){
      for(i in 1:length(index)){
        label_rename[pred_spot == index[i]] = paste0(k, letters[i])
      }
    }else{label_rename[pred_spot == index[1]] = as.character(k)}
  }
  
  colData(sce2)$BNPSpace = label_rename
  
  ##############################################
  
  label_mat = cbind(pred_spot, result_ARI$all)
  index = dim(result_ARI$ARI)[1] - which.max(result_ARI$ARI[, 4]) + 1
  
  #plot_list = list()
  for(k in 1:dim(label_mat)[2]){
    
    result_label = label_mat[, k]
    
    clusters = as.numeric(names(table(result_label)))
    pred_label = rep(NA, length(result_label))
    for(i in 1:length(result_label)){
      for(j in 1:length(clusters)){
        if(result_label[i] == clusters[j]){pred_label[i] = j}
      }
    }
    
    label_mat[, k] = pred_label
  }
  
  label_mat = label_mat[, dim(label_mat)[2]:1]
  
  for( k in index:(dim(label_mat)[2]-1)){
    p_table = table(label_mat[, k], label_mat[, k+1])
    for(i in 1:dim(p_table)[1]){
      if(p_table[i, i] != max(p_table[i, ])){
        
        ind1 = which(label_mat[, k+1] == i)
        
        label2 = which.max(p_table[i, ])
        
        ind2 = which(label_mat[, k+1] == label2)
        
        label_mat[, k+1][ind2] = i
        label_mat[, k+1][ind1] = label2
        
      }
      p_table = table(label_mat[, k], label_mat[, k+1])
      
    }
  }
  
  
  label_mat[, 13] = label_rename
  
  plot_list = list()
  
  
  for( k in (index):(dim(label_mat)[2])){
    
    colData(sce2)$BNPSpace = as.factor(label_mat[, k])
    
    p8 <- clusterPlot(sce2, label= colData(sce2)$BNPSpace, palette=NULL, size=0.05) +
      #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
      labs(title=paste0("K = ", length(unique(label_mat[, k]))), fill = "Domains") + 
      scale_fill_manual(values = color_pal[1:length(unique(label_mat[, k]))])+
      theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
            legend.key.height = unit(0.5, 'cm'), #change legend key height
            legend.key.width = unit(0.5, 'cm'), #change legend key width
            legend.title = element_text(size=tsi), #change legend title font size
            legend.text = element_text(size=si),#change legend text font size
            panel.border = element_rect(colour = "black", fill=NA, size=1) )
    
    plot_list = append(plot_list, list(p8))
  }
  return(plot_list)
}


plot_merge <- function(sampleID){
  
  plot_list = merge_list(sampleID)
  
  p = plot_list[[1]]
  for(i in 2:length(plot_list)){
    p = p + plot_list[[i]]
  }
  p = p + plot_layout(ncol = 4)
  
  p
}



library(SingleCellExperiment)
library(SC.MEB)
library(Matrix)
library(flexclust)
library(scater)
library(scran)
source("NSCFS/utils.R")


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
  #sampleID = 151507
  
  data_folder = paste0("real_data/DLPFC_data/", sampleID)
  
  file_name = paste0(data_folder, "/", sampleID, "_counts.RData")
  
  load(file_name)
  
  
  ### quality control
  sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
  ### construct spots network
  Adj = find_neighbors(sce2, "Visium", "lattice")
  neighbors = find_neighbor_index(Adj, "Visium")
  
  load(paste0(data_folder, "/", "reults_second_try.RData"))
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


plot_ARI = function(sampleID){
  library(ggplot2)
  library(ggrepel)
  data_folder = paste0("real_data/DLPFC_data/", sampleID)
  file_name = paste0(data_folder, "/", "diff_K_ARI.RData")
  load(file_name)
  ARI_K_results = ARI_K_results[, -c(5, 6)]
  colnames(ARI_K_results) = c("K", "SCMEB", "DRSC", "BayesSpace", "Louvain", "spaGCN")
  
  library(reshape)
  ARI_K_results_melt = melt(ARI_K_results, id.vars = "K", variable_name = "Method")
  colnames(ARI_K_results_melt) = c("K", "Method", "ARI")
  
  result = DLPFC_ARI(sampleID)
  result_ARI = as.data.frame(result$ARI[, c(1, 4)])
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
  
  ARI_K_results_melt$Method = factor(ARI_K_results_melt$Method, levels = c("BNPSpace", "BayesSpace", "spaGCN",
                                                                           "DRSC", "SCMEB", "Louvain"))
                                                                        
  
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


plot_PPI = function(sampleID){
  data_folder = paste0("real_data/DLPFC_data/", sampleID)
  file_name = paste0(data_folder, "/", sampleID, "_counts.RData")
  load(file_name)
  ### quality control
  sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
  
  load(paste0(data_folder, "/", "reults_second_try.RData"))
  result1 = res_summary$NSCFS
  index = which(rowData(sce2)$is.HVG)
  df = data.frame(gene = rep(1:length(result1$gamma[index])), PPI = result1$gamma[index])
  library(ggplot2)
  p <- ggplot() + geom_segment(aes(x = gene, y = 0, xend = gene, yend = PPI), data = df)
  p
  
  sum(result1$gamma > 0.5)
  sum(rowData(sce2)$is.HVG)
}

plot_exp <- function(sampleID){
  
  data_folder = paste0("real_data/DLPFC_data/", sampleID)
  file_name = paste0(data_folder, "/", sampleID, "_counts.RData")
  load(file_name)
  ### quality control
  sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
  ### construct spots network
  Adj = find_neighbors(sce2, "Visium", "lattice")
  neighbors = find_neighbor_index(Adj, "Visium")
  
  #construc Seurat object
  library(Seurat)
  seu = as.Seurat(sce2, counts = "counts", data = "logcounts", project = "sce_to_seurat")
  
  
  load(paste0(data_folder, "/", "reults_second_try.RData"))
  result1 = res_summary$NSCFS
  
  ground_truth = colData(sce2)$label
  pred_spot = result1$pred_label
  ARI_value = randIndex(table(pred_spot, ground_truth))
  temp = c(length(unique(pred_spot)), 0, 0, ARI_value)
  
  gamma = res_summary$NSCFS$gamma > 0.5
  
  result_ARI = DLPFC_ARI(sampleID)
  result_label = result_ARI$pred
  clusters = as.numeric(names(table(result_label)))
  pred_label = rep(NA, length(result_label))
  for(i in 1:length(result_label)){
    for(k in 1:length(clusters)){
      if(result_label[i] == clusters[k]){pred_label[i] = k}
    }
  }
  
  label_rename = pred_spot
  sum_table = table(pred_label, pred_spot)
  for(k in 1:length(unique(pred_label))){
    index = which(sum_table[k, ] != 0)
    if(length(index) > 1){
      for(i in 1:length(index)){
        label_rename[pred_spot == index[i]] = paste0(k, letters[i])
      }
    }else{label_rename[pred_spot == index[1]] = as.character(k)}
  }
  
  seu@meta.data$pred_label = pred_label
  seu@assays$originalexp@data@Dimnames[[1]] = seu@assays$originalexp@meta.features$gene_name
  
  seu = ScaleData(seu)
  Idents(seu) = as.factor(label_rename)
  names(seu@active.ident) = seu@assays$originalexp@counts@Dimnames[[2]]
  
  seu@assays$originalexp@scale.data =  seu@assays$originalexp@scale.data[gamma, ]
  seu@assays$originalexp@counts =  seu@assays$originalexp@counts[gamma, ]
  seu@assays$originalexp@data =  seu@assays$originalexp@data[gamma, ]
  seu@assays$originalexp@meta.features =  seu@assays$originalexp@meta.features[gamma, ]
  
  # Find markers differential each group
  
  seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  library(dplyr)
  top_gene <- seu.markers %>%
    group_by(cluster) %>% top_n(n=50, avg_log2FC)
  
  p = DoHeatmap(seu, features = top_gene$gene,
                group.bar = T, slot = "scale.data", size = 6,
                label = T, draw.lines = T, combine = T) + 
    theme(legend.text = element_text(size = 15),
          legend.title = element_text( size = 18, face='bold'),
          axis.text.y = element_blank()) +
    guides(colour = guide_colorbar("title")) +
    scale_fill_gradientn(colors = c("blue", "white", "red"))
  
  p
}


spatialPlotGene <- function(sampleID, platform = "Visium", features){
  
  data_folder = paste0("real_data/DLPFC_data/", sampleID)
  file_name = paste0(data_folder, "/", sampleID, "_counts.RData")
  load(file_name)
  ### quality control
  sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
  ### construct spots network
  Adj = find_neighbors(sce2, "Visium", "lattice")
  neighbors = find_neighbor_index(Adj, "Visium")
  library(patchwork)
  library(BayesSpace)
  metadata(sce1)$BayesSpace.data <- list()
  metadata(sce1)$BayesSpace.data$platform <- "Visium"
  metadata(sce1)$BayesSpace.data$is.enhanced <- FALSE
  
  metadata(sce2)$BayesSpace.data <- list()
  metadata(sce2)$BayesSpace.data$platform <- "Visium"
  metadata(sce2)$BayesSpace.data$is.enhanced <- FALSE
  rownames(sce2) = rowData(sce2)$gene_name
  
  mat = as.matrix(assay(sce2, "logcounts")) 
  scale_mat = matrix(0, nrow = dim(sce2)[1], ncol = dim(sce2)[2] )
  for(i in 1:dim(sce2)[1]){
    mu = mean(as.matrix(mat[i, ]))
    sd_i = sd(as.matrix(mat[i, ]))
    scale_mat[i, ] = (mat[i, ] - mu) / sd_i
  }
  
  rownames(scale_mat) = rownames(sce2)
  colnames(scale_mat) = colnames(sce2)
  assay(sce2, "scale.data") = scale_mat
  
  
  p1 = featurePlot(sce2, feature = features[1], assay.type = "scale.data", color = NA) + 
    labs(title = features[1]) + 
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 20, face = "bold"),
          legend.position = "none") + 
    scale_fill_stepsn(n.breaks = 20, colors = viridis::magma(20))
  
  
  p2 = featurePlot(sce2, feature = features[2], assay.type = "scale.data", color = NA) + 
    labs(title = features[2]) + 
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 20, face = "bold"),
          legend.position = "none") +
    scale_fill_stepsn(n.breaks = 20, colors = viridis::magma(20))
  
  p3 = featurePlot(sce2, feature = features[3], assay.type = "scale.data", color = NA) + 
    labs(title = features[3]) + 
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 20, face = "bold"), 
          legend.position = "none") + 
    scale_fill_stepsn(n.breaks = 20, colors = viridis::magma(20))
  
  p4 = featurePlot(sce2, feature = features[4], assay.type = "scale.data", color = NA) + 
    labs(title = features[4]) + 
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 20, face = "bold"),
          legend.position = "none") + 
    scale_fill_stepsn(n.breaks = 20, colors = viridis::magma(20))
  
  p5 = featurePlot(sce2, feature = features[5], assay.type = "scale.data", color = NA) + 
    labs(title = features[5]) + 
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 20, face = "bold"),
          legend.position = "none") + 
    scale_fill_stepsn(n.breaks = 20, colors = viridis::magma(20))
  
  p6 = featurePlot(sce2, feature = features[6], assay.type = "scale.data", color = NA) + 
    labs(title = features[6]) + 
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 20, face = "bold"),
          legend.position = "none") + 
    scale_fill_stepsn(n.breaks = 20, colors = viridis::magma(20))
  
  p7 = featurePlot(sce2, feature = features[7], assay.type = "scale.data", color = NA) + 
    labs(title = features[7]) + 
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 20, face = "bold"),
          legend.position = "none") + 
    scale_fill_stepsn(n.breaks = 20, colors = viridis::magma(20))
  
  p8 = featurePlot(sce2, feature = features[8], assay.type = "scale.data", color = NA) + 
    labs(title = features[8]) + 
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.title = element_text(size = 20, face = "bold"),
          legend.position = "none") + 
    scale_fill_stepsn(n.breaks = 20, colors = viridis::magma(20))
  
  
  plot_list = list(p1, p2, p3,
                   p4, p5, p6, 
                   p7, p8)
  
  p10 = featurePlot(sce2, feature = features[1], assay.type = "scale.data", color = NA) + 
    labs(title = features[1]) + 
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_blank(),
          legend.position = "right") + 
    scale_fill_stepsn(n.breaks = 100, colors = viridis::magma(20))
  
  library(cowplot)
  legend = get_legend(p10)
  
  p = plot_grid(plotlist = plot_list ,ncol = 4)
  p = plot_grid(p, legend, ncol = 2, rel_widths = c(4, 0.5))
  
  
  
  
  #### box plot
  load(paste0(data_folder, "/", "reults_second_try.RData"))
  result1 = res_summary$NSCFS
  
  ground_truth = colData(sce2)$label
  pred_spot = result1$pred_label
  ARI_value = randIndex(table(pred_spot, ground_truth))
  temp = c(length(unique(pred_spot)), 0, 0, ARI_value)
  
  gamma = res_summary$NSCFS$gamma > 0.5
  
  result_ARI = DLPFC_ARI(sampleID)
  result_label = result_ARI$pred
  clusters = as.numeric(names(table(result_label)))
  pred_label = rep(NA, length(result_label))
  for(i in 1:length(result_label)){
    for(k in 1:length(clusters)){
      if(result_label[i] == clusters[k]){pred_label[i] = k}
    }
  }
  label_rename = pred_spot
  sum_table = table(pred_label, pred_spot)
  for(k in 1:length(unique(pred_label))){
    index = which(sum_table[k, ] != 0)
    for(i in 1:length(index)){
      label_rename[label_rename == index[i]] = paste0(k, "-", i)
    }
  }
  
  df1 = data.frame(gene_name = rep(features[1], length(pred_label)),
                   Exp = assay(sce2, "scale.data")[which(rowData(sce2)$gene_name == features[1]), ], 
                   group = label_rename, 
                   Cluster = pred_label)
  df2 = data.frame(gene_name = rep(features[2], length(pred_label)),
                   Exp = assay(sce2, "scale.data")[which(rowData(sce2)$gene_name == features[2]), ], 
                   group = label_rename, 
                   Cluster = pred_label)
  
  df3 = data.frame(gene_name = rep(features[3], length(pred_label)),
                   Exp = assay(sce2, "scale.data")[which(rowData(sce2)$gene_name == features[3]), ], 
                   group = label_rename, 
                   Cluster = pred_label)
  df4 = data.frame(gene_name = rep(features[4], length(pred_label)),
                   Exp = assay(sce2, "scale.data")[which(rowData(sce2)$gene_name == features[4]), ], 
                   group = label_rename, 
                   Cluster = pred_label)
  df = rbind(df1, df2, df3, df4)
  df$gene_name = factor(df$gene_name, levels = features)
  df$Cluster = as.factor(df$Cluster)
  
  library(ggplot2)
  box_gene = ggplot(data = df, aes(x = factor(group), y = count, fill = Cluster)) + 
    geom_boxplot() + facet_wrap(~gene_name, scales = "free") + theme_bw()+
    theme(legend.title = element_text(size = 18, face = "bold"), 
          axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 90, size = 12, face = "bold"), 
          axis.title.y = element_text(size = 18, face = "bold"),
          strip.text = element_text(face = "bold", size = 18),
          legend.text = element_text(size = 15), 
          legend.position = "right")
  box_gene
  
  return(list(plot1 = p,
              plot2 = box_gene))
}


