
source("R/utils.R")
library(ggsci)
mypal = pal_d3(palette = "category20")(20)
library(scales)
new_pal = mypal
new_pal[3] = mypal[4]
new_pal[4] = mypal[3]
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


plot_cluster <- function(sampleID){
  
  color_pal = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF", "#2CA02CFF", "#9467BDFF",
                "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
                "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
                "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
  #show_col(color_pal)
  data_folder = "application/DLPFCdata"
  file_name = paste0(data_folder, "/", sampleID, "_counts.RData")
  load(file_name)
  
  ### quality control
  sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
  ### construct spots network
  Adj = find_neighbors(sce2, "Visium", "lattice")
  neighbors = find_neighbor_index(Adj, "Visium")
  
  load(paste0(data_folder, "/", "result_other_method.RData"))
  
  library(flexclust)
  spaGCN_ARI = randIndex(colData(sce2)$spaGCN, colData(sce2)$label)
  
  
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
  
  if(length(cluster_truth) == 7){
    colData(sce2)$label = factor(colData(sce2)$label, levels = c("Layer1", "Layer2", "Layer3", "Layer4", "Layer5", "Layer6", "WM"))
  }
  else{
    colData(sce2)$label = factor(colData(sce2)$label, levels = c("Layer3", "Layer4", "Layer5", "Layer6", "WM"))
  }
  
  
  p1 <- clusterPlot(sce2, label=colData(sce2)$label, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,1]) +
    labs(title="Ground Truth") + scale_fill_manual(values = color_pal[1:length(unique(colData(sce2)$label))]) +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1) )
  
  
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
          panel.border = element_rect(colour = "black", fill=NA, size=1) )
  
  
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
          panel.border = element_rect(colour = "black", fill=NA, size=1) )
  
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
          panel.border = element_rect(colour = "black", fill=NA, size=1) )
  
  
  colData(sce2)$kmeans = as.factor(res_summary$kmeans$label)
  p5 <- clusterPlot(sce2, label= colData(sce2)$kmeans, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("kmeans: ARI=", round(res_summary$kmeans$ARI, 3))) + 
    scale_fill_manual(values = color_assign(colData(sce2)$kmeans, colData(sce2)$label, color_pal)) +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1) )
  
  
  colData(sce2)$spaGCN = as.factor(as.numeric(colData(sce2)$spaGCN)+1)
  p6 <- clusterPlot(sce2, label= colData(sce2)$spaGCN, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("spaGCN: ARI=", round(spaGCN_ARI, 3))) + 
    scale_fill_manual(values = color_assign(colData(sce2)$spaGCN, colData(sce2)$label, color_pal)) +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1) )
  
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
          panel.border = element_rect(colour = "black", fill=NA, size=1) )
  
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
    labs(title=paste0("BNPSpace: ARI=", round(max(result_ARI$ARI[, 4]), 3))) + 
    scale_fill_manual(values = color_assign(colData(sce2)$BNPSpace, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1) )
  
  library(cowplot)
  p <- plot_grid(p1, p8, p2, p3, p4, p6, p7, p5, byrow = T, nrow = 2, ncol = 4)
  p
  
}



plot_ARI = function(sampleID){
  library(ggplot2)
  library(ggrepel)
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


plot_exp <- function(sampleID){
  
  data_folder = data_folder = "application/DLPFCdata"
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
  
  
  load(paste0(data_folder, "/", "BNPSpace_results.RData"))
  result1 = BNPSpace_res$ori
  
  ground_truth = colData(sce2)$label
  pred_spot = result1$pred_label
  ARI_value = randIndex(table(pred_spot, ground_truth))
  
  
  gamma = result1$gamma > 0.5
  
  result_ARI = BNPSpace_res$post
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
    scale_fill_gradientn(colors = c("blue", "white", "red"))
  
  p
}


cluster = plot_cluster(151509)
pARI = plot_ARI(151509)
pExp = plot_exp(151509)



