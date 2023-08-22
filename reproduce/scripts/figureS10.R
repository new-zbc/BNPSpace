merge_list <- function(sampleID){
  
  color_pal = c("#E377C2FF", "#9467BDFF", "#1F77B4FF", "#D62728FF", "#FF7F0EFF", "#2CA02CFF",
                "#7F7F7FFF", "#BCBD22FF", "#17BECFFF", "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF",
                "#FF9896FF", "#C5B0D5FF", "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF" ,
                "#9EDAE5FF")
  
  source("R/utils.R")
  
  data_folder = "application/DLPFCdata"
  file_name = paste0(data_folder, "/", sampleID, "_counts.RData")
  load(file_name)
  
  ### quality control
  sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
  ### construct spots network
  Adj = find_neighbors(sce2, "Visium", "lattice")
  neighbors = find_neighbor_index(Adj, "Visium")
  
  load(paste0(data_folder, "/", "BNPSpace_results.RData"))
  
  pred_spot = BNPSpace_res$ori$pred_label
  
  

  result_ARI = BNPSpace_res$post
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
  
  library(BayesSpace)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  
  
  metadata(sce2)$BayesSpace.data <- list()
  metadata(sce2)$BayesSpace.data$platform <- "Visium"
  metadata(sce2)$BayesSpace.data$is.enhanced <- FALSE
  
  si <- 8; tsi <- 10
  
  
  plot_list = list()
  for( k in index:(dim(label_mat)[2])){
    
    colData(sce2)$BNPSpace = as.factor(label_mat[, k])
    
    p8 <- clusterPlot(sce2, label= colData(sce2)$BNPSpace, palette=NULL, size=0.05) +
      #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
      labs(title=paste0("K = ", length(unique(label_mat[, k])))) + 
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

pMerge = plot_merge(151509)
ggsave(pMerge, filename = "reproduce/figures_and_tables/figureS10.png", width = 20, height = 15, units = "in")
