library(flexclust)
library(ggplot2)
plot_exp <- function(sampleID){
  
  source("R/utils.R")
  data_folder = "application/DLPFC"
  file_name = paste0(data_folder, "/data/", sampleID, "_counts.RData")
  load(file_name)
  ### quality control
  sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
  ### construct spots network
  Adj = find_neighbors(sce2, "Visium", "lattice")
  neighbors = find_neighbor_index(Adj, "Visium")
  
  #construct Seurat object
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
  #seu@assays$originalexp@data@Dimnames[[1]] = seu@assays$originalexp@meta.features$gene_name
  
  seu = ScaleData(seu)
  
  # scale_data = seu@assays$originalexp@scale.data
  # for(j in 1:dim(scale_data)[1]){
  #   scale_data[j, ] = (scale_data[j, ] - min(scale_data[j, ])) / (max(scale_data[j, ])- min(scale_data[j, ]))
  # }
  # seu@assays$originalexp@scale.data = scale_data
  
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
    guides(colour = guide_colourbar(title = "Cluster")) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text( size = 14, face='bold', angle = 90),
          #legend.position = "none",
          #legend.title = element_blank(),
          axis.text.y = element_blank()) +
    scale_fill_gradientn(colors = c(low = "blue", mid = "white", high="red")) +
    guides(fill = guide_colorbar(title = "Relative expression", title.position = "left"))
  
  p
}


pExp = plot_exp(151509)

ggsave(pExp, filename = "reproduce/figures_and_tables/figure4C.png", width = 11, height = 8, bg = "white")
