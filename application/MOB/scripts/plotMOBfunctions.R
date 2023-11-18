#### plot clusters
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


plot_cluster <- function(sce2){
  
  color_pal = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF", "#2CA02CFF", "#9467BDFF",
                "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
                "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
                "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
  
  library(flexclust)
  library(BayesSpace)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  
  metadata(sce2)$BayesSpace.data <- list()
  metadata(sce2)$BayesSpace.data$platform <- "ST"
  metadata(sce2)$BayesSpace.data$is.enhanced <- FALSE
  
  si <- 10; tsi <- 12
  
  p1 <- clusterPlot(sce2, label=colData(sce2)$label, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,1]) +
    labs(title="Ground Truth", fill = "Domains") + scale_fill_manual(values = color_pal[1:length(unique(colData(sce2)$label))]) +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1), 
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  
  
  
  colData(sce2)$SCMEB = as.factor(colData(sce2)$SCMEB)
  
  ARI = randIndex(table(colData(sce2)$label, colData(sce2)$SCMEB))
  
  
  p2 <- clusterPlot(sce2, label= colData(sce2)$SCMEB, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("SCMEB: ARI=", "0.390"), fill = "Domains") + 
    scale_fill_manual(values = color_assign(colData(sce2)$SCMEB, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  
  
  colData(sce2)$DRSC = as.factor(colData(sce2)$DRSC)
  
  ARI = randIndex(table(colData(sce2)$label, colData(sce2)$DRSC))
  
  
  p3 <- clusterPlot(sce2, label= colData(sce2)$DRSC, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("DRSC: ARI=", round(ARI, 3)), fill = "Domains") + 
    scale_fill_manual(values = color_assign(colData(sce2)$DRSC, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  
  colData(sce2)$BayesSpace = as.factor(colData(sce2)$BayesSpace)
  ARI = randIndex(table(colData(sce2)$label, colData(sce2)$BayesSpace))
  
  p4 <- clusterPlot(sce2, label= colData(sce2)$BayesSpace, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("BayesSpace: ARI=", round(ARI, 3)), fill = "Domains") + 
    scale_fill_manual(values = color_assign(colData(sce2)$BayesSpace, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  

  
  colData(sce2)$spaGCN = as.factor(colData(sce2)$spaGCN)
  ARI = randIndex(table(colData(sce2)$label, colData(sce2)$spaGCN))
  
  p6 <- clusterPlot(sce2, label= colData(sce2)$spaGCN, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("spaGCN: ARI=", round(ARI, 3)), fill = "Domains") + 
    scale_fill_manual(values = color_assign(colData(sce2)$spaGCN, colData(sce2)$label, color_pal)) +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  
  
  colData(sce2)$Louvain = as.numeric(colData(sce2)$Louvain) 
  colData(sce2)$Louvain = as.factor(colData(sce2)$Louvain)
  ARI = randIndex(table(colData(sce2)$label, colData(sce2)$Louvain))
  
  p7 <- clusterPlot(sce2, label= colData(sce2)$Louvain, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("Louvain: ARI=", round(ARI, 3)), fill = "Domains") + 
    scale_fill_manual(values = color_assign(colData(sce2)$Louvain, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  
  
  colData(sce2)$BNPSpace = as.factor(colData(sce2)$BNPSpace)
  ARI = 0.665
  
  p8 <- clusterPlot(sce2, label= colData(sce2)$BNPSpace, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("BNPSpace: ARI=", round(ARI, 3)), fill = "Domains") + 
    scale_fill_manual(values = color_assign(colData(sce2)$BNPSpace, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  
  library(cowplot)
  p <- plot_grid(p1, p8, p4, p6,  p2, p3, p7, byrow = T, nrow = 2, ncol = 4)
  p
}


plot_exp <- function(sce2){
  
  #construc Seurat object
  library(Seurat)
  seu = as.Seurat(sce2, counts = "counts", data = "logcounts", project = "sce_to_seurat")
  
  
  seu = ScaleData(seu)
  Idents(seu) = as.factor(seu$BNPSpace)
  #names(seu@active.ident) = seu@assays$originalexp@counts@Dimnames[[2]]
  
  # scale_data = seu@assays$originalexp@scale.data
  # for(j in 1:dim(scale_data)[1]){
  #   scale_data[j, ] = (scale_data[j, ] - min(scale_data[j, ])) / (max(scale_data[j, ])- min(scale_data[j, ]))
  # }
  # seu@assays$originalexp@scale.data = scale_data
  
  
  gamma = seu@assays$originalexp@meta.features$gamma > 0.8
  
  seu@assays$originalexp@scale.data =  seu@assays$originalexp@scale.data[gamma, ]
  seu@assays$originalexp@counts =  seu@assays$originalexp@counts[gamma, ]
  seu@assays$originalexp@data =  seu@assays$originalexp@data[gamma, ]
  seu@assays$originalexp@meta.features =  seu@assays$originalexp@meta.features[gamma, ]
  
  # Find markers differential each group
  
  seu.markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  #write.table(seu.markers, file = "real_data/MOB_data/results/DEG_analysis.txt", row.names = F)
  
  
  library(dplyr)
  top_gene <- seu.markers %>%
    group_by(cluster) %>% top_n(n=50, avg_log2FC)
  
  p = DoHeatmap(seu, features = top_gene$gene, group.by = "BNPSpace",
                group.bar = T, slot = "scale.data", size = 6, group.colors = c("#1F77B4FF", "#D62728FF", "#FF7F0EFF",  "#9467BDFF", "#2CA02CFF"),
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

