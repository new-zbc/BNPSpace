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

plot_cluster <- function(sce2, file){
  
  color_pal = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF", "#2CA02CFF", "#9467BDFF",
                "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
                "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
                "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
  
  load(file) 
  library(flexclust)
  #spaGCN_ARI = randIndex(colData(sce2)$spaGCN, colData(sce2)$label)
  
  library(BayesSpace)
  library(ggplot2)
  library(ggsci)
  library(patchwork)
  
  metadata(sce2)$BayesSpace.data <- list()
  metadata(sce2)$BayesSpace.data$platform <- "ST"
  metadata(sce2)$BayesSpace.data$is.enhanced <- FALSE
  
  si <- 8; tsi <- 10
  
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
  
  
  ARI_value = randIndex(colData(sce2)$BNPSpace, colData(sce2)$label)
  
  p8 <- clusterPlot(sce2, label= colData(sce2)$BNPSpace, palette=NULL, size=0.05) +
    #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
    labs(title=paste0("BNPSpace: ARI=", round(ARI_value, 3))) + 
    scale_fill_manual(values = color_assign(colData(sce2)$BNPSpace, colData(sce2)$label, color_pal))+
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.title = element_text(size=tsi), #change legend title font size
          legend.text = element_text(size=si),#change legend text font size
          panel.border = element_rect(colour = "black", fill=NA, size=1) )
  
  library(cowplot)
  p <- plot_grid(p1, p8, p2, p3, p4, p7, p5, byrow = T, nrow = 2, ncol = 4)
  p
}

plot_exp <- function(sce2){
  
  #construc Seurat object
  library(Seurat)
  seu = as.Seurat(sce2, counts = "counts", data = "logcounts", project = "sce_to_seurat")
  
  
  seu = ScaleData(seu)
  Idents(seu) = as.factor(seu$BNPSpace)
  #names(seu@active.ident) = seu@assays$originalexp@counts@Dimnames[[2]]
  
  gamma = seu@assays$originalexp@meta.features$gamma > 0.5
  
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
                group.bar = T, slot = "scale.data", size = 6,
                label = T, draw.lines = T, combine = T) + 
    theme(legend.text = element_text(size = 15),
          legend.title = element_text( size = 18, face='bold'),
          axis.text.y = element_blank()) +
    xlab("Clusters") +
    scale_fill_gradientn(colors = c("blue", "white", "red"))
  
  p
}


load("application/MOBdata/result_sce.RData")

p1 = plot_cluster(sce2, file = "application/MOBdata/result_other_method.RData")

p2 = plot_exp(sce2)

layout = "
AAAA
BBBB
"

p = p1 + p2 +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(0.8, 1))

#pdf("real_data/DLPFC_data/151509/results/summary.pdf")

ggsave(p, filename = "reproduce/figures_and_tables/figure4.png", width = 12, height = 9)
