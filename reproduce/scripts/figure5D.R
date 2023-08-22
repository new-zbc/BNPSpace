
source("R/utils.R")
spatialPlotGene <- function(sampleID, platform = "Visium", features){
  
  data_folder = data_folder = "application/DLPFCdata"
  file_name = paste0(data_folder, "/", sampleID, "_counts.RData")
  load(file_name)
  ### quality control
  sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
  ### construct spots network
  Adj = find_neighbors(sce2, "Visium", "lattice")
  neighbors = find_neighbor_index(Adj, "Visium")
  library(patchwork)
  library(BayesSpace)
  
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
  
  p
}

features = c("SCGB2A2", "SCGB1D2", "GFAP",
             "MBP", "COX6C", "FABP7", "MALAT1", "NEFL")


pGene = spatialPlotGene(151509, features = features)

ggsave(pGene, filename = "reproduce/figures_and_tables/figure5D.png", width = 14, height = 6, bg = "white")