
source("R/utils.R")
spatialPlotGene <- function(sampleID, platform = "Visium", features){
  
  data_folder = data_folder = "application/DLPFC"
  file_name = paste0(data_folder, "/data/", sampleID, "_counts.RData")
  load(file_name)
  ### quality control
  sce2 = spatialFilter(sce1, cutoff_sample = 100, cutoff_max = 5)
  ### construct spots network
  Adj = find_neighbors(sce2, "Visium", "lattice")
  neighbors = find_neighbor_index(Adj, "Visium")
  library(patchwork)
  library(BayesSpace)
  
  #construc Seurat object
  library(Seurat)
  seu = as.Seurat(sce2, counts = "counts", data = "logcounts", project = "sce_to_seurat")
  
  seu = ScaleData(seu)
  
  metadata(sce2)$BayesSpace.data <- list()
  metadata(sce2)$BayesSpace.data$platform <- "Visium"
  metadata(sce2)$BayesSpace.data$is.enhanced <- FALSE
  rownames(sce2) = rowData(sce2)$gene_name
  
  # mat = as.matrix(assay(sce2, "logcounts")) 
  # scale_mat = matrix(0, nrow = dim(sce2)[1], ncol = dim(sce2)[2] )
  # for(i in 1:dim(sce2)[1]){
  #   value_mean = mean(as.matrix(mat[i, ]))
  #   value_sd = sd(as.matrix(mat[i, ]))
  #   scale_mat[i, ] = (mat[i, ] - value_mean) / value_sd
  # }
  # 
  scale_mat =  seu@assays$originalexp@scale.data
  
  rownames(scale_mat) = rownames(sce2)
  colnames(scale_mat) = colnames(sce2)
  assay(sce2, "scale.data") = scale_mat
  
  tsi = 16
  
  p1 = featurePlot(sce2, feature = features[1], assay.type = "scale.data", color = NA) + 
    labs(title = features[1]) + 
    theme(plot.title = element_text(size = tsi, face = "bold.italic", hjust = 0.5),
          legend.title = element_text(size = 10, face = "bold"),
          legend.position = "none") + 
    scale_fill_gradientn(colors = c(low = "blue", mid = "white", high="red"))
  
  
  p2 = featurePlot(sce2, feature = features[2], assay.type = "scale.data", color = NA) + 
    labs(title = features[2]) + 
    theme(plot.title = element_text(size = tsi, face = "bold.italic", hjust = 0.5),
          legend.title = element_text(size = 10, face = "bold"),
          legend.position = "none") + 
    scale_fill_gradientn(colors = c(low = "blue", mid = "white", high="red"))
  
  p3 = featurePlot(sce2, feature = features[3], assay.type = "scale.data", color = NA) + 
    labs(title = features[3]) + 
    theme(plot.title = element_text(size = tsi, face = "bold.italic", hjust = 0.5),
          legend.title = element_text(size = 10, face = "bold"),
          legend.position = "none") + 
    scale_fill_gradientn(colors = c(low = "blue", mid = "white", high="red"))
  
  p4 = featurePlot(sce2, feature = features[4], assay.type = "scale.data", color = NA) + 
    labs(title = features[4]) + 
    theme(plot.title = element_text(size = tsi, face = "bold.italic", hjust = 0.5),
          legend.title = element_text(size = 10, face = "bold"),
          legend.position = "none") + 
    scale_fill_gradientn(colors = c(low = "blue", mid = "white", high="red"))
  
  p5 = featurePlot(sce2, feature = features[5], assay.type = "scale.data", color = NA) + 
    labs(title = features[5]) + 
    theme(plot.title = element_text(size = tsi, face = "bold.italic", hjust = 0.5),
          legend.title = element_text(size = 10, face = "bold"),
          legend.position = "none") + 
    scale_fill_gradientn(colors = c(low = "blue", mid = "white", high="red"))
  
  p6 = featurePlot(sce2, feature = features[6], assay.type = "scale.data", color = NA) + 
    labs(title = features[6]) + 
    theme(plot.title = element_text(size = tsi, face = "bold.italic", hjust = 0.5),
          legend.title = element_text(size = 10, face = "bold"),
          legend.position = "none") + 
    scale_fill_gradientn(colors = c(low = "blue", mid = "white", high="red"))
  
  # p7 = featurePlot(sce2, feature = features[7], assay.type = "scale.data", color = NA) + 
  #   labs(title = features[7]) + 
  #   theme(plot.title = element_text(size = 14, face = "bold.italic", hjust = 0.5),
  #         legend.title = element_text(size = 10, face = "bold"),
  #         legend.position = "none") + 
  #   scale_fill_gradientn(colors = c(low = "blue", mid = "white", high="red"))
  # 
  # p8 = featurePlot(sce2, feature = features[8], assay.type = "scale.data", color = NA) + 
  #   labs(title = features[8]) + 
  #   theme(plot.title = element_text(size = 14, face = "bold.italic", hjust = 0.5),
  #         legend.title = element_text(size = 10, face = "bold"),
  #         legend.position = "none") + 
  #   scale_fill_gradientn(colors = c(low = "blue", mid = "white", high="red"))
  
  
  plot_list = list(p1, p2, p3,
                   p4, p5, p6) 
  
  p10 = featurePlot(sce2, feature = features[3], assay.type = "scale.data", color = NA) + 
    labs(title = features[3]) + 
    theme(plot.title = element_text(size = 14, face = "bold"),
          legend.title = element_text(size = 14, face = "bold", angle = 90),
          legend.text = element_text(size = 14),
          legend.position = "right") + 
    scale_fill_gradientn(colors = c(low = "blue", mid = "white", high="red")) +
    guides(fill = guide_colorbar(title = "Relative expression", title.position = "left"))
  
  library(cowplot)
  legend = get_legend(p10)
  
  p = plot_grid(plotlist = plot_list ,ncol = 3)
  p = plot_grid(p, legend, ncol = 2, rel_widths = c(3, 0.5))
  
  p
}

features = c("SCGB2A2", "SCGB1D2", "COX6C",
             "MBP", "GFAP", "NEFL")


pGene = spatialPlotGene(151509, features = features)

ggsave(pGene, filename = "reproduce/figures_and_tables/figure4D.png", width = 10, height = 6, bg = "white")