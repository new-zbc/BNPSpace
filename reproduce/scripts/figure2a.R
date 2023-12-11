library(BayesSpace)
library(ggplot2)
library(ggsci)
library(patchwork)

color_pal = c("#1F77B4FF", "#FF7F0EFF", "#D62728FF", "#2CA02CFF", "#9467BDFF",
              "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF", "#17BECFFF",
              "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", "#C5B0D5FF",
              "#C49C94FF", "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")

##################################
## The first pattern 
#################################

load("simulation1/scenario1_1/data/1.RData")

metadata(sce)$BayesSpace.data <- list()
metadata(sce)$BayesSpace.data$platform <- "ST"
metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

p1 <- clusterPlot(sce, label=colData(sce)$label, palette=NULL, size=0.05) +
  #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,1]) +
  labs(title="Spatial pattern I", fill = "Domains") + 
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(face = "bold", size=14, hjust = 0.5), #change legend title font size
        legend.text = element_text(face = "bold", size=12),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) + scale_fill_manual(values = color_pal)
        #plot.title = element_text(face = "bold", size = 16, hjust = 0.5))



load("simulation1/scenario2_1/data/1.RData")
metadata(sce)$BayesSpace.data <- list()
metadata(sce)$BayesSpace.data$platform <- "ST"
metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

p2 <- clusterPlot(sce, label=colData(sce)$label, palette=NULL, size=0.05) +
  #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,1]) +
  labs(title="Spatial pattern II", fill = "Domains") + 
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(face = "bold", size=14, hjust = 0.5), #change legend title font size
        legend.text = element_text(face = "bold", size=12),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) + scale_fill_manual(values = color_pal)
        #plot.title = element_text(face = "bold", size = 16, hjust = 0.5))


load("simulation1/scenario3_1/data/1.RData")
metadata(sce)$BayesSpace.data <- list()
metadata(sce)$BayesSpace.data$platform <- "ST"
metadata(sce)$BayesSpace.data$is.enhanced <- FALSE

p3 <- clusterPlot(sce, label=colData(sce)$label, palette=NULL, size=0.05) +
  #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,1]) +
  labs(title="Spatial pattern III", fill = "Domains") + 
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(face = "bold", size=14, hjust = 0.5), #change legend title font size
        legend.text = element_text(face = "bold", size=12),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) + scale_fill_manual(values = color_pal)
        #plot.title = element_text(face = "bold", size = 16, hjust = 0.5))

#show_col(color_pal)
data_folder = "application/DLPFC"
file_name = paste0(data_folder, "/data/", 151509, "_counts.RData")
load(file_name)
metadata(sce1)$BayesSpace.data <- list()
metadata(sce1)$BayesSpace.data$platform <- "Visium"
metadata(sce1)$BayesSpace.data$is.enhanced <- FALSE

colData(sce1)$label = as.character(colData(sce1)$label)
colData(sce1)$label[colData(sce1)$label == "Layer1"] = "L1"
colData(sce1)$label[colData(sce1)$label == "Layer2"] = "L2"
colData(sce1)$label[colData(sce1)$label == "Layer3"] = "L3"
colData(sce1)$label[colData(sce1)$label == "Layer4"] = "L4"
colData(sce1)$label[colData(sce1)$label == "Layer5"] = "L5"
colData(sce1)$label[colData(sce1)$label == "Layer6"] = "L6"
colData(sce1)$label = factor(colData(sce1)$label, levels = c("L1", "L2", "L3", "L4", "L5", "L6", "WM"))


p4 <- clusterPlot(sce1, label=colData(sce1)$label, palette=NULL, size=0.05) +
  #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,1]) +
  labs(title="Spatial pattern IV", fill = "Domains") + 
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(face = "bold", size=14, hjust = 0.5), #change legend title font size
        legend.text = element_text(face = "bold", size=12),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) + scale_fill_manual(values = color_pal)
#plot.title = element_text(face = "bold", size = 16, hjust = 0.5))


plot_pattern = cowplot::plot_grid(p1, p2, p3, p4, nrow = 1, rel_widths = c(1,1,1,0.95))
ggsave(plot_pattern, filename = "reproduce/figures_and_tables/figure2a.png", width = 10, height = 4)

