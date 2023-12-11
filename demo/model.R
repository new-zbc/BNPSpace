library(Rcpp)
library(SingleCellExperiment)
library(RcppArmadillo)
library(flexclust)


source("R/utils.R")
source("R/main.R")

### construct spots network
Adj = find_neighbors(sce, "ST", "image")

neighbors = find_neighbor_index(Adj, "ST")
result = NSCFS(sce, neighbors, n_clusters = NULL, d = 1, n_iters = 300, seed = 1)
library(flexclust)
pred_label = result$pred_label
ARI_value = randIndex(table(pred_label, colData(sce)$label))


library(BayesSpace)
metadata(sce)$BayesSpace.data <- list()
metadata(sce)$BayesSpace.data$platform <- "ST"
metadata(sce)$BayesSpace.data$is.enhanced <- FALSE
p <- clusterPlot(sce, label= pred_label, palette=NULL, size=0.05) +
  #scale_fill_viridis_d(option = "A", labels = 1:clustNumMat[j,2]) +
  labs(title=paste0("BNPSpace: ARI=", round(ARI_value, 3)), fill = "Domains") + 
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10),#change legend text font size
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

data_gamma = data.frame(gene = sample(500), gamma_est = result$gamma, gamma_true = c(rep(TRUE, 20), rep(FALSE, 480)))

#data_gamma = data_gamma[sample(500), ]

library(ggplot2)
PPI_plot = ggplot(data = data_gamma) + geom_point(aes(x = gene, y = gamma_est, color = gamma_true))+
  ylab("Posterior probability of inclusion (PPI)") + xlab("Gene index") + 
  labs(title = "PPIs by BNPSpace", color = "DGs") +
  geom_hline(aes(yintercept = 0.5), color = "blue") +
  theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_text(size=12), #change legend title font size
        legend.text = element_text(size=10),#change legend text font size
        #panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

library(patchwork)
p + PPI_plot + plot_layout(nrow = 1)

