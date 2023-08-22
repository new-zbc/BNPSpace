dataRader = function(data_folder){
  data1 = read.table(file = paste0(data_folder, "/results/BNPSpace.txt"), header = TRUE, row.names = 1)
  
  data2 = read.table(file = paste0(data_folder, "/results/SCMEB.txt"), header = TRUE, row.names = 1)
  
  data3 = read.table(file = paste0(data_folder, "/results/DRSC.txt"), header = TRUE, row.names = 1)
  
  data4 = read.table(file = paste0(data_folder, "/results/BayesSpace.txt"), header = TRUE, row.names = 1)
  colnames(data4) = "BayesSpace"
  data5 = read.table(file = paste0(data_folder, "/results/kmeans.txt"), header = TRUE, row.names = 1)
  colnames(data5) = "kmean"
  
  data6 = read.table(file = paste0(data_folder, "/results/GMM.txt"), header = TRUE, row.names = 1)
  data7 = read.table(file = paste0(data_folder, "/results/Louvain.txt"), header = TRUE, row.names = 1)
  
  data8 = read.table(file = paste0(data_folder, "/results/BNPSpace_0.txt"), header = TRUE, row.names = 1)
  
  data_ARI = data.frame(BNPSpace = data1[,1], BNPSpace0 = data8[,1], SCMEB = data2[, 1],
                        DRSC = data3[, 1], BayesSpace = data4, Kmeans = data5,
                        GMM = data6[, 1], Louvain = data7[, 1])
  colnames(data_ARI)[2] = "BNPSpace w/o MRF"
  colnames(data_ARI)[4] = "DRSC"
  colnames(data_ARI)[6] = "Kmeans"
  
  library(reshape2)
  data = melt(data_ARI)
  colnames(data) = c("Method", "ARI")
  
  return(data)
}

###################################################################
###### show the ARI of simulation 2, figure 3A
###################################################################
data1 = dataRader("simulation2/scenario1_1")
data1 = cbind(data1, K = rep("K=7", dim(data1)[1]), pi = rep("0.1", dim(data1)[1]))

data2 = dataRader("simulation2/scenario1_2")
data2 = cbind(data2, K = rep("K=7", dim(data2)[1]), pi = rep("0.2", dim(data2)[1]))

data3 = dataRader("simulation2/scenario1_3")
data3 = cbind(data3, K = rep("K=7", dim(data3)[1]), pi = rep("0.3", dim(data3)[1]))

data_ARI = rbind(data1, data2, data3)

data_ARI$Method = factor(data_ARI$Method, levels = c("BNPSpace", "BNPSpace w/o MRF", "SCMEB", 
                                                     "DRSC", "BayesSpace", "Louvain", 
                                                     "Kmeans", "GMM"))
data_ARI$prob = factor(data_ARI$pi, levels = c("0.1", "0.2", "0.3"), 
                       labels = c("0.1" = expression(pi == "0.1"), "0.2" = expression(pi == "0.2"), "0.3" = expression(pi == "0.3")))

library(ggplot2)

box_ARI = ggplot(data = data_ARI, aes(x = Method, y = ARI, fill = Method)) + 
  stat_boxplot(geom ="errorbar", width=0.2, position=position_dodge(0.8)) +
  geom_boxplot() + facet_wrap(~prob, scales = "free", labeller = labeller(prob = label_parsed)) + 
  theme_bw()+
  theme(legend.title = element_text(size = 18, face = "bold"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 15), 
        legend.position = "bottom") 




dataRaderK = function(data_folder){
  data1 = read.table(file = paste0(data_folder, "/results/BNPSpace.txt"), header = TRUE, row.names = 1)
  
  data2 = read.table(file = paste0(data_folder, "/results/SCMEB.txt"), header = TRUE, row.names = 1)
  
  data3 = read.table(file = paste0(data_folder, "/results/DRSC.txt"), header = TRUE, row.names = 1)
  
  data4 = read.table(file = paste0(data_folder, "/results/BayesSpace.txt"), header = TRUE, row.names = 1)
  colnames(data4) = "BayesSpace"
  data5 = read.table(file = paste0(data_folder, "/results/kmeans.txt"), header = TRUE, row.names = 1)
  colnames(data5) = "kmean"
  
  data6 = read.table(file = paste0(data_folder, "/results/GMM.txt"), header = TRUE, row.names = 1)
  data7 = read.table(file = paste0(data_folder, "/results/Louvain.txt"), header = TRUE, row.names = 1)
  
  data8 = read.table(file = paste0(data_folder, "/results/BNPSpace_0.txt"), header = TRUE, row.names = 1)
  
  data_K = data.frame(BNPSpace = data1[,2], BNPSpace0 = data8[,2], SCMEB = data2[, 2],
                      DRSC = data3[, 2], GMM = data6[, 2], Louvain = data7[, 2])
  colnames(data_K)[2] = "BNPSpace w/o MRF"
  colnames(data_K)[4] = "DRSC"
  library(reshape2)
  data = melt(data_K)
  colnames(data) = c("Method", "K")
  return(data)
}


###################################################################
###### show the ARI of simulation 2, figure 3B
###################################################################
data1 = dataRaderK("simulation2/scenario1_1")
data1 = cbind(data1, K_true = rep("K_true = 7", dim(data1)[1]), pi = rep("0.1", dim(data1)[1]))

data2 = dataRaderK("simulation2/scenario1_2")
data2 = cbind(data2, K_true = rep("K_true = 7", dim(data2)[1]), pi = rep("0.2", dim(data2)[1]))

data3 = dataRaderK("simulation2/scenario1_3")
data3 = cbind(data3, K_true = rep("K_true = 7", dim(data3)[1]), pi = rep("0.3", dim(data3)[1]))

data_K = rbind(data1, data2, data3)

data_K$Method = factor(data_K$Method, levels = c("BNPSpace", "BNPSpace w/o MRF", "SCMEB", 
                                                 "DRSC", "Louvain", "GMM"))
data_K$pi = factor(data_K$pi, levels = c("0.1", "0.2", "0.3"), 
                   labels = c("0.1" = expression(pi == "0.1"), "0.2" = expression(pi == "0.2"), "0.3" = expression(pi == "0.3")))


####### first try for stacked barplot
data_K$K = as.character(data_K$K)
library(ggplot2)
data_K$value = rep(1, dim(data_K)[1])
data_K = aggregate(data_K$value, by=list(data_K$Method, data_K$pi, data_K$K), sum)
colnames(data_K) = c("Method", "pi", "K", "value")
data_K$K = factor(data_K$K, levels = c("2", "3", "4", "5", "6","7","8","9", "10"))

data_K$value = data_K$value / 50


barplot_K = ggplot(data = data_K, aes(x = Method, y = value, fill = K)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_grid(~ pi, scales = "free", labeller = labeller(pi = label_parsed)) + 
  theme_bw()+ ylab("Frequency") +
  theme(legend.title = element_text(size = 18, face = "bold"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 15), 
        legend.position = "bottom")

#ggsave(barplot_K, filename = "simulation2/simulation_2_K.png", width = 10, height = 5)

#### plot together
library(patchwork)
p = box_ARI + barplot_K + plot_annotation(tag_levels = "A") + plot_layout(ncol = 1)
ggsave(p, filename = "reproduce/figures_and_tables/figure3.png", width = 10, height = 10)



