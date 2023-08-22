
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
###### show the ARI of simulation 1, generate figure 1 
###################################################################
data1 = dataRader("simulation1/scenario1_1")
data1 = cbind(data1, K = rep("K=3", dim(data1)[1]), pi = rep("0.1", dim(data1)[1]))

data2 = dataRader("simulation1/scenario1_2")
data2 = cbind(data2, K = rep("K=3", dim(data2)[1]), pi = rep("0.2", dim(data2)[1]))

data3 = dataRader("simulation1/scenario1_3")
data3 = cbind(data3, K = rep("K=3", dim(data3)[1]), pi = rep("0.3", dim(data3)[1]))


data4 = dataRader("simulation1/scenario2_1")
data4 = cbind(data4, K = rep("K=5", dim(data4)[1]), pi = rep("0.1", dim(data4)[1]))

data5 = dataRader("simulation1/scenario2_2")
data5 = cbind(data5, K = rep("K=5", dim(data5)[1]), pi = rep("0.2", dim(data5)[1]))

data6 = dataRader("simulation1/scenario2_3")
data6 = cbind(data6, K = rep("K=5", dim(data6)[1]), pi = rep("0.3", dim(data6)[1]))


data7 = dataRader("simulation1/scenario3_1")
data7 = cbind(data7, K = rep("K=7", dim(data7)[1]), pi = rep("0.1", dim(data7)[1]))

data8 = dataRader("simulation1/scenario3_2")
data8 = cbind(data8, K = rep("K=7", dim(data8)[1]), pi = rep("0.2", dim(data8)[1]))

data9 = dataRader("simulation1/scenario3_3")
data9 = cbind(data9, K = rep("K=7", dim(data9)[1]), pi = rep("0.3", dim(data9)[1]))

data_ARI = rbind(data1, data2, data3, data4, data5, data6, data7, data8, data9)

data_ARI$Method = factor(data_ARI$Method, levels = c("BNPSpace", "BNPSpace w/o MRF","BayesSpace", "SCMEB", 
                                                     "Louvain", "Kmeans", "GMM","DRSC"))
data_ARI$pi = factor(data_ARI$pi, levels = c("0.1", "0.2", "0.3"), 
                     labels = c("0.1" = expression(pi == "0.1"), "0.2" = expression(pi == "0.2"), "0.3" = expression(pi == "0.3")))

library(ggplot2)

box_ARI = ggplot(data = data_ARI, aes(x = Method, y = ARI, fill = Method)) + 
  stat_boxplot(geom ="errorbar", width=0.15,position=position_dodge(0.8)) +
  geom_boxplot() + facet_grid(pi ~K, scales = "free", labeller = labeller(pi = label_parsed)) + 
  theme_bw()+
  theme(legend.title = element_text(size = 18, face = "bold"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 15), 
        legend.position = "bottom") 
#scale_fill_manual(values = c("#D62728FF", "#E377C2FF", "#9467BDFF", "#1F77B4FF", 
#  "#FF7F0EFF", "#2CA02CFF","#7F7F7FFF", "#BCBD22FF"))

ggsave(box_ARI, filename = "reproduce/figures_and_tables/figure1.png", width = 10, height = 10)


 


