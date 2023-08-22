
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
###### show the K estimation of simulation 1, figure 2 
###################################################################
data1 = dataRader("simulation1/scenario1_1")
data1 = cbind(data1, K_true = rep("K=3", dim(data1)[1]), pi = rep("0.1", dim(data1)[1]))

data2 = dataRader("simulation1/scenario1_2")
data2 = cbind(data2, K_true = rep("K=3", dim(data2)[1]), pi = rep("0.2", dim(data2)[1]))

data3 = dataRader("simulation1/scenario1_3")
data3 = cbind(data3, K_true = rep("K=3", dim(data3)[1]), pi = rep("0.3", dim(data3)[1]))


data4 = dataRader("simulation1/scenario2_1")
data4 = cbind(data4, K_true = rep("K=5", dim(data4)[1]), pi = rep("0.1", dim(data4)[1]))

data5 = dataRader("simulation1/scenario2_2")
data5 = cbind(data5, K_true = rep("K=5", dim(data5)[1]), pi = rep("0.2", dim(data5)[1]))

data6 = dataRader("simulation1/scenario2_3")
data6 = cbind(data6, K_true = rep("K=5", dim(data6)[1]), pi = rep("0.3", dim(data6)[1]))


data7 = dataRader("simulation1/scenario3_1")
data7 = cbind(data7, K_true = rep("K=7", dim(data7)[1]), pi = rep("0.1", dim(data7)[1]))

data8 = dataRader("simulation1/scenario3_2")
data8 = cbind(data8, K_true = rep("K=7", dim(data8)[1]), pi = rep("0.2", dim(data8)[1]))

data9 = dataRader("simulation1/scenario3_3")
data9 = cbind(data9, K_true = rep("K=7", dim(data9)[1]), pi = rep("0.3", dim(data9)[1]))

data_K = rbind(data1, data2, data3, data4, data5, data6, data7, data8, data9)

data_K$Method = factor(data_K$Method, levels = c("BNPSpace", "BNPSpace w/o MRF", "SCMEB", 
                                                 "DRSC", "Louvain", "GMM"))
data_K$pi = factor(data_K$pi, levels = c("0.1", "0.2", "0.3"), 
                   labels = c("0.1" = expression(pi == "0.1"), "0.2" = expression(pi == "0.2"), "0.3" = expression(pi == "0.3")))


####### first try for stacked barplot
data_K$K = as.character(data_K$K)
library(ggplot2)
data_K$value = rep(1, dim(data_K)[1])
data_K = aggregate(data_K$value, by=list(data_K$Method, data_K$K_true, data_K$pi, data_K$K), sum)
colnames(data_K) = c("Method", "K_true", "pi", "K", "value")
data_K$K = factor(data_K$K, levels = c("2", "3", "4", "5", "6","7","8","9", "10"))

data_K$value = data_K$value / 50


barplot_K = ggplot(data = data_K, aes(x = Method, y = value, fill = K)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  facet_grid(pi ~ K_true, scales = "free", labeller = labeller(pi = label_parsed)) + 
  theme_bw()+ ylab("Frequency") +
  theme(legend.title = element_text(size = 18, face = "bold"), 
        axis.title.x = element_blank(), 
        axis.text.x = element_text(face = "bold", size = 12, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 18, face = "bold"),
        strip.text = element_text(face = "bold", size = 18),
        legend.text = element_text(size = 15), 
        legend.position = "bottom") 
ggsave(barplot_K, filename = "reproduce/figures_and_tables/figure2.png", width = 10, height = 10)






