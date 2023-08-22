plot_PPI = function(result){
  gamma_est = rowData(result)$gamma
  gamma_est = gamma_est[order(gamma_est, decreasing = T)]
  df = data.frame(gene = rep(1:length(gamma_est)), PPI = gamma_est)
  library(ggplot2)
 
  p <- ggplot() + geom_point(aes(x = gene, y = PPI), data = df) + 
    geom_hline(aes(yintercept = 0.5), color = "blue", linetype = "dashed")+
    #geom_point(data = df[1:20,], aes(x = gene, y = PPI, color = "red")) + 
    theme_bw() + xlab("Gene index") + ylab("PPI")+ 
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 10,face = "bold"),
          legend.title = element_text(size = 10,face = "bold"),
          axis.text = element_text(size = 8, face = "bold"))
  p
}

load("application/MOBdata/result_sce.RData")
p5 = plot_PPI(sce2)


load("application/MOBdata/BIC_result.RData")

data_mDIC = data.frame(f = c(0, 0.5, 1, 1.5, 2, 2.5), mDIC = BIC_result[1:6])

p4 = ggplot(data_mDIC, aes(x = f, y = mDIC)) + 
  geom_line()+geom_point() + theme_bw() + ylab("pBIC") + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 10,face = "bold"),
        legend.title = element_text(size = 10,face = "bold"),
        axis.text = element_text(size = 8, face = "bold")) 

p = cowplot::plot_grid(p4, p5, ncol = 2, labels = c("a", "b"), rel_widths = c(1.1, 1) )
ggsave(p, filename = "reproduce/figures_and_tables/figureS3.png", width = 6.3, height = 3, units = "in")

