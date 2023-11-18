df = data.frame(fold = c(2.224, 1.312, 1.221), 
                method = c("BNPSpace", "SPARK", "SpatialDE"), 
                p_value = c(0.0001, 0.0912, 0.0315))
library(ggplot2)
p = ggplot(data = df) + 
  geom_col(aes(x = method, y = fold, alpha = 0.8)) + ylab("Fold enrichment (odd ratio)") +
  annotate(geom = "text", x = "BNPSpace", y = 2.4, label = "P < 0.0001", size = 5) + 
  annotate(geom = "text", x = "SPARK", y = 1.5, label = "P = 0.0912", size = 5) + 
  annotate(geom = "text", x = "SpatialDE", y = 1.4, label = "P = 0.0315", size = 5) + 
  theme_classic() + theme(axis.title.x = element_blank(), 
                          legend.position = "none", 
                          axis.title.y = element_text(face = "bold", size = 14),
                          axis.text.x = element_text(face = "bold", size = 14),
                          axis.text.y = element_text(face = "bold", size = 12))

ggsave(p, filename = "reproduce/figures_and_tables/figure4c.png", width = 4, height = 3, bg = "white")