
df = data.frame(fold = c(1.877, 1.368, 1.072), 
                method = c("BNPSpace", "SPARK", "SpatialDE"), 
                p_value = c(0.00162, 0.0228, 0.0728))
library(ggplot2)
p = ggplot(data = df) + 
  geom_col(aes(x = method, y = fold, alpha = 0.8)) + ylab("Fold enrichment (odd ratio)") +
  annotate(geom = "text", x = "BNPSpace", y = 2.0, label = "P = 0.00162", size = 4) + 
  annotate(geom = "text", x = "SPARK", y = 1.5, label = "P = 0.0228", size = 4) + 
  annotate(geom = "text", x = "SpatialDE", y = 1.2, label = "P = 0.0728", size = 4) + 
  theme_classic() + theme(axis.title.x = element_blank(), 
                     legend.position = "none", 
                     axis.title.y = element_text(face = "bold", size = 14),
                     axis.text.x = element_text(face = "bold", size = 14),
                     axis.text.y = element_text(face = "bold", size = 12))

ggsave(p, filename = "reproduce/figures_and_tables/figure5F.png", width = 4, height = 4, bg = "white")