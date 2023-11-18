source("reproduce/scripts/figure2a.R")
source("reproduce/scripts/figure2b.R")
library(patchwork)
p = plot_pattern + box_ARI + plot_layout(nrow = 2, heights = c(1, 3))

ggsave(p, filename = "reproduce/figures_and_tables/figure2.png", width = 13, height = 13)
ggsave(p, filename = "reproduce/figures_and_tables/figure2.pdf", width = 13, height = 13)
