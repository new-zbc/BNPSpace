source("application/MOB/scripts/plotMOBfunctions.R")
load("application/MOB/clustering_result_sce.RData")
sce = sce[, !is.na(colData(sce)$label)]

colData(sce)$BNPSpace[colData(sce)$BNPSpace == 4] = 0
colData(sce)$BNPSpace[colData(sce)$BNPSpace == 1] = 4
colData(sce)$BNPSpace[colData(sce)$BNPSpace == 0] = 1

colData(sce)$BNPSpace[colData(sce)$BNPSpace == 3] = 0
colData(sce)$BNPSpace[colData(sce)$BNPSpace == 1] = 3
colData(sce)$BNPSpace[colData(sce)$BNPSpace == 0] = 1

colData(sce)$BNPSpace[colData(sce)$BNPSpace == 2] = 0
colData(sce)$BNPSpace[colData(sce)$BNPSpace == 5] = 2
colData(sce)$BNPSpace[colData(sce)$BNPSpace == 0] = 5

colData(sce)$SCMEB[colData(sce)$SCMEB == 1] = 0
colData(sce)$SCMEB[colData(sce)$SCMEB == 4] = 1
colData(sce)$SCMEB[colData(sce)$SCMEB == 0] = 4

p1 = plot_cluster(sce)


p2 = plot_exp(sce)



layout = "
AAAA
BBBB
"

p = p1 + p2 + plot_layout(design = layout, heights = c(0.75, 1))


ggsave(p, filename = "reproduce/figures_and_tables/figure4ab.png", width = 16, height = 14)

ggsave(p, filename = "reproduce/figures_and_tables/figure4ab.pdf", width = 16, height = 14)


