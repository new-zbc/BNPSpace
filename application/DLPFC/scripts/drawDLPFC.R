library(SingleCellExperiment)
library(SC.MEB)
library(Matrix)
library(flexclust)
library(scater)
library(scran)
library(patchwork)
source("NSCFS/utils.R")
source("plot/DLPFC/plotDLPFCfunctions.R")

sampleID = 151509
#features = c("MOBP", "S100A11", "GFAP", "NEFL")
features = c("SCGB2A2", "SCGB1D2", "GFAP",
             "MBP", "COX6C", "FABP7", "MALAT1", "NEFL")

pARI = plot_ARI(sampleID)

pExp = plot_exp(sampleID)

save(pARI, pExp, file = "real_data/DLPFC_data/151509/results/ARI_Exp.RData")

pGene = spatialPlotGene(sampleID, features = features)

cluster = plot_cluster(sampleID)

save(cluster, file = "real_data/DLPFC_data/151509/results/cluster.RData")


ggsave(cluster, filename = "real_data/DLPFC_data/151509/results/cluster_result.pdf", width = 12, height = 6)
ggsave(cluster, filename = "real_data/DLPFC_data/151509/results/cluster_result.png", width = 12, height = 6)



x = merge_list(sampleID)

x[[10]] <- x[[10]] + theme(plot.title = element_blank(),
                           panel.border = element_blank(),
                           legend.title = element_text(face = "bold", size = 16),
                           legend.text = element_text(face = "bold", size = 12))

x[[10]] <- x[[10]] + scale_fill_manual(values = c("#E377C2FF", "#F7B6D2FF", "#8C564BFF", "#9467BDFF", 
                                                  "#2CA02CFF", "#1F77B4FF", "#17BECFFF", "#FF9896FF", 
                                                  "#FFBB78FF", "#AEC7E8FF", "#D62728FF", "#C49C94FF", "#C7C7C7FF"))
                                         
 

save(x, file = "real_data/DLPFC_data/151509/results/merge_list.RData")                                        
        
# layout = "
# AAAAAA
# BBBCCC
# DDDDEE
# "
# 
# p = cluster + pARI + pExp + pGene[[1]] +  x[[10]] + 
#   plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(1, 0.75, 0.75, 0.4, 0.5))

#pdf("real_data/DLPFC_data/151509/results/summary.pdf")


layout = "
AAABBB
CCCCDD
"

p = pARI + pExp + pGene[[1]] +  x[[10]] +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(0.75, 0.75, 0.4, 0.5))


ggsave(p, filename = "real_data/DLPFC_data/151509/results/postprocess.pdf", width = 20, height = 12.5, units = "in")
ggsave(p, filename = "real_data/DLPFC_data/151509/results/postprocess.png", width = 20, height = 12.5, units = "in")




pMerge = plot_merge(sampleID)
ggsave(pMerge, filename = "real_data/DLPFC_data/151509/results/merge.eps", width = 20, height = 15, units = "in")




###### another sample
sampleID = 151507
pARI = plot_ARI(sampleID)
pExp = plot_exp(sampleID)
#pExp
cluster = plot_cluster(sampleID)

layout = "
AAAAAA
BBBCCC
"
p = cluster + pARI + pExp +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(1, 0.75, 0.75))

ggsave(p, filename = "real_data/DLPFC_data/151507/results/summary.png", width = 20, height = 17.5)


###### another sample
sampleID = 151508
pARI = plot_ARI(sampleID)
pExp = plot_exp(sampleID)
#pExp
cluster = plot_cluster(sampleID)

layout = "
AAAAAA
BBBCCC
"
p = cluster + pARI + pExp +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(1, 0.75, 0.75))

ggsave(p, filename = "real_data/DLPFC_data/151508/results/summary.png", width = 20, height = 17.5)


###### another sample
sampleID = 151510
pARI = plot_ARI(sampleID)
pExp = plot_exp(sampleID)
#pExp
cluster = plot_cluster(sampleID)

layout = "
AAAAAA
BBBCCC
"
p = cluster + pARI + pExp +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(1, 0.75, 0.75))

ggsave(p, filename = "real_data/DLPFC_data/151510/results/summary.png", width = 20, height = 17.5)


###### another sample
sampleID = 151669
pARI = plot_ARI(sampleID)
pExp = plot_exp(sampleID)
#pExp
cluster = plot_cluster(sampleID)

layout = "
AAAAAA
BBBCCC
"
p = cluster + pARI + pExp +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(1, 0.75, 0.75))

ggsave(p, filename = "real_data/DLPFC_data/151669/results/summary.png", width = 20, height = 17.5)


###### another sample
sampleID = 151670
pARI = plot_ARI(sampleID)
pExp = plot_exp(sampleID)
#pExp
cluster = plot_cluster(sampleID)

layout = "
AAAAAA
BBBCCC
"
p = cluster + pARI + pExp +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(1, 0.75, 0.75))

ggsave(p, filename = "real_data/DLPFC_data/151670/results/summary.png", width = 20, height = 17.5)


###### another sample
sampleID = 151671
pARI = plot_ARI(sampleID)
pExp = plot_exp(sampleID)
#pExp
cluster = plot_cluster(sampleID)

layout = "
AAAAAA
BBBCCC
"
p = cluster + pARI + pExp +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(1, 0.75, 0.75))

ggsave(p, filename = "real_data/DLPFC_data/151671/results/summary.png", width = 20, height = 17.5)


###### another sample
sampleID = 151672
pARI = plot_ARI(sampleID)
pExp = plot_exp(sampleID)
#pExp
cluster = plot_cluster(sampleID)

layout = "
AAAAAA
BBBCCC
"
p = cluster + pARI + pExp +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(1, 0.75, 0.75))

ggsave(p, filename = "real_data/DLPFC_data/151672/results/summary.png", width = 20, height = 17.5)


###### another sample
sampleID = 151673
pARI = plot_ARI(sampleID)
pExp = plot_exp(sampleID)
#pExp
cluster = plot_cluster(sampleID)

layout = "
AAAAAA
BBBCCC
"
p = cluster + pARI + pExp +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(1, 0.75, 0.75))

ggsave(p, filename = "real_data/DLPFC_data/151673/results/summary.png", width = 20, height = 17.5)


###### another sample
sampleID = 151674
pARI = plot_ARI(sampleID)
pExp = plot_exp(sampleID)
#pExp
cluster = plot_cluster(sampleID)

layout = "
AAAAAA
BBBCCC
"
p = cluster + pARI + pExp +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(1, 0.75, 0.75))

ggsave(p, filename = "real_data/DLPFC_data/151674/results/summary.png", width = 20, height = 17.5)


###### another sample
sampleID = 151675
pARI = plot_ARI(sampleID)
pExp = plot_exp(sampleID)
#pExp
cluster = plot_cluster(sampleID)

layout = "
AAAAAA
BBBCCC
"
p = cluster + pARI + pExp +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(1, 0.75, 0.75))

ggsave(p, filename = "real_data/DLPFC_data/151675/results/summary.png", width = 20, height = 17.5)


###### another sample
sampleID = 151676
pARI = plot_ARI(sampleID)
pExp = plot_exp(sampleID)
#pExp
cluster = plot_cluster(sampleID)

layout = "
AAAAAA
BBBCCC
"
p = cluster + pARI + pExp +
  plot_annotation(tag_levels = "A") + plot_layout(design = layout, heights = c(1, 0.75, 0.75))

ggsave(p, filename = "real_data/DLPFC_data/151676/results/summary.png", width = 20, height = 17.5)


