#################################################################################
# Analyze E80 fetal pig data from Seurat object.
#   Compare clusters: Violin plots
#
#---------------------------------------------------------
# Ross Lampe
# Piedrahita Lab, College of Veterinary Medicine, North Carolina State University
# Last updated: 2022-10-06
#
#################################################################################
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(here)

source("../themes.R")

#################################################################################
# Violin plots
#   
#################################################################################
VlnPlot(obj_seurat, features=c("EPCAM","LGR5"),
        pt.size=1,
        sort=T, stack=T, flip=T) + 
  NoLegend() +
  labs(x = "Cluster")


#-------------------------------------------------------------------------------
# Plot separately, with dots
#-------------------------------------------------------------------------------
color_Mesc = "#198BFD"   #neutral blue
color_Epi = "#EC6464"    #neutral red
#color_EMT = "lightgreen"   #"lightgrey","darkgreen", "darkred"

col_VlnPlot <- c(color_Mesc,color_Mesc,color_Mesc,color_Mesc,color_Mesc,color_Mesc,color_Mesc,
                 color_Epi,color_Epi,color_Mesc,color_Epi,color_Mesc,color_Mesc,color_Epi)

p1 <- VlnPlot(obj_seurat, features = "EPCAM", pt.size = 0.2, cols=col_VlnPlot) +
  NoLegend() +
  theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) +
  theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.02, vjust = -1))

p2 <- VlnPlot(obj_seurat, features = "LGR5", pt.size = 0.5, cols=col_VlnPlot) +
    NoLegend() +
    theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank()) +
    theme(plot.title = element_text(size = 14, face = "bold",hjust = 0.02, vjust = -1))

plot_grid(p1, p2,
          nrow=2, ncol=1)


