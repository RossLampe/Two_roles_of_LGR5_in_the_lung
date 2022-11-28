
#################################################################################
# Analyze E80 fetal pig data from Seurat object.
#   DotPlots, FeaturePlots of transcript expression.
#
#---------------------------------------------------------
# Ross Lampe
# Piedrahita Lab, College of Veterinary Medicine, North Carolina State University
# Last updated: 2022-09-08
#
#################################################################################
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(here)

source("../themes.R")


#-------------------------------------------------------------------------------
# Collect the markers
#-------------------------------------------------------------------------------
df_Markers_Ref <- read.csv(here("../Markers","KEGG_pathways","KEGG_TGFbeta_pathway_shortlist_ssc04350.csv"), header = FALSE)
vec_TGFbeta_Markers <- df_Markers_Ref$V1

df_Markers_Ref <- read.csv(here("../Markers","KEGG_pathways","KEGG_WNT_pathway_Adultlist_ssc04310.csv"), header = FALSE)
vec_WNT_Markers <- df_Markers_Ref$V1

df_Markers_Ref <- read.csv(here("../Markers","Epithelial_Markers.csv"), header = FALSE)
vec_Budtip_Epi_Markers <- df_Markers_Ref$V1

df_Markers_Ref <- read.csv(here("../Markers","Neural_markers_Gerber.csv"), header = FALSE)
vec_Neural_Markers <- df_Markers_Ref$V1

df_Markers_Ref <- read.csv(here("../Markers","Neural_markers_shortlist_2.csv"), header = FALSE)
vec_Neural_Markers_Shortlist <- df_Markers_Ref$V1

df_Markers_Ref <- read.csv(here("../Markers","Neural_markers_shortlist_2_SC_PnC.csv"), header = FALSE)
vec_Neural_Markers_Shortlist_SC_PnC <- df_Markers_Ref$V1

df_Markers_Ref <- read.csv(here("../Markers","Neural_markers_shortlist_2_EpC_EnC.csv"), header = FALSE)
vec_Neural_Markers_Shortlist_EpC_EnC <- df_Markers_Ref$V1

vec_LGR5_other <- c("LGR5","EPCAM")


#-------------------------------------------------------------------------------
# Compute the log2fold change of obj_seurat_CombinedClusters
#-------------------------------------------------------------------------------
fc_DE_CombinedClusters <- FoldChange(obj_seurat_CombinedClusters,
                                     ident.1 = "Epithelial",
                                     ident.2 = "Mesenchymal")

#Save the rownames (genes) to a column.
fc_DE_CombinedClusters$gene <- rownames(fc_DE_CombinedClusters)

#Remove undefined genes
fc_DE_CombinedClusters[!grepl('ENSSSC',fc_DE_CombinedClusters$gene),]
head(fc_DE_CombinedClusters)



################################################################################
# Epithelial Markers: Epithelial cluster.
#
################################################################################
vec_geneSymbols <- c(vec_Budtip_Epi_Markers)
vec_geneSymbols <- toupper(unique(c(vec_geneSymbols)))
vec_geneSymbols <- sort(vec_geneSymbols)

p1 <- DotPlot(obj_seurat, features=as.character(vec_geneSymbols),
              idents = vec_Epi_Clusters, 
              assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("Epithelial Markers; scRNAseq GD80")
p1 <- theme_DotPlot(p1)
p1

#-------------------------------------------------------------------------------
# Write the dotpot log2fold change data to disk.
#-------------------------------------------------------------------------------
#Filter all_markers to just the vec_geneSymbols and sort
fc_DE_CombinedClusters_filtered <- fc_DE_CombinedClusters[fc_DE_CombinedClusters$gene %in% vec_geneSymbols,]
fc_DE_CombinedClusters_sorted <- fc_DE_CombinedClusters_filtered[order(fc_DE_CombinedClusters_filtered$avg_log2FC, decreasing = TRUE), ]
head(fc_DE_CombinedClusters_sorted,20)

write.csv(fc_DE_CombinedClusters_sorted, here("Results","scRNAseq GD80_Epithelial_Markers_Epi_Cluster.csv"))



################################################################################
# Epithelial Markers: Epi vs Mesc clusters.
#
################################################################################
vec_geneSymbols <- c(vec_Budtip_Epi_Markers)
vec_geneSymbols <- toupper(unique(c(vec_geneSymbols)))
vec_geneSymbols <- sort(vec_geneSymbols)

p1 <- DotPlot(obj_seurat_CombinedClusters, features=as.character(vec_geneSymbols),
              #idents = vec_Epi_Clusters, 
              assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("Epithelial Markers; scRNAseq GD80")
p1 <- theme_DotPlot(p1)
p1

#-------------------------------------------------------------------------------
# Write the dotpot log2fold change data to disk.
#-------------------------------------------------------------------------------
#Filter all_markers to just the vec_geneSymbols and sort
fc_DE_CombinedClusters_filtered <- fc_DE_CombinedClusters[fc_DE_CombinedClusters$gene %in% vec_geneSymbols,]
fc_DE_CombinedClusters_sorted <- fc_DE_CombinedClusters_filtered[order(fc_DE_CombinedClusters_filtered$avg_log2FC, decreasing = TRUE), ]
head(fc_DE_CombinedClusters_sorted,20)

write.csv(fc_DE_CombinedClusters_sorted, here("Results","scRNAseq GD80_Epithelial_Markers_EpiVsMesc_Clusters.csv"))




################################################################################
# WNT markers 
#
################################################################################
vec_geneSymbols <- c(vec_LGR5_other,
                     vec_WNT_Markers)
vec_geneSymbols <- toupper(unique(vec_geneSymbols))  

p1 <- DotPlot(obj_seurat_CombinedClusters, features=as.character(vec_geneSymbols),
              assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("KEGG_WNT_pathway_ssc04310; scRNAseq GD80")
p1 <- theme_DotPlot(p1)
p1

#-------------------------------------------------------------------------------
# Write the dotpot log2fold change data to disk.
#-------------------------------------------------------------------------------
#Filter all_markers to just the vec_geneSymbols and sort
fc_DE_CombinedClusters_filtered <- fc_DE_CombinedClusters[fc_DE_CombinedClusters$gene %in% vec_geneSymbols,]
fc_DE_CombinedClusters_sorted <- fc_DE_CombinedClusters_filtered[order(fc_DE_CombinedClusters_filtered$avg_log2FC, decreasing = TRUE), ]
head(fc_DE_CombinedClusters_sorted)

write.csv(fc_DE_CombinedClusters_sorted, here("Results","scRNAseq GD80_KEGG_WNT_pathway.csv"))



################################################################################
# TGFbeta Markers
#
################################################################################
vec_geneSymbols <- c(vec_LGR5_other,
                     vec_TGFbeta_Markers)

vec_geneSymbols <- toupper(unique(vec_geneSymbols)) 

p1 <- DotPlot(obj_seurat_CombinedClusters, features=as.character(vec_geneSymbols),
              assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("KEGG_TGFβ_ssc04350; scRNAseq GD80")
p1 <- theme_DotPlot(p1)
p1

#-------------------------------------------------------------------------------
# Write the dotpot log2fold change data to disk.
#-------------------------------------------------------------------------------
#Filter all_markers to just the vec_geneSymbols and sort
fc_DE_CombinedClusters_filtered <- fc_DE_CombinedClusters[fc_DE_CombinedClusters$gene %in% vec_geneSymbols,]
fc_DE_CombinedClusters_sorted <- fc_DE_CombinedClusters_filtered[order(fc_DE_CombinedClusters_filtered$avg_log2FC, decreasing = TRUE), ]
head(fc_DE_CombinedClusters_sorted)

write.csv(fc_DE_CombinedClusters_sorted, here("Results","scRNAseq GD80_KEGG_TGFβ.csv"))



################################################################################
# Combined WNT, TGFbeta and budtip markers
#
################################################################################
vec_geneSymbols <- c(vec_LGR5_other,
                     vec_Budtip_Epi_Markers,
                     vec_WNT_Markers,
                     vec_TGFbeta_Markers)

vec_geneSymbols <- toupper(unique(vec_geneSymbols))  

p1 <- DotPlot(obj_seurat_CombinedClusters, features=as.character(vec_geneSymbols),
              assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("Budtip, WNT and TGF-β markers; scRNAseq GD80")
p1 <- theme_DotPlot(p1)
p1

#-------------------------------------------------------------------------------
# Write the dotpot log2fold change data to disk.
#-------------------------------------------------------------------------------
#Filter all_markers to just the vec_geneSymbols and sort
fc_DE_CombinedClusters_filtered <- fc_DE_CombinedClusters[fc_DE_CombinedClusters$gene %in% vec_geneSymbols,]
fc_DE_CombinedClusters_sorted <- fc_DE_CombinedClusters_filtered[order(fc_DE_CombinedClusters_filtered$avg_log2FC, decreasing = TRUE), ]
head(fc_DE_CombinedClusters_sorted)

write.csv(fc_DE_CombinedClusters_sorted, here("Results","scRNAseq GD80_Budtip, WNT and TGF-β markers.csv"))





################################################################################
# Volcano plot, Neural markers: SC+PnC, EpC,EnC
#
################################################################################
#-------------------------------------------------------------------------------
# SC+PnC
#-------------------------------------------------------------------------------
vec_geneSymbols <- unique(toupper(vec_Neural_Markers_Shortlist_SC_PnC))
vec_geneSymbols <- sort(vec_geneSymbols)

p1 <- DotPlot(obj_seurat_CombinedClusters, features=as.character(vec_geneSymbols),
              assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("SC, PnC Neural Markers; scRNAseq  GD80")
p1 <- theme_DotPlot(p1)
p1

#-------------------------------------------------------------------------------
# EpC_EnC
#-------------------------------------------------------------------------------
vec_geneSymbols <- unique(toupper(vec_Neural_Markers_Shortlist_EpC_EnC))
vec_geneSymbols <- sort(vec_geneSymbols)

p1 <- DotPlot(obj_seurat_CombinedClusters, features=as.character(vec_geneSymbols),
              assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("EpC, EnC Neural Markers; scRNAseq  GD80")
p1 <- theme_DotPlot(p1)
p1



################################################################################
# Neural markers
# 
################################################################################
#vec_geneSymbols <- toupper(unique(vec_Neural_Markers)) 
vec_geneSymbols <- toupper(unique(vec_Neural_Markers_Shortlist))
vec_geneSymbols <- sort(vec_geneSymbols)

p1 <- DotPlot(obj_seurat_CombinedClusters, features=as.character(vec_geneSymbols),
              assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("Neural Markers; scRNAseq  GD80")
p1 <- theme_DotPlot(p1)
p1


#-------------------------------------------------------------------------------
# p1 <- DotPlot(obj_seurat_CombinedClusters, features=as.character(vec_geneSymbols),
#               #idents = vec_Mesc_Clusters, 
#               assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
#   ggtitle("Neural Markers; Pig Lung, Fetal P111, Epithelial Clusters; Piedrahita")
#   #ggtitle("Neural Markers; Fetal Pig Lung, Fetal P111, Mesenchymal Clusters; Piedrahita")
# p1 <- theme_DotPlot(p1)
# 
# p2 <- DotPlot(obj_seurat, features=as.character(vec_geneSymbols),
#               #idents = vec_Epi_Clusters, 
#               idents = vec_Mesc_Clusters, 
#               assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
#   #ggtitle("Neural Markers; Fetal Pig Lung, Fetal P111, Epithelial Clusters; Piedrahita")
#   ggtitle("Neural Markers; Fetal Pig Lung, Fetal P111, Mesenchymal Clusters; Piedrahita")
# p2 <- theme_DotPlot(p2)
# 
# plot_grid(p1 + NoLegend(),
#           p2 + NoLegend(),
#           nrow=, ncol=1, 
#           labels=c("", "", ""))


#-------------------------------------------------------------------------------
# Neural markers - long list
#-------------------------------------------------------------------------------
# df_Neural_Markers_Ref <- read.csv(here("../Markers","Neural_markers_shared.csv"), header = FALSE)
# 
# vec_LGR5_other <- c("LGR5","EPCAM")
# vec_Neural <- sort(df_Neural_Markers_Ref$V1)
# vec_geneSymbols <- c(vec_LGR5_other,
#                      vec_Neural)
# 
# #Convert uppercase gene symbols to mouse symbols (like COL1A1 to Col1a1)
# vec_geneSymbols <- toupper(unique(vec_geneSymbols))  
# length(vec_geneSymbols)
# #[1] 214
# 
# vec_geneSymbols_A = vec_geneSymbols[1:72]
# vec_geneSymbols_B = vec_geneSymbols[73:146]
# vec_geneSymbols_C = vec_geneSymbols[147:215]
# 
# #DotPlot
# p1 <- DotPlot(obj_seurat, features=as.character(vec_geneSymbols_A),
#               #idents = vec_Epi_Clusters, 
#               assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
#   ggtitle("Neural Markers, Series A; Fetal Pig Lung; Piedrahita")
# p1 <- theme_DotPlot(p1)
# 
# p2 <- DotPlot(obj_seurat, features=as.character(vec_geneSymbols_B),
#              # idents = vec_Epi_Clusters, 
#               assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
#   ggtitle("Neural Markers, Series B; Fetal Pig Lung; Piedrahita")
# p2 <- theme_DotPlot(p2)
# 
# p3 <- DotPlot(obj_seurat, features=as.character(vec_geneSymbols_C),
#              # idents = vec_Epi_Clusters, 
#               assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
#   ggtitle("Neural Markers, Series C; Fetal Pig Lung; Piedrahita")
# p3 <- theme_DotPlot(p3)
# 
# plot_grid(p1 + NoLegend(),
#           p2 + NoLegend(),
#           p3 + NoLegend(),  
#           nrow=3, ncol=1, 
#           labels=c("", "", ""))

#-------------------------------------------------------------------------------
# log2fold change between Developing airway and bud tip clusters
#-------------------------------------------------------------------------------
# fc_DE <- FoldChange(obj_seuratv)
# 
# #Save the rownames (genes) to a column.
# fc_DE$gene <- rownames(fc_DE)
# 
# #Filter all_markers to just the vec_geneSymbols. Sort
# fc_DE_filtered <- fc_DE[fc_DE$gene %in% vec_geneSymbols,]
# fc_DE_sorted <- fc_DE_filtered[order(fc_DE_filtered$avg_log2FC, decreasing = TRUE), ]  
# 
# #Remove undefined genes
# fc_DE_sorted <- fc_DE_sorted[!grepl('ENSSSC',fc_DE_sorted$gene),] 
# 
# write.csv(fc_DE_sorted, here("Results","fc_Neural_cluster_LGR5_vs_EPCAM.csv"))




################################################################################
# Collagen, WNT, Neural related markers - shortlist
#   Mesenchymal clusters
# 
################################################################################
# p_Collagen_Wnt_neural_Markers_Ref <- read.csv(here("../Markers","Collagen_Neural_WNT_shortlist.csv"), header = FALSE)
# 
# vec_LGR5_other <- c("LGR5","VIM","EPCAM")
# vec_Common_Markers <- p_Collagen_Wnt_neural_Markers_Ref$V1
# # vec_Fibro_Neural <- df_Collagen_Markers_Ref$V1[df_shortlist_common_markers$V2 == "Fibro; Neural"]
# # vec_Pericyte_Endothelial_Advential_VSMC<- df_Collagen_Markers_Ref$V1[df_shortlist_common_markers$V2 == "Pericyte; Endothelial; Advential; VSMC"]
# # vec_Wnt_ECM <- df_Collagen_Markers_Ref$V1[df_Collagen_Markers_Ref$V2 == "Wnt; ECM"]
# 
# vec_geneSymbols <- c(vec_Common_Markers)
# vec_geneSymbols <- sort(vec_geneSymbols)
# 
# #Convert to uppercase gene symbols.
# vec_geneSymbols <- toupper(unique(vec_geneSymbols))
# #length(vec_geneSymbols)
# 
# #DotPlot
# p_Collagen_Wnt_neural <- DotPlot(obj_seurat_CombinedClusters, features=as.character(vec_geneSymbols),
#                                 # idents = vec_Mesc_Clusters,
#                                  assay="RNA", cols = colors_DotPlot) + 
#   RotatedAxis() +
#   ggtitle("Fetal Pig Lung, E80; Collagen, Neural, WNT Markers")
# p_Collagen_Wnt_neural <- theme_DotPlot(p_Collagen_Wnt_neural)
# p_Collagen_Wnt_neural
# 
# 
# #-------------------------------------------------------------------------------
# # FeaturePlot: pct.exp, avg.exp.scaled by cluster
# #-------------------------------------------------------------------------------
# # % Expression
# p_Collagen_Wnt_neural_FeaturePlot_PC_Exp <- ggplot(p_Collagen_Wnt_neural$data,aes(x=id, y=pct.exp, fill=id)) +
#   geom_col() +
#   labs(fill='') +   #legend title
#   xlab('') + ylab('% Expression') +     #% Expression
#   facet_wrap(~features.plot) +
#   ggtitle("CWN: Percent Expression; Fetal Pig Lung; Piedrahita")  +
#   theme(axis.text.x=element_blank())    #cluster names are long. Identify by colors.
# 
# # Average Expression
# p_Collagen_Wnt_neural_FeaturePlot_AvgExp<- ggplot(p_Collagen_Wnt_neural$data,aes(x=id, y=avg.exp.scaled, fill=id)) +
#   geom_col() +
#   labs(fill='') +   #legend title
#   xlab('') + ylab('Average Expression') +    #
#   facet_wrap(~features.plot) +
#   ggtitle("CWN: Scaled Average Expression; Fetal Pig Lung; Piedrahita")  +
#   theme(axis.text.x=element_blank())    #cluster names are long.
# 
# plot_grid(p_Collagen_Wnt_neural_FeaturePlot_PC_Exp + NoLegend(),
#           p_Collagen_Wnt_neural_FeaturePlot_AvgExp,
#           nrow=1, ncol=2,
#           labels=c("A", ""),
#           rel_widths = c(1, 1.4))

#-------------------------------------------------------------------------------
# log2fold change between Developing airway and bud tip clusters
#-------------------------------------------------------------------------------
# fc_DE <- FoldChange(obj_seurat,
#                     ident.1 = "LGR5+ Fibroblast",  #"Developing Airway 2"
#                     ident.2 = "EPCAM+")
# 
# #Save the rownames (genes) to a column.
# fc_DE$gene <- rownames(fc_DE)
# 
# #Filter all_markers to just the vec_geneSymbols. Sort
# fc_DE_filtered <- fc_DE[fc_DE$gene %in% vec_geneSymbols,]
# fc_DE_sorted <- fc_DE_filtered[order(fc_DE_filtered$avg_log2FC, decreasing = TRUE), ]  
# 
# #Remove undefined genes
# fc_DE_sorted <- fc_DE_sorted[!grepl('ENSSSC',fc_DE_sorted$gene),] 
# 
# write.csv(fc_DE_sorted, here("Results","fc_Collagen_Neural_WNT_cluster_LGR5_vs_EPCAM.csv"))



################################################################################
# Collagen-WNT-Neural transcripts upregulated in the adult pig lung (D111) 
#   Mesenchymal clusters
#
# ################################################################################
# df_Markers_Ref <- read.csv(here("../Markers","Collagen_Neural_WNT_adult_pig_upregulated.csv"), header = TRUE)
# head(df_Markers_Ref,2)
# 
# vec_LGR5_other <- c("LGR5","VIM","EPCAM")
# vec_Common_Markers <- df_Markers_Ref$Gene
# vec_geneSymbols <- c(vec_Common_Markers)
# vec_geneSymbols <- sort(vec_geneSymbols)
# 
# #Convert to uppercase gene symbols
# vec_geneSymbols <- toupper(unique(vec_geneSymbols))  
# 
# #DotPlot
# p_Collagen_Wnt_neural <- DotPlot(obj_seurat_mesc, features=as.character(vec_geneSymbols),
#                                  assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
#   ggtitle("Collagen, Neural, WNT Markers Upregulated in the Adult Pig; Fetal Pig Lung; Mesenchymal Clusters; Piedrahita")
# p_Collagen_Wnt_neural <- theme_DotPlot(p_Collagen_Wnt_neural)
# p_Collagen_Wnt_neural
# 
# 
# #-------------------------------------------------------------------------------
# # FeaturePlot: pct.exp, avg.exp.scaled by cluster
# #-------------------------------------------------------------------------------
# # % Expression
# p_Collagen_Wnt_neural_FeaturePlot_PC_Exp <- ggplot(p_Collagen_Wnt_neural$data,aes(x=id, y=pct.exp, fill=id)) +
#   geom_col() +
#   labs(fill='') +   #legend title
#   xlab('') + ylab('% Expression') +     #% Expression
#   facet_wrap(~features.plot) + 
#   ggtitle("CWN: Percent Expression; Fetal Pig Lung; Mesc Clusters")  +
#   theme(axis.text.x=element_blank())    #cluster names are long. Identify by colors.
# 
# # Average Expression
# p_Collagen_Wnt_neural_FeaturePlot_AvgExp<- ggplot(p_Collagen_Wnt_neural$data,aes(x=id, y=avg.exp.scaled, fill=id)) +
#   geom_col() +
#   labs(fill='') +   #legend title
#   xlab('') + ylab('Average Expression') +    #
#   facet_wrap(~features.plot) + 
#   ggtitle("CWN: Scaled Average Expression; Fetal Pig Lung; Mesc Clusters")  +
#   theme(axis.text.x=element_blank())    #cluster names are long.
# 
# plot_grid(p_Collagen_Wnt_neural_FeaturePlot_PC_Exp + NoLegend(),
#           p_Collagen_Wnt_neural_FeaturePlot_AvgExp,  
#           nrow=1, ncol=2, 
#           labels=c("A", ""), 
#           rel_widths = c(1, 1.2))



################################################################################
# WNT markers 
#
################################################################################
df_Markers_Ref <- read.csv(here("../Markers","KEGG_pathways","KEGG_WNT_pathway_ssc04310.csv"), header = FALSE)
#df_Markers_Ref <- read.csv(here("../Markers","Wnt_signaling_shortlist.csv"), header = FALSE)

vec_LGR5_other <- c("LGR5","EPCAM")
vec_Common_Markers <- sort(df_Markers_Ref$V1)
vec_geneSymbols <- c(vec_LGR5_other,
                     vec_Common_Markers)

vec_geneSymbols <- toupper(unique(vec_geneSymbols))  

p1 <- DotPlot(obj_seurat, features=as.character(vec_geneSymbols),
              
              assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("KEGG_WNT_pathway_ssc04310; Pig Lung, Fetal P111; Piedrahita")
p1 <- theme_DotPlot(p1)
p1


#-------------------------------------------------------------------------------
# WNT markers, shortlist
#-------------------------------------------------------------------------------
#df_Markers_Ref <- read.csv(here("../Markers","KEGG_pathways","KEGG_WNT_pathway_ssc04310.csv"), header = FALSE)
df_Markers_Ref <- read.csv(here("../Markers","KEGG_pathways","KEGG_WNT_pathway_Shortlist_ssc04310.csv"), header = FALSE)

vec_LGR5_other <- c("LGR5","EPCAM")
vec_Common_Markers <- sort(df_Markers_Ref$V1)
vec_geneSymbols <- c(vec_LGR5_other,
                     vec_Common_Markers)

vec_geneSymbols <- toupper(unique(vec_geneSymbols))  

p1 <- DotPlot(obj_seurat_CombinedClusters, features=as.character(vec_geneSymbols),
              assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("KEGG_WNT_pathway_ssc04310; Pig Lung, Fetal P111; Piedrahita")
p1 <- theme_DotPlot(p1)
p1


#-------------------------------------------------------------------------------
# log2fold change between Developing airway and bud tip clusters
#-------------------------------------------------------------------------------
# fc_DE <- FoldChange(obj_seurat,
#                     ident.1 = "LGR5+ Fibroblast",  #"Developing Airway 2"
#                     ident.2 = "EPCAM+")
# 
# #Save the rownames (genes) to a column.
# fc_DE$gene <- rownames(fc_DE)
# 
# #Filter all_markers to just the vec_geneSymbols. Sort
# fc_DE_filtered <- fc_DE[fc_DE$gene %in% vec_geneSymbols,]
# fc_DE_sorted <- fc_DE_filtered[order(fc_DE_filtered$avg_log2FC, decreasing = TRUE), ]  
# 
# #Remove undefined genes
# fc_DE_sorted <- fc_DE_sorted[!grepl('ENSSSC',fc_DE_sorted$gene),] 
# 
# write.csv(fc_DE_sorted, here("Results","fc_WNT_cluster_LGR5_vs_EPCAM.csv"))



# ################################################################################
# # Differential gene expression - Collagen
# #
# ################################################################################
# df_Collagen_Markers_Ref <- read.csv(here("../Markers","Collagen_related_markers.csv"), header = FALSE)
# 
# vec_LGR5_other <- c("LGR5","VIM","EPCAM")
# vec_Collagen <- sort(df_Collagen_Markers_Ref$V1)
# vec_geneSymbols <- c(vec_Collagen)  
# 
# #Convert to uppercase gene symbols 
# vec_geneSymbols <- toupper(unique(vec_geneSymbols))  
# 
# pCollagen_1 <- DotPlot(obj_seurat, features=as.character(vec_geneSymbols),
#                        assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
#   ggtitle("Collagen; Fetal Pig Lung; Piedrahita")
# pCollagen_1 <- theme_DotPlot(pCollagen_1)
# pCollagen_1
# 
# # save_plot(filename=file.path(here("Results","Collagen_fetal_pig_lung_Pietrahita.pdf")),
# #           pCollagen_1, base_asp = 5)   
# 
# #-------------------------------------------------------------------------------
# # Gene count & percent by cluster
# #-------------------------------------------------------------------------------
# pCollagen_1 <- DotPlot(object = obj_seurat, features = as.character(vec_geneSymbols))
# head(pCollagen_1$data,3)
# #            avg.exp   pct.exp features.plot id avg.exp.scaled
# # COL1A1 17.74319066 99.987840        COL1A1  0      0.7923403
# # COL2A1  0.01264591  1.009241        COL2A1  0     -0.5014515
# # COL3A1 18.34982977 99.963521        COL3A1  0      0.7259203
# 
# #write.csv(pCollagen_1$data, here("Gene_percent_by_Collagen_Cluster.csv"))
# 
# #-------------------------------------------------------------------------------
# # FeaturePlot: pct.exp, avg.exp.scaled by cluster
# #-------------------------------------------------------------------------------
# # % Expression
# p_Collagen_FeaturePlot_PC_Exp <- ggplot(pCollagen_1$data,aes(x=id, y=pct.exp, fill=id)) +
#   geom_col() +
#   labs(fill='') +   #legend title
#   xlab('') + ylab('PE') +     #% Expression
#   facet_wrap(~features.plot) + 
#   ggtitle("Percent Expression; Fetal Pig Lung; Piedrahita")  +
#   theme(axis.text.x=element_blank())    #cluster names are long. Identify by colors.
# #p_Collagen_FeaturePlot_PC_Exp
# 
# # Average Expression
# p_Collagen_FeaturePlot_AvgExp<- ggplot(pCollagen_1$data,aes(x=id, y=avg.exp.scaled, fill=id)) +
#   geom_col() +
#   labs(fill='') +   #legend title
#   xlab('') + ylab('AE') +    #Average Expression
#   facet_wrap(~features.plot) + 
#   ggtitle("Scaled Average Expression; Fetal Pig Lung; Piedrahita")  +
#   theme(axis.text.x=element_blank())    #cluster names are long.
# #p_Collagen_FeaturePlot_AvgExp
# 
# plot_grid(p_Collagen_FeaturePlot_PC_Exp + NoLegend(),
#           p_Collagen_FeaturePlot_AvgExp,  
#           nrow=1, ncol=2, 
#           labels=c("B", ""), 
#           rel_widths = c(1, 1.4))
# 
# #-------------------------------------------------------------------------------
# # log2fold change between Developing airway and bud tip clusters
# #-------------------------------------------------------------------------------
# #FoldChange returns avg_logFC: natural log fold-change of the average expression between the two groups. 
# #Positive values indicate that the gene is more highly expressed in the first group
# fc_DE <- FoldChange(obj_seurat,
#                     ident.1 = "LGR5+ Fibroblast",  #"Developing Airway 2"
#                     ident.2 = "EPCAM+")
# #head(fc_DE,3)
# #                      avg_log2FC pct.1 pct.2               gene
# # ENSSSCG00000000002 -0.004238966 0.002 0.005 ENSSSCG00000000002
# # TTC38              -0.090828417 0.043 0.111              TTC38
# # CDPF1               0.064873773 0.169 0.128              CDPF1
# 
# #Save the rownames (genes) to a column.
# fc_DE$gene <- rownames(fc_DE)
# 
# #Filter all_markers to just the vec_geneSymbols. Sort
# fc_DE_filtered <- fc_DE[fc_DE$gene %in% vec_geneSymbols,]
# fc_DE_sorted <- fc_DE_filtered[order(fc_DE_filtered$avg_log2FC, decreasing = TRUE), ]  
# 
# #Remove undefined genes
# fc_DE_sorted <- fc_DE_sorted[!grepl('ENSSSC',fc_DE_sorted$gene),] 
# #view(fc_Collagen_sorted)
# 
# #write.csv(fc_DE_sorted, here("Results","fc_Collagen_cluster_LGR5_vs_EPCAM.csv"))
# 
# 
# 
# 
# ################################################################################
# # Differential gene expression - RSPONDIN
# #
# ################################################################################
# vec_LGR5_other <- c("LGR5","VIM","EPCAM", "ETV5","SFTPB","SFTPC","NGFR")
# vec_RSPO <- c("RSPO1","RSPO2","RSPO3","RSPO4")
# vec_geneSymbols <- c(vec_LGR5_other, vec_RSPO)
# 
# #Convert to uppercase gene symbols 
# vec_geneSymbols <- toupper(unique(vec_geneSymbols))  
# 
# p1 <- DotPlot(obj_seurat, features=as.character(vec_geneSymbols),
#                        assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
#   ggtitle("RSPO; Fetal Pig Lung; Piedrahita")
# p1 <- theme_DotPlot(p1)
# p1
# 
# #-------------------------------------------------------------------------------
# # Inspect the clusters 
# #-------------------------------------------------------------------------------
# p1 <- DimPlot(obj_seurat, reduction ="umap", pt.size = 0.5, 
#               label=TRUE, label.size=5, repel=TRUE) + NoLegend() + NoAxes()
# p2 <- FeaturePlot(obj_seurat,"VIM",order=T, cols = colors_FeaturePlot) + 
#   NoLegend() + NoAxes()
# p3 <- FeaturePlot(obj_seurat,"EPCAM",order=T, cols = colors_FeaturePlot) + 
#   ylab(NULL) + NoLegend() + NoAxes()
# p4 <- FeaturePlot(obj_seurat,"LGR5",order=T, cols = colors_FeaturePlot) + 
#   ylab(NULL) + NoAxes()
# 
# # p5 <- FeaturePlot(obj_seurat,"RSPO1",order=T, cols = colors_FeaturePlot) + 
# #   NoLegend() + NoAxes()
# p5 <- FeaturePlot(obj_seurat,"NGFR",order=T, cols = colors_FeaturePlot) + 
#    NoLegend() + NoAxes()
# p6 <- FeaturePlot(obj_seurat,"RSPO2",order=T, cols = colors_FeaturePlot) + 
#   NoLegend() + NoAxes()
# p7 <- FeaturePlot(obj_seurat,"RSPO3",order=T, cols = colors_FeaturePlot) + 
#   ylab(NULL) + NoLegend() + NoAxes()
# p8 <- FeaturePlot(obj_seurat,"ETV5",order=T, cols = colors_FeaturePlot) +
#   ylab(NULL) + NoLegend() + NoAxes()
# 
# p9 <- FeaturePlot(obj_seurat,"SFTPB",order=T, cols = colors_FeaturePlot) + 
#   NoLegend() + NoAxes()
# p10 <- FeaturePlot(obj_seurat,"SFTPC",order=T, cols = colors_FeaturePlot) + 
#   NoLegend() + NoAxes()
# p11 <- FeaturePlot(obj_seurat,"NGFR",order=T, cols = colors_FeaturePlot) + 
#   NoLegend() + NoAxes()
# p12 <- FeaturePlot(obj_seurat,"WNT5A",order=T, cols = colors_FeaturePlot) + 
#   NoLegend() + NoAxes()
# 
# plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,p11,p12,
#           nrow=3, ncol=4, 
#           rel_widths = c(1, 1, 1, 1.2, 1, 1, 1, 1.2, 1, 1, 1, 1.2))  
# 
# 
# 
# 

################################################################################
# Markers of fibroblasts and peripheral nerve associated fibroblasts
# 
#   “Single-cell analysis uncovers fibroblast heterogeneity and criteria for fibroblast and mural cell identification and discrimination” 
#     (2020) https://www.nature.com/articles/s41467-020-17740-1
################################################################################
#-------------------------------------------------------------------------------
# Fibroblast markers, general organ
#-------------------------------------------------------------------------------
# df_Fibro_Markers_Ref <- read.csv(here("../Markers","Fibroblast_markers_general_organs.csv"), header = FALSE)
# 
# vec_LGR5_other <- c("LGR5","EPCAM")  #"UCHL1" = "PGP9.5"
# vec_Fibro <- sort(df_Fibro_Markers_Ref$V1)
# vec_geneSymbols <- c(vec_LGR5_other,
#                      vec_Fibro)
# 
# vec_geneSymbols <- toupper(unique(vec_geneSymbols))  
# length(vec_geneSymbols)
# 
# p1 <- DotPlot(obj_seurat_CombinedClusters, features=as.character(vec_geneSymbols),
#               assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
#   ggtitle("Fibroblast Markers; Pig Lung, Fetal P111; Piedrahita")
# p1 <- theme_DotPlot(p1)
# p1
# #The following requested variables were not found: ABCA8A, CELF2, DPEP1, ENTPD2, HEG1, IGFBP6, LY6A, PCOLCE2, S100A16, SCARA5


#-------------------------------------------------------------------------------
# Fibroblast markers, lung
#   > Figure 
#-------------------------------------------------------------------------------
# df_Fibro_Markers_Ref <- read.csv(here("../Markers","Fibroblast_markers_lung.csv"), header = FALSE)
# 
# vec_LGR5_other <- c("LGR5","EPCAM")
# vec_Fibro <- sort(df_Fibro_Markers_Ref$V1)
# vec_geneSymbols <- c(vec_LGR5_other,
#                      vec_Fibro)
# 
# vec_geneSymbols <- toupper(unique(vec_geneSymbols))  
# length(vec_geneSymbols)
# 
# p1 <- DotPlot(obj_seurat_CombinedClusters, features=as.character(vec_geneSymbols),
#               #idents = c(vec_Epi_Clusters, vec_Mesc_Clusters, 13),
#               assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
#   ggtitle("Fibroblast Markers; Pig Lung, Fetal E111")
# p1 <- theme_DotPlot(p1)
# p1
# #The following requested variables were not found: FGF1, TRM1


