#################################################################################
# KEGG pathway DotPlots and Heatmaps of Transcript Expression.
#
#---------------------------------------------------------
# Ross Lampe
# Piedrahita Lab, College of Veterinary Medicine, North Carolina State University
# Last updated: 2022-09-05
#
#################################################################################
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(here)

library(edgeR)

source("../themes.R")

#Read the complete set of markers for all clusters.
all_markers <- read.csv(here("../Markers","Fetal_E80_all_markers.csv")) 
colnames(all_markers)[1] <- "Gene"

#Remove undefined genes
res_signif_AvsB <- all_markers[!grepl('ENSSSC',all_markers$Gene),]  #
nrow(res_signif_AvsB)
head(res_signif_AvsB,2)
#                 Gene avg_log2FC pct.1 pct.2    p_val_adj cluster   gene
# DCN     0.000000e+00  1.5847308 0.979 0.638  0.00000e+00       0    DCN
# CRABP1 2.974545e-274  0.9959538 0.994 0.803 4.71614e-270       0 CRABP1



################################################################################
# Heatmap, KEGG_WNT_pathway_ssc04310 
#     (LGR5+,VIM+) vs EPCAM+ clusters
################################################################################
df_Markers_Ref <- read.csv(here("../Markers","KEGG_pathways","KEGG_WNT_pathway_ssc04310.csv"), header = FALSE)

vec_LGR5_other <- c("LGR5","VIM","EPCAM")
vec_Common_Markers <- df_Markers_Ref$V1
vec_geneSymbols <- c(vec_Common_Markers)
vec_geneSymbols <- sort(vec_geneSymbols)

#Convert to uppercase gene symbols
vec_geneSymbols <- toupper(unique(vec_geneSymbols))  
length(vec_geneSymbols)

#-------------------------------------------------------------------------------
# Dotplot
#-------------------------------------------------------------------------------
p1 <- DotPlot(obj_seurat, features=as.character(vec_geneSymbols),
              assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("KEGG_WNT_pathway_ssc04310; Fetal Pig")
p1 <- theme_DotPlot(p1)
p1

#-------------------------------------------------------------------------------
# Heatmap, Reduce the assay to only the marker genes
#-------------------------------------------------------------------------------
# #Filter all_markers to just the vec_geneSymbols.
# sigMat <- res_signif_AvsB[res_signif_AvsB$Gene %in% vec_geneSymbols,]
# nrow(sigMat)
# #view(sigMat)
# 
# 
# #Match to marker list, group the markers by cluster and rank by avg_log2FC  
# df_topGenes <- sigMat %>% group_by(cluster) %>% top_n(100, avg_log2FC)
# 
# #Remove undefined genes
# df_topGenes <- df_topGenes[!grepl('ENSSSC',df_topGenes$gene),] 
# head(df_topGenes,2)
# #   Gene        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene   
# #   <chr>       <dbl>      <dbl> <dbl> <dbl>     <dbl>   <int> <chr>  
# # 1 VIM     4.75e-140     0.302  1     0.961 8.90e-136       0 VIM    
# # 2 COL14A1 6.88e- 83     0.207  0.914 0.899 1.29e- 78       0 COL14A1
# 
# pHeatMap <- DoHeatmap(obj_seurat, as.character(df_topGenes$gene, 
#                                                               group.bar = TRUE, 
#                                                               label=FALSE,
#                                                               slim.col.label = TRUE,
#                                                               remove.key = TRUE,
#                                                               obj_seurat)) +
#                     labs(title = "KEGG_WNT_pathway_ssc04310")
# pHeatMap <- theme_DoHeatmap(pHeatMap)
# pHeatMap



################################################################################
# KEGG_TGFbeta_pathway_ssc04350
#
################################################################################
df_Markers_Ref <- read.csv(here("../Markers","KEGG_pathways","KEGG_TGFbeta_pathway_ssc04350.csv"), header = FALSE)

vec_LGR5_other <- c("LGR5","VIM","EPCAM")
vec_Common_Markers <- df_Markers_Ref$V1
vec_geneSymbols <- c(vec_Common_Markers)
vec_geneSymbols <- sort(vec_geneSymbols)

#Convert to uppercase gene symbols
vec_geneSymbols <- toupper(unique(vec_geneSymbols))  
length(vec_geneSymbols)

#-------------------------------------------------------------------------------
# Dotplot
#-------------------------------------------------------------------------------
p1 <- DotPlot(obj_seurat_mesc, features=as.character(vec_geneSymbols),
                       assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("KEGG_TGFbeta_pathway_ssc04350; Fetal Pig")
p1 <- theme_DotPlot(p1)
p1

#-------------------------------------------------------------------------------
# Heatmap, Reduce the assay to only the marker genes
#-------------------------------------------------------------------------------
# #Filter all_markers to just the vec_geneSymbols.
# sigMat <- res_signif_AvsB[res_signif_AvsB$Gene %in% vec_geneSymbols,]
# nrow(sigMat)
# #view(sigMat)
# 
# #Match to marker list, group the markers by cluster and rank by avg_log2FC  
# df_topGenes <- sigMat %>% group_by(cluster) %>% top_n(100, avg_log2FC)
# nrow(df_topGenes)
# 
# #Remove undefined genes
# df_topGenes <- df_topGenes[!grepl('ENSSSC',df_topGenes$gene),] 
# head(df_topGenes,2)
# #   Gene        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene   
# #   <chr>       <dbl>      <dbl> <dbl> <dbl>     <dbl>   <int> <chr>  
# # 1 VIM     4.75e-140     0.302  1     0.961 8.90e-136       0 VIM    
# # 2 COL14A1 6.88e- 83     0.207  0.914 0.899 1.29e- 78       0 COL14A1
# 
# #write.csv(df_topGenes, file=here("Results", "FetalP80_KEGG TGF-β pathway, ssc04310.csv"), row.names = FALSE)
# 
# 
# pHeatMap <- DoHeatmap(obj_seurat, as.character(df_topGenes$gene, 
#                                                group.bar = TRUE, 
#                                                label=FALSE,
#                                                slim.col.label = TRUE,
#                                                remove.key = TRUE,
#                                                obj_seurat)) +
#   labs(title = "KEGG_TGFbeta_pathyway_ssc04350")
# pHeatMap <- theme_DoHeatmap(pHeatMap)
# pHeatMap



