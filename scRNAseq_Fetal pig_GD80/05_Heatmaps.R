#################################################################################
# Analyze E80 fetal pig data from Seurat object.
#   Heatmaps of transcript expression.
#
#---------------------------------------------------------
# Ross Lampe
# Piedrahita Lab, College of Veterinary Medicine, North Carolina State University
# Last updated: 2022-08-26
#
#################################################################################
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(here)

source("../themes.R")


################################################################################
# Heatmap of top DE expressed genes
#
################################################################################
#Group the markers by cluster and rank by avg_log2FC   
df_topGenes <- all_markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
#head(df_topGenes,2)

#Remove undefined genes
#df_topGenes <- df_topGenes[!grepl('ENSSSC',df_topGenes$gene),] 

#-------------------------------------------------------------------------------
# Heatmap of top DE expressed genes
#-------------------------------------------------------------------------------
pHeatMap <- DoHeatmap(obj_seurat, as.character(df_topGenes$gene, 
                                   group.bar = TRUE, 
                                   label=FALSE,
                                   slim.col.label = TRUE,
                                   remove.key = TRUE)) +
                    labs(title = "Top Differentially Expressed Transcripts") 
pHeatMap <- theme_DoHeatmap (pHeatMap)
pHeatMap
    



################################################################################
# Heatmap, WNT markers
#
################################################################################
#df_WNT_Markers_Ref <- read.csv(here("../Markers","Wnt_signaling_shortlist.csv"), header = FALSE)
df_WNT_Markers_Ref <- read.csv(here("../Markers","KEGG_pathways","KEGG_WNT_pathway_Shortlist_ssc04310.csv"), header = FALSE)

vec_LGR5_other <- c("LGR5","EPCAM")
vec_Common_Markers <- df_WNT_Markers_Ref$V1
vec_geneSymbols <- c(vec_LGR5_other,
                     vec_Common_Markers)

#Convert uppercase gene symbols to mouse symbols.
vec_geneSymbols <- toupper(unique(vec_geneSymbols))  
length(vec_geneSymbols)

#Filter all_markers to just the vec_geneSymbols.
sigMat <- all_markers[all_markers$gene %in% vec_geneSymbols,]
#nrow(sigMat)

#Match to marker list, group the markers by cluster and rank by avg_log2FC  
df_topGenes <- sigMat %>% group_by(cluster) %>% top_n(100, avg_log2FC)
#head(df_topGenes)

#Remove undefined genes
df_topGenes <- df_topGenes[!grepl('ENSSSC',df_topGenes$gene),] 

#Heatmap 
pHeatMap <- DoHeatmap(obj_seurat_CombinedClusters, as.character(df_topGenes$gene, 
                                               group.bar = TRUE, 
                                               label=FALSE,
                                               slim.col.label = TRUE,
                                               remove.key = TRUE)) +
                    labs(title = "WNT Expression")
pHeatMap <- theme_DoHeatmap(pHeatMap)
pHeatMap







