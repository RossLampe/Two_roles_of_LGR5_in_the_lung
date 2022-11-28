#################################################################################
# Find Markers
#
#---------------------------------------------------------
# Ross Lampe
# Piedrahita Lab, College of Veterinary Medicine, North Carolina State University
# Last updated: 2022-09-06
#
#################################################################################
library(Seurat)
library(presto)
library(dplyr)
library(ggplot2)
library(cowplot)
library(readr)
library(here)

source("../themes.R")


################################################################################
# Find markers (differentially expressed genes) for all clusters
# 
################################################################################
obj_seurat.AllMarkers <- FindAllMarkers(obj_seurat, 
                                    # only.pos = T, 
                                     min.pct = 0.5, 
                                     logfc.threshold = 0.05)

all_markers <- obj_seurat.AllMarkers

# Drop ENSSSCG names from list
all_markers <- all_markers[!grepl('ENSSSC',all_markers$gene),] 

#Drop ribosomal genes
all_markers <- all_markers[-grep("^RP[LS]", all_markers$gene),]

all_markers <- all_markers %>%
              group_by(cluster)

#top300Genes <- as.data.frame(all_markers %>% group_by(cluster) %>% top_n(n = 300, wt = avg_log2FC))

all_markers <- as.data.frame(all_markers)
head(all_markers,2)
#           p_val avg_log2FC pct.1 pct.2     p_val_adj cluster  gene
# 1 3.844175e-214  0.8678391  0.99 0.905 6.094940e-210       0   MGP
# 2 9.114374e-169  0.3743659  1.00 0.987 1.445084e-164       0 RACK1

write.csv(all_markers, here("Results","scRNAseq_GD80_AllMarkers.csv"))

#-------------------------------------------------------------------------------
# Dotplot markers, all clusters
#-------------------------------------------------------------------------------
# p1 <- DotPlot(obj_seurat, features=all_markers,
#               assay="RNA", cols = colors_DotPlot) + RotatedAxis() + 
#   ggtitle("GD80; Top Markers")
# p1 <- theme_DotPlot(p1)
# p1




################################################################################
# Create a refined set of markers, all clusters: Presto
################################################################################
# Use Presto to identify top marker genes of each cluster, quickly
vargenes<- presto::wilcoxauc(obj_seurat, 'seurat_clusters', seurat_assay = 'RNA')
head(vargenes)

vargenes <- presto::wilcoxauc(obj_seurat,
                              #groups_use = c("Mesenchymal", "Epithelial"),
                              seurat_assay = 'RNA')
head(vargenes)

# Drop ENSSSCG names from list
vargenes_filtered <- vargenes[!grepl('ENSSSC',vargenes$feature),] 

#colnames(vargenes_filtered)[1] <- "Gene"
head(vargenes_filtered,1)

#-------------------------------------------------------------------------------
# Inspect the top markers
#-------------------------------------------------------------------------------
top_vargenes = top_markers(vargenes_filtered, n = 12, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
head(top_vargenes, 2)
nrow(top_vargenes)

# write.table(top_vargenes,"Piedrahita_10X_alevin_ensembl_seurat_clustered_prestoMarkers.txt",row.names=F,quote=F,sep="\t")
write.csv(top_vargenes, here("Results","scRNAseq_GD80_All_Clusters_Labeled_Presto_Markers.csv"))

all_markers_top_vargenes <- top_vargenes %>%
  select(-rank) %>%
  unclass() %>%
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

length(all_markers_top_vargenes)

p1 <- DotPlot(obj_seurat, features=rev(all_markers_top_vargenes),
              assay="RNA",cols = colors_DotPlot) + RotatedAxis() +
  ggtitle("GD80; Top Markers")
p1 <- theme_DotPlot(p1)
p1



################################################################################
################################################################################
# Find ALL markers between (LGR5+) and (EPCAM+) clusters
# 
################################################################################
################################################################################
# group_Mesc <- c(0:6,9,11,12)
# group_Epi <- c(7,8,10,13)
# 
# obj_seurat.mesc_vs_epi_markers = FindMarkers(obj_seurat, 
#                                              ident.1 = group_Mesc, 
#                                              ident.2=group_Epi, 
#                                              assay="RNA", 
#                                              slot="scale.data")

obj_seurat_CombinedClusters.AllMarkers = FindAllMarkers(obj_seurat_CombinedClusters, 
                                                              # ident.1 = "Mesenchymal", 
                                                              # ident.2= "Epithelial", 
                                                              assay="RNA", 
                                                              slot="scale.data")

#Drop ENSSSCG names from list
combined_AllMarkers_markers <- obj_seurat_CombinedClusters.AllMarkers
combined_AllMarkers_markers <- combined_AllMarkers_markers[!grepl('ENSSSC',rownames(combined_AllMarkers_markers)),] 

#Drop ribosomal genes
combined_AllMarkers_markers <- combined_AllMarkers_markers[-grep("^RP[LS]", rownames(combined_AllMarkers_markers)),]
head(arrange(combined_AllMarkers_markers,desc(avg_diff)), 20)

write.csv(combined_AllMarkers_markers, here("Results","scRNAseq_GD80_combined_AllMarkers_markers.csv"))


#-------------------------------------------------------------------------------
# Inspect using a DotPlot
#-------------------------------------------------------------------------------
p_combined_AllMarkers_markers<- rownames(head(arrange(combined_AllMarkers_markers, desc(avg_diff)), n = 100))

p1 <- DotPlot(obj_seurat_CombinedClusters, 
              features=as.character(p_combined_AllMarkers_markers),
              assay="RNA", 
              cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("GD80; Mesc vs Epi Clusters, All Markers")
p1 <- theme_DotPlot(p1)
p1



################################################################################
################################################################################
# Find markers between (LGR5+) and (EPCAM+) clusters
# 
################################################################################
################################################################################
# group_Mesc <- c(0:6,9,11,12)
# group_Epi <- c(7,8,10,13)
# 
# obj_seurat.mesc_vs_epi_markers = FindMarkers(obj_seurat, 
#                                              ident.1 = group_Mesc, 
#                                              ident.2=group_Epi, 
#                                              assay="RNA", 
#                                              slot="scale.data")

obj_seurat_CombinedClusters.mesc_vs_epi_markers = FindMarkers(obj_seurat_CombinedClusters, 
                                             ident.1 = "Mesenchymal", 
                                             ident.2= "Epithelial", 
                                             assay="RNA", 
                                             slot="scale.data")

#Drop ENSSSCG names from list
mesc_vs_epi_markers <- obj_seurat_CombinedClusters.mesc_vs_epi_markers
mesc_vs_epi_markers <- mesc_vs_epi_markers[!grepl('ENSSSC',rownames(mesc_vs_epi_markers)),] 

#Drop ribosomal genes
mesc_vs_epi_markers <- mesc_vs_epi_markers[-grep("^RP[LS]", rownames(mesc_vs_epi_markers)),]
head(arrange(mesc_vs_epi_markers,desc(avg_diff)), 5)

#write.csv(mesc_vs_epi_markers, here("../Markers","scRNAseq_GD80_mesc_vs_epi_markers.csv"))
write.csv(mesc_vs_epi_markers, here("Results","scRNAseq_GD80_Markers_Mesc_vs_Epi.csv"))


#-------------------------------------------------------------------------------
# Inspect using a DotPlot
#-------------------------------------------------------------------------------
p_mesc_vs_epi_markers <- rownames(head(arrange(mesc_vs_epi_markers, desc(avg_diff)), n = 100))

p1 <- DotPlot(obj_seurat_CombinedClusters, 
              features=as.character(p_mesc_vs_epi_markers),
              assay="RNA", 
              cols = colors_DotPlot) + RotatedAxis() + 
  ggtitle("GD80; Mesc vs Epi Clusters Top Markers")
p1 <- theme_DotPlot(p1)
p1



################################################################################
# Create a refined set of markers: FindMarkers(obj_seurat_CombinedClusters,...)
#
################################################################################
#Sort by avg_diff
top_mesc_vs_epi_markers <- head(arrange(mesc_vs_epi_markers, desc(avg_diff)), n = 100)
head(top_mesc_vs_epi_markers, 2)
#        p_val avg_diff pct.1 pct.2 p_val_adj
# COL1A2     0 2.586972 0.821 0.024         0
# COL3A1     0 2.526602 0.767 0.021         0

write.csv(top_mesc_vs_epi_markers, here("Results","scRNAseq_GD80_Markers_Mesc_vs_Epi.csv"))

#Inspect
vec_top_mesc_vs_epi_markers <- rownames(top_mesc_vs_epi_markers)

p1 <- DotPlot(obj_seurat_CombinedClusters, features=vec_top_mesc_vs_epi_markers,
        assay="RNA",cols = colors_DotPlot) + RotatedAxis() +
  ggtitle("GD80; Mesc vs Epi Clusters Top Markers")
p1 <- theme_DotPlot(p1)
p1



################################################################################
# Create a refined set of markers, Epi vs Mesc clusters: Presto
################################################################################
# Use Presto to identify top marker genes of each cluster, quickly
vargenesTest<- presto::wilcoxauc(obj_seurat_CombinedClusters, 'seurat_clusters', seurat_assay = 'RNA')
head(vargenesTest)
#              feature group     avgExpr        logFC statistic       auc         pval         padj     pct_in   pct_out
# 1 ENSSSCG00000000002     0 0.001773825 -0.022891089   8688470 0.4878856 4.054623e-12 3.313714e-11  0.1875293  2.611091
# 2              TTC38     0 0.040773833 -0.014787509   8705840 0.4888610 9.548368e-05 2.296214e-04  4.0318800  6.348066

vargenes <- presto::wilcoxauc(obj_seurat_CombinedClusters,
                             groups_use = c("Mesenchymal", "Epithelial"),
                             seurat_assay = 'RNA')
head(vargenes)
#              feature       group    avgExpr       logFC statistic       auc         pval         padj    pct_in    pct_out
# 1 ENSSSCG00000000002 Mesenchymal 0.02205134  0.01606539   6187504 0.5070298 8.640020e-04 1.448530e-03  2.295584  0.8995502
# 2              TTC38 Mesenchymal 0.04226746 -0.08081312   5386413 0.4413851 8.748790e-65 8.331055e-64  4.339746 16.4167916


# Drop ENSSSCG names from list
vargenes_filtered <- vargenes[!grepl('ENSSSC',vargenes$feature),] 

#colnames(vargenes_filtered)[1] <- "Gene"
head(vargenes_filtered,1)

#-------------------------------------------------------------------------------
# Inspect the top markers
#-------------------------------------------------------------------------------
top_vargenes = top_markers(vargenes_filtered, n = 50, auc_min = 0.5, pct_in_min = 20, pct_out_max = 20)
head(top_vargenes, 2)
nrow(top_vargenes)

# write.table(top_vargenes,"Piedrahita_10X_alevin_ensembl_seurat_clustered_prestoMarkers.txt",row.names=F,quote=F,sep="\t")
write.csv(top_vargenes, here("Results","scRNAseq_GD80_Presto_Markers_Mesc_vs_Epi.csv"))

all_markers_top_vargenes <- top_vargenes %>%
  select(-rank) %>%
  unclass() %>%
  stack() %>%
  pull(values) %>%
  unique() %>%
  .[!is.na(.)]

length(all_markers_top_vargenes)

p1 <- DotPlot(obj_seurat_CombinedClusters, features=rev(all_markers_top_vargenes),
        assay="RNA",cols = colors_DotPlot) + RotatedAxis() +
  ggtitle("GD80; Mesc vs Epi Clusters Top Markers")
p1 <- theme_DotPlot(p1)
p1


