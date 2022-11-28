################################################################################
# Create a Seurat object from GD80 fetal pig data.
#   The cells were sorted to collect only GFP+ LGR5 cells for scRNAseq.
#   5284 cells were collected.
#
#---------------------------------------------------------
# Ross Lampe
# Piedrahita Lab, College of Veterinary Medicine, North Carolina State University
# Last updated: 2022-09-12
#
################################################################################
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(readr)
library(here)

source("../themes.R")


################################################################################
# Create a Seurat object
#
################################################################################
sl1 <- readRDS(file = here("Seurat_Data","Piedrahita_pig_10X_alevin_quants_mat_sl1.rds"))
sl1.annot = sl1$counts
colnames(sl1.annot) = paste0("Piedrahita_10X_",colnames(sl1.annot))


#-------------------------------------------------------------------------------
#Create Seurat object and remove low quality cells.
# Perform initial filtering in order to exclude genes that are expressed in fewer than min.cells, 
# and to exclude cells that contain fewer than min.features expressed genes
#
#-------------------------------------------------------------------------------
# Create and inspect the Seurat object
#-------------------------------------------------------------------------------
obj_seurat = CreateSeuratObject(counts = sl1.annot, 
                                min.cells=100, 
                                min.features=100, 
                                project = "Piedrahita_10X") 

# Check raw dimensions. Returns: genes, cells
dim(obj_seurat)
#[1] 15855 72851

head(obj_seurat,2)
#                                 orig.ident nCount_RNA nFeature_RNA
# Piedrahita_10X_TGAACGTGTTTACGAC Piedrahita    26667.0         6342
# Piedrahita_10X_TCAATTCTCTACGCAA Piedrahita    24388.5         5626



################################################################################
# QC - Filter out low quality reads using selected thresholds
#   
################################################################################
obj_seurat<- subset(obj_seurat, subset = nCount_RNA > 2000 & nFeature_RNA > 1000)
# dim(obj_seurat)


################################################################################
# Seurat standard workflow
#   
################################################################################
#---------------------------------------------------------
# Normalize and scale data
#---------------------------------------------------------
obj_seurat <- NormalizeData(obj_seurat)
obj_seurat <- FindVariableFeatures(obj_seurat, selection.method = "vst", nfeatures = 2000)
obj_seurat <- ScaleData(obj_seurat)


#---------------------------------------------------------
# Determine the dimensionality of the dataset
# by finding the "elbow" of an elbow plot.
#---------------------------------------------------------
obj_seurat<- RunPCA(obj_seurat, verbose = FALSE, npcs = 50)

# Run an ElbowPlot to assess the number of principle components (PCs) needed.
ElbowPlot(obj_seurat)

# Perform dimensionality reduction
obj_seurat <- RunUMAP(obj_seurat, dims = 1:50, verbose = FALSE)


#-------------------------------------------------------------------------------
# Visually explore the similarity in gene expression between the clusters 
# by plotting the clusters in PCA space.
#-------------------------------------------------------------------------------
DimPlot(obj_seurat, reduction = "pca", label=TRUE, label.size=6)


#-------------------------------------------------------------------------------
# Cluster the cells
#   https://hbctraining.github.io/scRNA-seq/lessons/07_SC_clustering_cells_SCT.html
#-------------------------------------------------------------------------------
obj_seurat <- FindNeighbors(obj_seurat, dims = 1:50, verbose = FALSE)
obj_seurat <- FindClusters(obj_seurat, verbose = FALSE, resolution = 0.75, algorithm=2)

#levels(obj_seurat)   # Enumerate the clusters



################################################################################
# Inspect the Seurat object
#   
################################################################################
#Inspect the genes
# nrow(x = obj_seurat)  #genes
# rownames(obj_seurat@assays$RNA@counts)

#-------------------------------------------------------------------------------
# Inspect the clusters 
#-------------------------------------------------------------------------------
p1 <- DimPlot(obj_seurat, reduction ="umap", pt.size = 0.5, 
              label=TRUE, label.size=5, repel=TRUE) + NoLegend() + NoAxes()
p2 <- FeaturePlot(obj_seurat,"VIM",order=T, cols = colors_FeaturePlot) + 
  NoLegend() + NoAxes()
p3 <- FeaturePlot(obj_seurat,"EPCAM",order=T, cols = colors_FeaturePlot) + 
  ylab(NULL) + NoLegend() + NoAxes()
p4 <- FeaturePlot(obj_seurat,"LGR5",order=T, cols = colors_FeaturePlot) + 
  ylab(NULL) + NoAxes()

plot_grid(p1, p4,
          nrow=1, ncol=2, 
          rel_widths = c(1, 1.2)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

plot_grid(p1, p2, p3, p4,
          nrow=1, ncol=4, 
         #labels=c("A", "B", "C", "D"),  #,labels="AUTO"
          rel_widths = c(1, 1, 1, 1.2)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


#-------------------------------------------------------------------------------
# Explore the similarity in gene expression between the clusters by plotting the clusters in PCA space.
#-------------------------------------------------------------------------------
DimPlot(obj_seurat, reduction = "pca", label=TRUE, label.size=6)


#-------------------------------------------------------------------------------
# Save the Seurat object as an RDS file.
#-------------------------------------------------------------------------------
#saveRDS(obj_seurat, file = here("Seurat_Data","Seurat_14clusters.rds"))

#save.image("Seurat_Data","Seurat_14clusters.Rdata")



################################################################################
# Define the subsets clusters: epithelial, mesenchymal
#   
################################################################################
vec_Epi_Clusters <- c(7,8,10,13)
vec_Mesc_Clusters <- c(0:6,9,11,12)


################################################################################
# Combine Epi and Mesc clusters
#   
################################################################################
obj_seurat_CombinedClusters <- obj_seurat
levels(obj_seurat_CombinedClusters)   # Enumerate the clusters
# [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13"

s_Mes <- "Mesenchymal"
s_Epi <- "Epithelial"
new.cluster.ids <- c(s_Mes, s_Mes, s_Mes, s_Mes, s_Mes, s_Mes, s_Mes, 
                     s_Epi,  s_Epi,  s_Mes,  s_Epi, s_Mes, s_Mes,  s_Epi)

names(new.cluster.ids) <- levels(obj_seurat_CombinedClusters)
obj_seurat_CombinedClusters <- RenameIdents(obj_seurat_CombinedClusters, new.cluster.ids)


p1 <- DimPlot(obj_seurat, reduction ="umap", pt.size = 0.5, 
              label=TRUE, label.size=5, repel=TRUE) + NoLegend() + NoAxes()
p2 <- DimPlot(obj_seurat_CombinedClusters, reduction = "umap", pt.size = 0.5, 
              label=TRUE, label.size=5, repel=TRUE) + NoLegend() + NoAxes()
plot_grid(p1, p2,
          nrow=1, ncol=2,
          rel_widths = c(1, 1)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


