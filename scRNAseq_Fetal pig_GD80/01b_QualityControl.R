#################################################################################
# Quality Control Assessment
#   to prepare for cell-level filtering
#   https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
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
library(cowplot)
library(readr)
library(here)

source("../themes.R")


################################################################################
#
################################################################################
# Create metadata dataframe
metadata <- obj_seurat@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)



################################################################################
# QC: mitochondrial genes
# 
################################################################################
# grep("^MT-",rownames(obj_seurat@assays$RNA@counts),value = TRUE)
# 
# obj_seurat$percent.mito <- PercentageFeatureSet(obj_seurat,pattern="^MT-[LS]|mt-[ls]")
# head(obj_seurat$percent.mito)
# 
# metadata <- obj_seurat@meta.data
# 
# # Make diagnostic plots to help with filtering.
# p <- ggplot(obj_seurat@meta.data, aes(nCount_RNA, nFeature_RNA, color=percent.mito))
# p <- p + geom_point(size=0.2)
# p 
# 
# p <- metadata %>%
#   ggplot(aes(percent.mito)) + 
#   geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
#   ggtitle("Distribution of Percentage Mito") +
#   geom_vline(xintercept = 10)
# p 
# 
# # QC: Basic filtering - remove low quality reads - mito percentages.
# selected_mito <- WhichCells(obj_seurat, expression = percent.mito < 0.2)
# obj_seurat <- subset(obj_seurat, cells = selected_mito)



################################################################################
# QC: View the ribosomal genes
# The ribosomal proteins are highly expressed they will make up a larger proportion of the transcriptional landscape 
# when fewer of the lowly expressed genes are detected. And we can plot the different QC-measures as scatter plots.
################################################################################
grep("^RP[LS]",rownames(obj_seurat@assays$RNA@counts),value = TRUE)
# # [1] "RPS18"          "RPL10A"         "RPS6KL1"        "RPS6KA5"        "RPS11"          "RPL22"          "RPS6KA1"        "RPS8"
# # [9] "RPS5"           "RPS6KA2"        "RPS12"          "RPL4"           "RPS29"          "RPL36AL"        "RPS6"           "RPL5"
# # [17] "RPS3A"          "RPL34"          "RPSA"           "RPL24"          "RPS6KA3"        "RPS4X"          "RPS6KA6"        "RPLP2"
# # [25] "RPS6KB2"        "RPS13"          "RPS28"          "RPS23"          "RPS3"           "RPS25"          "RPS6KC1"        "RPL26L1"
# # [33] "RPL11"          "RPL26"          "RPS16"          "RPL27"          "RPS6KA4"        "RPS21"          "RPL7L1"         "RPL29"
# # [41] "RPL17-C18orf32" "RPS27L"         "RPS9"           "RPS17"          "RPL8"           "RPL28"          "RPL10"          "RPL35"
# # [49] "RPLP0"          "RPS20"          "RPS15"          "RPL15"          "RPL32"          "RPL3L"          "RPL38"          "RPL37"
# # [57] "RPL36A-HNRNPH2" "RPL6"           "RPL23"          "RPL7"           "RPS15A"         "RPL22L1"        "RPL3"

obj_seurat$percent.ribo <- PercentageFeatureSet(obj_seurat,pattern="^RP[LS]|rp[ls]")
head(obj_seurat$percent.ribo)
# Lampe_10X_TGAACGTGTTTACGAC Lampe_10X_TCAATTCTCTACGCAA Lampe_10X_CCTCTCCCAGTTAAAG Lampe_10X_CCACTTGGTTTCGCTC Lampe_10X_GTAATCGTCTCAGTCC
# 16.22600                   17.25823                   14.55269                   12.57743                   15.90369

metadata <- obj_seurat@meta.data

# Make diagnostic plots to help with filtering.
p <- ggplot(metadata, aes(nCount_RNA, nFeature_RNA, color=percent.ribo))
p <- p + geom_point(size=0.2)
p

p <- metadata %>%
  ggplot(aes(percent.ribo)) +
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Ribosomal") +
  geom_vline(xintercept = 10)
p

# QC: Basic filtering - remove low quality reads - ribo percentages.
# selected_ribo <- WhichCells(obj_seurat, expression = percent.ribo > 0.05)
# obj_seurat <- subset(obj_seurat, cells = selected_ribo)


#-------------------------------------------------------------------------------
# QC: Visualize the nCount_RNA feature-feature relationships.
#
# Set thresholds and filter out low quality reads.
#-------------------------------------------------------------------------------
p1 <- FeatureScatter(obj_seurat, feature1 = "nCount_RNA", feature2 = "percent.ribo")
p2 <- FeatureScatter(obj_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p1 + p2

#Define some QC parameters from inspection of this scatterplot
min_features_RNA = 1000 # Filter out cells with fewer genes to remove dead cells
min_Count_RNA = 2000
#max_Count_RNA = 18000  # Filter out cells with more UMIs to catch a few remaining doublets



################################################################################
# QC: Filter out low quality reads using selected thresholds
################################################################################
# obj_seurat <- subset(obj_seurat,
#                      nFeature_RNA > min_features_RNA &
#                      nCount_RNA > min_Count_RNA)
#  & nCount_RNA < max_Count_RNA)



################################################################################
# Clean up
################################################################################
rm(metadata)


