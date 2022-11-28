################################################################################
# Heatmap differential gene expression from bulk RNAseq.
#
#-------------------------------------------------------------------------------
# Ross Lampe
# Piedrahita Lab, College of Veterinary Medicine, North Carolina State University
# Last updated: 2022-08-28
# 
################################################################################
library(dplyr)
library(cowplot)
library(EnhancedVolcano)
library(pheatmap)
library(here)

source("../themes.R")


################################################################################
# Filter and shape the rlogNorm data
#
################################################################################
res_AvsB <- read.csv(here("bulkRNA_D111 Adult Pigs", "rlogNorm__D111_AdultPig_GFPpos_vs_GFPneg_Annotated11.1.105.csv"), header = TRUE)

df_res_AvsB_shaped <- res_AvsB


#-------------------------------------------------------------------------------
# Remove low expressing genes
#-------------------------------------------------------------------------------
n_minTotalCount = 10
df_res_AvsB_shaped <- subset(df_res_AvsB_shaped,
                             ((df_res_AvsB_shaped$GFPneg_Rep1_Cnt + df_res_AvsB_shaped$GFPneg_Rep2_Cnt) >=n_minTotalCount) &
                               ((df_res_AvsB_shaped$GFPpos_Rep1_Cnt + df_res_AvsB_shaped$GFPpos_Rep2_Cnt) >=n_minTotalCount))


#-------------------------------------------------------------------------------
# 
#-------------------------------------------------------------------------------
#Set the row names.
rownames(df_res_AvsB_shaped) <- df_res_AvsB_shaped[,"Gene"]
#head(res_AvsB,2)

df_res_AvsB_shaped <- df_res_AvsB_shaped

#Sort summary list by p-value and log2FoldChange
df_res_AvsB_shaped <- df_res_AvsB_shaped %>% arrange(desc(log2FoldChange),desc(padj))  #Sorting a Dataframe in R Dplyr
head(df_res_AvsB_shaped,2)
#                                  Gene    baseMean log2FoldChange       pvalue         padj GFPneg_Rep1 GFPpos_Rep1 GFPneg_Rep2 GFPpos_Rep2 GFPneg_Rep1.1 GFPpos_Rep1.1 GFPneg_Rep2.1 GFPpos_Rep2.1
# SLC30A8                       SLC30A8    41.08775       8.065317 4.560000e-06 1.017490e-04     0.00000    68.44416    0.000000    95.90683             0           119             0           152
# LGR5                             LGR5   954.84920       7.039042 5.520000e-29 6.670000e-25    19.67646  1850.86808    9.893774  1938.95851            10          3218             7          3073

#-------------------------------------------------------------------------------
# Filter by padj and log2FoldChange to find significant genes
#-------------------------------------------------------------------------------
padj.cutoff <- 0.05 # False Discovery Rate (FDR) cutoff

#Filter by padj.cutoff & abs(log2FoldChange) to report top DE genes
# res_signif_AvsB <- df_res_AvsB_shaped %>%
#   mutate(threshold_OE = (padj < padj.cutoff) & ((log2FoldChange) > 3.85)|(log2FoldChange < -5.5))
res_signif_AvsB <- df_res_AvsB_shaped %>%
  mutate(threshold_OE =  (padj < padj.cutoff) & ((log2FoldChange) > 1.5)|(log2FoldChange < -3))


#-------------------------------------------------------------------------------
#Clean up df
#-------------------------------------------------------------------------------
res_signif_AvsB <- res_signif_AvsB %>% tidyr::drop_na()                    #Use tidyr function to remove all rows in which NA appears

res_signif_AvsB <- res_signif_AvsB[res_signif_AvsB$threshold_OE == TRUE,]  #remove rows that fail the threshold tests
res_signif_AvsB <- res_signif_AvsB[,!(names(res_signif_AvsB) %in% c("threshold_OE"))]  #Remove the threshold_OE column

#Remove undefined genes
res_signif_AvsB <- res_signif_AvsB[!grepl('ENSSSC',res_signif_AvsB$Gene),]  #
nrow(res_signif_AvsB)
head(res_signif_AvsB,2)
#        Gene    baseMean log2FoldChange   pvalue     padj Pig1_GFP_Neg Pig1_GFP_Pos Pig2_GFP_Neg Pig2_GFP_Pos Pig1_GFP_NegCnt Pig1_GFP_PosCnt Pig2_GFP_NegCnt Pig2_GFP_PosCnt
# 2      LGR5   954.84920       7.039042 5.52e-29 6.67e-25    19.676455    1850.8681     9.893774    1938.9585              10            3218               7            3073
# 5    CRABP1 25950.47840       6.742196 3.16e-28 2.55e-24   623.743624   25480.7823   337.801700   77359.5859             317           44302             239          122605




################################################################################
# Heatmap top differential gene expression using counts data.
#   LGR5+ vs LGR5- by log2FoldChange and padj
#
################################################################################
#Get list of significant genes, sorted by padj
vec_geneSymbols <- rownames(res_signif_AvsB[1:120, ])

#Convert case
vec_geneSymbols <- toupper(unique(vec_geneSymbols))  
vec_geneSymbols <- sort(vec_geneSymbols)

#Reduce the size of the assay to only the significant genes.
signif_AvsB_counts <- res_signif_AvsB    #res_signif_AvsB[,c(1,6:9)]
nrow(signif_AvsB_counts)     

sigMat <- signif_AvsB_counts[signif_AvsB_counts$Gene %in% vec_geneSymbols,]
nrow(sigMat)  

sigGenes <- sigMat$Gene

#Reduce the dds_vst_AvsB to include only significant genes
dds_vst_AvsB_signif <- dds_vst_AvsB[dds_vst_AvsB$Gene %in% sigGenes,]
head(dds_vst_AvsB_signif,2)

#Sort in order
dds_vst_AvsB_signif_sorted <- dds_vst_AvsB_signif[order(match(dds_vst_AvsB_signif$Gene, sigGenes)),]  #sigGenes
head(dds_vst_AvsB_signif_sorted,2)

#Drop the "Gene" column
dds_vst_AvsB_signif_sorted$Gene <- NULL
nrow(dds_vst_AvsB_signif_sorted)

#-------------------------------------------------------------------------------
# Heatmap of top DE expressed genes
#-------------------------------------------------------------------------------
pheatmap(dds_vst_AvsB_signif_sorted,     # Flip rows and columns by transpose the input matrix (with t())
         scale = "row",
         legend = TRUE, 
         main = 'Top DE Genes; Pig Lung, PD111',
         show_rownames=TRUE,
         fontsize_row=8,
         fontsize_col=10,
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         drop_levels=TRUE)   



################################################################################
# Heatmap, Gerber Neural markers
#
################################################################################
df_Markers_Ref <- read.csv(here("../Markers","Neural_markers_Gerber.csv"), header = FALSE)

vec_Common_Markers <- df_Markers_Ref$V1
vec_geneSymbols <- c(vec_Common_Markers)

#Convert case
vec_geneSymbols <- toupper(unique(vec_geneSymbols))  
vec_geneSymbols <- sort(vec_geneSymbols)

#Reduce the size of the assay to only the significant genes.
signif_AvsB_counts <- res_signif_AvsB    #res_signif_AvsB[,c(1,6:9)]
nrow(signif_AvsB_counts)     

sigMat <- signif_AvsB_counts[signif_AvsB_counts$Gene %in% vec_geneSymbols,]
nrow(sigMat)  

sigGenes <- sigMat$Gene

#Reduce the dds_vst_AvsB to include only significant genes
dds_vst_AvsB_signif <- dds_vst_AvsB[dds_vst_AvsB$Gene %in% sigGenes,]
head(dds_vst_AvsB_signif,2)

#Sort in order
dds_vst_AvsB_signif_sorted <- dds_vst_AvsB_signif[order(match(dds_vst_AvsB_signif$Gene, vec_geneSymbols)),]  #sigGenes
head(dds_vst_AvsB_signif_sorted,2)

#Drop the "Gene" column
dds_vst_AvsB_signif_sorted$Gene <- NULL

nrow(dds_vst_AvsB_signif_sorted)

pheatmap(t(dds_vst_AvsB_signif_sorted),     # Flip rows and columns by transpose the input matrix (with t())
         scale = "row",
         legend = TRUE, 
         main = "Neural Markers; Pig Lung, PD111",
         show_rownames=TRUE,
         fontsize_row=8,
         fontsize_col=10,
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         drop_levels=TRUE)   

