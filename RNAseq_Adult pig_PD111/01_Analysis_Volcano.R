################################################################################
# Volcano plot differential gene expression from bulk RNAseq.
#
#-------------------------------------------------------------------------------
# Ross Lampe
# Piedrahita Lab, College of Veterinary Medicine, North Carolina State University
# Last updated: 2022-08-28
# 
################################################################################
library(dplyr)
library(cowplot)
library(EnhancedVolcano)   #BiocManager::install("EnhancedVolcano")
library(here)

source("../themes.R")


################################################################################
# Read the DE results file.
# Note: the data are VST normalized which includes a log2 transformation.
#
################################################################################
res_AvsB <- read.csv(here("bulkRNA_D111 Adult Pigs", "rlogNorm__D111_AdultPig_GFPpos_vs_GFPneg_Annotated11.1.105.csv"), header = TRUE)
df_res_AvsB_shaped <- res_AvsB
head(df_res_AvsB_shaped,2)


#-------------------------------------------------------------------------------
# Inspect the results
#-------------------------------------------------------------------------------
#head(df_res_AvsB_shaped)
#summary(df_res_AvsB_shaped)

# How many adjusted p-values were less than 0.01?
sum(df_res_AvsB_shaped$padj < 0.01, na.rm=TRUE)

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



################################################################################
# Filter and shape the DE table.
#
################################################################################
#-------------------------------------------------------------------------------
# Filter by padj and log2FoldChange to find significant genes
#-------------------------------------------------------------------------------
# Set thresholds
padj.cutoff <- 0.05 # False Discovery Rate (FDR) cutoff
lfc.cutoff <- 0

#Filter by padj.cutoff & abs(log2FoldChange)
res_signif_AvsB <- df_res_AvsB_shaped %>%
    mutate(threshold_OE =  padj < padj.cutoff & abs(log2FoldChange) != lfc.cutoff)

#Clean up df
res_signif_AvsB <- res_signif_AvsB %>% tidyr::drop_na()                    #Use tidyr function to remove all rows in which NA appears
res_signif_AvsB <- res_signif_AvsB[res_signif_AvsB$threshold_OE == TRUE,]  #remove rows that fail the threshold tests
head(res_signif_AvsB,3)

#Drop the threshold_OE column
res_signif_AvsB <- res_signif_AvsB[ ,-c(15)]

#Save the gene names to row names.
rownames(res_signif_AvsB) <- res_signif_AvsB[, 1] 

#Sort summary list by p-value and log2FoldChange
res_signif_AvsB <- res_signif_AvsB %>% arrange(desc(log2FoldChange),desc(padj))  #Sorting a Dataframe in R Dplyr

#-------------------------------------------------------------------------------
# Remove low expressing genes
#-------------------------------------------------------------------------------
n_minTotalCount = 10
res_signif_AvsB_filtered <- subset(res_signif_AvsB,
                             ((res_signif_AvsB$GFPneg_Rep1_Cnt + res_signif_AvsB$GFPneg_Rep2_Cnt) >= n_minTotalCount) &
                             ((res_signif_AvsB$GFPpos_Rep1_Cnt + res_signif_AvsB$GFPpos_Rep2_Cnt) >= n_minTotalCount))
head(res_signif_AvsB_filtered,2)



################################################################################
# Volcano plot, Neural markers
#
################################################################################
vec_geneSymbols <- unique(toupper(vec_Neural_Markers_Shortlist))

#Filter for only dataframe rows with vec_geneSymbols
res_filtered <- res_signif_AvsB_filtered[res_signif_AvsB_filtered$Gene %in% vec_geneSymbols, ] 
head(res_filtered)

#-------------------------------------------------------------------------------
# Enhanced volcano plot (log2FoldChange vs -log10(padj))
#-------------------------------------------------------------------------------
p_V_AvsB <- EnhancedVolcano(res_filtered,
                            lab = res_filtered$Gene,
                            x = 'log2FoldChange',
                            y = 'padj',
                            ylab = bquote(~-Log[10] ~ pAdj),
                            title = "Neural Markers; Pig Lung, PD111; GFP+ vs GFP-",
                            titleLabSize = 14,
                            pCutoff = 10e-2,
                            FCcutoff = 0.58,
                            pointSize = 2.0,
                            labSize = 3.5,
                            colAlpha = 1,
                            legendPosition = 'none',
                            max.overlaps = 30,     # max.overlaps = 15
                            drawConnectors = TRUE,
                            maxoverlapsConnectors = NULL,
                            widthConnectors = 0.2,
                            colConnectors = "grey30")

p1 = p_V_AvsB + coord_cartesian(xlim=c(-5,5), ylim=c(-0.1,18.1))
p_V_AvsB

#write.csv(res_filtered, here("Results","PD111_Neural_Markers.csv"))




################################################################################
# Volcano plot, Neural markers: SC+PnC, EpC,EnC
#
################################################################################
#-------------------------------------------------------------------------------
# SC+PnC
#-------------------------------------------------------------------------------
vec_geneSymbols <- unique(toupper(vec_Neural_Markers_Shortlist_SC_PnC))

res_filtered <- res_signif_AvsB_filtered[res_signif_AvsB_filtered$Gene %in% vec_geneSymbols, ] 
head(res_filtered)

p_V_AvsB <- EnhancedVolcano(res_filtered,
                            lab = res_filtered$Gene,
                            x = 'log2FoldChange',
                            y = 'padj',
                            ylab = bquote(~-Log[10] ~ pAdj),
                            title = "SC, PnC Neural Markers; Pig Lung, PD111; GFP+ vs GFP-",
                            titleLabSize = 14,
                            pCutoff = 10e-2,
                            FCcutoff = 0.58,
                            pointSize = 2.0,
                            labSize = 3.5,
                            colAlpha = 1,
                            legendPosition = 'none',
                            max.overlaps = 30,     # max.overlaps = 15
                            drawConnectors = TRUE,
                            maxoverlapsConnectors = NULL,
                            widthConnectors = 0.2,
                            colConnectors = "grey30")

p1 = p_V_AvsB + coord_cartesian(xlim=c(-4,4), ylim=c(-0.1,5))
p_V_AvsB


#-------------------------------------------------------------------------------
# EpC_EnC
#-------------------------------------------------------------------------------
vec_geneSymbols <- unique(toupper(vec_Neural_Markers_Shortlist_EpC_EnC))

res_filtered <- res_signif_AvsB_filtered[res_signif_AvsB_filtered$Gene %in% vec_geneSymbols, ] 
head(res_filtered)

p_V_AvsB <- EnhancedVolcano(res_filtered,
                            lab = res_filtered$Gene,
                            x = 'log2FoldChange',
                            y = 'padj',
                            ylab = bquote(~-Log[10] ~ pAdj),
                            title = "EpC_EnC Neural Markers; Pig Lung, PD111; GFP+ vs GFP-",
                            titleLabSize = 14,
                            pCutoff = 10e-2,
                            FCcutoff = 0.58,
                            pointSize = 2.0,
                            labSize = 3.5,
                            colAlpha = 1,
                            legendPosition = 'none',
                            max.overlaps = 30,     # max.overlaps = 15
                            drawConnectors = TRUE,
                            maxoverlapsConnectors = NULL,
                            widthConnectors = 0.2,
                            colConnectors = "grey30")

p1 = p_V_AvsB + coord_cartesian(xlim=c(-5,5), ylim=c(-0.1,18.1))
p_V_AvsB



################################################################################
# Volcano plot, KEGG_WNT_pathway_ssc04310
#
################################################################################
vec_geneSymbols <- c(vec_LGR5_other,
                     vec_WNT_Markers)

vec_geneSymbols <- unique(toupper(vec_geneSymbols))  

res_filtered <- res_signif_AvsB_filtered[res_signif_AvsB_filtered$Gene %in% vec_geneSymbols, ] 
head(res_filtered)

#-------------------------------------------------------------------------------
# Enhanced volcano plot (log2FoldChange vs -log10(padj))
#-------------------------------------------------------------------------------
p_V_AvsB <- EnhancedVolcano(res_filtered,
                            lab = res_filtered$Gene,
                            x = 'log2FoldChange',
                            y = 'padj',
                            ylab = bquote(~-Log[10] ~ pAdj),
                            title = "KEGG_WNT_ssc04310; PD111; GFP+vsGFP-",
                            titleLabSize = 14,
                            pCutoff = 10e-2,
                            FCcutoff = 0.58,    #2,
                            pointSize = 2.0,
                            labSize = 3.5,
                            colAlpha = 1,
                            legendPosition = 'none',
                            max.overlaps = 30,     # max.overlaps = 15
                            drawConnectors = TRUE,
                            maxoverlapsConnectors = NULL,
                            widthConnectors = 0.2,
                            colConnectors = "grey30")

p1 = p_V_AvsB + coord_cartesian(ylim=c(-0.1,25))
p1

write.csv(res_filtered, here("Results","PD111_WNT_markers.csv"))



################################################################################
# Volcano plot, TGFbeta_pathway
#
################################################################################
vec_geneSymbols <- c(vec_TGFbeta_Markers)

vec_geneSymbols <- toupper(unique(vec_geneSymbols)) 

#Filter for only dataframe rows with vec_geneSymbols
res_filtered <- res_signif_AvsB_filtered[res_signif_AvsB_filtered$Gene %in% vec_geneSymbols, ] 
head(res_filtered)

#-------------------------------------------------------------------------------
# Enhanced volcano plot (log2FoldChange vs -log10(padj))
#-------------------------------------------------------------------------------
p_V_AvsB <- EnhancedVolcano(res_filtered,
                            lab = res_filtered$Gene,
                            x = 'log2FoldChange',
                            y = 'padj',
                            ylab = bquote(~-Log[10] ~ pAdj),
                            title = "KEGG_TGFβ_ssc04350; PD111; GFP+vsGFP-",
                            titleLabSize = 14,
                            pCutoff = 10e-2,
                            FCcutoff = 0.58,
                            pointSize = 2.0,
                            labSize = 3.5,
                            colAlpha = 1,
                            legendPosition = 'none',
                            max.overlaps = 30,     # max.overlaps = 15
                            drawConnectors = TRUE,
                            maxoverlapsConnectors = NULL,
                            widthConnectors = 0.2,
                            colConnectors = "grey30")

p1 = p_V_AvsB + coord_cartesian(ylim=c(-0.1,10))
p1

# Write the log2fold data to disk.
write.csv(res_filtered, here("Results","PD111_TGF-β markers.csv"))



################################################################################
# Volcano plot, Epithelial markers
#
################################################################################
vec_geneSymbols <- c(vec_Budtip_Epi_Markers)

vec_geneSymbols <- unique(toupper(vec_geneSymbols)) 

#Filter for only dataframe rows with vec_geneSymbols
res_filtered <- res_signif_AvsB_filtered[res_signif_AvsB_filtered$Gene %in% vec_geneSymbols, ] 
head(res_filtered)

#-------------------------------------------------------------------------------
# Enhanced volcano plot (log2FoldChange vs -log10(padj))
#-------------------------------------------------------------------------------
p_V_AvsB <- EnhancedVolcano(res_filtered,
                            lab = res_filtered$Gene,
                            x = 'log2FoldChange',
                            y = 'padj',
                            ylab = bquote(~-Log[10] ~ pAdj),
                            title = "Epithelial Markers; Pig Lung, PD111; GFP+ vs GFP-",
                            titleLabSize = 14,
                            pCutoff = 0.05,
                            FCcutoff = 0.58,
                            pointSize = 2.0,
                            labSize = 3.5,
                            colAlpha = 1,
                            legendPosition = 'none',
                            max.overlaps = 30,     # max.overlaps = 15
                            drawConnectors = TRUE,
                            maxoverlapsConnectors = NULL,
                            widthConnectors = 0.2,
                            colConnectors = "grey30")

p1 = p_V_AvsB + coord_cartesian(ylim=c(-0.1,12))
p1

write.csv(res_filtered, here("Results","PD111_Budtip_markers.csv"))
