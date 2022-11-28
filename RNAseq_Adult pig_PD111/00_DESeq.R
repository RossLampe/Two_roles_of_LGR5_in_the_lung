################################################################################
# Create differential gene expression files from counts,
# for each technical replicate.
#
# Ross Lampe
# Last updated: 2022-09-02
#
#-------------------------------------------------------------------------------
#"Piedrahita_pig_bulkRNA_salmon_transcriptLevel_quants_summed_by_ENSSSCG"
# (GFPpos-vs-GFPneg)
#   GFPneg_Rep1 
#   GFPneg_Rep1
#   GFPneg_Rep2
#   GFPneg_Rep2
# 
################################################################################
library(dplyr)
library(DESeq2)
library(ggplot2)
library(vsn)  #BiocManager::install("vsn")
library(here)


################################################################################
# Read counts files into dataframes.
# 
################################################################################
# Read the data into a dataframe
df_countData <- read.table(here("Data", "Piedrahita_pig_bulkRNA_salmon_transcriptLevel_quants_summed_by_ENSSSCG.csv"),
                                header=TRUE, sep=',', stringsAsFactors=FALSE) 

#-------------------------------------------------------------------------------
# Save the gene names to row names.
# Drop the Gene column to leave only count columns.
#-------------------------------------------------------------------------------
#Eliminate the few duplicated gene names, like ENSSSCG00000010931
df_countData <- df_countData[!duplicated(df_countData$Gene),]

nrow(df_countData)
#[1] 29263
head(df_countData,2)  
#            Gene GFPneg_Rep1 GFPpos_Rep1 GFPneg_Rep2 GFPpos_Rep2
# 5S_rRNA 5S_rRNA           0           0           0           0
# A1CF       A1CF           3           0          30           0

#Save the gene names to row names.
rownames(df_countData) <- df_countData[, 1] 

#Save the rownames
rownames  <- row.names(df_countData)

#Drop the Gene column to leave only count columns.
df_countData <- df_countData[ ,-c(1)]
head(df_countData,2)
#         GFPneg_Rep1 GFPpos_Rep1 GFPneg_Rep2 GFPpos_Rep2
# 5S_rRNA           0           0           0           0
# A1CF              3           0          30           0

#-------------------------------------------------------------------------------
# Convert the dataframe to a matrix.
#-------------------------------------------------------------------------------
mtrx_countData <- as.matrix(df_countData)
#head(mtrx_countData)

# Assign conditions 
(condition <- factor(c("neg","pos","neg","pos")))



################################################################################
# Construct DESeqDataSet objects
#
################################################################################
#Create a colData dataframe and instantiate the DESeqDataSet.
(df_countData <- data.frame(row.names=colnames(mtrx_countData), condition))
head(df_countData)

# Prepare dataset for DESeq.
dds_AvsB <- DESeqDataSetFromMatrix(countData=mtrx_countData, colData=df_countData, design=~condition)

#Run the DESeq pipeline for differential expression steps.
dds_AvsB <- DESeq(dds_AvsB, quiet=TRUE)

#Standard analysis
res_AvsB <- results(dds_AvsB)
head(res_AvsB,2)
#            baseMean log2FoldChange     lfcSE      stat      pvalue       padj
#           <numeric>      <numeric> <numeric> <numeric>   <numeric>  <numeric>
# 5S_rRNA     0.00000             NA        NA        NA          NA         NA
# A1CF       12.07621      -7.775581  2.318061 -3.354347 0.000795525 0.00501675


################################################################################
# Shape and annotate the data, then write to a file.
#
################################################################################
#-------------------------------------------------------------------------------
# Combine columns for writing to file
#-------------------------------------------------------------------------------
df_res_AvsB_shaped <- cbind(Gene <- rownames, res_AvsB$baseMean, res_AvsB$log2FoldChange, res_AvsB$stat, res_AvsB$pvalue, res_AvsB$padj, 
                            counts(dds_AvsB, normalized=TRUE), counts(dds_AvsB, normalized=FALSE))
head(df_res_AvsB_shaped)
#view(df_output_lfc_AvsB)

colnames(df_res_AvsB_shaped)[1] <- "Gene"
colnames(df_res_AvsB_shaped)[2] <- "baseMean"
colnames(df_res_AvsB_shaped)[3] <- "log2FoldChange"
colnames(df_res_AvsB_shaped)[4] <- "stat"
colnames(df_res_AvsB_shaped)[5] <- "pvalue"
colnames(df_res_AvsB_shaped)[6] <- "padj"

#Replace all NA with o
df_res_AvsB_shaped[is.na(df_res_AvsB_shaped)] <- 0
head(df_res_AvsB_shaped)


#-------------------------------------------------------------------------------
# Write to file
#------------------------------------------------------------------------------_
write.csv(df_res_AvsB_shaped,
          file=here("bulkRNA_D111 Adult Pigs", 
                    "rlogNorm__D111_AdultPig_GFPpos_vs_GFPneg_Annotated11.1.105.csv"), row.names = FALSE)



################################################################################
# Plot the counts of some of the top differentially expressed genes as a sanity check
#     https://biohpc.cornell.edu/doc/RNA-Seq-2020-exercise2.html
#     http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#model-matrix-not-full-rank
################################################################################
#gene=which.min(res_AvsB$padj)
# d <- plotCounts(dds_AvsB, 
#                 gene=which.min(res_AvsB$padj), 
#                 intgroup="condition", 
#                 returnData=TRUE)
# ggplot(d, aes(x=condition, y=count)) +
#   geom_point(position=position_jitter(w=0.1,h=0)) +
#   scale_y_log10(breaks=c(10,100,1000,4000)) + 
#   ggtitle("min(res_dds_AvsB$padj)") +
#   xlab('Condition') + ylab('Log10(count)')   #axis labels 

#gene="LGR5"
d <- plotCounts(dds_AvsB, 
                gene="LGR5", 
                intgroup="condition", 
                returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) +
  geom_point(position=position_jitter(w=0.12,h=0)) +
  scale_y_log10(breaks=c(10,100,1000,2000)) + 
  ggtitle("LGR5; PD111")
  labs(x='Condition', 
       y='Log10(count)',
       title = "LGR5") 



################################################################################
# Count data transformations - VST.
#-------------------------------------------------------------------------------
# The point of these two transformations, the VST and the rlog, 
# is to remove the dependence of the variance on the mean, 
# particularly the high variance of the logarithm of count data when the mean is low. 
#     http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#de
################################################################################
# variance stabilizing transformation.
# The transformed data is on the log2 scale for large counts.
vst_AvsB <- vst(dds_AvsB, blind=FALSE)
head(vst_AvsB)

# Convert the DESeq transformed object to a data frame
dds_vst_AvsB <- assay(vst_AvsB)
dds_vst_AvsB <- as.data.frame(dds_vst_AvsB)
dds_vst_AvsB$Gene <- rownames(dds_vst_AvsB)
head(dds_vst_AvsB,2)
# GFPneg_Rep1 GFPpos_Rep1 GFPneg_Rep2 GFPpos_Rep2    Gene
# 5S_rRNA    4.896495    4.896495    4.896495    4.896495 5S_rRNA

#Rename the columns
colnames(dds_vst_AvsB)[1] <- "Pig1_GFP_Neg"
colnames(dds_vst_AvsB)[2] <- "Pig1_GFP_Pos"
colnames(dds_vst_AvsB)[3] <- "Pig2_GFP_Neg"
colnames(dds_vst_AvsB)[4] <- "Pig2_GFP_Pos"


#Reorder the columns to put Gene in the first position.
dds_vst_AvsB <- select(dds_vst_AvsB,5,1,2,3,4)
#dds_vst_AvsB <- select(dds_vst_AvsB,1,3,5,2,4)
head(dds_vst_AvsB,2)
#            Gene Pig1_GFP_Pos Pig2_GFP_Pos Pig1_GFP_Neg Pig2_GFP_Neg
# 5S_rRNA 5S_rRNA     4.896495     4.896495     4.896495     4.896495
# A1CF       A1CF     4.896495     4.896495     5.533569     6.529301


#-------------------------------------------------------------------------------
#Plot to confirm that the variance stabilized data has a standard deviation 
#that is roughly constant along the whole dynamic range.
#-------------------------------------------------------------------------------
#meanSdPlot(assay(vst_AvsB))

#-------------------------------------------------------------------------------
# PCA plot of vst data.
#-------------------------------------------------------------------------------
plotPCA(vst_AvsB, intgroup="condition")
