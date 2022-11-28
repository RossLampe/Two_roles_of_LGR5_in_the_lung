################################################################################
# Gene Set Enrichment Analysis (GSEA)
# Molecular Signatures Database analysis (MSigDb)
#
#-------------------------------------------------------------------------------
# Ross Lampe
# Piedrahita Lab, College of Veterinary Medicine, North Carolina State University
# Last updated: 2022-10-08
# 
#-------------------------------------------------------------------------------
# Molecular Signatures Database analysis (MSigDb)
#   https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm
#   https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
#   https://www.biostars.org/tag/msigdb/
#
# Collections
#   H: hallmark gene sets
#   C1: positional gene sets
#   C2: curated gene sets
#   C3: motif gene sets
#   C4: computational gene sets
#   C5: GO gene sets
#   C6: oncogenic signatures
#   C7: immunologic signatures
#   C8: Cell type signature gene sets that contain curated cluster markers for cell types identified in single-cell sequencing studies of human tissue.
#
# The pig and human onthologies are very similar, but the human is much more studied.
# Use the human.
################################################################################
library(dplyr)
library(stringr)
library(ggplot2)
library(clusterProfiler)
library(fgsea)
library(msigdbr)
library(here)

source("../themes.R")

# packageVersion("msigdbr")   # [1] ‘7.5.1’
# msigdbr_species()
# 
# df_MSigDb_all_gene_sets = msigdbr(species = "Sus scrofa")
# # head(df_MSigDb_all_gene_sets)
# 
# df_MSigDb_all_gene_sets %>%
#   dplyr::filter(gs_cat == "H") %>%
#   head()


################################################################################
# Use the ranked list of genes generate using GO_KEGG.R
################################################################################
# df_Markers_ranked_Top_ordered is generated in "05_GO_KEGG.R"
head(df_Markers_ranked_Top_ordered,20)

vec_geneSymbols_Top <- df_Markers_ranked_Top_ordered$Gene



################################################################################
# H: hallmark gene sets
################################################################################
df_MSigDb_H <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol, entrez_gene, ensembl_gene)

# df_MSigDb_Ss <- msigdbr(species = "Sus scrofa", category = "H") %>% 
#   dplyr::select(gs_name, gene_symbol, entrez_gene, ensembl_gene)
# head(df_MSigDb_Ss)

em <- enricher(vec_geneSymbols_Top, TERM2GENE=df_MSigDb_H)
em <- as.data.frame(em)
rownames(em) <- NULL

#Inspect
head(em)



################################################################################
# MSigDb gene set enrichment analysis
# C3: motif gene sets
################################################################################
df_MSigDb_C3 <- msigdbr(species = "Homo sapiens", category = "C3") %>%
  dplyr::select(gs_name, gene_symbol, entrez_gene, ensembl_gene)

# df_MSigDb_C3 <- msigdbr(species = "Sus scrofa", category = "C3") %>% 
#   dplyr::select(gs_name, gene_symbol, entrez_gene, ensembl_gene)
# head(df_MSigDb_C3)

em <- enricher(vec_geneSymbols_Top, TERM2GENE=df_MSigDb_C3)
head(em)



################################################################################
# MSigDb cluster markers analysis
# C8: Gene sets that contain curated cluster markers for cell types identified in single-cell sequencing studies of human tissue
################################################################################
df_MSigDb_C8 <- msigdbr(species = "Homo sapiens", category = "C8") %>%
  dplyr::select(gs_name, gene_symbol, entrez_gene, ensembl_gene)

# df_MSigDb_C8 <- msigdbr(species = "Sus scrofa", category = "C8") %>% 
#   dplyr::select(gs_name, gene_symbol, entrez_gene, ensembl_gene)
# head(df_MSigDb_C8)

em <- enricher(vec_geneSymbols_Top, TERM2GENE=df_MSigDb_C8)

em <- as.data.frame(em)
rownames(em) <- NULL
em <- em %>% 
  select_if(!names(.) %in% c('ID','geneID')) 

head(em,2)

#Convert the GeneRatio (= count/setSize) to a decimal
em$x_GeneRatio <- sapply(strsplit(em$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

#Order by enrichment (GeneRatio)
em_ordered <- em[order(-em$x_GeneRatio), ]
head(em_ordered,10)

# VISUALIZATION
p1 <- ggplot(head(em_ordered,15), 
       aes(x = x_GeneRatio, y = reorder(Description, x_GeneRatio))) +   
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_colour_gradient(low="red", high = "blue") +
  # scale_x_log10(name="p.adjust") +
  theme_bw() +
  xlab("Enrichment score") +
  ylab(NULL) +
  #coord_flip() +
  ggtitle("GSEA C8 cell type enrichment\nGFP+ vs GFP-; RNAseq_PD111 \nranked by L2F, pAdj, n=500")
p1 <- theme_GO_dotplot(p1)
p1

