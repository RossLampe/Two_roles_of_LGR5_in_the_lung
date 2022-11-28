###############################################################################
# Gene Set Enrichment Analysis (GSEA)
# Molecular Signatures Database analysis (MSigDb)
#   https://www.gsea-msigdb.org/gsea/msigdb/index.jsp
#   http://www.gsea-msigdb.org/gsea/msigdb/mouse/annotate.jsp
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
library(clusterProfiler)
library(ggplot2)
library(fgsea)
library(msigdbr)

packageVersion("msigdbr")   # [1] ‘7.5.1’
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
head(df_Markers_Mesc_vs_Epi_ranked_Top_ordered,20)
vec_geneSymbols_Mesc_vs_Epi_Top <- df_Markers_Mesc_vs_Epi_ranked_Top_ordered$Gene

#-------------------------------------------------------------------------------
# Epi top markers
#-------------------------------------------------------------------------------
head(df_Markers_Epi_ranked_Top,10)
vec_geneSymbols_Epi_Top = df_Markers_Epi_ranked_Top$Gene

#-------------------------------------------------------------------------------
# Mesc top  markers
#-------------------------------------------------------------------------------
head(df_Markers_Mesc_ranked_Top,10)
vec_geneSymbols_Mesc_Top = df_Markers_Mesc_ranked_Top$Gene



################################################################################
# H: hallmark gene sets
################################################################################
df_MSigDb_H <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol, entrez_gene, ensembl_gene)

# df_MSigDb_Ss <- msigdbr(species = "Sus scrofa", category = "H") %>% 
#   dplyr::select(gs_name, gene_symbol, entrez_gene, ensembl_gene)
# head(df_MSigDb_Ss)

#-------------------------------------------------------------------------------
# Epi top markers
#-------------------------------------------------------------------------------
em <- enricher(vec_geneSymbols_Epi_Top, TERM2GENE=df_MSigDb_H)
em <- as.data.frame(em)
rownames(em) <- NULL

#Inspect
head(em)

em <- as.data.frame(em)
rownames(em) <- NULL

#Convert the GeneRatio (= count/setSize) to a decimal
em$dec_GeneRatio <- sapply(strsplit(em$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

ggplot(head(em,15), aes(x = dec_GeneRatio, y = reorder(Description, dec_GeneRatio))) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_colour_gradient(low="red", high = "blue") +
  theme_bw() +
  xlab("Enrichment score") +
  ylab(NULL) +
  ggtitle("GSEA H Hallmark gene sets \nEpithelial Clusters; scRNAseq GD80") +
  theme(plot.title = element_text(hjust = 0.5))


#-------------------------------------------------------------------------------
# Mesc top  markers
#-------------------------------------------------------------------------------
em <- enricher(vec_geneSymbols_Mesc_Top, TERM2GENE=df_MSigDb_H)
em <- as.data.frame(em)
rownames(em) <- NULL

#Convert the GeneRatio (= count/setSize) to a decimal
em$dec_GeneRatio <- sapply(strsplit(em$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

ggplot(head(em,15), aes(x = dec_GeneRatio, y = reorder(Description, dec_GeneRatio))) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_colour_gradient(low="red", high = "blue") +
  theme_bw() +
  xlab("Enrichment score") +
  ylab(NULL) +
  ggtitle("GSEA H Hallmark gene sets  \nMesenchymal Clusters; scRNAseq GD80") +
  theme(plot.title = element_text(hjust = 0.5))



################################################################################
# MSigDb gene set enrichment analysis
# C3: motif gene sets
################################################################################
df_MSigDb_C3 <- msigdbr(species = "Homo sapiens", category = "C3") %>%
  dplyr::select(gs_name, gene_symbol, entrez_gene, ensembl_gene)

# df_MSigDb_C3 <- msigdbr(species = "Sus scrofa", category = "C3") %>% 
#   dplyr::select(gs_name, gene_symbol, entrez_gene, ensembl_gene)
# head(df_MSigDb_C3)

#-------------------------------------------------------------------------------
# Epi top markers
#-------------------------------------------------------------------------------
em <- enricher(vec_geneSymbols_Epi_Top, TERM2GENE=df_MSigDb_C3)
head(em)

em <- as.data.frame(em)
rownames(em) <- NULL

#Convert the GeneRatio (= count/setSize) to a decimal
em$dec_GeneRatio <- sapply(strsplit(em$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

ggplot(head(em,15), aes(x = dec_GeneRatio, y = reorder(Description, dec_GeneRatio))) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_colour_gradient(low="red", high = "blue") +
  theme_bw() +
  xlab("Enrichment score") +
  ylab(NULL) +
  ggtitle("GSEA C3 motif gene sets \nEpithelial Clusters; scRNAseq GD80") +
  theme(plot.title = element_text(hjust = 0.5))


#-------------------------------------------------------------------------------
# Mesc top  markers
#-------------------------------------------------------------------------------
em <- enricher(vec_geneSymbols_Mesc_Top, TERM2GENE=df_MSigDb_C3)
em <- as.data.frame(em)
rownames(em) <- NULL

#Convert the GeneRatio (= count/setSize) to a decimal
em$dec_GeneRatio <- sapply(strsplit(em$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

ggplot(head(em,15), aes(x = dec_GeneRatio, y = reorder(Description, dec_GeneRatio))) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_colour_gradient(low="red", high = "blue") +
  theme_bw() +
  xlab("Enrichment score") +
  ylab(NULL) +
  ggtitle("GSEA C3 motif gene sets  \nMesenchymal Clusters; scRNAseq GD80") +
  theme(plot.title = element_text(hjust = 0.5))



##############################################################################
# MSigDb cluster markers analysis
# C8: Gene sets that contain curated cluster markers for cell types identified in single-cell sequencing studies of human tissue
################################################################################
df_MSigDb_C8 <- msigdbr(species = "Homo sapiens", category = "C8") %>%
  dplyr::select(gs_name, gene_symbol, entrez_gene, ensembl_gene)

# df_MSigDb_C8 <- msigdbr(species = "Sus scrofa", category = "C8") %>% 
#   dplyr::select(gs_name, gene_symbol, entrez_gene, ensembl_gene)
# head(df_MSigDb_C8)

head(df_MSigDb_C8)

#-------------------------------------------------------------------------------
# Epi top markers
#-------------------------------------------------------------------------------
em <- enricher(vec_geneSymbols_Epi_Top, TERM2GENE=df_MSigDb_C8)
em <- as.data.frame(em)
rownames(em) <- NULL
head(em)

#Convert the GeneRatio (= count/setSize) to a decimal
em$dec_GeneRatio <- sapply(strsplit(em$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

# VISUALIZATION
# #png(filename = "GSEA_cell_type-enrichment-2.png", width = 1200, height = 600)
ggplot(head(em,15), aes(x = dec_GeneRatio, y = reorder(Description, dec_GeneRatio))) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_colour_gradient(low="red", high = "blue") +
  theme_bw() +
  xlab("Enrichment score") +
  ylab(NULL) +
  ggtitle("GSEA C8 ccell type enrichment \nEpithelial Clusters; scRNAseq GD80 \nranked by L2F, pAdj<0.05; n=500") +
  theme(plot.title = element_text(hjust = 0.5))


#-------------------------------------------------------------------------------
# Mesc top  markers
#-------------------------------------------------------------------------------
em <- enricher(vec_geneSymbols_Mesc_Top, TERM2GENE=df_MSigDb_C8)
em <- as.data.frame(em)
rownames(em) <- NULL

#Convert the GeneRatio (= count/setSize) to a decimal
em$dec_GeneRatio <- sapply(strsplit(em$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

# VISUALIZATION
# #png(filename = "GSEA_cell_type-enrichment-2.png", width = 1200, height = 600)
ggplot(head(em,15), aes(x = dec_GeneRatio, y = reorder(Description, dec_GeneRatio))) + 
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_colour_gradient(low="red", high = "blue") +
  theme_bw() +
  xlab("Enrichment score") +
  ylab(NULL) +
  ggtitle("GSEA C8 cell type enrichment \nMesenchymal Clusters; scRNAseq GD80") +
  theme(plot.title = element_text(hjust = 0.5))

