#################################################################################
# Label clusters
#
#---------------------------------------------------------
# Ross Lampe
# Piedrahita Lab, College of Veterinary Medicine, North Carolina State University
# Last updated: 2022-10-03
#
#################################################################################

#Order of vector: 0,1,2,3,4,...
cluster.ids <- c("Mes 1",
                 "Mes 2",
                 "Mes 3",
                 "Mes 4",
                 "Mes 5",
                 "Mes 6",
                 "Mes 7",
                 "Epi 1",
                 "Epi 2",
                 "Mes 8",
                 "Epi 3",
                 "Mes 9",
                 "Mes 10",
                 "Epi 4")

names(cluster.ids) <- levels(obj_seurat)

obj_seurat <- RenameIdents(obj_seurat, cluster.ids)
levels(obj_seurat)

vec_Epi_Clusters <- c("Epi 1","Epi 2","Epi 3","Epi 4")
vec_Mesc_Clusters <- c("Mes 1",
                       "Mes 2",
                       "Mes 3",
                       "Mes 4",
                       "Mes 5",
                       "Mes 6",
                       "Mes 7","Mes 8","Mes 9","Mes 10")



DimPlot(obj_seurat, reduction ="umap", pt.size = 0.5, 
        label=TRUE, label.size=5, repel=TRUE) + NoLegend() + NoAxes()
