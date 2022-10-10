library(Seurat)
library(URD)

# Load the dataset
# Convert to sparse matrices for efficiency
E35.data <- as.matrix(read.table("~/Documents/Lab/Single_Cell_Seq/Embryo_SCS_2018_matrices/sc.counts_E35.matrix", sep = "\t"))
E50.data <- as.matrix(read.table("~/Documents/Lab/Single_Cell_Seq/Embryo_SCS_2018_matrices/sc.counts_E50.matrix", sep = "\t"))
E65.data <- as.matrix(read.table("~/Documents/Lab/Single_Cell_Seq/Embryo_SCS_2018_matrices/sc.counts_E65.matrix", sep = "\t"))
E80.data <- as.matrix(read.table("~/Documents/Lab/Single_Cell_Seq/Embryo_SCS_2018_matrices/sc.counts_E80.matrix", sep = "\t"))
E95.data <- as.matrix(read.table("~/Documents/Lab/Single_Cell_Seq/Embryo_SCS_2018_matrices/sc.counts_E95.matrix", sep = "\t"))
E110.data <- as.matrix(read.table("~/Documents/Lab/Single_Cell_Seq/Embryo_SCS_2018_matrices/E110_adjusted.txt", sep = "\t"))
E125.data <- as.matrix(read.table("~/Documents/Lab/Single_Cell_Seq/Embryo_SCS_2018_matrices/E125_adjusted.txt", sep = "\t"))
E145.data <- as.matrix(read.table("~/Documents/Lab/Single_Cell_Seq/Embryo_SCS_2018_matrices/E145_adjusted.txt", sep = "\t"))
juv.data <- as.matrix(read.table("~/Documents/Lab/Single_Cell_Seq/Embryo_SCS_2018_matrices/sc.counts_J.matrix", sep = "\t"))


#Merge the datasets
list_embryos <- list(E35.data,E50.data,E65.data,E80.data,E95.data,E110.data,E125.data,E145.data,juv.data)
int_embryo_matrix <- dsCombineDGE(list_embryos)
write.table(int_embryo_matrix, "~/Documents/Lab/Single_Cell_Seq/UMAP_clustering/embryo_juv_adjust_3_6_20/count_matrix_adjusted_3_6_20.txt", sep = "\t")

int_embryo_matrix <- as.matrix(read.table("~/Documents/Lab/Single_Cell_Seq/UMAP_clustering/embryo_juv_adjust_3_6_20/count_matrix_adjusted_3_6_20.txt", sep = "\t"))

# Initialize the Seurat object with the raw (non-normalized data).
int_embryo <- CreateSeuratObject(counts = int_embryo_matrix, min.cells = 5, min.features =200)

# Filter the data
VlnPlot(int_embryo, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
#eliminate all cells that have more than 4000 genes
int_embryo<- subset(int_embryo, subset = nFeature_RNA > 500 & nFeature_RNA < 4000)

# Normalizing the data
int_embryo <- NormalizeData(object = int_embryo, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
int_embryo <- FindVariableFeatures(object = int_embryo, do.plot = F)

# Scaling the data and removing unwanted sources of variation
int_embryo <- ScaleData(object = int_embryo, display.progress = F)

# Perform linear dimensional reduction
int_embryo <- RunPCA(object = int_embryo, pc.genes = int_embryo@var.genes, do.print = TRUE)


# Run heatmap (optional)
DimHeatmap(object=int_embryo, dims = 1:20, cells = 100, balanced = TRUE)

#Elbowplot
ElbowPlot(object = int_embryo)


#jackstraw
int_embryo <- JackStraw(object = int_embryo, num.replicate = 100)
int_embryo <- ScoreJackStraw(object = int_embryo, dims = 1:20)
JackStrawPlot(object = int_embryo, dims = 1:20)


####################Run non-linear dimensional reduction (UMAP)
library(RColorBrewer)
library(colorRamps)
display.brewer.all()
colors <- brewer.pal(n = 9, name = "Set1")
colors <- rainbow(9)
colors <- primary.colors(9)

#rename the stages
current.stage.ids <- c("E351","E352","E353","E501", "E502","E503","E651","E652","E653","E801","E802","E803", "E951","E952","E953" ,"E1101","E1102","E1103","E1251","E1252", "E1253","E1451", "E1452", "E1453", "J1", "J2", "J3", "J4")
new.stage.ids <- c("E35","E35","E35","E50","E50","E50", "E65","E65","E65","E80","E80","E80", "E95","E95","E95", "E110","E110","E110","E125","E125","E125","E145", "E145","E145","J", "J","J", "J")
int_embryo@meta.data$orig.ident<- plyr::mapvalues(x = int_embryo@meta.data$orig.ident, from = current.stage.ids, to = new.stage.ids)

int_embryo <- FindNeighbors(object = int_embryo, dims = 1:20)
int_embryo <- FindClusters(int_embryo, resolution = 1.5, print.output = 0, save.SNN = T)
int_embryo <-  RunUMAP(int_embryo, reduction.use = "pca", dims= 1:20, min.dist= 0.3, n.neighbors =78 )
DimPlot(object = int_embryo, reduction = "umap", group.by= "orig.ident", cols = colors)
DimPlot(object = int_embryo, reduction = "umap", group.by= "orig.ident")
DimPlot(object = int_embryo, reduction = "umap", group.by= "orig.ident", split.by = "orig.ident", cols = colors)
DimPlot(object = int_embryo, reduction = "umap",label = T)

FeaturePlot(int_embryo, c("98012200-nec2"),pt.size = 0.5)


#cluster 9 is still weird and throwing things off. cluster 9, which comes from E110, is very similar to that of E145
#removing here
int_embryo_2.0 <- subset(int_embryo, idents = c(0:8,10:46))
########check to make sure that cluster0 is indeed gone
DimPlot(object = int_embryo_2.0, reduction = "umap", label = T)
DimPlot(object = int_embryo_2.0, reduction = "umap", group.by= "orig.ident", cols = colors)
########extract the new count matrix for downstream use
int_embryo_2.0_counts <- int_embryo_2.0@assays$RNA@data


###########################################################################################
#analysis of embryo and juvenile dataset without the funky clusters

# Initialize the Seurat object with the raw (non-normalized data).
int_embryo <- CreateSeuratObject(counts = int_embryo_2.0_counts, min.cells = 5, min.features =200)

# Filter the data
VlnPlot(int_embryo, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
#eliminate all cells that have more than 4000 genes
int_embryo<- subset(int_embryo, subset = nFeature_RNA > 500 & nFeature_RNA < 4000)

# Normalizing the data
int_embryo <- NormalizeData(object = int_embryo, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
int_embryo <- FindVariableFeatures(object = int_embryo, do.plot = F)

# Scaling the data and removing unwanted sources of variation
int_embryo <- ScaleData(object = int_embryo, display.progress = F)

# Perform linear dimensional reduction
int_embryo <- RunPCA(object = int_embryo, pc.genes = int_embryo@var.genes, do.print = TRUE)


# Run heatmap (optional)
DimHeatmap(object=int_embryo, dims = 1:20, cells = 100, balanced = TRUE)

#Elbowplot
ElbowPlot(object = int_embryo)


#jackstraw
int_embryo <- JackStraw(object = int_embryo, num.replicate = 100)
int_embryo <- ScoreJackStraw(object = int_embryo, dims = 1:20)
JackStrawPlot(object = int_embryo, dims = 1:20)


####################Run non-linear dimensional reduction (UMAP)
library(RColorBrewer)
library(colorRamps)
display.brewer.all()
colors <- brewer.pal(n = 9, name = "Set1")
colors <- rainbow(9)
colors <- primary.colors(9)

#rename the stages
current.stage.ids <- c("E351","E352","E353","E501", "E502","E503","E651","E652","E653","E801","E802","E803", "E951","E952","E953" ,"E1101","E1102","E1103","E1251","E1252", "E1253","E1451", "E1452", "E1453", "J1", "J2", "J3", "J4")
new.stage.ids <- c("E35","E35","E35","E50","E50","E50", "E65","E65","E65","E80","E80","E80", "E95","E95","E95", "E110","E110","E110","E125","E125","E125","E145", "E145","E145","J", "J","J", "J")
int_embryo@meta.data$orig.ident<- plyr::mapvalues(x = int_embryo@meta.data$orig.ident, from = current.stage.ids, to = new.stage.ids)

int_embryo <- FindNeighbors(object = int_embryo, dims = 1:20)
int_embryo <- FindClusters(int_embryo, resolution = 3.5, print.output = 0, save.SNN = T)
int_embryo <-  RunUMAP(int_embryo, reduction.use = "pca", dims= 1:20, min.dist= 0.3, n.neighbors =78)
DimPlot(object = int_embryo, reduction = "umap", group.by= "orig.ident", cols = c("E35" = "#F8766D", "E50" = "#D39200", "E65" = "#83AA00", "E80" = "#00BA38", "E95" = "#00C19F", "E110" = "#00B9E3", "E125" = "#619CFF", "E145" = "#DB72FB", "J" = "#FF61C3"))




DimPlot(object = int_embryo, reduction = "umap", group.by= "orig.ident", split.by = "orig.ident", cols = colors)

DimPlot(object = int_embryo, reduction = "umap", group.by= "orig.ident")

DimPlot(object = int_embryo, reduction = "umap", label = T)

current.stage.ids <- c("A-E35","B-E50", "C-E65","D-E80", "E-E95","F-E110","G-E125","H-E145", "J")
new.stage.ids <- c("1-E35","2-E50", "3-E65","4-E80", "5-E95","6-E110","7-E125","8-E145", "9-J")
int_embryo@meta.data$orig.ident<- plyr::mapvalues(x = int_embryo@meta.data$orig.ident, from = current.stage.ids, to = new.stage.ids)


VlnPlot(int_embryo, "98028414-ctcf", group.by = "orig.ident")

counts <- int_embryo@assays$RNA@data
metadata <- int_embryo@meta.data
write.table(metadata, "~/Documents/Lab/Single_Cell_Seq/Scanpy/PAGA_5_21_20/og_metadata.txt", sep = '\t')
embedding <- int_embryo@reductions$umap@cell.embeddings
write.table(embedding, "~/Documents/Lab/Single_Cell_Seq/Scanpy/PAGA_5_21_20/cell_embedding.txt", sep = '\t')


clust30 <- FindMarkers(int_embryo, ident.1 = 30, min.diff.pct = 0.25, only.pos = TRUE)
write.table(clust30, "~/Documents/Lab/Single_Cell_Seq/UMAP_clustering/embryo_juv_adjust_3_6_20/clust30_markers.txt", sep = '\t')

#now I need to isolate the neoblasts from other cell types present here
plot <- DimPlot(object = int_embryo, reduction = "umap", label = T)
select.cells <- CellSelector(plot = plot)
Idents(int_embryo, cells = select.cells) <- "NewCells"

subset_x<- SubsetData(int_embryo, ident.use = "NewCells")
DimPlot(object = subset_x, reduction = "umap", label = T)
DimPlot(object = subset_x, reduction = "umap", label = T, group.by= "orig.ident")
subset_x<- SubsetData(int_embryo, ident.use = "NewCells")


#subset out only the E145 cells
subset_x_145 <- subset(subset_x, subset = orig.ident == "E145")
DimPlot(object = subset_x_145, reduction = "umap", label = T, group.by= "orig.ident")
meta_E145_neoblast<- subset_x_145@meta.data
write.table(meta_E145_neoblast, "~/Documents/Lab/Single_Cell_Seq/UMAP_clustering/embryo_juv_adjust_3_6_20/meta_E145_neoblast.txt", sep = '\t')

#subset out only the E125 cells

subset_x_125 <- subset(subset_x, subset = orig.ident == "E125")
DimPlot(object = subset_x_125, reduction = "umap", label = T, group.by= "orig.ident")
meta_E125_neoblast<- subset_x_125@meta.data
write.table(meta_E125_neoblast, "~/Documents/Lab/Single_Cell_Seq/UMAP_clustering/embryo_juv_adjust_3_6_20/meta_E125_neoblast.txt", sep = '\t')




manual_neoblast_metadata<- subset_x@meta.data
manual_neoblast_metadata<- manual_neoblast_metadata[,-4]
write.table(manual_neoblast_metadata, "~/Documents/Lab/Single_Cell_Seq/Scanpy/PAGA_5_21_20/manual_neoblast_metadata.txt", sep = '\t')

E125.markers <- FindMarkers(int_embryo, ident.1 = "E125", group.by = 'orig.ident', min.diff.pct = 0.3, 
                                only.pos = TRUE)


FeaturePlot(int_embryo, c("98043523-piwl1-2"),pt.size = 0.5)

FeaturePlot(subset_x, c("98010043-actb-2"),pt.size = 0.5)


#adding the new metadata that is annotated with juvenile clusters
new_metadata <- as.matrix(read.table("~/Documents/Lab/Single_Cell_Seq/UMAP_clustering/embryo_juv_adjust_3_6_20/juv_neoblast_annot/outputfile.txt", sep = "\t"))

res <- new_metadata[,4]
int_embryo@meta.data <- cbind(int_embryo@meta.data, res)

#creating new umap plot with the annotated juvenile clusters
DimPlot(object = int_embryo, reduction = "umap", group.by= "res", label = T)

#found that a good proportion of the juvenile neoblasts have muscle-like cells. here we are isolating this population of cells by setting resolution to 3, and finding its markers. 
cluster28.markers <- FindMarkers(object = int_embryo, ident.1 =28, min.pct = 0.25)

FeaturePlot(int_embryo, c("98105511-grn-6"),pt.size = 0.5)

subset_int<- SubsetData(int_embryo, subset.name = "res", accept.value = "28")
subset_x<- SubsetData(int_embryo, ident.use = "28")


DimPlot(object = subset_int, reduction = "umap", group.by = "orig.ident", label = T)

DimPlot(object = subset_x, reduction = "umap",label = T, group.by = "orig.ident")


####################cell lineage markers

#endodermal
FeaturePlot(int_embryo, c("98051971-fabp5"),pt.size = 0.5)

#epidermal
FeaturePlot(int_embryo, c("98008272-cah6"),pt.size = 0.5)

#neural


#neoblast
FeaturePlot(int_embryo, c("98046269-ryr1"), pt.size = 0.5)

#muscle
FeaturePlot(int_embryo, c("98058687-adsv"), pt.size = 0.5)
FeaturePlot(int_embryo, c("98120717-gli3"), pt.size = 0.5)




##################for subsampling for downstream use in URD
set.seed(111)
int_embryo.subset <- subset(int_embryo, downsample = 2200)

DimPlot(object = int_embryo.subset, reduction = "umap", group.by= "orig.ident")

DimPlot(object = int_embryo.subset, reduction = "umap", label = T)

write.table(int_embryo.subset@assays$RNA@data, "~/Desktop/45000_subsample_counts.txt", sep = "\t")
write.table(int_embryo.subset@meta.data, "~/Desktop/45000_subsample_metadata.txt", sep = "\t" )


