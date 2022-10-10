facs_matrix <- as.matrix(read.table("/Users/JulianKimura/Desktop/neoblast_smartseq/P1_P2_P3_quantfiles/combined_counts_filtered.txt", sep = "\t", header = T, row.names = 1))


facs_smart <- CreateSeuratObject(counts = facs_matrix, min.cells = 5, min.features =200)

VlnPlot(facs_smart, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
#facs_smart<- subset(facs_smart, subset = nFeature_RNA > 3000)


# Normalizing the data
facs_smart <- NormalizeData(object = facs_smart, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
facs_smart <- FindVariableFeatures(object = facs_smart, do.plot = F)

# Scaling the data and removing unwanted sources of variation
facs_smart <- ScaleData(object = facs_smart, display.progress = F)

# Perform linear dimensional reduction
facs_smart <- RunPCA(object = facs_smart, pc.genes = facs_smart@var.genes, do.print = TRUE)


#jackstraw
facs_smart <- JackStraw(object = facs_smart, num.replicate = 100)
facs_smart <- ScoreJackStraw(object = facs_smart, dims = 1:20)
JackStrawPlot(object = facs_smart, dims = 1:20)


#umapc(1:12,14,17,19)
facs_smart <- FindNeighbors(object = facs_smart, dims = 1:20)
facs_smart <- FindClusters(facs_smart, resolution = 1.0, print.output = 0, save.SNN = T)
facs_smart <-  RunUMAP(facs_smart, reduction.use = "pca", dims= 1:20, min.dist=1, n.neighbors =49 )
DimPlot(object = facs_smart, reduction = "umap",label = T, pt.size = 2)
plot+ggtitle("myplot")

metadata <- facs_smart@meta.data

write.table(metadata, "~/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/combined_metadata.txt", sep = "\t" )

new_metadata <- as.matrix(read.table("/Users/JulianKimura/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/combined_metadata.txt", sep = "\t", header = T, row.names = 1))


#changing adding a new column called stages in the metadata
stages <- new_metadata[,6]
facs_smart@meta.data <- cbind(facs_smart@meta.data, stages)
dim_plot <- DimPlot(object = facs_smart, reduction = "umap", group.by = "facs", pt.size = 2)
dim_plot$data$stages <- factor(x = dim_plot$data$stages, levels = c("R", "G"))
dim_plot



vlnplot <- VlnPlot(facs_smart, features ="facs")
vlnplot+geom_boxplot(width=0.3, fill = "white")


RidgePlot(facs_smart, "facs")
VlnPlot(facs_smart, "facs")


plot <- DimPlot(object = facs_smart, reduction = "umap", cols = c("0" = "#4CAF50","1" = "#64B5F6", "2"= "#BA68C8", "3" = "#3F51B5", "4" = "#FF6F00", "5" ="#FF69F1", "6" = "#F0CE24", "7" = "#D82339"),pt.size = 2, label = T)
plot+ggtitle("myplot")

red_markers <- FindMarkers(object = facs_smart, ident.1 = "R" , group.by = "stages",min.pct = 0.25)



green.markers <- FindMarkers(object = facs_smart, group.by = "stages", ident.1 = "G" ,min.pct = 0.25)
red.markers <- FindMarkers(object = facs_smart, group.by = "stages", ident.1 = "R" ,min.pct = 0.25)
R_G_cluster0.markers <- FindMarkers(object = facs_smart, ident.1 = "0" ,min.pct = 0.25)
R_G_cluster1.markers <- FindMarkers(object = facs_smart, ident.1 = "1" ,min.pct = 0.25)
R_G_cluster2.markers <- FindMarkers(object = facs_smart, ident.1 = "2" ,min.pct = 0.25)
R_G_cluster3.markers <- FindMarkers(object = facs_smart, ident.1 = "3" ,min.pct = 0.25)
R_G_cluster4.markers <- FindMarkers(object = facs_smart, ident.1 = "4" ,min.pct = 0.25)
R_G_cluster5.markers <- FindMarkers(object = facs_smart, ident.1 = "5" ,min.pct = 0.25)
R_G_cluster6.markers <- FindMarkers(object = facs_smart, ident.1 = "6" ,min.pct = 0.25)
R_G_cluster7.markers <- FindMarkers(object = facs_smart, ident.1 = "7" ,min.pct = 0.25)

write.table(R_G_cluster0.markers, "~/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/clust_0.txt", sep = "\t" )
write.table(R_G_cluster1.markers, "~/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/clust_1.txt", sep = "\t" )
write.table(R_G_cluster2.markers, "~/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/clust_2.txt", sep = "\t" )
write.table(R_G_cluster3.markers, "~/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/clust_3.txt", sep = "\t" )
write.table(R_G_cluster4.markers, "~/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/clust_4.txt", sep = "\t" )
write.table(R_G_cluster5.markers, "~/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/clust_5.txt", sep = "\t" )
write.table(R_G_cluster6.markers, "~/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/clust_6.txt", sep = "\t" )
write.table(R_G_cluster7.markers, "~/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/clust_7.txt", sep = "\t" )


FeaturePlot(object=facs_smart, c("98050217-msx1-2"), pt.size = 2)



facs_smart@meta.data$stages <- factor(x = facs_smart@meta.data$stages, levels = c("R", "G"))

VlnPlot(facs_smart, features = c("98043523-piwl1-2", "98132150-h10-6", "98029236-pcna"), ncol = 3, group.by = "stages")
RidgePlot(facs_smart, features = "98043523-piwl1-2", group.by = "stages")


####subset the red only cells
red_only <- subset(facs_smart, subset = stages =="R")
DimPlot(red_only, reduction = "umap", label = T)
RidgePlot(red_only, features = "98043523-piwl1-2")
VlnPlot(red_only, features = "98043523-piwl1-2")
FeaturePlot(object=red_only, c("98049"), pt.size = 2)


green_only <- subset(facs_smart, subset = stages =="G")
DimPlot(green_only, reduction = "umap", label = T)
RidgePlot(green_only, features = "98043523-piwl1-2")
VlnPlot(green_only, features = "98043523-piwl1-2")
FeaturePlot(object=green_only, c("98049625-nfic"), pt.size = 2)



VlnPlot(facs_smart, features = c("98043523-piwl1-2", "98132150-h10-6", "98029236-pcna"), ncol = 3, group.by = "stages")


vlnplot+geom_boxplot(width=0.3, aes(x=Identity, y= Expression Level,  fill = "white"))

box_plot_df = data.frame(piwi1= facs_smart[["RNA"]]@data["98043523-piwl1-2",], identity= facs_smart@meta.data$stages)

t.test(piwi1 ~ identity, data = box_plot_df)

ggplot(box_plot_df, mapping = aes(x=identity, y=piwi1)) + geom_boxplot(mapping = aes(fill=identity)) +geom_signif(comparisons = c("G", "R")) + geom_signif(comparisons = "R", "G", map_signif_level = T)


####isolating the neoblast cluster then making a pie chart
neoblast_only <- subset(facs_smart, idents = "1")
DimPlot(neoblast_only, reduction = "umap", label = T, group.by = "stages")


neoblast_table$Var1 <- as.character(neoblast_table$Var1)

bp <- ggplot(neoblast_table, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = c( "#F8766D","#00BFC4"))
theme(legend.title = element_blank())
bp + coord_polar("y", start=0)



#Making a stacked bar graph

pt <- table(Idents(facs_smart), facs_smart$stages)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

green_pie_matrix <- pt[-c(1, 2, 3, 4, 5,6,7,8), ]
red_pie_matrix <- pt[-c(9,10,11,12,13,14,15,16), ]

hue_pal()(2)                             # Identify hex codes





#for making the pie chart for green cells
bp <- ggplot(green_pie_matrix, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = c( "#4CAF50","#64B5F6", "#BA68C8", "#3F51B5", "#FF6F00", "#FF69F1", "#F0CE24", "#D82339"))
  theme(legend.title = element_blank())
bp + coord_polar("y", start=0)


#for making the pie chart for red cells
bp <- ggplot(red_pie_matrix, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = c("#4CAF50","#64B5F6", "#BA68C8", "#3F51B5", "#FF6F00", "#FF69F1", "#F0CE24", "#D82339")) +
  theme(legend.title = element_blank())
bp + coord_polar("y", start=0)


#performing a chi squared test on the frequency of clusters
chi_squared_table <- as.matrix(read.table("/Users/JulianKimura/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/chi_squared_matrix.txt", sep = "\t", header = T, row.names = 1))

chisq <- chisq.test(chi_squared_table)
chisq

chisq$observed
round(chisq$expected,2)
round(chisq$residuals, 3)

library(corrplot)
pal <- colorRampPalette(c("#F9C449", "#F04393"))
corrplot(chisq$residuals, is.cor = FALSE, col = pal(400))

pal <- colorRampPalette(c("#031B88", "#FB7BBE"))
corrplot(chisq$residuals, is.cor = FALSE, col = pal(400))

pal <- colorRampPalette(c("#D39902", "#A60A3D"))

pal <- colorRampPalette(c("#7FACD6", "#A5678E"))

pal <- colorRampPalette(c("#FFFFFF", "#000000"))



#####################now it is time to subset out only the green and only the red to see if we can identify the heterogeneity of these cells

Idents(object = facs_smart) <- facs_smart@meta.data$stages
green_only <- subset(facs_smart, idents ="G")


green_only <- NormalizeData(object = green_only, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
green_only <- FindVariableFeatures(object = green_only, do.plot = F)

# Scaling the data and removing unwanted sources of variation
green_only <- ScaleData(object = green_only, display.progress = F)

# Perform linear dimensional reduction
green_only <- RunPCA(object = green_only, pc.genes = green_only@var.genes, do.print = TRUE)


#jackstraw
green_only <- JackStraw(object = green_only, num.replicate = 100)
green_only <- ScoreJackStraw(object = green_only, dims = 1:20)
JackStrawPlot(object = green_only, dims = 1:20)


#umapc(1:12,14,17,19)c(1:11,14,15)
green_only <- FindNeighbors(object = green_only, dims = 1:20)
green_only <- FindClusters(green_only, resolution = 1.0, print.output = 0, save.SNN = T)
green_only <-  RunUMAP(green_only, reduction.use = "pca", dims= 1:20, min.dist= 0.1, n.neighbors =49 )
DimPlot(object = green_only, reduction = "umap",label = T)

green_cluster0.markers <- FindMarkers(object = green_only, ident.1 = "0" ,min.pct = 0.25)
green_cluster1.markers <- FindMarkers(object = green_only, ident.1 = "1" ,min.pct = 0.25)
green_cluster2.markers <- FindMarkers(object = green_only, ident.1 = "2" ,min.pct = 0.25)
green_cluster3.markers <- FindMarkers(object = green_only, ident.1 = "3" ,min.pct = 0.25)


FeaturePlot(object=green_only, c("98009983-tbx2-2"))



############now we will subset only the red cells

red_only <- subset(facs_smart, idents ="R")

red_only <- NormalizeData(object = red_only, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
red_only <- FindVariableFeatures(object = red_only, do.plot = F)

# Scaling the data and removing unwanted sources of variation
red_only <- ScaleData(object = red_only, display.progress = F)

# Perform linear dimensional reduction
red_only <- RunPCA(object = red_only, pc.genes = red_only@var.genes, do.print = TRUE)


#jackstraw
red_only <- JackStraw(object = red_only, num.replicate = 100)
red_only <- ScoreJackStraw(object = red_only, dims = 1:20)
JackStrawPlot(object = red_only, dims = 1:20)


#umap
red_only <- FindNeighbors(object = red_only, dims = 1:20)
red_only <- FindClusters(red_only, resolution = 1.0, print.output = 0, save.SNN = T)
red_only <-  RunUMAP(red_only, reduction.use = "pca", dims= 1:20, min.dist= 0.1, n.neighbors =49 )
DimPlot(object = red_only, reduction = "umap",label = T)

red_cluster0.markers <- FindMarkers(object = red_only, ident.1 = "0" ,min.pct = 0.25)
red_cluster1.markers <- FindMarkers(object = red_only, ident.1 = "1" ,min.pct = 0.25)
red_cluster2.markers <- FindMarkers(object = red_only, ident.1 = "2" ,min.pct = 0.25)


FeaturePlot(object=red_only, c("98021231-dyh5"))

#############subsetting the dataset based off of the level of red fluorescence

####adding facs fluorescence levels onto the umap
facs_metadata <- as.matrix(read.table("/Users/JulianKimura/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/new_metadata.txt", sep = "\t", header = T, row.names = 1))

#changing adding a new column called facs in the metadata
facs <- facs_metadata[,7]
facs_smart@meta.data$nCount_RNA <- facs

facs <- as.numeric(facs_metadata[,7])
facs_smart@meta.data<- cbind(facs_smart@meta.data, facs)

FeaturePlot(facs_smart, features = "facs")

VlnPlot(facs_smart, "facs")

#subset out the facs fluorescence greater than 0.5
facs_red<- subset(facs_smart, subset = facs > 0.5)

VlnPlot(facs_red, "facs", group.by = "stages")

DimPlot(object = facs_red, reduction = "umap",label = T, pt.size = 2)


#subset out the green cells
VlnPlot(facs_smart, "facs", group.by = "stages")

facs_green<- subset(facs_smart, subset = facs < 0.2)

VlnPlot(facs_green, "facs", group.by = "stages")

DimPlot(object = facs_green, reduction = "umap",label = T, pt.size = 2)
plot+ggtitle("myplot")

#extract the matrix
green_subset_matrix <- facs_green@assays$RNA@counts


#merge the two matrices together
library(URD)
list_matrix <- list(green_subset_matrix, red_subset_matrix)
gated_matrix <- dsCombineDGE(list_matrix)

facs_gated <- CreateSeuratObject(counts = gated_matrix, min.cells = 5, min.features =200)

VlnPlot(facs_smart, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
#facs_smart<- subset(facs_smart, subset = nFeature_RNA > 3000)


# Normalizing the data
facs_gated <- NormalizeData(object = facs_gated, normalization.method = "LogNormalize", scale.factor = 10000)

# Detection of variable genes across the single cells
facs_gated <- FindVariableFeatures(object = facs_gated, do.plot = F)

# Scaling the data and removing unwanted sources of variation
facs_gated <- ScaleData(object = facs_gated, display.progress = F)

# Perform linear dimensional reduction
facs_gated <- RunPCA(object = facs_gated, pc.genes = facs_gated@var.genes, do.print = TRUE)


#jackstraw
facs_gated <- JackStraw(object = facs_gated, num.replicate = 100)
facs_gated <- ScoreJackStraw(object = facs_gated, dims = 1:20)
JackStrawPlot(object = facs_gated, dims = 1:20)


#umapc(1:12,14,17,19)
facs_gated <- FindNeighbors(object = facs_gated, dims = 1:20)
facs_gated <- FindClusters(facs_gated, resolution = 1.0, print.output = 0, save.SNN = T)
facs_gated <-  RunUMAP(facs_gated, reduction.use = "pca", dims= 1:20, min.dist=1, n.neighbors =49 )
DimPlot(object = facs_gated, reduction = "umap",label = T, pt.size = 2)

gated_metadata <- facs_gated@meta.data

write.table(gated_metadata, "~/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/gated_metadata.txt", sep = "\t" )

gated_metadata <- as.matrix(read.table("/Users/JulianKimura/Desktop/neoblast_smartseq/P1_P2_P3_Seurat/gated_metadata.txt", sep = "\t", header = T, row.names = 1))


#changing adding a new column called stages in the metadata
stages <- gated_metadata[,6]
facs_gated@meta.data <- cbind(facs_gated@meta.data, stages)
dim_plot <- DimPlot(object = facs_gated, reduction = "umap", group.by = "stages", pt.size = 2)
dim_plot$data$stages <- factor(x = dim_plot$data$stages, levels = c("R", "G"))
dim_plot

facs_gated@meta.data$stages <- factor(x = facs_gated@meta.data$stages, levels = c("R", "G"))

VlnPlot(facs_gated, features = c("98043523-piwl1-2", "98132150-h10-6", "98029236-pcna"), ncol = 3, group.by = "stages")
RidgePlot(facs_gated, features = "98043523-piwl1-2", group.by = "stages")
