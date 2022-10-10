###########
library(Seurat)
library(magrittr)
library(dplyr)
require(scales)
library(URD)
library(cowplot)
####data was subsetted from the integrated UMAP dataset down to 33000 cells. The count matrix was then transposed in python to extract the list of cells. Next, a subset of the combined, full metadata from the integrated umap dataset was pulled out using a python script. 

#create large matrix of everything
list_embryos <- list(E35_matrix, E50_matrix, E65_matrix, E80_matrix, E95_matrix, E110_matrix,E125_matrix,E145_matrix, juv_matrix)
embryo_juv_matrix <- dsCombineDGE(list_embryos)

# count.data must be provided as a matrix, not a data.frame
library(Matrix)
int_matrix <- read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/URD/subsamped_cleaned_updated_tips/original/45000_subsample_counts.txt")
int_matrix <- as(as.matrix(int_matrix), "dgCMatrix")
embryo_metadata <- read.table("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/URD/subsamped_cleaned_updated_tips/original/updated_tips_metadata.txt")

#editing the metadata so that it only has individual stage names inputted.
current.stage.ids <- c("E351","E352","E353","E501", "E502","E503","E651","E652","E653","E801","E802","E803", "E951","E952","E953" ,"E1101","E1102","E1103","E1251","E1252", "E1253","E1451", "E1452", "E1453", "J1", "J2", "J3", "J4")
new.stage.ids <- c("E35","E35","E35","E50","E50","E50", "E65","E65","E65","E80","E80","E80", "E95","E95","E95", "E110","E110","E110","E125","E125","E125","E145", "E145","E145","J", "J","J", "J")
embryo_metadata$orig.ident<- plyr::mapvalues(x = embryo_metadata$orig.ident, from = current.stage.ids, to = new.stage.ids)

embryo <- createURD(count.data = int_matrix, meta = embryo_metadata, min.cells=5, min.counts=200)


embryo@group.ids$stage <- as.character(embryo@meta[rownames(embryo@group.ids),"orig.ident"])
embryo@group.ids

stages <- sort(unique(embryo@group.ids$stage))
stages

#you are taking the variable genes per stage here. 
var.by.stage <- lapply(seq(1,9,1), function(n) {findVariableGenes(embryo, cells.fit=cellsInCluster(embryo, "stage", stages), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, main.use=paste0("Stages ", stages[n-2], " to ", stages[n]), do.plot=T)})


# Combine the results from each group of stages into a single list of variable genes and load into the URD object
var.genes <- sort(unique(unlist(var.by.stage)))
embryo@var.genes <- var.genes



save.image("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/URD/subsamped_cleaned_updated_tips/original/bro_juv_subsample_pca.RData")


# Calculate PCA and consider those PCs that with standard deviation 2x expected by noise as significant
embryo <- calcPCA(embryo, mp.factor = 2)
embryo <- calcTsne(object = embryo)


embryo <- calcDM(embryo, knn = 181, sigma=NULL)

plotDimArray(embryo, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "Diffusion Map(Sigma NULL, nn=213): Stage", label="stage", plot.title="", legend=F)


root.cells <- cellsInCluster(embryo, "stage", "E35")


embryo.floods <- floodPseudotime(embryo, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)

embryo <- floodPseudotimeProcess(embryo, embryo.floods, floods.name="pseudotime")

pseudotimePlotStabilityOverall(embryo)


embryo.juv <- urdSubset(embryo, cells.keep=cellsInCluster(embryo, "stage", "J"))


embryo.juv@var.genes <- var.by.stage[[9]]

set.seed(20)
embryo.juv <- calcTsne(embryo.juv, which.dims = 1:10, perplexity = 50)


embryo@group.ids[rownames(embryo.juv@group.ids), "tip.clusters"] <- embryo.juv@meta$`res.1.5`


# Determine the parameters of the logistic used to bias the transition probabilities. The procedure
# is relatively robust to this parameter, but the cell numbers may need to be modified for larger
# or smaller data sets.
embryo.ptlogistic <- pseudotimeDetermineLogistic(embryo, "pseudotime", optimal.cells.forward=80, max.cells.back=100, do.plot = T)

# Bias the transition matrix acording to pseudotime
embryo.biased.tm <-as.matrix(pseudotimeWeightTransitionMatrix(embryo, "pseudotime", logistic.params=embryo.ptlogistic))

# Simulate the biased random walks from each tip
embryo.walks <- simulateRandomWalksFromTips(embryo, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = embryo.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000, verbose = F)

# Process the biased random walks into visitation frequencies
embryo <- processRandomWalksFromTips(embryo, embryo.walks, verbose = F)

# Load the cells used for each tip into the URD object
embryo.tree <- loadTipCells(embryo, "tip.clusters")

# Build the tree
embryo.tree <- buildTree(embryo.tree, pseudotime = "pseudotime", tips.use=c(0:19), divergence.method = "preference", cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)

#plot the tree

#rename the stages
current.stage.ids <- c("E35","E50","E65","E80","E95","E110","E125","E145", "J")
new.stage.ids <- c("A-E35","B-E50","C-E65","D-E80","E-E95", "F-E110","G-E125","H-E145", "I-J")
embryo.tree@group.ids$stage<- plyr::mapvalues(x = embryo.tree@group.ids$stage, from = current.stage.ids, to = new.stage.ids)



library(RColorBrewer)
library(colorRamps)
colors <- brewer.pal(n = 9, name = "Set1")


plotTree(embryo.tree, "segment", title="URD Tree Segment", label.segments = TRUE)
plotTree(embryo.tree, "stage", title="URD Tree Stages")


embryo.tree <- treeForceDirectedLayout(embryo.tree, num.nn=213, cut.unconnected.segments=2, verbose=T)

plotTreeForce(embryo.tree, "stage", title = "stage", title.cex = 2, title.line=2.5)

#calling differential expression

#decide what tips to run DE for
tips.to.run <- c("11")
#genes.use <- NULL # Calculate for all genes
genes.use <- embryo.tree@var.genes
# Calculate the markers of each other population.
gene.markers <- list()
for (tipn in 1:length(tips.to.run)) {
  tip <- tips.to.run[tipn]
  print(paste0(Sys.time(), ": ", tip))
  markers <- aucprTestAlongTree(embryo.tree, pseudotime="pseudotime", tips=tip, log.effect.size=0.4, auc.factor = 1.25, max.auc.threshold = 0.1, frac.must.express = 0.1, frac.min.diff = 0, genes.use=genes.use, root="25", only.return.global=F, must.beat.sibs=0.8, report.debug=T, segs.to.skip = NULL)
  saveRDS(markers, file=paste0("/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/URD/subsamped_cleaned_updated_tips/original/", tip, "tip11_root25.rds"))
  gene.markers[[tip]] <- markers
}
write.table(gene.markers[[1]]$diff.exp, file="/Users/JulianKimura/Documents/Lab/Single_Cell_Seq/URD/subsamped_cleaned_updated_tips/original/tip11_root25.txt", sep="\t")

#plot genes on tree
plotTree(embryo.tree, "98048776-hes1", title = "")
plotTree(embryo.tree, "98048063-e2f6", title = "")
plotTree(embryo.tree, "98009983-tbx2-2", title = "")


plotTree(embryo.tree, "98003853-dlx1", title = "")


#########subsetting the tree based off of particular lineages
plotTree(embryo.tree, "segment", title="URD Tree Segment", label.segments = TRUE)
plotTree(embryo.tree, "stage", title="URD Tree Segment", label.segments = TRUE)


branch30_neoblast_subset <- urdSubset(embryo.tree, cells.keep=cellsInCluster(embryo.tree, clustering = "segment", cluster = c("30")))

late_branch30_subset <- urdSubset(branch30_neoblast_subset, cells.keep=cellsInCluster(embryo.tree, clustering = "stage", cluster = c("E110", "E125", "E145", "J")))

#, "F-E110", "G-E125", "H-E145"

plotTree(late_branch30_subset, "stage", title="URD Tree Stages", plot.tree =T)


late_branch30_subset <- as.matrix(late_branch30_subset@count.data)

write.table(late_branch30_subset, "~/Documents/Lab/embryo_singlecell_paper/late_branch30_neoblast_matrix.txt", sep = "\t")




