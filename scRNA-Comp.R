#!/usr/bin/RScript
basedir <- "" #as.character(args[1]) # this is the folder with the subfolder "Counts"

# Take 'all' htseq-count results and melt them in to one big dataframe
library(Seurat)
library(SingleCellExperiment)
library(org.Mm.eg.db)
library(scater)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
library(ggplot2)
library(cowplot)
library(scran)

setwd(basedir)

# Get list of count files
cntdir <- paste(basedir, "Counts", sep="/")
pat <- ".count" # count file extension
myfiles <- list.files(path = cntdir,
                      pattern = pat,
                      all.files = TRUE,
                      recursive = FALSE,
                      ignore.case = FALSE,
                      include.dirs = FALSE)

# read each file as array element of DT and rename the last 2 cols
# we created a list of single sample tables
DT <- list()
for (i in 1:length(myfiles) ) {
  infile = paste(cntdir, myfiles[i], sep = "/")
  DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
  cnts <- gsub("(.*).count", "\\1", myfiles[i])
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}

# merge all elements based on first ID columns
data <- DT[[myfiles[1]]]

# we now add each other table with the ID column as key
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(data, y, by = c("ID"))
  data <- z
}

# ID column becomes rownames
rownames(data) <- data$ID
data <- data[,-1]

## add total counts per sample
data <- rbind(data, tot.counts=colSums(data))

# transpose table for readability
data.all.summary <- data[grep("^ENS", rownames(data), perl=TRUE, invert=TRUE), ]

# write summary to file
write.csv(data.all.summary, file = "htseq_counts_all-summary.csv")

data.allgenes <- data[grep("^ENS", rownames(data), perl=TRUE, invert=FALSE), ]

# write data to file
write.csv(data.allgenes, file = "htseq_counts_all.csv")

# setting up the single cell object
sce <- SingleCellExperiment(list(counts=as.matrix(data.allgenes)))

# Mapping ensembl id's to symbol
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL, rowData(sce)$SYMBOL)

# Identifying ID's located on ChrM
location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$ENSEMBL, column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="chrM")

# Read metadata for all experiments and attach to single cell data
metadata <- read.delim("SingleCellBarcodes-Git.txt", check.names=FALSE, header=TRUE)
m <- match(colnames(sce), metadata[["FileName"]])
stopifnot(all(!is.na(m)))
metadata <- metadata[m,]
colData(sce)$Cond <- factor(metadata[["Prep-Cell"]])
colData(sce)$Well <- factor(metadata[["Well"]])
colData(sce)$Prep <- factor(metadata[["Prep"]])
colData(sce)$plate_position <- factor(metadata[["Well2"]])
colData(sce)$Cond2 <- factor(metadata[["Cond2"]])

# QC metrics on mitocondrial genes
mito <- which(rowData(sce)$CHR=="chrM")
sce <- addPerCellQC(sce, subsets=list(Mt=mito))

# QC plots
plotColData(sce, y="sum", x="Prep", show_median=TRUE)
plotColData(sce, y="detected", x="Prep", show_median=TRUE)
plotColData(sce, y="total", x="Prep", show_median=TRUE)
plotColData(sce, y="percent_top_500", x="Prep", show_median=TRUE)
plotColData(sce, y="subsets_Mt_percent", x="Prep", show_median=TRUE)
plotColData(sce, x = "sum", y="detected", colour_by="Prep") 

# Split into separate single cell objects
sce <- computeSumFactors(sce)
sce <- logNormCounts(sce) # see below.

# Normalize all kits separately, then combine and visualize all together
sce.seurat <- as.Seurat(sce)
sce.seurat.list <- SplitObject(sce.seurat, split.by = "Prep")
sce.seurat.list <- sce.seurat.list[c("KIT1","KIT2", "KIT3", "KIT4", "PUBLISHED")]
for (i in 1:length(sce.seurat.list)) {
  sce.seurat.list[[i]] <- NormalizeData(sce.seurat.list[[i]], verbose = FALSE)
  sce.seurat.list[[i]] <- FindVariableFeatures(sce.seurat.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}
sce.seurat.anchors <- FindIntegrationAnchors(object.list = sce.seurat.list, dims = 1:10, k.filter = 25)
sce.seurat.integrated <- IntegrateData(anchorset = sce.seurat.anchors, dims = 1:10)

DefaultAssay(sce.seurat.integrated) <- "integrated"
sce.seurat.integrated <- ScaleData(sce.seurat.integrated, verbose = FALSE)
sce.seurat.integrated <- RunPCA(sce.seurat.integrated, npcs = 10, verbose = FALSE)
sce.seurat.integrated <- RunUMAP(sce.seurat.integrated, reduction = "pca", dims = 1:10)

DimPlot(sce.seurat.integrated, reduction = "umap", group.by = c("Prep", "Cond2"))

sce.seurat <- as.Seurat(sce)
sce.seurat.list <- SplitObject(sce.seurat, split.by = "Prep")
sce.seurat.list <- sce.seurat.list[c("KIT1","KIT2", "KIT3", "KIT4", "PUBLISHED")]

for (i in 1:length(sce.seurat.list)) {
  sce.seurat.list[[i]] <- SCTransform(sce.seurat.list[[i]], verbose = FALSE)
}

sce.seurat.features <- SelectIntegrationFeatures(object.list = sce.seurat.list, nfeatures = 3000)
sce.seurat.list <- PrepSCTIntegration(object.list = sce.seurat.list, anchor.features = sce.seurat.features, verbose = FALSE)
sce.seurat.anchors <- FindIntegrationAnchors(object.list = sce.seurat.list, normalization.method = "SCT", anchor.features = sce.seurat.features, verbose = FALSE, k.filter = 25, dims = 1:10)
sce.seurat.integrated <- IntegrateData(anchorset = sce.seurat.anchors, normalization.method = "SCT", verbose = FALSE)
sce.seurat.integrated <- RunPCA(sce.seurat.integrated, verbose = FALSE)
sce.seurat.integrated <- RunUMAP(sce.seurat.integrated, dims = 1:10)

DimPlot(sce.seurat.integrated, group.by = c("Prep", "Cond2"), combine = FALSE, reduction = "pca")
DimPlot(sce.seurat.integrated, group.by = c("Prep", "Cond2"))
