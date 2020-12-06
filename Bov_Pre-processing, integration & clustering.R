#Load packages
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)

#Load Bovine fetal ovary datasets, sample 1 and sample 2 
Bov1 <- Read10X(data.dir = "~/BovSample1/filtered_feature_bc_matrix")
Bov2 <- Read10X(data.dir = "~/BovSample2/filtered_feature_bc_matrix")

#Crate Seurat object for each sample
Bov1 <- CreateSeuratObject(counts = Bov1, min.cells = 5, min.features = 200)
Bov2 <- CreateSeuratObject(counts = Bov2, min.cells = 5, min.features = 200)

#Add a column to the meta.data to identify each sample
Bov1@meta.data[, "sample"] <- "Sample 1"
Bov2@meta.data[, "sample"] <- "Sample 2"

#Add a column to the meta.data to store percentage of mitochondrial genes
Bov1[["percent.mt"]] <- PercentageFeatureSet(Bov1, pattern = "^MT-")
Bov2[["percent.mt"]] <- PercentageFeatureSet(Bov2, pattern = "^MT-")

#Remove unwanted cells
Bov1 <- subset(Bov1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
Bov2 <- subset(Bov2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)

# Visualize QC metrics
VlnPlot(Bov1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Bov2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Normalize each dataset
Bov1 <- NormalizeData(Bov1, verbose = FALSE)
Bov2 <- NormalizeData(Bov2, verbose = FALSE)

#Identificate highly variable genes
Bov1 <- FindVariableFeatures(Bov1, selection.method = "vst", nfeatures = 2000)
Bov2 <- FindVariableFeatures(Bov2, selection.method = "vst", nfeatures = 2000)

#Find anchors and integrate datasets
anchors <- FindIntegrationAnchors(object.list = list(Bov1, Bov2), dims = 1:30)
Bov <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(Bov) <- "integrated"

#Clustering
Bov <- ScaleData(Bov, verbose = FALSE)
Bov <- RunPCA(Bov, npcs = 30, verbose = FALSE)
Bov <- FindNeighbors(Bov, reduction = "pca", dims = 1:20)
Bov <- FindClusters(Bov, resolution = 0.15)
Bov <- RunUMAP(Bov, reduction = "pca", dims = 1:20)

#Store integrated dataset
saveRDS(Bov, file = "~/BovineIntegrated.rds")