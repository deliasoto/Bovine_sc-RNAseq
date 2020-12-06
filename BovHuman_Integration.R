#Load packages
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(readxl)
library(tidyverse)

#Load Bovine and Human dataset from Li et al., 2017
Bov1 <- Read10X(data.dir = "~/BovSample1/filtered_feature_bc_matrix")
Bov2 <- Read10X(data.dir = "~/BovSample2/filtered_feature_bc_matrix")
Human <- read.delim("~/HumanDataset/TPM_and_UMI_counts/FGC_umi_counts.xls", row.names = 1)
Human <- as(as.matrix(Human), "dgCMatrix")

#Crate Seurat object for each dataset
Bov1 <- CreateSeuratObject(counts = Bov1, min.cells = 5, min.features = 200)
Bov2 <- CreateSeuratObject(counts = Bov2, min.cells = 5, min.features = 200)
Human<- CreateSeuratObject(counts = Human, min.cells = 5, min.features = 200)

#Add a column to the meta.data to identify the species
Bov1@meta.data[, "species"] <- "Bovine"
Bov2@meta.data[, "species"] <- "Bovine"
Human@meta.data[, "species"] <- "Human"

#Add a column to the Bovine meta.data to store percentage of mitochondrial genes (Human dataset does not have mito genes)
Bov1[["percent.mt"]] <- PercentageFeatureSet(Bov1, pattern = "^MT-")
Bov2[["percent.mt"]] <- PercentageFeatureSet(Bov2, pattern = "^MT-")

#Remove unwanted cells, including male cells
Bov1 <- subset(Bov1, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
Bov2 <- subset(Bov2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
Human <- subset(Human, subset = nFeature_RNA > 2000 & nFeature_RNA < 10000 & nCount_RNA > 100000 & nCount_RNA < 1100000)
Human <- subset(Human, orig.ident=="F")

#Store cell labels from Li et al., 2017
clusters <- read.delim("~/Documents/Analysis/Datasets/HumanDatasets/TPM_and_UMI_counts/Clustering Information.txt", header = 1, sep = "\t", check.names = F)
Human@meta.data <- merge(Human@meta.data, clusters, by.x = 0, by.y = "Cell", sort = F, all = T) %>% mutate(cluster = fct_explicit_na(Cluster, na_level = "Female_NA")) %>% column_to_rownames("Row.names")

#Normalize & identify highly variable genes
Bov1 <- NormalizeData(Bov1, verbose = FALSE)
Bov2 <- NormalizeData(Bov2, verbose = FALSE)
Human <- NormalizeData(Human, verbose = FALSE)
Bov1 <- FindVariableFeatures(Bov1, selection.method = "vst", nfeatures = 2000)
Bov2 <- FindVariableFeatures(Bov2, selection.method = "vst", nfeatures = 2000)
Human <- FindVariableFeatures(Human, selection.method = "vst", nfeatures = 2000)

#Find anchors and integrate datasets
anchors <- FindIntegrationAnchors(object.list = list(Human, Bov1, Bov2), dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

#Cell-cycle scoring and regression
s_cellcyclegenes <- read_excel("~/s_cellcyclegenes.xlsx")
g2m_cellcyclegenes <- read_excel("~/g2m_cellcyclegenes.xlsx")
s.genes <- s_cellcyclegenes$s
g2m.genes <- g2m_cellcyclegenes$g2_m
integrated <- CellCycleScoring(integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
integrated <- ScaleData(integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(integrated))

#PCA & clustering
ElbowPlot(integrated)
integrated <- RunPCA(integrated, npcs = 20, verbose = FALSE)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.5)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20)
integrated <- RunPCA(integrated, npcs = 20, verbose = FALSE)

#Subsetting Human and Bovine germ cells
Germ <- subset(integrated, idents = c(4,8))

