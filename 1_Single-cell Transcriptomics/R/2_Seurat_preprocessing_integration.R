
#start by reading in the data. The Read10X() function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix. 
#The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).
#We next use the count matrix to create a Seurat object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset. 
#For a technical discussion of the Seurat object structure, check out our GitHub Wiki. For example, the count matrix is stored in pbmc[["RNA"]]@counts.
setwd("~/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/Seurat_analysis")

library(dplyr)
library(patchwork)
library(tidyverse)
library(cowplot)
library(viridis)
library(ggplot2)
library(dittoSeq)
library(Seurat)
library(Azimuth)
library(SeuratDisk)
library(SeuratData)
library(viridis)
library(RColorBrewer)
library(DoubletFinder)
library(DropletUtils)
library(dsb)
library(harmony)

#Methods for the Seurat class can be found with the following:
utils::methods(class = 'Seurat')
#Methods for the Assay class can be found with the following:
utils::methods(class = 'Assay')
#Methods for the DimReduc class can be found with the following:
utils::methods(class = 'DimReduc')

# 1- Loading data from 10X multi-modal experiments
BM1.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2586_Test_run_PBMC_BM_lonza_exp1/Files/Sample_BM-1/outs/filtered_feature_bc_matrix/")
BM2.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2586_Test_run_PBMC_BM_lonza_exp1/Files/Sample_BM-2/outs/filtered_feature_bc_matrix/")
BM3.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2586_Test_run_PBMC_BM_lonza_exp1/Files/Sample_BM-3/outs/filtered_feature_bc_matrix/")
BM4.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2586_Test_run_PBMC_BM_lonza_exp1/Files/Sample_BM-4/outs/filtered_feature_bc_matrix/")
BM5.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2622_Test_run_PBMC_BM_lonza_exp2/Files/Sample_BM-5/outs/filtered_feature_bc_matrix/")
BM6.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2622_Test_run_PBMC_BM_lonza_exp2/Files/Sample_BM-6/outs/filtered_feature_bc_matrix/")
BM7.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2622_Test_run_PBMC_BM_lonza_exp2/Files/Sample_BM-7/outs/filtered_feature_bc_matrix/")
BM8.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2622_Test_run_PBMC_BM_lonza_exp2/Files/Sample_BM-8/outs/filtered_feature_bc_matrix/")

# BM-9 = BM1.3
BM34.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective_Study_1/BM_lonza/scRNAseq/2862_Test_run_PBMC_BM_lonza_exp3/Files/Sample_BM-9/outs/filtered_feature_bc_matrix/")

# BM-10 = BM2.3
BM35.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective_Study_1/BM_lonza/scRNAseq/2862_Test_run_PBMC_BM_lonza_exp3/Files/Sample_BM-10/outs/filtered_feature_bc_matrix/")

# BM12 = BM4.3
BM36.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective_Study_1/BM_lonza/scRNAseq/2862_Test_run_PBMC_BM_lonza_exp3/Files/Sample_BM-12/outs/filtered_feature_bc_matrix/")

# BM-13 = BM5.1
BM37.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective_Study_1/BM_lonza/scRNAseq/2862_Test_run_PBMC_BM_lonza_exp3/Files/Sample_BM-13/outs/filtered_feature_bc_matrix/")

# BM15 = BM3 + BM4
BM38.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective_Study_1/BM_lonza/scRNAseq/2862_Test_run_PBMC_BM_lonza_exp3/Files/Sample_BM-15/outs/filtered_feature_bc_matrix/")

# BM16 = BM1= BM4 + BM5
BM39.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective_Study_1/BM_lonza/scRNAseq/2862_Test_run_PBMC_BM_lonza_exp3/Files/Sample_BM-16/outs/filtered_feature_bc_matrix/")


# 1.1- Creates a Seurat object based on the scRNA-seq data

# Initialize the Seurat object with the raw (non-normalized data).
# https://www.biostars.org/p/9468066/
BM1 <- CreateSeuratObject(counts = BM1.data$`Gene Expression`, project = "HD_BM1.1", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM1

BM2 <- CreateSeuratObject(counts = BM2.data$`Gene Expression`, project = "HD_BM2.1", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM2

BM3 <- CreateSeuratObject(counts = BM3.data$`Gene Expression`, project = "HD_BM3.1", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM3

BM4 <- CreateSeuratObject(counts = BM4.data$`Gene Expression`, project = "HD_BM4.1", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM4

BM5 <- CreateSeuratObject(counts = BM5.data$`Gene Expression`, project = "HD_BM1.2", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM5

BM6 <- CreateSeuratObject(counts = BM6.data$`Gene Expression`, project = "HD_BM2.2", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM6

BM7 <- CreateSeuratObject(counts = BM7.data$`Gene Expression`, project = "HD_BM3.2", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM7

BM8 <- CreateSeuratObject(counts = BM8.data$`Gene Expression`, project = "HD_BM4.2", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM8

# Samples from 3rd 10X run
BM34 <- CreateSeuratObject(counts = BM34.data$`Gene Expression`, project = "HD_BM1.3", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM34

BM35 <- CreateSeuratObject(counts = BM35.data$`Gene Expression`, project = "HD_BM2.3", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM35

BM36 <- CreateSeuratObject(counts = BM36.data$`Gene Expression`, project = "HD_BM4.3", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM36

BM37 <- CreateSeuratObject(counts = BM37.data$`Gene Expression`, project = "HD_BM5.1", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM37

BM38 <- CreateSeuratObject(counts = BM38.data$`Gene Expression`, project = "HD_BM3/4", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM38

BM39 <- CreateSeuratObject(counts = BM39.data$`Gene Expression`, project = "HD_BM1/4/5", min.cells = 3, min.features = 0) # min.cells = 3 takes only genes that are expressed in at least 3 cells. 
BM39

# What does data in a count matrix look like?
# Lets examine a few genes in the first thirty cells
# The . values in the matrix represent 0s (no molecules detected). Since most values in an scRNA-seq matrix are 0, 
# Seurat uses a sparse-matrix representation whenever possible. 
# This results in significant memory and speed savings for Drop-seq/inDrop/10x data.
BM1.data$`Gene Expression`[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# 1.2- Create a new assay to store ADT information
rownames(x = BM1.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                     x = rownames(x = BM1.data[["Antibody Capture"]]))
rownames(x = BM2.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                     x = rownames(x = BM2.data[["Antibody Capture"]]))
rownames(x = BM3.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                     x = rownames(x = BM3.data[["Antibody Capture"]]))
rownames(x = BM4.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                     x = rownames(x = BM4.data[["Antibody Capture"]]))
rownames(x = BM5.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                     x = rownames(x = BM5.data[["Antibody Capture"]]))
rownames(x = BM6.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                     x = rownames(x = BM6.data[["Antibody Capture"]]))
rownames(x = BM7.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                     x = rownames(x = BM7.data[["Antibody Capture"]]))
rownames(x = BM8.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                     x = rownames(x = BM8.data[["Antibody Capture"]]))

rownames(x = BM34.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                      x = rownames(x = BM34.data[["Antibody Capture"]]))
rownames(x = BM35.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                      x = rownames(x = BM35.data[["Antibody Capture"]]))
rownames(x = BM36.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                      x = rownames(x = BM36.data[["Antibody Capture"]]))
rownames(x = BM37.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                      x = rownames(x = BM37.data[["Antibody Capture"]]))
rownames(x = BM38.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                      x = rownames(x = BM38.data[["Antibody Capture"]]))
rownames(x = BM39.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "",
                                                      x = rownames(x = BM39.data[["Antibody Capture"]]))

adt_assay1 <- CreateAssayObject(counts = BM1.data$`Antibody Capture`, project = "HD_BM1.1")
adt_assay2 <- CreateAssayObject(counts = BM2.data$`Antibody Capture`, project = "HD_BM2.1")
adt_assay3 <- CreateAssayObject(counts = BM3.data$`Antibody Capture`, project = "HD_BM3.1")
adt_assay4 <- CreateAssayObject(counts = BM4.data$`Antibody Capture`, project = "HD_BM4.1")
adt_assay5 <- CreateAssayObject(counts = BM5.data$`Antibody Capture`, project = "HD_BM1.2")
adt_assay6 <- CreateAssayObject(counts = BM6.data$`Antibody Capture`, project = "HD_BM2.2")
adt_assay7 <- CreateAssayObject(counts = BM7.data$`Antibody Capture`, project = "HD_BM3.2")
adt_assay8 <- CreateAssayObject(counts = BM8.data$`Antibody Capture`, project = "HD_BM4.2")

adt_assay34 <- CreateAssayObject(counts = BM34.data$`Antibody Capture`, project = "HD_BM1.3")
adt_assay35 <- CreateAssayObject(counts = BM35.data$`Antibody Capture`, project = "HD_BM2.3")
adt_assay36 <- CreateAssayObject(counts = BM36.data$`Antibody Capture`, project = "HD_BM4.3")
adt_assay37 <- CreateAssayObject(counts = BM37.data$`Antibody Capture`, project = "HD_BM5.1")
adt_assay38 <- CreateAssayObject(counts = BM38.data$`Antibody Capture`, project = "HD_BM3/4")
adt_assay39 <- CreateAssayObject(counts = BM39.data$`Antibody Capture`, project = "HD_BM1/4/5")


# 1.2.1- Edit ADT assay to remove the HTOs
adt_assay34_ADT <- LayerData( adt_assay34, assay = "counts")
adt_assay34_ADT <- adt_assay34_ADT[-(which(rownames(adt_assay34) %in% c("HTO1", "HTO2", "HTO3", "HTO4", "HTO5"))),]
adt_assay34_ADT <- CreateAssayObject(counts = adt_assay34_ADT, project = "HD_BM1.3") #use this to add in the seurat object as ADT assay

adt_assay35_ADT <- LayerData( adt_assay35, assay = "counts")
adt_assay35_ADT <- adt_assay35_ADT[-(which(rownames(adt_assay35) %in% c("HTO1", "HTO2", "HTO3", "HTO4", "HTO5"))),]
adt_assay35_ADT <- CreateAssayObject(counts = adt_assay35_ADT, project = "HD_BM2.3") #use this to add in the seurat object as ADT assay

adt_assay36_ADT <- LayerData( adt_assay36, assay = "counts")
adt_assay36_ADT <- adt_assay36_ADT[-(which(rownames(adt_assay36) %in% c("HTO1", "HTO2", "HTO3", "HTO4", "HTO5"))),]
adt_assay36_ADT <- CreateAssayObject(counts = adt_assay36_ADT, project = "HD_BM4.3") #use this to add in the seurat object as ADT assay

adt_assay37_ADT <- LayerData( adt_assay37, assay = "counts")
adt_assay37_ADT <- adt_assay37_ADT[-(which(rownames(adt_assay37) %in% c("HTO1", "HTO2", "HTO3", "HTO4", "HTO5"))),]
adt_assay37_ADT <- CreateAssayObject(counts = adt_assay37_ADT, project = "HD_BM5.1") #use this to add in the seurat object as ADT assay

adt_assay38_ADT <- LayerData( adt_assay38, assay = "counts")
adt_assay38_ADT <- adt_assay38_ADT[-(which(rownames(adt_assay38) %in% c("HTO1", "HTO2", "HTO3", "HTO4", "HTO5"))),]
adt_assay38_ADT <- CreateAssayObject(counts = adt_assay38_ADT, project = "HD_BM3/4") #use this to add in the seurat object as ADT assay

adt_assay39_ADT <- LayerData( adt_assay39, assay = "counts")
adt_assay39_ADT <- adt_assay39_ADT[-(which(rownames(adt_assay39) %in% c("HTO1", "HTO2", "HTO3", "HTO4", "HTO5"))),]
adt_assay39_ADT <- CreateAssayObject(counts = adt_assay39_ADT, project = "HD_BM1/4/5") #use this to add in the seurat object as ADT assay

# 1.2.2- Create a HTO assay
adt_assay34_HTO <- LayerData( adt_assay34, assay = "counts")
adt_assay34_HTO <- adt_assay34_HTO[(which(rownames(adt_assay34) %in% c("HTO1", "HTO2", "HTO3", "HTO4", "HTO5"))),]
adt_assay34_HTO <- CreateAssayObject(counts = adt_assay34_HTO, project = "HD_BM1.3") #use this to add in the seurat object as ADT assay

adt_assay35_HTO <- LayerData( adt_assay35, assay = "counts")
adt_assay35_HTO <- adt_assay35_HTO[(which(rownames(adt_assay35) %in% c("HTO1", "HTO2", "HTO3", "HTO4", "HTO5"))),]
adt_assay35_HTO <- CreateAssayObject(counts = adt_assay35_HTO, project = "HD_BM2.3") #use this to add in the seurat object as ADT assay

adt_assay36_HTO <- LayerData( adt_assay36, assay = "counts")
adt_assay36_HTO <- adt_assay36_HTO[(which(rownames(adt_assay36) %in% c("HTO1", "HTO2", "HTO3", "HTO4", "HTO5"))),]
adt_assay36_HTO <- CreateAssayObject(counts = adt_assay36_HTO, project = "HD_BM4.3") #use this to add in the seurat object as ADT assay

adt_assay37_HTO <- LayerData( adt_assay37, assay = "counts")
adt_assay37_HTO <- adt_assay37_HTO[(which(rownames(adt_assay37) %in% c("HTO1", "HTO2", "HTO3", "HTO4", "HTO5"))),]
adt_assay37_HTO <- CreateAssayObject(counts = adt_assay37_HTO, project = "HD_BM5.1") #use this to add in the seurat object as ADT assay

adt_assay38_HTO <- LayerData( adt_assay38, assay = "counts")
adt_assay38_HTO <- adt_assay38_HTO[(which(rownames(adt_assay38) %in% c("HTO3", "HTO4"))),] #remove hashtags that should mot be stained cells here
adt_assay38_HTO <- CreateAssayObject(counts = adt_assay38_HTO, project = "HD_BM3/4") #use this to add in the seurat object as ADT assay

adt_assay39_HTO <- LayerData( adt_assay39, assay = "counts")
adt_assay39_HTO <- adt_assay39_HTO[(which(rownames(adt_assay39) %in% c("HTO1", "HTO4", "HTO5"))),] #remove hashtags that should mot be stained cells here
adt_assay39_HTO <- CreateAssayObject(counts = adt_assay39_HTO, project = "HD_BM1/4/5") #use this to add in the seurat object as ADT assay


# add this assay to the previously created Seurat object
BM1[["ADT"]] <- adt_assay1
BM2[["ADT"]] <- adt_assay2
BM3[["ADT"]] <- adt_assay3
BM4[["ADT"]] <- adt_assay4
BM5[["ADT"]] <- adt_assay5
BM6[["ADT"]] <- adt_assay6
BM7[["ADT"]] <- adt_assay7
BM8[["ADT"]] <- adt_assay8

BM34[["ADT"]] <- adt_assay34_ADT
BM35[["ADT"]] <- adt_assay35_ADT
BM36[["ADT"]] <- adt_assay36_ADT
BM37[["ADT"]] <- adt_assay37_ADT
BM38[["ADT"]] <- adt_assay38_ADT
BM39[["ADT"]] <- adt_assay39_ADT

# add HTO assay
BM34[["HTO"]] <- adt_assay34_HTO
BM35[["HTO"]] <- adt_assay35_HTO
BM36[["HTO"]] <- adt_assay36_HTO
BM37[["HTO"]] <- adt_assay37_HTO
BM38[["HTO"]] <- adt_assay38_HTO
BM39[["HTO"]] <- adt_assay39_HTO


# Validate that the object now contains multiple assays
Assays(BM1)
Assays(BM2)
Assays(BM3)
Assays(BM4)
Assays(BM5)
Assays(BM6)
Assays(BM7)
Assays(BM8)

Assays(BM34)
Assays(BM35)
Assays(BM36)
Assays(BM37)
Assays(BM38)
Assays(BM39)


# 1.3- Metadata addition ----

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats - calculate and add the percentage of MT genes
BM1[["percent.mt"]] <- PercentageFeatureSet(BM1, pattern = "^MT-")
BM2[["percent.mt"]] <- PercentageFeatureSet(BM2, pattern = "^MT-")
BM3[["percent.mt"]] <- PercentageFeatureSet(BM3, pattern = "^MT-")
BM4[["percent.mt"]] <- PercentageFeatureSet(BM4, pattern = "^MT-")
BM5[["percent.mt"]] <- PercentageFeatureSet(BM5, pattern = "^MT-")
BM6[["percent.mt"]] <- PercentageFeatureSet(BM6, pattern = "^MT-")
BM7[["percent.mt"]] <- PercentageFeatureSet(BM7, pattern = "^MT-")
BM8[["percent.mt"]] <- PercentageFeatureSet(BM8, pattern = "^MT-")

BM34[["percent.mt"]] <- PercentageFeatureSet(BM34, pattern = "^MT-")
BM35[["percent.mt"]] <- PercentageFeatureSet(BM35, pattern = "^MT-")
BM36[["percent.mt"]] <- PercentageFeatureSet(BM36, pattern = "^MT-")
BM37[["percent.mt"]] <- PercentageFeatureSet(BM37, pattern = "^MT-")
BM38[["percent.mt"]] <- PercentageFeatureSet(BM38, pattern = "^MT-")
BM39[["percent.mt"]] <- PercentageFeatureSet(BM39, pattern = "^MT-")


# Calculation of the percentage and absolute numbers of human transcripts/UMIs
BM1[["UMI"]] <- BM1$nCount_RNA
BM2[["UMI"]] <- BM2$nCount_RNA
BM3[["UMI"]] <- BM3$nCount_RNA
BM4[["UMI"]] <- BM4$nCount_RNA
BM5[["UMI"]] <- BM5$nCount_RNA
BM6[["UMI"]] <- BM6$nCount_RNA
BM7[["UMI"]] <- BM7$nCount_RNA
BM8[["UMI"]] <- BM8$nCount_RNA

BM34[["UMI"]] <- BM34$nCount_RNA
BM35[["UMI"]] <- BM35$nCount_RNA
BM36[["UMI"]] <- BM36$nCount_RNA
BM37[["UMI"]] <- BM37$nCount_RNA
BM38[["UMI"]] <- BM38$nCount_RNA
BM39[["UMI"]] <- BM39$nCount_RNA


# Calculation of the number of  genes detected per cell
BM1[["Genes"]] <- BM1$nFeature_RNA 
BM2[["Genes"]] <- BM2$nFeature_RNA 
BM3[["Genes"]] <- BM3$nFeature_RNA 
BM4[["Genes"]] <- BM4$nFeature_RNA 
BM5[["Genes"]] <- BM5$nFeature_RNA 
BM6[["Genes"]] <- BM6$nFeature_RNA 
BM7[["Genes"]] <- BM7$nFeature_RNA 
BM8[["Genes"]] <- BM8$nFeature_RNA 

BM34[["Genes"]] <- BM34$nFeature_RNA 
BM35[["Genes"]] <- BM35$nFeature_RNA 
BM36[["Genes"]] <- BM36$nFeature_RNA 
BM37[["Genes"]] <- BM37$nFeature_RNA 
BM38[["Genes"]] <- BM38$nFeature_RNA 
BM39[["Genes"]] <- BM39$nFeature_RNA 

# add metadata - sample_id and experiment_date
BM1$sample_id = "HD_BM1.1"
BM2$sample_id = "HD_BM2.1"
BM3$sample_id = "HD_BM3.1"
BM4$sample_id = "HD_BM4.1"
BM5$sample_id = "HD_BM1.2"
BM6$sample_id = "HD_BM2.2"
BM7$sample_id = "HD_BM3.2"
BM8$sample_id = "HD_BM4.2"

BM34$sample_id = "HD_BM1.3"
BM35$sample_id = "HD_BM2.3"
BM36$sample_id = "HD_BM4.3"
BM37$sample_id = "HD_BM5.1"
BM38$sample_id = "HD_BM3/4"
BM39$sample_id = "HD_BM1/4/5"

BM1$experiment_date = "09_22"
BM2$experiment_date = "09_22"
BM3$experiment_date = "09_22"
BM4$experiment_date = "09_22"
BM5$experiment_date = "10_22"
BM6$experiment_date = "10_22"
BM7$experiment_date = "10_22"
BM8$experiment_date = "10_22"

BM34$experiment_date = "06_23"
BM35$experiment_date = "06_23"
BM36$experiment_date = "06_23"
BM37$experiment_date = "06_23"
BM38$experiment_date = "06_23"
BM39$experiment_date = "06_23"

# 1.4- Note that all operations below are performed on the RNA assay Set and verify that the default assay is RNA
DefaultAssay(BM1) <- "RNA"

# Data information ----
# Different ways to look at the data. Doesn't change it
head(BM1)
dim(BM1) # Gives you the numbers of genes (first number) and cells (second number) in your seurat object

saveRDS(BM1, "BM1.rds")
saveRDS(BM2, "BM2.rds")
saveRDS(BM3, "BM3.rds")
saveRDS(BM4, "BM4.rds")
saveRDS(BM5, "BM5.rds")
saveRDS(BM6, "BM6.rds")
saveRDS(BM7, "BM7.rds")
saveRDS(BM8, "BM8.rds")

saveRDS(BM34, "BM34.rds")
saveRDS(BM35, "BM35.rds")
saveRDS(BM36, "BM36.rds")
saveRDS(BM37, "BM37.rds")
saveRDS(BM38, "BM38.rds")
saveRDS(BM39, "BM39.rds")

BM1 <- readRDS("./RDS_files/Lonza_dataset/BM1.rds")
BM2 <- readRDS("./RDS_files/Lonza_dataset/BM2.rds")
BM3 <- readRDS("./RDS_files/Lonza_dataset/BM3.rds")
BM4 <- readRDS("./RDS_files/Lonza_dataset/BM4.rds")
BM5 <- readRDS("./RDS_files/Lonza_dataset/BM5.rds")
BM6 <- readRDS("./RDS_files/Lonza_dataset/BM6.rds")
BM7 <- readRDS("./RDS_files/Lonza_dataset/BM7.rds")
BM8 <- readRDS("./RDS_files/Lonza_dataset/BM8.rds")

BM34 <- readRDS("./RDS_files/Lonza_dataset/BM34.rds")
BM35 <- readRDS("./RDS_files/Lonza_dataset/BM35.rds")
BM36 <- readRDS("./RDS_files/Lonza_dataset/BM36.rds")
BM37 <- readRDS("./RDS_files/Lonza_dataset/BM37.rds")
BM38 <- readRDS("./RDS_files/Lonza_dataset/BM38.rds")
BM39 <- readRDS("./RDS_files/Lonza_dataset/BM39.rds")

##########################################################################################################################################################################################################################################################################################################################################

# 2- DEMULTIPLEX POOLED SAMPLES with HTODemux and log-ratio (CLR) transformation

#Reference: https://satijalab.org/seurat/articles/hashing_vignette.html

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
BM38 <- readRDS("BM38.rds")
BM38 <- NormalizeData(BM38, assay = "HTO", normalization.method = "CLR") #test dsb normalisation before demultiiplexing
BM38 <- HTODemux(BM38, assay = "HTO", kfunc = "clara", positive.quantile = 0.65)

BM39 <- NormalizeData(BM39, assay = "HTO", normalization.method = "CLR") #test dsb normalisation before demultiiplexing
BM39 <- HTODemux(BM39, assay = "HTO", kfunc = "clara", positive.quantile = 0.9)

# 2.1- Visualize demultiplexing results
# Output from running HTODemux() is saved in the object metadata. 
# We can visualize how many cells are classified as singlets, doublets and negative/ambiguous cells.
# Global classification results
table(BM38$HTO_classification.global) #only 33& of cells were singlets - most cells were negative to these HTOs with quantile = 0.99 >> I used a lower quantile to improve numbers. 62% singlets
table(BM39$HTO_classification.global) # reducing the positive quantile not necessarily increases the % of singlets. 84% singlets

# 2.2- Group cells based on the max HTO signal
Idents(BM38) <- "HTO_maxID"
RidgePlot(BM38, assay = "HTO", features = rownames(BM38[["HTO"]])[1:2], ncol = 2)

Idents(BM39) <- "HTO_maxID"
RidgePlot(BM39, assay = "HTO", features = rownames(BM39[["HTO"]])[1:3], ncol = 3)

# Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
FeatureScatter(BM38, feature1 = "HTO3", feature2 = "HTO4")

# Compare number of UMIs for singlets, doublets and negative cells
Idents(BM38) <- "HTO_classification.global"
VlnPlot(BM38, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# 2.3- Extract the singlets
BM38_singlet <- subset(BM38, idents = "Singlet")

Idents(BM39) <- "HTO_classification.global"
BM39_singlet <- subset(BM39, idents = "Singlet")

# 2.4- Split the samples
Idents(BM38_singlet) <- "hash.ID"
BM3.3 <- subset(BM38_singlet, idents = "HTO3") #1603 cells
BM4.3 <- subset(BM38_singlet, idents = "HTO4") #788 cells

saveRDS(BM3.3, "BM3.3.rds")
saveRDS(BM4.3, "BM4.3.rds")

Idents(BM39_singlet) <- "hash.ID"
BM1.3.2 <- subset(BM39_singlet, idents = "HTO1") #153 cells
BM4.3.2 <- subset(BM39_singlet, idents = "HTO4") #49 cells
BM5.1.2 <- subset(BM39_singlet, idents = "HTO5") #122 cells

saveRDS(BM1.3.2, "BM1.3.2.rds")
saveRDS(BM4.3.2, "BM4.3.2.rds") 
saveRDS(BM5.1.2, "BM5.1.2.rds") 

######################################################################################################################################################################################################################################################################################################################################################################################################

# 3- CITE-seq Normalisation using dsb package - per each sample

# References: https://github.com/niaid/dsb
# https://cran.r-project.org/web/packages/dsb/vignettes/end_to_end_workflow.html
# Normalizing and denoising protein expression data from droplet-based single cell proﬁling
# https://github.com/niaid/dsb_manuscript/#instructions

library(dsb)
DefaultAssay(BM1) <- "ADT"
DefaultAssay(BM2) <- "ADT"
DefaultAssay(BM3) <- "ADT"
DefaultAssay(BM4) <- "ADT"
DefaultAssay(BM5) <- "ADT"
DefaultAssay(BM6) <- "ADT"
DefaultAssay(BM7) <- "ADT"
DefaultAssay(BM8) <- "ADT"

DefaultAssay(BM34) <- "ADT"
DefaultAssay(BM35) <- "ADT"
DefaultAssay(BM36) <- "ADT"
DefaultAssay(BM37) <- "ADT"
DefaultAssay(BM38) <- "ADT"
DefaultAssay(BM39) <- "ADT"

# 3.1- Read raw data using the Seurat function "Read10X" 
BM1_raw = Seurat::Read10X("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2586_Test_run_PBMC_BM_lonza_exp1/Files/Sample_BM-1/outs/raw_feature_bc_matrix")
BM2_raw = Seurat::Read10X("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2586_Test_run_PBMC_BM_lonza_exp1/Files/Sample_BM-2/outs/raw_feature_bc_matrix")
BM3_raw = Seurat::Read10X("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2586_Test_run_PBMC_BM_lonza_exp1/Files/Sample_BM-3/outs/raw_feature_bc_matrix")
BM4_raw = Seurat::Read10X("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2586_Test_run_PBMC_BM_lonza_exp1/Files/Sample_BM-4/outs/raw_feature_bc_matrix")
BM5_raw = Seurat::Read10X("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2622_Test_run_PBMC_BM_lonza_exp2/Files/Sample_BM-5/outs/raw_feature_bc_matrix")
BM6_raw = Seurat::Read10X("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2622_Test_run_PBMC_BM_lonza_exp2/Files/Sample_BM-6/outs/raw_feature_bc_matrix")
BM7_raw = Seurat::Read10X("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2622_Test_run_PBMC_BM_lonza_exp2/Files/Sample_BM-7/outs/raw_feature_bc_matrix")
BM8_raw = Seurat::Read10X("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/2622_Test_run_PBMC_BM_lonza_exp2/Files/Sample_BM-8/outs/raw_feature_bc_matrix")


BM34_raw <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective_Study_1/BM_lonza/scRNAseq/2862_Test_run_PBMC_BM_lonza_exp3/Files/Sample_BM-9/outs/raw_feature_bc_matrix")
BM35_raw <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective_Study_1/BM_lonza/scRNAseq/2862_Test_run_PBMC_BM_lonza_exp3/Files/Sample_BM-10/outs/raw_feature_bc_matrix")
BM36_raw <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective_Study_1/BM_lonza/scRNAseq/2862_Test_run_PBMC_BM_lonza_exp3/Files/Sample_BM-12/outs/raw_feature_bc_matrix")
BM37_raw <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective_Study_1/BM_lonza/scRNAseq/2862_Test_run_PBMC_BM_lonza_exp3/Files/Sample_BM-13/outs/raw_feature_bc_matrix")
BM38_raw <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective_Study_1/BM_lonza/scRNAseq/2862_Test_run_PBMC_BM_lonza_exp3/Files/Sample_BM-15/outs/raw_feature_bc_matrix")
BM39_raw <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective_Study_1/BM_lonza/scRNAseq/2862_Test_run_PBMC_BM_lonza_exp3/Files/Sample_BM-16/outs/raw_feature_bc_matrix")

#cells = Seurat::Read10X("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/scRNAseq/2586_Test_run_PBMC_BM_lonza_exp1/Files/Sample_BM-1/outs/filtered_feature_bc_matrix/")

# 3.2- Define cell-containing barcodes and separate cells and empty drops - defining background signal
BM1_stained_cells = colnames(BM1)
BM1_background = setdiff(colnames(BM1_raw$`Gene Expression`), BM1_stained_cells)

BM2_stained_cells = colnames(BM2)
BM2_background = setdiff(colnames(BM2_raw$`Gene Expression`), BM2_stained_cells)

BM3_stained_cells = colnames(BM3)
BM3_background = setdiff(colnames(BM3_raw$`Gene Expression`), BM3_stained_cells)

BM4_stained_cells = colnames(BM4)
BM4_background = setdiff(colnames(BM4_raw$`Gene Expression`), BM4_stained_cells)

BM5_stained_cells = colnames(BM5)
BM5_background = setdiff(colnames(BM5_raw$`Gene Expression`), BM5_stained_cells)

BM6_stained_cells = colnames(BM6)
BM6_background = setdiff(colnames(BM6_raw$`Gene Expression`), BM6_stained_cells)

BM7_stained_cells = colnames(BM7)
BM7_background = setdiff(colnames(BM7_raw$`Gene Expression`), BM7_stained_cells)

BM8_stained_cells = colnames(BM8)
BM8_background = setdiff(colnames(BM8_raw$`Gene Expression`), BM8_stained_cells)



BM34_stained_cells = colnames(BM34)
BM34_background = setdiff(colnames(BM34_raw$`Gene Expression`), BM34_stained_cells)

BM35_stained_cells = colnames(BM35)
BM35_background = setdiff(colnames(BM35_raw$`Gene Expression`), BM35_stained_cells)

BM36_stained_cells = colnames(BM36)
BM36_background = setdiff(colnames(BM36_raw$`Gene Expression`), BM36_stained_cells)

BM37_stained_cells = colnames(BM37)
BM37_background = setdiff(colnames(BM37_raw$`Gene Expression`), BM37_stained_cells)

BM38_stained_cells = colnames(BM38)
BM38_background = setdiff(colnames(BM38_raw$`Gene Expression`), BM38_stained_cells)

BM39_stained_cells = colnames(BM39)
BM39_background = setdiff(colnames(BM39_raw$`Gene Expression`), BM39_stained_cells)


# 3.3- Split the raw data into separate matrices for RNA and ADT
BM1_prot = BM1_raw$`Antibody Capture`
BM1_rna = BM1_raw$`Gene Expression`

BM2_prot = BM2_raw$`Antibody Capture`
BM2_rna = BM2_raw$`Gene Expression`

BM3_prot = BM3_raw$`Antibody Capture`
BM3_rna = BM3_raw$`Gene Expression`

BM4_prot = BM4_raw$`Antibody Capture`
BM4_rna = BM4_raw$`Gene Expression`

BM5_prot = BM5_raw$`Antibody Capture`
BM5_rna = BM5_raw$`Gene Expression`

BM6_prot = BM6_raw$`Antibody Capture`
BM6_rna = BM6_raw$`Gene Expression`

BM7_prot = BM7_raw$`Antibody Capture`
BM7_rna = BM7_raw$`Gene Expression`

BM8_prot = BM8_raw$`Antibody Capture`
BM8_rna = BM8_raw$`Gene Expression`

BM34_prot = BM34_raw$`Antibody Capture`
BM34_rna = BM34_raw$`Gene Expression`

BM35_prot = BM35_raw$`Antibody Capture`
BM35_rna = BM35_raw$`Gene Expression`

BM36_prot = BM36_raw$`Antibody Capture`
BM36_rna = BM36_raw$`Gene Expression`

BM37_prot = BM37_raw$`Antibody Capture`
BM37_rna = BM37_raw$`Gene Expression`

BM38_prot = BM38_raw$`Antibody Capture`
BM38_rna = BM38_raw$`Gene Expression`

BM39_prot = BM39_raw$`Antibody Capture`
BM39_rna = BM39_raw$`Gene Expression`


# 3.4- Create background adt matrix
BM1_md = data.frame(
  rna.size = log10(Matrix::colSums(BM1_rna)), 
  prot.size = log10(Matrix::colSums(BM1_prot)), 
  n.gene = Matrix::colSums(BM1_rna > 0))

BM2_md = data.frame(
  rna.size = log10(Matrix::colSums(BM2_rna)), 
  prot.size = log10(Matrix::colSums(BM2_prot)), 
  n.gene = Matrix::colSums(BM2_rna > 0))

BM3_md = data.frame(
  rna.size = log10(Matrix::colSums(BM3_rna)), 
  prot.size = log10(Matrix::colSums(BM3_prot)), 
  n.gene = Matrix::colSums(BM3_rna > 0))

BM4_md = data.frame(
  rna.size = log10(Matrix::colSums(BM4_rna)), 
  prot.size = log10(Matrix::colSums(BM4_prot)), 
  n.gene = Matrix::colSums(BM4_rna > 0))

BM5_md = data.frame(
  rna.size = log10(Matrix::colSums(BM5_rna)), 
  prot.size = log10(Matrix::colSums(BM5_prot)), 
  n.gene = Matrix::colSums(BM5_rna > 0))

BM6_md = data.frame(
  rna.size = log10(Matrix::colSums(BM6_rna)), 
  prot.size = log10(Matrix::colSums(BM6_prot)), 
  n.gene = Matrix::colSums(BM6_rna > 0))

BM7_md = data.frame(
  rna.size = log10(Matrix::colSums(BM7_rna)), 
  prot.size = log10(Matrix::colSums(BM7_prot)), 
  n.gene = Matrix::colSums(BM7_rna > 0))

BM8_md = data.frame(
  rna.size = log10(Matrix::colSums(BM8_rna)), 
  prot.size = log10(Matrix::colSums(BM8_prot)), 
  n.gene = Matrix::colSums(BM8_rna > 0))

BM34_md = data.frame(
  rna.size = log10(Matrix::colSums(BM34_rna)), 
  prot.size = log10(Matrix::colSums(BM34_prot)), 
  n.gene = Matrix::colSums(BM34_rna > 0))

BM35_md = data.frame(
  rna.size = log10(Matrix::colSums(BM35_rna)), 
  prot.size = log10(Matrix::colSums(BM35_prot)), 
  n.gene = Matrix::colSums(BM35_rna > 0))

BM36_md = data.frame(
  rna.size = log10(Matrix::colSums(BM36_rna)), 
  prot.size = log10(Matrix::colSums(BM36_prot)), 
  n.gene = Matrix::colSums(BM36_rna > 0))

BM37_md = data.frame(
  rna.size = log10(Matrix::colSums(BM37_rna)), 
  prot.size = log10(Matrix::colSums(BM37_prot)), 
  n.gene = Matrix::colSums(BM37_rna > 0))

BM38_md = data.frame(
  rna.size = log10(Matrix::colSums(BM38_rna)), 
  prot.size = log10(Matrix::colSums(BM38_prot)), 
  n.gene = Matrix::colSums(BM38_rna > 0))

BM39_md = data.frame(
  rna.size = log10(Matrix::colSums(BM39_rna)), 
  prot.size = log10(Matrix::colSums(BM39_prot)), 
  n.gene = Matrix::colSums(BM39_rna > 0))


# 3.5- Add indicator for barcodes Cell Ranger called as cells

BM1_md$drop.class = ifelse(rownames(BM1_md) %in% BM1_stained_cells, 'cell', 'background')
BM2_md$drop.class = ifelse(rownames(BM2_md) %in% BM2_stained_cells, 'cell', 'background')
BM3_md$drop.class = ifelse(rownames(BM3_md) %in% BM3_stained_cells, 'cell', 'background')
BM4_md$drop.class = ifelse(rownames(BM4_md) %in% BM4_stained_cells, 'cell', 'background')
BM5_md$drop.class = ifelse(rownames(BM5_md) %in% BM5_stained_cells, 'cell', 'background')
BM6_md$drop.class = ifelse(rownames(BM6_md) %in% BM6_stained_cells, 'cell', 'background')
BM7_md$drop.class = ifelse(rownames(BM7_md) %in% BM7_stained_cells, 'cell', 'background')
BM8_md$drop.class = ifelse(rownames(BM8_md) %in% BM8_stained_cells, 'cell', 'background')

BM34_md$drop.class = ifelse(rownames(BM34_md) %in% BM34_stained_cells, 'cell', 'background')
BM35_md$drop.class = ifelse(rownames(BM35_md) %in% BM35_stained_cells, 'cell', 'background')
BM36_md$drop.class = ifelse(rownames(BM36_md) %in% BM36_stained_cells, 'cell', 'background')
BM37_md$drop.class = ifelse(rownames(BM37_md) %in% BM37_stained_cells, 'cell', 'background')
BM38_md$drop.class = ifelse(rownames(BM38_md) %in% BM38_stained_cells, 'cell', 'background')
BM39_md$drop.class = ifelse(rownames(BM39_md) %in% BM39_stained_cells, 'cell', 'background')

# 3.6- Remove barcodes with no evidence of capture in the experiment

BM1_md = BM1_md[BM1_md$rna.size > 0 & BM1_md$prot.size > 0, ]
BM2_md = BM2_md[BM2_md$rna.size > 0 & BM2_md$prot.size > 0, ]
BM3_md = BM3_md[BM3_md$rna.size > 0 & BM3_md$prot.size > 0, ]
BM4_md = BM4_md[BM4_md$rna.size > 0 & BM4_md$prot.size > 0, ]
BM5_md = BM5_md[BM5_md$rna.size > 0 & BM5_md$prot.size > 0, ]
BM6_md = BM6_md[BM6_md$rna.size > 0 & BM6_md$prot.size > 0, ]
BM7_md = BM7_md[BM7_md$rna.size > 0 & BM7_md$prot.size > 0, ]
BM8_md = BM8_md[BM8_md$rna.size > 0 & BM8_md$prot.size > 0, ]

BM34_md = BM34_md[BM34_md$rna.size > 0 & BM34_md$prot.size > 0, ]
BM35_md = BM35_md[BM35_md$rna.size > 0 & BM35_md$prot.size > 0, ]
BM36_md = BM36_md[BM36_md$rna.size > 0 & BM36_md$prot.size > 0, ]
BM37_md = BM37_md[BM37_md$rna.size > 0 & BM37_md$prot.size > 0, ]
BM38_md = BM38_md[BM38_md$rna.size > 0 & BM38_md$prot.size > 0, ]
BM39_md = BM39_md[BM39_md$rna.size > 0 & BM39_md$prot.size > 0, ]

library(ggplot2)
# Basic scatter plot
ggplot(BM39_md, aes(x=prot.size, y=rna.size, color=drop.class)) + geom_point() + 
  theme_bw() + scale_color_brewer(palette="Dark2")

BM1_background_drops = rownames(
  BM1_md[ BM1_md$prot.size > 1.5 & 
        BM1_md$prot.size < 3 & 
        BM1_md$rna.size < 2.5, ]
) 
length(BM1_background_drops)

BM2_background_drops = rownames(
  BM2_md[ BM2_md$prot.size > 1.5 & 
            BM2_md$prot.size < 3 & 
            BM2_md$rna.size < 2.5, ]
) 
length(BM2_background_drops)

BM3_background_drops = rownames(
  BM3_md[ BM3_md$prot.size > 1.5 & 
            BM3_md$prot.size < 3 & 
            BM3_md$rna.size < 2.5, ]
) 
length(BM3_background_drops)

BM4_background_drops = rownames(
  BM4_md[ BM4_md$prot.size > 1.5 & 
            BM4_md$prot.size < 3 & 
            BM4_md$rna.size < 2.5, ]
) 
length(BM4_background_drops)

BM5_background_drops = rownames(
  BM5_md[ BM5_md$prot.size > 1.5 & 
            BM5_md$prot.size < 3 & 
            BM5_md$rna.size < 2.5, ]
) 

BM6_background_drops = rownames(
  BM6_md[ BM6_md$prot.size > 1.5 & 
            BM6_md$prot.size < 3 & 
            BM6_md$rna.size < 2.5, ]
) 

BM7_background_drops = rownames(
  BM7_md[ BM7_md$prot.size > 1.5 & 
            BM7_md$prot.size < 3 & 
            BM7_md$rna.size < 2.5, ]
) 

BM8_background_drops = rownames(
  BM8_md[ BM8_md$prot.size > 1.5 & 
            BM8_md$prot.size < 3 & 
            BM8_md$rna.size < 2.5, ]
) 

BM34_background_drops = rownames(
  BM34_md[ BM34_md$prot.size > 1.5 & 
            BM34_md$prot.size < 3 & 
            BM34_md$rna.size < 2.5, ]
) 
length(BM34_background_drops)

BM35_background_drops = rownames(
  BM35_md[ BM35_md$prot.size > 1.5 & 
            BM35_md$prot.size < 3 & 
            BM35_md$rna.size < 2.5, ]
) 
length(BM35_background_drops)

BM36_background_drops = rownames(
  BM36_md[ BM36_md$prot.size > 1.5 & 
            BM36_md$prot.size < 3 & 
            BM36_md$rna.size < 2.5, ]
) 
length(BM36_background_drops)

BM37_background_drops = rownames(
  BM37_md[ BM37_md$prot.size > 1.5 & 
            BM37_md$prot.size < 3 & 
            BM37_md$rna.size < 2.5, ]
) 
length(BM37_background_drops)

BM38_background_drops = rownames(
  BM38_md[ BM38_md$prot.size > 1.5 & 
            BM38_md$prot.size < 3 & 
            BM38_md$rna.size < 2.5, ]
) 
length(BM38_background_drops)

BM39_background_drops = rownames(
  BM39_md[ BM39_md$prot.size > 1.5 & 
            BM39_md$prot.size < 3 & 
            BM39_md$rna.size < 2.5, ]
) 
length(BM39_background_drops)


BM1_background.adt.mtx = as.matrix(BM1_prot[ , BM1_background_drops])
BM2_background.adt.mtx = as.matrix(BM2_prot[ , BM2_background_drops])
BM3_background.adt.mtx = as.matrix(BM3_prot[ , BM3_background_drops])
BM4_background.adt.mtx = as.matrix(BM4_prot[ , BM4_background_drops])
BM5_background.adt.mtx = as.matrix(BM5_prot[ , BM5_background_drops])
BM6_background.adt.mtx = as.matrix(BM6_prot[ , BM6_background_drops])
BM7_background.adt.mtx = as.matrix(BM7_prot[ , BM7_background_drops])
BM8_background.adt.mtx = as.matrix(BM8_prot[ , BM8_background_drops])

BM34_background.adt.mtx = as.matrix(BM34_prot[ , BM34_background_drops])
BM35_background.adt.mtx = as.matrix(BM35_prot[ , BM35_background_drops])
BM36_background.adt.mtx = as.matrix(BM36_prot[ , BM36_background_drops])
BM37_background.adt.mtx = as.matrix(BM37_prot[ , BM37_background_drops])
BM38_background.adt.mtx = as.matrix(BM38_prot[ , BM38_background_drops])
BM39_background.adt.mtx = as.matrix(BM39_prot[ , BM39_background_drops])

# 3.7- Integrating dsb with Seurat object - add dsb normalized matrix "cell.adt.dsb" to the "CITE" data (not counts!) slot 

# renaming rows of the backgroung matrix to match the rownames of the ADT counts matrix of Seurat object
new_names <- c("CD3","CD11b","CD20","CD31","CD34","CD44",
               "CD71", "CD49d", "IgG1", "IgG2a", "HTO1", "HTO2", "HTO3", "HTO4", "HTO5")

# define isotype controls 
isotype.controls = c("IgG1", "IgG2a")

rownames(BM2_background.adt.mtx) <- new_names
BM2_background.adt.mtx <- as.matrix(BM2_background.adt.mtx)

rownames(BM3_background.adt.mtx) <- new_names
BM3_background.adt.mtx <- as.matrix(BM3_background.adt.mtx)

rownames(BM4_background.adt.mtx) <- new_names
BM4_background.adt.mtx <- as.matrix(BM4_background.adt.mtx)

rownames(BM5_background.adt.mtx) <- new_names
BM5_background.adt.mtx <- as.matrix(BM5_background.adt.mtx)

rownames(BM6_background.adt.mtx) <- new_names
BM6_background.adt.mtx <- as.matrix(BM6_background.adt.mtx)

rownames(BM7_background.adt.mtx) <- new_names
BM7_background.adt.mtx <- as.matrix(BM7_background.adt.mtx)

rownames(BM8_background.adt.mtx) <- new_names
BM8_background.adt.mtx <- as.matrix(BM8_background.adt.mtx)

rownames(BM34_background.adt.mtx) <- new_names
BM34_background.adt.mtx <- as.matrix(BM34_background.adt.mtx)

rownames(BM35_background.adt.mtx) <- new_names
BM35_background.adt.mtx <- as.matrix(BM35_background.adt.mtx)

rownames(BM36_background.adt.mtx) <- new_names
BM36_background.adt.mtx <- as.matrix(BM36_background.adt.mtx)

rownames(BM37_background.adt.mtx) <- new_names
BM37_background.adt.mtx <- as.matrix(BM37_background.adt.mtx)

rownames(BM38_background.adt.mtx) <- new_names
BM38_background.adt.mtx <- as.matrix(BM38_background.adt.mtx)

rownames(BM39_background.adt.mtx) <- new_names
BM39_background.adt.mtx <- as.matrix(BM39_background.adt.mtx)

# Generate background matrixes for ADT and HTO- skip this and do as below...
BM34_background.hto.mtx <- BM34_background.adt.mtx[11:15,]
BM34_background.hto.mtx <- as.matrix(BM34_background.hto.mtx)

BM34_background.adt.mtx2 <- BM34_background.adt.mtx[1:10,]
BM34_background.adt.mtx2 <- as.matrix(BM34_background.adt.mtx2)


# normalize and denoise ADT data with dsb with 
BM1_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM1@assays[["ADT"]]@counts, 
  empty_drop_matrix = BM1_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

library("Matrix")
BM1_cells.dsb.norm <- as(BM1_cells.dsb.norm, "dgCMatrix")
#cells.dsb.norm
BM1[["dsb"]] = Seurat::CreateAssayObject(data = BM1_cells.dsb.norm)

# normalize and denoise with dsb with 
BM2_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM2@assays[["ADT"]]@counts, 
  empty_drop_matrix = BM2_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)
 
library("Matrix")
BM2_cells.dsb.norm <- as(BM2_cells.dsb.norm, "dgCMatrix")
BM2[["dsb"]] = Seurat::CreateAssayObject(data = BM2_cells.dsb.norm)

# normalize and denoise with dsb with 
BM3_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM3@assays[["ADT"]]@counts, 
  empty_drop_matrix = BM3_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

library("Matrix")
BM3_cells.dsb.norm <- as(BM3_cells.dsb.norm, "dgCMatrix")
BM3[["dsb"]] = Seurat::CreateAssayObject(data = BM3_cells.dsb.norm)

# normalize and denoise with dsb with 
BM4_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM4@assays[["ADT"]]@counts, 
  empty_drop_matrix = BM4_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

library("Matrix")
BM4_cells.dsb.norm <- as(BM4_cells.dsb.norm, "dgCMatrix")
BM4[["dsb"]] = Seurat::CreateAssayObject(data = BM4_cells.dsb.norm)

# normalize and denoise with dsb with 
BM5_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM5@assays[["ADT"]]@counts, 
  empty_drop_matrix = BM5_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

library("Matrix")
BM5_cells.dsb.norm <- as(BM5_cells.dsb.norm, "dgCMatrix")
BM5[["dsb"]] = Seurat::CreateAssayObject(data = BM5_cells.dsb.norm)

# normalize and denoise with dsb with 
BM6_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM6@assays[["ADT"]]@counts, 
  empty_drop_matrix = BM6_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

library("Matrix")
BM6_cells.dsb.norm <- as(BM6_cells.dsb.norm, "dgCMatrix")
BM6[["dsb"]] = Seurat::CreateAssayObject(data = BM6_cells.dsb.norm)


# normalize and denoise with dsb with 
BM7_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM7@assays[["ADT"]]@counts, 
  empty_drop_matrix = BM7_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

library("Matrix")
BM7_cells.dsb.norm <- as(BM7_cells.dsb.norm, "dgCMatrix")
BM7[["dsb"]] = Seurat::CreateAssayObject(data = BM7_cells.dsb.norm)


# normalize and denoise with dsb with 
BM8_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM8@assays[["ADT"]]@counts, 
  empty_drop_matrix = BM8_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

library("Matrix")
BM8_cells.dsb.norm <- as(BM8_cells.dsb.norm, "dgCMatrix")
BM8[["dsb"]] = Seurat::CreateAssayObject(data = BM8_cells.dsb.norm)



# normalize and denoise ADT data with dsb with 
BM34_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM34@assays[["ADT"]]@counts, 
  empty_drop_matrix = BM34_background.adt.mtx2, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

# Add normalised data to Seurat object. 
# Make sure to add the dsb normalized matrix cell.adt.dsb to the data slot, not the counts slot.
library("Matrix")
BM34_cells.dsb.norm <- as(BM34_cells.dsb.norm, "dgCMatrix")
#cells.dsb.norm
BM34[["ADT_dsb"]] = Seurat::CreateAssayObject(data = BM34_cells.dsb.norm)

# normalize and denoise HTO data with dsb with 
# add Isotype control back for normalisation then remove it afterwards.
adt_assay34_ADT_HTO <- LayerData(adt_assay34, assay = "counts")
adt_assay34_ADT_HTO <- CreateAssayObject(counts = adt_assay34_ADT_HTO, project = "HD_BM1.3") #use this to add in the seurat object as ADT assay
BM34[["ADT_HTO"]] <- adt_assay34_ADT_HTO

BM34_cells.hto.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM34@assays[["ADT_HTO"]]@counts, 
  empty_drop_matrix = BM34_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

# remove ADT data and leave only HTO normalised data
BM34_cells.hto.norm = BM34_cells.hto.norm[11:15,]

library("Matrix")
BM34_cells.hto.norm <- as(BM34_cells.hto.norm, "dgCMatrix")
#cells.dsb.norm
BM34[["HTO_dsb"]] = Seurat::CreateAssayObject(data = BM34_cells.hto.norm)
BM34[["ADT_HTO"]] <- NULL


# normalize and denoise ADT data with dsb with 
adt_assay35_ADT_HTO <- LayerData(adt_assay35, assay = "counts")
adt_assay35_ADT_HTO <- CreateAssayObject(counts = adt_assay35_ADT_HTO, project = "HD_BM2.3") #use this to add in the seurat object as ADT assay
BM35[["ADT_HTO"]] <- adt_assay35_ADT_HTO

BM35_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM35@assays[["ADT_HTO"]]@counts, 
  empty_drop_matrix = BM35_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

# Add normalised data to Seurat object. 
# Make sure to add the dsb normalized matrix cell.adt.dsb to the data slot, not the counts slot.
BM35_cells.adt.norm = BM35_cells.dsb.norm[1:10,]
BM35_cells.hto.norm = BM35_cells.dsb.norm[11:15,]

library("Matrix")
BM35_cells.adt.norm <- as(BM35_cells.adt.norm, "dgCMatrix")
BM35_cells.hto.norm <- as(BM35_cells.hto.norm, "dgCMatrix")
#cells.dsb.norm
BM35[["ADT_dsb"]] = Seurat::CreateAssayObject(data = BM35_cells.adt.norm)
BM35[["HTO_dsb"]] = Seurat::CreateAssayObject(data = BM35_cells.hto.norm)
BM35[["ADT_HTO"]] <- NULL

# normalize and denoise ADT data with dsb with 
adt_assay36_ADT_HTO <- LayerData(adt_assay36, assay = "counts")
adt_assay36_ADT_HTO <- CreateAssayObject(counts = adt_assay36_ADT_HTO, project = "HD_BM4.3") #use this to add in the seurat object as ADT assay
BM36[["ADT_HTO"]] <- adt_assay36_ADT_HTO

BM36_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM36@assays[["ADT_HTO"]]@counts, 
  empty_drop_matrix = BM36_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

# Add normalised data to Seurat object. 
# Make sure to add the dsb normalized matrix cell.adt.dsb to the data slot, not the counts slot.
BM36_cells.adt.norm = BM36_cells.dsb.norm[1:10,]
BM36_cells.hto.norm = BM36_cells.dsb.norm[11:15,]

library("Matrix")
BM36_cells.adt.norm <- as(BM36_cells.adt.norm, "dgCMatrix")
BM36_cells.hto.norm <- as(BM36_cells.hto.norm, "dgCMatrix")
#cells.dsb.norm
BM36[["ADT_dsb"]] = Seurat::CreateAssayObject(data = BM36_cells.adt.norm)
BM36[["HTO_dsb"]] = Seurat::CreateAssayObject(data = BM36_cells.hto.norm)
BM36[["ADT_HTO"]] <- NULL

# normalize and denoise ADT data with dsb with 
adt_assay37_ADT_HTO <- LayerData(adt_assay37, assay = "counts")
adt_assay37_ADT_HTO <- CreateAssayObject(counts = adt_assay37_ADT_HTO, project = "HD_BM5.1") #use this to add in the seurat object as ADT assay
BM37[["ADT_HTO"]] <- adt_assay37_ADT_HTO

BM37_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM37@assays[["ADT_HTO"]]@counts, 
  empty_drop_matrix = BM37_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

# Add normalised data to Seurat object. 
# Make sure to add the dsb normalized matrix cell.adt.dsb to the data slot, not the counts slot.
BM37_cells.adt.norm = BM37_cells.dsb.norm[1:10,]
BM37_cells.hto.norm = BM37_cells.dsb.norm[11:15,]

library("Matrix")
BM37_cells.adt.norm <- as(BM37_cells.adt.norm, "dgCMatrix")
BM37_cells.hto.norm <- as(BM37_cells.hto.norm, "dgCMatrix")
#cells.dsb.norm
BM37[["ADT_dsb"]] = Seurat::CreateAssayObject(data = BM37_cells.adt.norm)
BM37[["HTO_dsb"]] = Seurat::CreateAssayObject(data = BM37_cells.hto.norm)
BM37[["ADT_HTO"]] <- NULL

# normalize and denoise ADT data with dsb with 
adt_assay38_ADT_HTO <- LayerData(adt_assay38, assay = "counts")
adt_assay38_ADT_HTO <- CreateAssayObject(counts = adt_assay38_ADT_HTO, project = "HD_BM3/4") #use this to add in the seurat object as ADT assay
BM38[["ADT_HTO"]] <- adt_assay38_ADT_HTO

BM38_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM38@assays[["ADT_HTO"]]@counts, 
  empty_drop_matrix = BM38_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

# Add normalised data to Seurat object. 
# Make sure to add the dsb normalized matrix cell.adt.dsb to the data slot, not the counts slot.
BM38_cells.adt.norm = BM38_cells.dsb.norm[1:10,]
BM38_cells.hto.norm = BM38_cells.dsb.norm[11:15,]

library("Matrix")
BM38_cells.adt.norm <- as(BM38_cells.adt.norm, "dgCMatrix")
BM38_cells.hto.norm <- as(BM38_cells.hto.norm, "dgCMatrix")
#cells.dsb.norm
BM38[["ADT_dsb"]] = Seurat::CreateAssayObject(data = BM38_cells.adt.norm)
BM38[["HTO_dsb"]] = Seurat::CreateAssayObject(data = BM38_cells.hto.norm)
BM38[["HTO_dsb"]] = Seurat::CreateAssayObject(counts = BM38_cells.hto.norm)
BM38[["ADT_HTO"]] <- NULL

# normalize and denoise ADT data with dsb with 
adt_assay39_ADT_HTO <- LayerData(adt_assay39, assay = "counts")
adt_assay39_ADT_HTO <- CreateAssayObject(counts = adt_assay39_ADT_HTO, project = "HD_BM1/4/5") #use this to add in the seurat object as ADT assay
BM39[["ADT_HTO"]] <- adt_assay39_ADT_HTO

BM39_cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM39@assays[["ADT_HTO"]]@counts, 
  empty_drop_matrix = BM39_background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

# Add normalised data to Seurat object. 
# Make sure to add the dsb normalized matrix cell.adt.dsb to the data slot, not the counts slot.
BM39_cells.adt.norm = BM39_cells.dsb.norm[1:10,]
BM39_cells.hto.norm = BM39_cells.dsb.norm[11:15,]

library("Matrix")
BM39_cells.adt.norm <- as(BM39_cells.adt.norm, "dgCMatrix")
BM39_cells.hto.norm <- as(BM39_cells.hto.norm, "dgCMatrix")
#cells.dsb.norm
BM39[["ADT_dsb"]] = Seurat::CreateAssayObject(data = BM39_cells.adt.norm)
BM39[["HTO_dsb"]] = Seurat::CreateAssayObject(data = BM39_cells.hto.norm)
BM39[["ADT_HTO"]] <- NULL


saveRDS(BM1, "BM1_dsb_normalised.rds")
saveRDS(BM2, "BM2_dsb_normalised.rds")
saveRDS(BM3, "BM3_dsb_normalised.rds")
saveRDS(BM4, "BM4_dsb_normalised.rds")
saveRDS(BM5, "BM5_dsb_normalised.rds")
saveRDS(BM6, "BM6_dsb_normalised.rds")
saveRDS(BM7, "BM7_dsb_normalised.rds")
saveRDS(BM8, "BM8_dsb_normalised.rds")

saveRDS(BM34, "BM34_dsb_normalised.rds")
saveRDS(BM35, "BM35_dsb_normalised.rds")
saveRDS(BM36, "BM36_dsb_normalised.rds")
saveRDS(BM37, "BM37_dsb_normalised.rds")
saveRDS(BM38, "BM38_dsb_normalised.rds")
saveRDS(BM39, "BM39_dsb_normalised.rds")

BM1 <- readRDS("./RDS_files/BM1_dsb_normalised..rds")
BM2 <- readRDS("./RDS_files/BM2_dsb_normalised..rds")
BM3 <- readRDS("./RDS_files/BM3_dsb_normalised..rds")
BM4 <- readRDS("./RDS_files/BM4_dsb_normalised..rds")
BM5 <- readRDS("./RDS_files/BM5_dsb_normalised..rds")
BM6 <- readRDS("./RDS_files/BM6_dsb_normalised..rds")
BM7 <- readRDS("./RDS_files/BM7_dsb_normalised..rds")
BM8 <- readRDS("./RDS_files/BM8_dsb_normalised..rds")

##########################################################################################################################################################################################################################################################################################################################################

# 4- DEMULTIPLEX POOLED SAMPLES with HTODemux and dsb transformation

#Reference: https://satijalab.org/seurat/articles/hashing_vignette.html

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
#BM38 <- readRDS("BM38.rds")
BM38 <- NormalizeData(BM38, assay = "HTO", normalization.method = "CLR") #test dsb normalisation before demultiiplexing
BM38 <- HTODemux(BM38, assay = "HTO", kfunc = "clara", positive.quantile = 0.65)

BM39 <- NormalizeData(BM39, assay = "HTO", normalization.method = "CLR") #test dsb normalisation before demultiiplexing
BM39 <- HTODemux(BM39, assay = "HTO", kfunc = "clara", positive.quantile = 0.9)

# 2.1- Visualize demultiplexing results
# Output from running HTODemux() is saved in the object metadata. 
# We can visualize how many cells are classified as singlets, doublets and negative/ambiguous cells.
# Global classification results
table(BM38$HTO_classification.global) #only 33& of cells were singlets - most cells were negative to these HTOs with quantile = 0.99 >> I used a lower quantile to improve numbers. 62% singlets
table(BM39$HTO_classification.global) # reducing the positive quantile not necessarily increases the % of singlets. 84% singlets

# 2.2- Group cells based on the max HTO signal
Idents(BM38) <- "HTO_maxID"
RidgePlot(BM38, assay = "HTO", features = rownames(BM38[["HTO"]])[1:2], ncol = 2)

Idents(BM39) <- "HTO_maxID"
RidgePlot(BM39, assay = "HTO", features = rownames(BM39[["HTO"]])[1:3], ncol = 3)

# Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
FeatureScatter(BM38, feature1 = "HTO3", feature2 = "HTO4")

# Compare number of UMIs for singlets, doublets and negative cells
Idents(BM38) <- "HTO_classification.global"
VlnPlot(BM38, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# 2.3- Extract the singlets
BM38_singlet <- subset(BM38, idents = "Singlet")

Idents(BM39) <- "HTO_classification.global"
BM39_singlet <- subset(BM39, idents = "Singlet")

# 2.4- Split the samples
Idents(BM38_singlet) <- "hash.ID"
BM3.3 <- subset(BM38_singlet, idents = "HTO3") #1603 cells
BM4.3 <- subset(BM38_singlet, idents = "HTO4") #788 cells

saveRDS(BM3.3, "BM3.3_dsb_normalised.rds")
saveRDS(BM4.3, "BM4.3_dsb_normalised.rds")

Idents(BM39_singlet) <- "hash.ID"
BM1.3.2 <- subset(BM39_singlet, idents = "HTO1") #153 cells
BM4.3.2 <- subset(BM39_singlet, idents = "HTO4") #49 cells
BM5.1.2 <- subset(BM39_singlet, idents = "HTO5") #122 cells

saveRDS(BM1.3.2, "BM1.3.2_dsb_normalised.rds")
saveRDS(BM4.3.2, "BM4.3.2_dsb_normalised.rds") 
saveRDS(BM5.1.2, "BM5.1.2_dsb_normalised.rds") 

######################################################################################################################################################################################################################################################################################################################################################################################################


# 5- Create one merged Seurat object before pre-processing and QC as we will use the same parameters across all samples

# Reference 1: https://satijalab.org/seurat/articles/merge_vignette.html
# Reference 2: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets

# Add dataset labels as cell.ids just in case you have overlapping barcodes between the datasets. 
# Merge datasets into one single seurat object
BM_data <- merge(BM1, c(BM2, BM3, BM4, BM5, BM6, BM7, BM8), add.cell.ids = c("HD_BM1.1", "HD_BM2.1", "HD_BM3.1", "HD_BM4.1",
                                                                             "HD_BM1.2", "HD_BM2.2", "HD_BM3.2", "HD_BM4.2"), project = "HD_BM")

saveRDS(BM_data, "BM_data.rds")
BM_data <- readRDS("./RDS_files/Lonza_dataset/BM_data.rds")

# Once you have created the merged Seurat object, the count matrices and individual count matrices and objects are not needed anymore. 
# It is a good idea to remove them and run garbage collect to free up some memory.
# remove all objects that will not be used.
rm(counts, Patient.object)
# run garbage collect to free up memory
gc()

# Here it is how the count matrix and the metatada look like for every cell.
as.data.frame(BM_data@assays$RNA@counts[1:8, 1:2])
head(BM_data@meta.data, 8)
unique(sapply(X = strsplit(colnames(BM_data), split = "_"), FUN = "[", 1))
table(BM_data$orig.ident)

# 5- Standard pre-processing workflow

# The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. 
# These represent the selection and filtration of cells based on QC metrics, 
# data normalization and scaling, and the detection of highly variable features.

# 5.1- Calculate QC

# Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. 
# A few QC metrics commonly used by the community include:
# The number of unique genes detected in each cell.
# Low-quality cells or empty droplets will often have very few genes
# Cell doublets or multiplets may exhibit an aberrantly high gene count

# Calculate and add the percentage of ribossomal genes
BM_data <- PercentageFeatureSet(BM_data, "^RP[SL]", col.name = "percent_ribo")

# Calculate proportion hemoglobin genes
BM_data <- PercentageFeatureSet(BM_data, "^HB[^(P)]", col.name = "percent_hb")

# Visualize QC metrics as a violin plot before filtering
feats <- c("Genes", "UMI", "percent.mt", "percent_ribo", "percent_hb")
VlnPlot(BM_data, features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(BM_data, feature1 = "UMI", feature2 = "percent.mt")
plot2 <- FeatureScatter(BM_data, feature1 = "UMI", feature2 = "Genes")
plot1 + plot2

# 5.2- Filtering

# 5.2.1- Detection-based filtering and removal of cells with high % of mitochondrial genes (optional: high % of ribossomal genes)

# A standard approach is to filter cells with low amount of reads as well as genes that are present in at least a certain amount of cells. 
# Here we will only consider cells with at least 200 detected genes. 
# Please note that those values are highly dependent on the library preparation method used.
BM_subset <- subset(BM_data, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
head(BM_subset)
dim(BM_subset) # Gives you the numbers of genes (first number) and cells (second number) in your seurat object
head(BM_subset) 
unique(sapply(X = strsplit(colnames(BM_subset), split = "_"), FUN = "[", 1))
table(BM_subset$orig.ident)

# 5.2.2- Filter genes

# First see which genes contribute the most to such reads. We can for instance plot the percentage of counts per gene.
# Compute the relative expression of each gene per cell Use sparse matrix
# operations, if your dataset is large, doing matrix devisions the regular way will take a very long time.
par(mar = c(4, 8, 2, 1))
C <- BM_subset@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE) #As the level of expression of mitochondrial and MALAT1 genes are judged as mainly technical, it can be wise to remove them from the dataset bofore any further analysis.

# Filter MALAT1
# BM_subset <- BM_subset[!grepl("MALAT1", rownames(BM_subset)), ] #skip this, for whatever reason it removes rhe ADT data

# Filter Mitocondrial
# BM_subset2 <- BM_subset2[!grepl("^MT-", rownames(BM_subset)), ]
# Filter Ribossomal gene (optional if that is a problem on your data) data.filt
# <- data.filt[ ! grepl('^RP[SL]', rownames(data.filt)), ]
# Filter Hemoglobin gene (optional if that is a problem on your data)
#data.filt <- data.filt[!grepl("^HB[^(P)]", rownames(data.filt)), ]
dim(BM_subset)

# Plot filtered QC metrics

feats <- c("Genes", "UMI", "percent.mt", "percent_ribo", "percent_hb")
VlnPlot(BM_subset, features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
plot1 <- FeatureScatter(BM_subset, feature1 = "UMI", feature2 = "percent.mt")
plot2 <- FeatureScatter(BM_subset, feature1 = "UMI", feature2 = "Genes")
plot1 + plot2

saveRDS(BM_subset, "BM_subset.rds")

# 5.2.3- Sample sex

genes.file = "./genes.table.csv"

if (!file.exists(genes.file)) {
  suppressMessages(require(biomaRt))
  
  # initialize connection to mart, may take some time if the sites are
  # unresponsive.
  mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  
  # fetch chromosome info plus some other annotations
  genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name",
                                                   "description", "gene_biotype", "chromosome_name", "start_position"), mart = mart,
                                    useCache = F))
  
  if (!dir.exists("data/results")) {
    dir.create("data/results")
  }
  if (is.data.frame(genes.table)) {
    write.csv(genes.table, file = genes.file)
  }
  
  if (!file.exists(genes.file)) {
    download.file("https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/master/labs/misc/genes.table.csv",
                  destfile = "data/results/genes.table.csv")
    genes.table = read.csv(genes.file)
  }
  
} else {
  genes.table = read.csv(genes.file)
}

genes.table <- genes.table[genes.table$external_gene_name %in% rownames(BM_subset),
]

# Now that we have the chromosome information, we can calculate per cell the proportion of reads that comes from chromosome Y.

chrY.gene = genes.table$external_gene_name[genes.table$chromosome_name == "Y"]
BM_subset$pct_chrY = colSums(BM_subset@assays$RNA@counts[chrY.gene, ])/colSums(BM_subset@assays$RNA@counts)
# Then plot XIST expression vs chrY proportion. 
# As you can see, the samples are clearly on either side, even if some cells do not have detection of either.
FeatureScatter(BM_subset, feature1 = "XIST", feature2 = "pct_chrY")
VlnPlot(BM_subset, features = c("XIST", "pct_chrY"))

# 6- Normalize the data

# After removing unwanted cells from the dataset, the next step is to normalize the data. 
# By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
# Normalized values are stored in BM1_subset[["RNA"]]@data.

BM_subset <- NormalizeData(BM_subset, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

# 6.1- Identification of highly variable features (feature selection)

# We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others).
# We and others have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.
# By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.

BM_subset <- FindVariableFeatures(BM_subset, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(BM_subset), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(BM_subset)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# 6.2- Scaling the data

# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(BM_subset)
BM_subset <- ScaleData(BM_subset, features = all.genes)

# 7- Predict and remove doublets with DoubletFinder - should be done after UMAP and before clustering and run this per sample.

# Reference: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets
# Github page https://github.com/chris-mcginnis-ucsf/DoubletFinder

# Most doublet detectors simulates doublets by merging cell counts and predicts doublets as cells that have similar embeddings as the simulated doublets. 
# Most such packages need an assumption about the number/proportion of expected doublets in the dataset. 
# The data you are using is subsampled, but the orignial datasets contained about 5 000 cells per sample, hence we can assume that they loaded about 9 000 cells and should have a doublet rate at about 4%.

#OBS: Ideally doublet prediction should be run on each sample separately, especially if your different samples have different proportions of celltypes. 
# Here, we will use DoubletFinder to predict doublet cells. 
# But before doing doublet detection we need to run scaling, variable gene selection and pca, as well as UMAP for visualization.

# run doubletFinder, selecting first 10 PCs and a pK value of 0.9. 
# To optimize the parameters, you can run the paramSweep function in the package.

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
BM_subset <- RunPCA(BM_subset)
BM_subset <- RunUMAP(BM_subset, dims = 1:30)

sweep.res.list_BM <- paramSweep_v3(BM_subset, PCs = 1:10, sct = FALSE)
sweep.stats_BM <- summarizeSweep(sweep.res.list_BM, GT = FALSE)
bcmvn_BM <- find.pK(sweep.stats_BM)
barplot(bcmvn_BM$BCmetric, names.arg = bcmvn_BM$pK, las=2) # pK = 0.01

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
sweep.res.list_BM <- paramSweep_v3(BM_subset, PCs = 1:10, sct = FALSE)
gt.calls <- BM_subset@meta.data[rownames(sweep.res.list_BM[[1]]), "GT"]   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
sweep.stats_BM <- summarizeSweep(sweep.res.list_BM, GT = TRUE, GT.calls = gt.calls)
bcmvn_BM <- find.pK(sweep.stats_BM)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
#homotypic.prop <- modelHomotypic(BM_subset@meta.data$predicted.celltype.l2)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.004*nrow(BM_subset@meta.data))  ## Assuming 4% doublet formation rate - tailor for your dataset
#nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
BM_subset <- doubletFinder_v3(BM_subset, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi) #Run 1
#BM_query <- doubletFinder_v3(BM_query, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.001_25") #Run 2

# Plot results - name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(BM_subset@meta.data)[grepl("DF.classification", colnames(BM_subset@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(BM_subset, group.by = "orig.ident") + NoAxes(),
                   DimPlot(BM_subset, group.by = DF.name) + NoAxes())

VlnPlot(BM_subset, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

# Remove all predicted doublets from our data.
dim(BM_subset[, BM_subset@meta.data[, DF.name] == "Doublet"])
BM_subset = BM_subset[, BM_subset@meta.data[, DF.name] == "Singlet"]
dim(BM_subset)

# 8- Detecting and removing empty droplets with DoubletUtils - run this per sample before merging or integration?

# Reference: https://bioconductor.org/packages/devel/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html#detecting-empty-droplets
# https://rockefelleruniversity.github.io/scRNA-seq/exercises/answers/exercise2_answers.html

# Empty droplets often contain RNA from the ambient solution, resulting in non-zero counts after debarcoding. 
# The emptyDrops function is designed to distinguish between empty droplets and cells. It does so by testing each barcode’s expression profile for significant deviation from the ambient profile. 
# Given a matrix my.counts containing UMI counts for all barcodes, we call:

# To use DropletUtils convert Seurat object to SCE object
BM.sce <- as.SingleCellExperiment(BM_subset, assay = "RNA")

# Draw a knee plot and identify inflection point and knee point.
bcrank <- barcodeRanks(counts(BM.sce))
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)

legend("bottomleft", legend=c("Inflection", "Knee"), 
       col=c("darkgreen", "dodgerblue"), lty=2, cex=1.2)

# Identify non-empty droplets and compare to the results with hard cut-off
# set limit as “100” (filter out droplets with UMI counts less than 100)
# FDR cut-off for non-empty droplet is 0.01 - using a UMI threshold of 100 and FDR of 1%.

set.seed(100)
limit <- 100   
e.out <- emptyDrops(counts(BM.sce),lower=limit, test.ambient=TRUE) # it did not identified empty droplets or ambient RNA

e.out
summary(e.out$FDR <= 0.01)
is.cell <- e.out$FDR <= 0.01
sum(is.cell, na.rm=TRUE) #number of true droplets with cells

# The p-values are calculated by permutation testing, hence the need to set a seed. 
# The Limited field indicates whether a lower p-value could be obtained by increasing the number of permutations. 
# If there are any entries with FDR above the desired threshold and Limited==TRUE, it indicates that npts should be increased in the emptyDrops call.
table(Limited=e.out$Limited, Significant=is.cell)

# Concordance by testing with FDR and limited
table(Sig=e.out$FDR <= 0.001, Limited=e.out$Limited)

# Removing empty droplets
BM.sce2 <- BM.sce[,which(e.out$FDR <= 0.01)]

# 9- Calculate and regress out cell-cycle scores
# To score a gene list, the algorithm calculates the difference of mean expression of the given list and the mean expression of reference genes. 
# To build the reference, the function randomly chooses a bunch of genes matching the distribution of the expression of the given list. 
# Cell cycle scoring adds three slots in data, a score for S phase, a score for G2M phase and the predicted cell cycle phase.
# Mitigate the effects of cell cycle heterogeneity in scRNA-seq data by calculating cell cycle phase scores based on canonical markers, and regressing these out of the data during pre-processing. 

BM_subset <- CellCycleScoring(object = BM_subset, g2m.features = cc.genes$g2m.genes,
                              s.features = cc.genes$s.genes)

VlnPlot(BM_subset, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
        ncol = 2, pt.size = 0.1)

# Visualize the distribution of cell cycle markers
RidgePlot(BM_subset, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), group.by = "orig.ident", ncol = 2)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
BM_subset<- RunPCA(BM_subset, features = c(cc.genes$s.genes, cc.genes$g2m.genes))
DimPlot(BM_subset, reduction = "pca", 
        group.by = "Phase")

# Quantifications
# Code to quantify the numbers of cells in each cell cycle and the proportion of cells 
table(BM_subset$Phase)
prop.table(table(BM_subset$Phase))

# Optional: we suggest regressing out the difference between the G2M and S phase scores. 
# This means that signals separating non-cycling cells and cycling cells will be maintained, 
# but differences in cell cycle phase among proliferating cells (which are often uninteresting), will be regressed out of the data
BM_subset$CC.Difference <- BM_subset$S.Score - BM_subset$G2M.Score
BM_subset <- ScaleData(BM_subset, vars.to.regress = "CC.Difference", features = rownames(BM_subset))


# 10- SCTransform

# Different algorithm to normalise and scale the data, the results are stored in a different assay, so our previous normalisation is not overwritten
# Biological heterogeneity in single-cell RNA-seq data is often confounded by technical factors including sequencing depth.
# This procedure omits the need for heuristic steps including pseudocount addition or log-transformation and improves common downstream analytical tasks 
# such as variable gene selection, dimensional reduction, and differential expression.
# Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
# During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage
# The latest version of sctransform also supports using glmGamPoi package which substantially improves the speed of the learning procedure. It can be invoked by specifying method="glmGamPoi".

library(sctransform)
BM_subset <- SCTransform(BM_subset, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

saveRDS(BM_subset, "BM_subset.rds")
BM_subset <- readRDS("./RDS_files/BM_subset.rds")

# 10.1- Perform linear dimensional reduction

# Calculate PCs using variable features determined by SCTransform (3000 by default)
BM_subset <- RunPCA(BM_subset, assay = "SCT", npcs = 50)

# Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()
VizDimLoadings(BM_subset, dims = 1:2, reduction = "pca")
DimPlot(BM_subset, reduction = "pca")

# In particular DimHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. 
# Both cells and features are ordered according to their PCA scores. 
# Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. 
# Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated feature sets.
DimHeatmap(BM_subset, dims = 1, cells = 1000, balanced = TRUE)
DimHeatmap(BM_subset, dims = 1:15, cells = 500, balanced = TRUE)

# 10.2- Determine the ‘dimensionality’ of the dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
BM_subset <- JackStraw(BM_subset, num.replicate = 100) #cannot be run on SCTransform-normalized data.
BM_subset <- ScoreJackStraw(BM_subset, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(BM_subset, ndims = 50)

# 10.3- Run non-linear dimensional reduction (UMAP/tSNE) before integration

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
BM_subset <- RunUMAP(BM_subset, dims = 1:30)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual samples
DimPlot(BM_subset, reduction = "umap", label = FALSE) #try other ploting packages and explore the function

dittoDimPlot(BM_subset, var = "sample_id",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

# individual 10X runs
dittoDimPlot(BM_subset, var = "experiment_date",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("10X runs on UMAP")

saveRDS(BM_subset, "BM_subset.rds")

# 11- Harmony integration
BM_integrated <- RunHarmony(BM_subset, 
                                group.by.vars = c("sample_id", "experiment_date"), 
                                reduction = "pca", assay.use = "SCT", reduction.save = "harmony")

BM_integrated <- RunUMAP(BM_integrated, reduction = "harmony", assay = "SCT", dims = 1:40)

# individual samples
DimPlot(BM_integrated, reduction = "umap", label = FALSE) #try other ploting packages and explore the function

dittoDimPlot(BM_integrated, var = "sample_id",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

# individual 10X runs
dittoDimPlot(BM_integrated, var = "experiment_date",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("10X runs on UMAP")

saveRDS(BM_integrated, "BM_integrated.rds")
save.image(file = "BM_Lonza_Preprocessing_Merging__Integration.RData")




