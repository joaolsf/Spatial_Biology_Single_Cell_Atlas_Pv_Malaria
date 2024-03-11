
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

# add ADT assay to the previously created Seurat object
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

# 1.4- Note that all operations below are performed on the RNA assay Set and verify that the
# default assay is RNA
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

BM1 <- readRDS("BM1.rds")
BM2 <- readRDS("BM2.rds")
BM3 <- readRDS("BM3.rds")
BM4 <- readRDS("BM4.rds")
BM5 <- readRDS("BM5.rds")
BM6 <- readRDS("BM6.rds")
BM7 <- readRDS("BM7.rds")
BM8 <- readRDS("BM8.rds")

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

# 3- Create one merged Seurat object before pre-processing and QC as we will use the same parameters across all samples

# Reference 1: https://satijalab.org/seurat/articles/merge_vignette.html
# Reference 2: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets

# Add dataset labels as cell.ids just in case you have overlapping barcodes between the datasets. 
# Merge datasets into one single seurat object
BM_data <- merge(BM1, c(BM2, BM3, BM4, BM5, BM6, BM7, BM8), add.cell.ids = c("HD_BM1.1", "HD_BM2.1", "HD_BM3.1", "HD_BM4.1",
                                                                             "HD_BM1.2", "HD_BM2.2", "HD_BM3.2", "HD_BM4.2"), project = "HD_BM")

saveRDS(BM_data, "BM_data.rds")

# Once you have created the merged Seurat object, the count matrices and individual count matrices and objects are not needed anymore. 
# It is a good idea to remove them and run garbage collect to free up some memory.
# remove all objects that will not be used.
rm(BM1, BM1.data, BM2, BM2.data, BM3, BM3.data, BM4, BM4.data, BM5, BM5.data, BM6, BM6.data, BM7, BM7.data, BM8, BM8.data, adt_assay1, adt_assay2,
   adt_assay3, adt_assay4, adt_assay5, adt_assay6, adt_assay7, adt_assay8)
# run garbage collect to free up memory
gc()

# Here it is how the count matrix and the metatada look like for every cell.
as.data.frame(BM_data@assays$RNA@counts[1:8, 1:2])
head(BM_data@meta.data, 8)
unique(sapply(X = strsplit(colnames(BM_data), split = "_"), FUN = "[", 1))
table(BM_data$orig.ident)


# 4- Standard pre-processing workflow

# The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. 
# These represent the selection and filtration of cells based on QC metrics, 
# data normalization and scaling, and the detection of highly variable features.

# 4.1- Calculate QC

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

# 4.2- Filtering

# 4.2.1- Detection-based filtering and removal of cells with high % of mitochondrial genes (optional: high % of ribossomal genes)

# A standard approach is to filter cells with low amount of reads as well as genes that are present in at least a certain amount of cells. 
# Here we will only consider cells with at least 200 detected genes. 
# Please note that those values are highly dependent on the library preparation method used.
BM_subset <- subset(BM_data, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
head(BM_subset)
dim(BM_subset) # Gives you the numbers of genes (first number) and cells (second number) in your seurat object
head(BM_subset) 
unique(sapply(X = strsplit(colnames(BM_subset), split = "_"), FUN = "[", 1))
table(BM_subset$orig.ident)

# 4.2.2- Filter genes

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
BM_subset <- BM_subset[!grepl("MALAT1", rownames(BM_subset)), ]

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

# 4.2.3- Sample sex

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

# 5- Normalize the data

# After removing unwanted cells from the dataset, the next step is to normalize the data. 
# By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
# Normalized values are stored in BM1_subset[["RNA"]]@data.

BM_subset <- NormalizeData(BM_subset, normalization.method = "LogNormalize", scale.factor = 10000)

# 5.1- Identification of highly variable features (feature selection)

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

# 5.2- Scaling the data

# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in pbmc[["RNA"]]@scale.data

all.genes <- rownames(BM_subset)
BM_subset <- ScaleData(BM_subset, features = all.genes)

# 5.3- SCTransform

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

# 5.4- Calculate cell-cycle scores
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

# 6- Perform linear dimensional reduction

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

# 6.1- Determine the ‘dimensionality’ of the dataset

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
BM_subset <- JackStraw(BM_subset, num.replicate = 100) #cannot be run on SCTransform-normalized data.
BM_subset <- ScoreJackStraw(BM_subset, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

ElbowPlot(BM_subset, ndims = 50)

# 6.2- Run non-linear dimensional reduction (UMAP/tSNE) before integration

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

# 7- Harmony integration
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
save.image(file = "BM_Lonza_Merging_Preprocessing_Integration.RData")

# 8- Clustering analysis

BM_integrated <- FindNeighbors(BM_integrated, reduction = "harmony")
BM_integrated <- FindClusters(BM_integrated, resolution = 1)

# Look at cluster IDs of the first 5 cells
head(Idents(BM_subset), 5)

# UMAP of the identified clusters
DimPlot(BM_integrated, reduction = "umap", label = FALSE) #try other ploting packages and explore the function

dittoDimPlot(BM_integrated, var = "seurat_clusters",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Seurat clusters on UMAP")

saveRDS(BM_integrated, "BM_integrated.rds")

# 8.1- Finding differentially expressed features (cluster biomarkers) - marker genes

# find markers that define clusters via differential expression. 
# By default, it identifies positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. 
# FindAllMarkers() automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, 
# and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. 
# You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are unlikely to be highly discriminatory. 
# As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. 
# While there is generally going to be a loss in power, the speed increases can be significant and the most highly differentially expressed features will likely still rise to the top.

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(BM1_subset) <- "SCT" 
BM1.markers <- FindAllMarkers(BM1_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
clipr::write_clip(BM1_subset.marker) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

BM1.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

# find all markers of cluster 2
cluster2.markers <- FindMarkers(BM1_subset, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 12 from clusters 16 and 3
cluster12.markers <- FindMarkers(BM1_subset, ident.1 = 12, ident.2 = c(16, 3), min.pct = 0.25)
head(cluster12.markers, n = 5)

# Seurat has several tests for differential expression which can be set with the test.use parameter (see our DE vignette for details). 
# For example, the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect).
cluster0.markers <- FindMarkers(BM1_subset, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# We include several tools for visualizing marker expression. 
# VlnPlot() (shows expression probability distributions across clusters), and 
# FeaturePlot() (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations.
# We also suggest exploring RidgePlot(), CellScatter(), and DotPlot() as additional methods to view your dataset.

VlnPlot(BM1_subset, features = c("TFRC", "GYPA"))
VlnPlot(BM1_subset, features = c("nCount_RNA"))

# you can plot raw counts as well
VlnPlot(BM1_subset, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(BM1_subset, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
BM1.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
# changing the default color
DoHeatmap(BM1_subset, features = top10$gene, size = 2) + NoLegend() + theme(axis.text.y = element_text(size = 3))
# changing the default color
DoHeatmap(BM1_subset, features = top10$gene, size = 2) + NoLegend() + theme(axis.text.y = element_text(size = 3)) +
    scale_fill_viridis()

# Plots with dittoSeq
# UMAP colored by cell type and expression - dittoDimPlot
p1 <- dittoDimPlot(BM1_subset, var = "seurat_clusters",
                   reduction.use = "umap", size = 0.5,
                   do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")
p1

# Heatmap visualization - DittoHeatmap
dittoHeatmap(BM1_subset, genes = top10$gene,
             assay = "SCT", order.by = c("seurat_clusters"), cluster_rows = FALSE,
             cluster_cols = FALSE, scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("seurat_clusters"), fontsize_row = 4)

# Quantifications
# Code to quantify the numbers of cells in each cluster and the proportion of cells 

table(BM1_subset$seurat_clusters)
prop.table(table(BM1_subset$seurat_clusters))
tapply(BM1_subset$nCount_RNA, BM1_subset$seurat_clusters,  median) #median number of RNA molecules per cluster
# Explore more Dittoseq plotting functions, such as bar plots for cell numbers


# 8.2.1- Marker gene Identification using RNA slot 
# Calculates the genes upregulated in each cluster. It's important that you switch to the RNA assay for that   

DefaultAssay(BM1_subset) <- "RNA" 
BM1_subset.marker2 <- FindAllMarkers(object = BM1_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
clipr::write_clip(BM1_subset.marker2) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# Generate a list of marker genes that you would like to visualise
Marker <- c("TFRC", "GYPA")  

# Different plotting options. Feel free to play around
DotPlot(BM1_subset, features = Marker)
VlnPlot(BM1_subset, features = Marker)
FeaturePlot(BM1_subset, Marker, cols = viridis(100, direction = -1)) # cols changes the colors of the plot. I really like viridis for those heatmaps


# 8- Mapping and annotating query datasets

# 8.1- Using the output from the Azimuth app - https://azimuth.hubmapconsortium.org

BM1_subset[["cell_type"]] <- BM1_azimuth[["pred.df"]][["predicted.celltype.l2"]]
p1 <- dittoDimPlot(BM1_subset, var = "cell_type",
                   reduction.use = "umap", size = 0.5,
                   do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")
p1

p2 <- dittoBarPlot(BM1_subset, var = "cell_type", group.by = "orig.ident")
p2

# 8.2- Run Azimuth Locally

# 5.2.1- Download the Azimuth reference and extract the archive
# Load the reference - it has to be a list object!
# Change the file path based on where the reference is located on your system.
reference <- LoadReference(path = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/scRNAseq/2586_Test_run_PBMC_BM_lonza_exp1/Files/Sample_BM-1/Analysis/reference") 

# 5.2.2- Load the query object for mapping
# Change the file path based on where the query file is located on your system.
query <- LoadFileInput(path = "BM1.rds") # it is already loaded
query <- BM1_subset

query <- ConvertGeneNames(
  object = query,
  reference.names = rownames(x = reference$map),
  homolog.table = 'https://seurat.nygenome.org/azimuth/references/homologs.rds'
)

# 5.2.3- Calculate nCount_RNA and nFeature_RNA if the query does not contain them already
if (!all(c("nCount_RNA", "nFeature_RNA") %in% c(colnames(x = query[[]])))) {
  calcn <- as.data.frame(x = Seurat:::CalcN(object = query))
  colnames(x = calcn) <- paste(
    colnames(x = calcn),
    "RNA",
    sep = '_'
  )
  query <- AddMetaData(
    object = query,
    metadata = calcn
  )
  rm(calcn)
}

# 5.2.4- Calculate percent mitochondrial genes if the query contains genes
# matching the regular expression "^MT-"
if (any(grepl(pattern = '^MT-', x = rownames(x = query)))) {
  query <- PercentageFeatureSet(
    object = query,
    pattern = '^MT-',
    col.name = 'percent.mt',
    assay = "RNA"
  )
}

# 5.2.5- Filter cells based on the thresholds for nCount_RNA and nFeature_RNA
# you set in the app
cells.use <- query[["nCount_RNA", drop = TRUE]] <= 206031 &
  query[["nCount_RNA", drop = TRUE]] >= 500 &
  query[["nFeature_RNA", drop = TRUE]] <= 7500 &
  query[["nFeature_RNA", drop = TRUE]] >= 200

# If the query contains mitochondrial genes, filter cells based on the thresholds for percent.mt you set in the app
if ("percent.mt" %in% c(colnames(x = query[[]]))) {
  cells.use <- cells.use & (query[["percent.mt", drop = TRUE]] <= 10 &
                              query[["percent.mt", drop = TRUE]] >= 0)
}

# 5.2.6- Remove filtered cells from the query - DO NOT SKIP THIS STEP
query <- query[, cells.use]

# 5.2.7- Preprocess with SCTransform
query <- SCTransform(
  object = query,
  assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = reference$map),
  reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

# 5.2.8- Find anchors between query and reference
anchors <- FindTransferAnchors(
  reference = reference$map,
  query = query,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)

# 5.2.9- Transfer cell type labels and impute protein expression
#
# Transferred labels are in metadata columns named "predicted.*"
# The maximum prediction score is in a metadata column named "predicted.*.score"
# The prediction scores for each class are in an assay named "prediction.score.*"
# The imputed assay is named "impADT" if computed

refdata <- lapply(X = "celltype.l2", function(x) {
  reference$map[[x, drop = TRUE]]
})
names(x = refdata) <- "celltype.l2"
if (FALSE) {
  refdata[["impADT"]] <- GetAssayData(
    object = reference$map[['ADT']], #looks like the Azimuth ref does not contain ADT data but this can be imputed from refs that contain ADT data.
    slot = 'data'
  )
}
query <- TransferData(
  reference = reference$map,
  query = query,
  dims = 1:50,
  anchorset = anchors,
  refdata = refdata,
  n.trees = 20,
  store.weights = TRUE
)

# Calculate the embeddings of the query data on the reference SPCA
query <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = reference$map,
  query = query,
  reductions = "pcaproject",
  reuse.weights.matrix = TRUE
)

# Calculate the query neighbors in the reference with respect to the integrated embeddings
query[["query_ref.nn"]] <- FindNeighbors(
  object = Embeddings(reference$map[["refDR"]]),
  query = Embeddings(query[["integrated_dr"]]),
  return.neighbor = TRUE,
  l2.norm = TRUE
)

# The reference used in the app is downsampled compared to the reference on which
# the UMAP model was computed. This step, using the helper function NNTransform,
# corrects the Neighbors to account for the downsampling.
query <- NNTransform(
  object = query,
  meta.data = reference$map[[]]
)

# Project the query to the reference UMAP.
query[["proj.umap"]] <- RunUMAP(
  object = query[["query_ref.nn"]],
  reduction.model = reference$map[["refUMAP"]],
  reduction.key = 'UMAP_'
)


# Calculate mapping score and add to metadata
query <- AddMetaData(
  object = query,
  metadata = MappingScore(anchors = anchors),
  col.name = "mapping.score"
)

# VISUALIZATIONS

# First predicted metadata field, change to visualize other predicted metadata
id <- "celltype.l2"[1]
predicted.id <- paste0("predicted.", id)

p2 <- DimPlot(query, group.by = "orig.ident")
p2

# DimPlot of the reference
DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE) + NoLegend()

p1 <- dittoDimPlot(reference$plot, var = "celltype.l2",
                   reduction.use = "refUMAP", size = 0.5,
                   do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")
p1

# DimPlot of the query, colored by predicted cell type
DimPlot(object = query, reduction = "proj.umap", group.by = predicted.id, label = TRUE) + NoLegend()

p1 <- dittoDimPlot(query, var = "predicted.celltype.l2",
                   reduction.use = "proj.umap", size = 0.5,
                   do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")
p1

p2 <- dittoBarPlot(query, var = "predicted.celltype.l2", group.by = "orig.ident")
p2

# Plot the score for the predicted cell type of the query
FeaturePlot(object = query, features = paste0(predicted.id, ".score"), reduction = "proj.umap")
VlnPlot(object = query, features = paste0(predicted.id, ".score"), group.by = predicted.id) + NoLegend()

# Plot the mapping score
FeaturePlot(object = query, features = "mapping.score", reduction = "proj.umap")
VlnPlot(object = query, features = "mapping.score", group.by = predicted.id) + NoLegend()

# Plot the prediction score for the class CD16 Mono
FeaturePlot(object = query, features = "CD16 Mono", reduction = "proj.umap")
VlnPlot(object = query, features = "CD16 Mono", group.by = predicted.id) + NoLegend()

# Plot an RNA feature
FeaturePlot(object = query, features = "TFRC", reduction = "proj.umap")
VlnPlot(object = query, features = "TFRC", group.by = predicted.id) + NoLegend()

FeaturePlot(object = query, features = "GYPA", reduction = "proj.umap")
VlnPlot(object = query, features = "GYPA", group.by = predicted.id) + NoLegend()

# Plot an imputed protein feature - if the ref dataset contains ADT data and it was tranferred to the query dataset
if (FALSE) {
  FeaturePlot(object = query, features = "CD3-TotalSeqB", reduction = "proj.umap")
  VlnPlot(object = query, features = "CD3-TotalSeqB", group.by = predicted.id) + NoLegend()
}

saveRDS(query, "BM1_Azimuth_local.rds")

# 5.3- Map to any published dataset 

# 5.3.1- Multimodal reference mapping - Option 1
# Reference: https://satijalab.org/seurat/articles/multimodal_reference_mapping.html

# This vignette introduces the process of mapping query datasets to annotated references in Seurat.
# Supervised analysis guided by a reference dataset can help to enumerate cell states that would be challenging to find with unsupervised analysis. 
# We have previously demonstrated how to use reference-mapping approach to annotate cell labels in a query dataset. 
# In Seurat v4, we have substantially improved the speed and memory requirements for integrative tasks including reference mapping, and also include new functionality to project query cells onto a previously computed UMAP visualization.

# In this vignette, we demonstrate how to use a previously established reference to interpret an scRNA-seq query:
  
# Annotate each query cell based on a set of reference-defined cell states
# Project each query cell onto a previously computed UMAP visualization
# Impute the predicted levels of surface proteins that were measured in the CITE-seq reference

# Load reference dataset available through SeuratData
library(SeuratData)
#load reference data
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
#load query data - it is the BM1_query

# The reference dataset contains a WNN graph, reflecting a weighted combination of the RNA and protein data in this CITE-seq experiment.
# We can compute a UMAP visualization based on this graph. We set return.model = TRUE, which will enable us to project query datasets onto this visualization.
bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
              reduction.key = "wnnUMAP_", return.model = TRUE)
DimPlot(bm, group.by = "celltype.l2", reduction = "wnn.umap") 

# Computing an sPCA transformation

# As described in our manuscript, we first compute a ‘supervised’ PCA. 
# This identifies the transformation of the transcriptome data that best encapsulates the structure of the WNN graph. 
# This allows a weighted combination of the protein and RNA measurements to ‘supervise’ the PCA, and highlight the most relevant sources of variation. 
# After computing this transformation, we can project it onto a query dataset. 
# We can also compute and project a PCA projection, but recommend the use of sPCA when working with multimodal references that have been constructed with WNN analysis.
# The sPCA calculation is performed once, and then can be rapidly projected onto each query dataset.

bm <- ScaleData(bm, assay = 'RNA')
bm <- RunSPCA(bm, assay = 'RNA', graph = 'wsnn')

# Computing a cached neighbor index

# Since we will be mapping multiple query samples to the same reference, we can cache particular steps that only involve the reference. 
# This step is optional but will improve speed when mapping multiple samples.
# We compute the first 50 neighbors in the sPCA space of the reference. 
# We store this information in the spca.annoy.neighbors object within the reference Seurat object and also cache the annoy index data structure (via cache.index = TRUE).

bm <- FindNeighbors(
  object = bm,
  reduction = "spca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors", 
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

# IMPORTANT: If you want to save and load a cached index for a Neighbor object generated with method = "annoy" and cache.index = TRUE, use the SaveAnnoyIndex()/LoadAnnoyIndex() functions. 
# Importantly, this index cannot be saved normally to an RDS or RDA file, so it will not persist correctly across R session restarts or saveRDS/readRDS for the Seurat object containing it. 
# Instead, use LoadAnnoyIndex() to add the Annoy index to the Neighbor object every time R restarts or you load the reference Seurat object from RDS. 
# The file created by SaveAnnoyIndex() can be distributed along with a reference Seurat object, and added to the Neighbor object in the reference.
bm[["spca.annoy.neighbors"]]
SaveAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./reftmp.idx")
#bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = bm[["spca.annoy.neighbors"]], file = "./Multimodal_Seurat_reference/reftmp.idx")

# Query dataset preprocessing

# Here we will demonstrate mapping multiple donor bone marrow samples to the multimodal bone marrow reference. 
# These query datasets are derived from the Human Cell Atlas (HCA) Immune Cell Atlas Bone marrow dataset and are available through SeuratData. 
# This dataset is provided as a single merged object with 8 donors. 
# We first split the data back into 8 separate Seurat objects, one for each original donor to map individually.sing
# if running this in the merged/integrated samples dataset - if running in individual sample, skip it!

library(dplyr)
library(SeuratData)
InstallData('hcabm40k')
hcabm40k.batches <- SplitObject(hcabm40k, split.by = "orig.ident")

# We then normalize the query in the same manner as the reference. Here, the reference was normalized using log-normalization via NormalizeData(). 
# If the reference had been normalized using SCTransform(), the query must be normalized with SCTransform() as well.
BM1_query2 <- BM1_query #this is not necessary!
DefaultAssay(BM1_query2) <- "RNA" 
BM1_query2 <- NormalizeData(BM1_query2, verbose = FALSE) #maybe use the original seurat object, not sure this extra normalisation step would affect the data.

# Mapping

# We then find anchors between each donor query dataset and the multimodal reference. 
# This command is optimized to minimize mapping time, by passing in a pre-computed set of reference neighbors, and turning off anchor filtration.
DefaultAssay(BM1_query) <- "RNA" # make sure to change the default assat to match the reference default assay

# Command below if mapping several samples in a loop
anchors <- list()
for (i in 1:length(hcabm40k.batches)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = bm,
    query = hcabm40k.batches[[i]],
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
}
# Adapt command to run each sample individually
anchors <- FindTransferAnchors(
    reference = bm,
    query = BM1_query,
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )

# We then individually map each of the datasets.
# Command below if mapping several samples in a loop
for (i in 1:length(hcabm40k.batches)) {
  hcabm40k.batches[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = hcabm40k.batches[[i]],
    reference = bm, 
    refdata = list(
      celltype = "celltype.l2", 
      predicted_ADT = "ADT"),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
}
# Adapt command to run each sample individually
BM1_query <- MapQuery(
    anchorset = anchors, 
    query = BM1_query,
    reference = bm, 
    refdata = list(
      celltype = "celltype.l2", 
      predicted_ADT = "ADT"),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )

rownames(BM1_query@assays[["predicted_ADT"]]) #full list of ADT features in the imputed ADT SLOT
 
# Explore the mapping results
p1 <- DimPlot(BM1_query, reduction = 'ref.umap', group.by = 'predicted.celltype', label = TRUE) + NoLegend()
p1

# If mapping several samples at once: We can also merge all the objects into one dataset. Note that they have all been integrated into a common space, defined by the reference. 
# We can then visualize the results together.
# Merge the batches 
hcabm40k <- merge(hcabm40k.batches[[1]], hcabm40k.batches[2:length(hcabm40k.batches)], merge.dr = "ref.umap")
DimPlot(hcabm40k, reduction = "ref.umap", group.by =  "predicted.celltype", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()

# We can visualize gene expression, cluster prediction scores, and (imputed) surface protein levels in the query cells:
DefaultAssay(BM1_query) <- 'SCT'
p3 <- FeaturePlot(BM1_query, features = c("GYPA", "TFRC", "CD34", "CD3E"), reduction = 'ref.umap', 
                  max.cutoff = 3, ncol = 2)
p3
# cell type prediction scores
DefaultAssay(BM1_query) <- 'prediction.score.celltype'
p4 <- FeaturePlot(BM1_query, features = c("CD14 Mono", "CD4 Naive", "CD8 Naive", "Memory B"), reduction = 'ref.umap', ncol = 2, 
                  cols = c("lightgrey", "darkred"))
p4

p5 <- FeaturePlot(BM1_query, features = c("MAIT", "NK", "Prog-RBC", "HSC"), reduction = 'ref.umap', ncol =2, 
                  cols = c("lightgrey", "darkred"))
p5

# imputed protein levels
DefaultAssay(BM1_query) <- 'predicted_ADT'
p6 <- FeaturePlot(BM1_query, features = c("CD45RA","CD45RO", "CD14", "CD16"), reduction = 'ref.umap',
                  min.cutoff = 'q10', max.cutoff = 'q99', cols = c("lightgrey", "darkgreen") ,
                  ncol = 2)
p6

p7 <- FeaturePlot(BM1_query, features = c("CD19", "CD3", "CD8a", "CD127-IL7Ra"), reduction = 'ref.umap',
                  min.cutoff = 'q10', max.cutoff = 'q99', cols = c("lightgrey", "darkgreen") ,
                  ncol = 2)
p7

p8 <- FeaturePlot(BM1_query, features = c("CD56", "CD57", "CD161", "CD11c"), reduction = 'ref.umap',
                  min.cutoff = 'q10', max.cutoff = 'q99', cols = c("lightgrey", "darkgreen") ,
                  ncol = 2)
p8

# Explore the mapping results - to compare with the Azimuth ref mapping, use the Azimuth ref reduction to generate the UMAPs

p1 <- DimPlot(BM1_query, reduction = 'proj.umap', group.by = 'predicted.celltype', label = TRUE) + NoLegend()
p1

p1 <- dittoDimPlot(BM1_query, var = "predicted.celltype",
                   reduction.use = "proj.umap", size = 0.75,
                   do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")
p1

p2 <- dittoBarPlot(BM1_query, var = "predicted.celltype", group.by = "orig.ident")
p2

# We can visualize gene expression, cluster prediction scores, and (imputed) surface protein levels in the query cells:
DefaultAssay(BM1_query) <- 'SCT'
p3 <- FeaturePlot(BM1_query, features = c("GYPA", "TFRC", "CD34", "CD3E"), reduction = 'proj.umap', 
                  max.cutoff = 3, ncol = 2)
p3
# cell type prediction scores
DefaultAssay(BM1_query) <- 'prediction.score.celltype'
p4 <- FeaturePlot(BM1_query, features = c("CD14 Mono", "CD4 Naive", "CD8 Naive", "Memory B"), reduction = 'proj.umap', ncol = 2, 
                  cols = c("lightgrey", "darkred"))
p4

p5 <- FeaturePlot(BM1_query, features = c("MAIT", "NK", "Prog-RBC", "HSC"), reduction = 'proj.umap', ncol =2, 
                  cols = c("lightgrey", "darkred"))
p5

# imputed protein levels
DefaultAssay(BM1_query) <- 'predicted_ADT'
p6 <- FeaturePlot(BM1_query, features = c("CD45RA","CD45RO", "CD14", "CD16"), reduction = 'proj.umap',
                  min.cutoff = 'q10', max.cutoff = 'q99', cols = c("lightgrey", "darkgreen") ,
                  ncol = 2)
p6

p7 <- FeaturePlot(BM1_query, features = c("CD19", "CD3", "CD8a", "CD127-IL7Ra"), reduction = 'proj.umap',
                  min.cutoff = 'q10', max.cutoff = 'q99', cols = c("lightgrey", "darkgreen") ,
                  ncol = 2)
p7

p8 <- FeaturePlot(BM1_query, features = c("CD56", "CD57", "CD161", "CD11c"), reduction = 'proj.umap',
                  min.cutoff = 'q10', max.cutoff = 'q99', cols = c("lightgrey", "darkgreen") ,
                  ncol = 2)
p8

# Computing a new UMAP visualization

# In the previous examples, we visualize the query cells after mapping to the reference-derived UMAP. 
# Keeping a consistent visualization can assist with the interpretation of new datasets. 
# However, if there are cell states that are present in the query dataset that are not represented in the reference, they will project to the most similar cell in the reference. 
# This is the expected behavior and functionality as established by the UMAP package, but can potentially mask the presence of new cell types in the query which may be of interest.
# In our manuscript, we map a query dataset containing developing and differentiated neutrophils, which are not included in our reference. 
# We find that computing a new UMAP (‘de novo visualization’) after merging the reference and query can help to identify these populations, as demonstrated in Supplementary Figure 8. 
# In the ‘de novo’ visualization, unique cell states in the query remain separated. 
# We emphasize that if users are attempting to map datasets where the underlying samples are not BM or PBMC, or contain cell types that are not present in the reference, computing a ‘de novo’ visualization is an important step in interpreting their dataset.

# merge reference and query
bm$id <- 'reference'
BM1_query$id <- 'query'
refquery <- merge(bm, BM1_query)
refquery[["spca"]] <- merge(bm[["spca"]], BM1_query[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)
DimPlot(refquery, group.by = 'id', shuffle = TRUE)
DimPlot(refquery, group.by = 'predicted.celltype', label = TRUE) + NoLegend()
DimPlot(refquery, group.by = 'predicted.celltype.l2', label = TRUE) + NoLegend()

table(BM1_query$predicted.celltype)
table(BM1_query$predicted.celltype.l2)

# 6- Differential expression testing

# 6.1- Perform default differential expression tests

# The bulk of Seurat’s differential expression features can be accessed through the FindMarkers() function. 
# As a default, Seurat performs differential expression based on the non-parametric Wilcoxon rank sum test. 
# This replaces the previous default test (‘bimod’). 
# To test for differential expression between two specific groups of cells, specify the ident.1 and ident.2 parameters. 

# which assay should be used for DGE comparisons? 
# Reference: Current best practices in single-cell RNA-seq analysis: a tutorial" described to use the raw pre-normalised measured values.

# list options for groups to perform differential expression on
levels(BM1_query)

# 6.2- Prefilter features or cells to increase the speed of DE testing

# To increase the speed of marker discovery, particularly for large datasets, Seurat allows for pre-filtering of features or cells. 
# For example, features that are very infrequently detected in either group of cells, or features that are expressed at similar average levels, are unlikely to be differentially expressed. 
# Example use cases of the min.pct, logfc.threshold, min.diff.pct, and max.cells.per.ident parameters are demonstrated below.

# Use min.pct to pre-filter features that are detected at <X% frequency;
# Use logfc.threshold to pre-filter features that have less than a X-fold change between the average expression;
# Use min.diff.pct to pre-filter features whose detection percentages across the two groups are similar;
# Use max.cells.per.ident to subsample each group to a maximum of 200 cells. Can be very useful for large clusters, or computationally-intensive DE tests;
# Increasing min.pct, logfc.threshold, and min.diff.pct, will increase the speed of DE testing, but could also miss features that are prefiltered.

#  Find differentially expressed features between Late Eryth and all other cells, only search for positive markers
Late_Eryth_markers <- FindMarkers(BM1_query, assay = "RNA", slot = "data", ident.1 = "Late Eryth", only.pos = TRUE, logfc.threshold = 0.1)
print((Late_Eryth_markers))
clipr::write_clip(Late_Eryth_markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# Find all markers distinguishing early vs late erythrocytes
Eryth_markers <- FindMarkers(BM1_query, assay = "RNA", slot = "data", ident.1 = "Early Eryth", ident.2 = "Late Eryth", min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.1)
print(Eryth_markers)
clipr::write_clip(Eryth_markers)

# Exploring the expression of canonical marker genes
Idents(BM1_query) <- 'predicted.celltype.l2'
VlnPlot(BM1_query, features = c("TFRC", "ITGA4", "GYPA"), assay = "RNA", slot = "data", sort = TRUE, ncol=1) + NoLegend()
VlnPlot(BM1_query, features = c("TFRC", "ITGA4", "GYPA"), assay = "SCT", slot = "data", sort = TRUE, ncol=1) + NoLegend()

# Using the RNA normalised slot
BM1_RNA_data_markers <- FindAllMarkers(BM1_query, assay = "RNA", slot = "data", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
clipr::write_clip(BM1_RNA_data_markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste
saveRDS(BM1_RNA_data_markers, "BM1_RNA_data_markers.rds")

# top 5 upregulated genes per cluster
BM1_RNA_data_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  print(n=100)

# Using the SCT normalised slot

BM1_SCT_data_markers <- FindAllMarkers(BM1_query, assay = "SCT", slot = "data", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
clipr::write_clip(BM1_SCT_data_markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste
saveRDS(BM1_SCT_data_markers, "BM1_SCT_data_markers.rds")

BM1_SCT_data_markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  print(n=100)

# 7- Visualisations in Seurat

# Reference: https://satijalab.org/seurat/articles/visualization_vignette.html
# Five visualizations of marker feature expression

# 7.1- # Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
features <- c("TFRC", "ITGA4", "GYPA", "CD34")
features2 <- c("MS4A3", "CD247", "CCL5", "CD69")
features3 <- c("GZMB", "PF4", "SDC1", "CSF3R")
features4 <- c("PECAM1", "ICAM1", "ANGPT1", "SELE", "KITLG", "SPP1")

p1 <- RidgePlot(BM1_query, features = features, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
#+ theme(plot.title = element_text(size = rel(0.5)), axis.text = element_text(size = rel(0.5)))
p1

p1.5 <- RidgePlot(BM1_query, features = c("TFRC", "ITGA4", "GYPA"), assay = "RNA", slot = "data", idents = c("Late Eryth", "Early Eryth"), ncol = 2, sort = TRUE, combine = TRUE) 
p1.5

p2 <- RidgePlot(BM1_query, features = features2, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
p2

p3 <- RidgePlot(BM1_query, features = features3, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
p3

p4 <- RidgePlot(BM1_query, features = features4, assay = "RNA", slot = "data", ncol = 3, sort = TRUE, combine = TRUE) 
p4

# 7.2- Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(BM1_query, features = features, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE)
VlnPlot(BM1_query, features = features2, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE)
VlnPlot(BM1_query, features = features3, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE)
VlnPlot(BM1_query, features = features4, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE)

VlnPlot(BM1_query, features = c("nCount_RNA"), sort = TRUE, combine = TRUE)
VlnPlot(BM1_query, features = c("nFeature_RNA"), sort = TRUE, combine = TRUE)
VlnPlot(BM1_query, features = c("percent.mt"), sort = TRUE, combine = TRUE)

# 7.3- Feature plot - visualize feature expression in low-dimensional space
DefaultAssay(BM1_query) <- 'RNA'
FeaturePlot(BM1_query, features = features, ncol = 2)
FeaturePlot(BM1_query, features = features2, ncol = 2)
FeaturePlot(BM1_query, features = features3, ncol = 2)
FeaturePlot(BM1_query, features = features4, ncol = 2)

# 7.4- Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
features5 <- c("TFRC", "ITGA4", "GYPA", "CD34", "MS4A3", "CD247", "CCL5", "CD69", "GZMB", "PF4", "SDC1", "CSF3R", "PECAM1", "ICAM1", "ANGPT1", "SELE", "KITLG", "SPP1")
DotPlot(BM1_query, features = features5, assay = "RNA", cluster.idents = TRUE) + RotatedAxis()
DotPlot(BM1_query, features = features) + RotatedAxis()
DotPlot(BM1_query, features = features2) + RotatedAxis()
DotPlot(BM1_query, features = features3) + RotatedAxis()
DotPlot(BM1_query, features = features4) + RotatedAxis()

# 7.5- Single cell heatmap of feature expression
BM1_RNA_data_markers <- readRDS("./5_DGE/BM1_RNA_data_markers.rds")
write.csv(BM1_RNA_data_markers, "BM1_RNA_data_markers.csv")

#DoHeatmap(subset(BM1_query, downsample = 1000), features = features5, size = 2, combine = TRUE)
BM1_RNA_data_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
# Plot the top10 marker genes per cluster
DoHeatmap(BM1_query, features = top10$gene, size = 2, assay = "RNA", cells = 1:500) + NoLegend() + theme(axis.text.y = element_text(size = 3))
# changing the default color
DoHeatmap(BM1_query, features = top10$gene, size = 2, assay = "RNA", cells = 1:500) + NoLegend() + theme(axis.text.y = element_text(size = 3)) +
  scale_fill_viridis()

# Subset Seurat object top plot heatmap of specific clusters
DoHeatmap(subset(BM1_query, downsample = 1000, idents = c("Late Eryth", "Early Eryth")), features = top10$gene, size = 2, assay = "RNA", combine = TRUE) + theme(axis.text.y = element_text(size = 3))

# 8- Visualisations in DittoSeq

# Reference: https://bioconductor.org/packages/devel/bioc/vignettes/dittoSeq/inst/doc/dittoSeq.html
# All plotting functions use these colors, stored in dittoColors(), by default.

# 8.1- Functions:

# DimPlot/ (I)FeaturePlot / UMAPPlot / etc.	dittoDimPlot / multi_dittoDimPlot
# VlnPlot / RidgePlot	dittoPlot / multi_dittoPlot
# DotPlot	dittoDotPlot
# FeatureScatter / GenePlot	dittoScatterPlot
# DoHeatmap	dittoHeatmap*
# [No Seurat Equivalent]	dittoBarPlot / dittoFreqPlot
# [No Seurat Equivalent]	dittoDimHex / dittoScatterHex
# [No Seurat Equivalent]	dittoPlotVarsAcrossGroups
# SpatialDimPlot, SpatialFeaturePlot, etc.	dittoSpatial (coming soon!)

# 8.2- Input:

# Seurat has had inconsistency in input names from version to version. dittoSeq drew some of its parameter names from previous Seurat-equivalents to ease cross-conversion, 
# but continuing to blindly copy their parameter standards will break people’s already existing code. Instead, dittoSeq input names are guaranteed to remain consistent across versions, 
# unless a change is required for useful feature additions.

# Seurat Viz Input(s)	dittoSeq Equivalent(s)
# object	SAME
# features	var / vars (generally the 2nd input, so name not needed!) OR genes & metas for dittoHeatmap()
# cells (cell subsetting is not always available)	cells.use (consistently available)
# reduction & dims	reduction.use & dim.1, dim.2
# pt.size	size (or jitter.size)
# group.by	SAME
# split.by	SAME
# shape.by	SAME and also available in dittoPlot()
# fill.by	color.by (can be used to subset group.by further!)
# assay / slot	SAME
# order = logical	order but = “unordered” (default), “increasing”, or “decreasing”
# cols	color.panel for discrete OR min.color, max.color for continuous
# label & label.size & repel	do.label & labels.size & labels.repel
# interactive	do.hover = via plotly conversion
# [Not in Seurat]	data.out, do.raster, do.letter, do.ellipse, add.trajectory.lineages and others!

# 8.3- Helper Functions

# dittoSeq’s helper functions make it easy to determine the metadata, gene, and dimensionality reduction options for plotting.
getMetas(BM1_query)
# Query for the presence of a metadata slot
isMeta("nCount_RNA", BM1_query)
# Retrieve metadata values:
meta("predicted.celltype.l2", BM1_query)[1:10]
# Retrieve unique values of a metadata
metaLevels("predicted.celltype.l2", BM1_query)
metaLevels("predicted.celltype", BM1_query)
metaLevels("seurat_clusters", BM1_query)
# Retrieve all gene names
DefaultAssay(BM1_query) <- "SCT"
getGenes(BM1_query)[1:10]
# Query for the presence of a gene(s)
isGene("CD3E", BM1_query)
isGene(c("CD3E","ENO1","INS","non-gene"), BM1_query, return.values = TRUE)
# Retrieve gene expression values:
gene("CD3E", BM1_query)[1:10]
# Retrieve all dimensionality reductions
getReductions(BM1_query)

# 8.4- Visualisations

#There are many different types of dittoSeq visualizations. 
# Each has intuitive defaults which allow creation of immediately usable plots. 
# Each also has many additional tweaks available through discrete inputs that can help ensure you can create precisely-tuned, deliberately-labeled, publication-quality plots out-of-the-box.

# 8.4.1- dittoDimPlot & dittoScatterPlot
# These show cells/samples data overlaid on a scatter plot, with the axes of dittoScatterPlot() being gene expression or metadata data and with the axes of dittoDimPlot() being dimensionality reductions like tsne, pca, umap or similar.

# Visualisation of specific genes/features in UMAP
dittoDimPlot(BM1_query, "CD34", reduction.use = "proj.umap", size = 0.5,
             do.label = TRUE, assay = "RNA",
             slot = "data", do.contour = FALSE,
             contour.color = "lightblue", # Optional, black by default
             contour.linetype = "dashed", max = 1) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("CD34 expression")

dittoScatterPlot(
  object = BM1_query,
  x.var = "GYPA", y.var = "TFRC",
  color.var = "predicted.celltype.l2")

# Visualisation of QC in UMAP
dittoDimPlot(BM1_query, "nCount_RNA", reduction.use = "proj.umap", size = 0.5,
             do.label = TRUE, assay = "RNA",
             slot = "data", do.contour = FALSE,
             contour.color = "lightblue", # Optional, black by default
             contour.linetype = "dashed", max = 50000) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("UMI counts")

dittoDimPlot(BM1_query, "nFeature_RNA", reduction.use = "proj.umap", size = 0.5,
             do.label = TRUE, assay = "RNA",
             slot = "data", do.contour = FALSE,
             contour.color = "lightblue", # Optional, black by default
             contour.linetype = "dashed") + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Number of Genes")

dittoDimPlot(BM1_query, "percent.mt", reduction.use = "proj.umap", size = 0.5,
             do.label = TRUE, assay = "RNA",
             slot = "data", do.contour = FALSE,
             contour.color = "lightblue", # Optional, black by default
             contour.linetype = "dashed") + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Percentage of mitochondrial genes")

dittoScatterPlot(
  object = BM1_query,
  x.var = "nCount_RNA", y.var = "nFeature_RNA",
  color.var = "percent.mt")

# Visualisation of prediction scores in UMAP
dittoDimPlot(BM1_query, "predicted.celltype.score")
dittoDimPlot(BM1_query, "predicted.celltype.l2.score")

# Visualisation of cell types in UMAP
dittoDimPlot(BM1_query, var = "predicted.celltype",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")

dittoDimPlot(BM1_query, var = "predicted.celltype.l2",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")

# 8.4.2- dittoPlot (and dittoRidgePlot + dittoBoxPlot wrappers)

# These display continuous cells/samples’ data on a y-axis (or x-axis for ridgeplots) grouped on the x-axis by sample, age, condition, or any discrete grouping metadata. 
# Data can be represented with violin plots, box plots, individual points for each cell/sample, and/or ridge plots. 
# The plots input controls which data representations are used. 
# The group.by input controls how the data are grouped in the x-axis. And the color.by input controls the colors that fill in violin, box, and ridge plots.

# dittoPlot() is the main function, but dittoRidgePlot() and dittoBoxPlot() are wrappers which essentially just adjust the default for the plots input from c(“jitter”, “vlnplot”) to c(“ridgeplot”) or c(“boxplot”,“jitter”), respectively.

dittoPlot(BM1_query, "CD34", assay = "RNA", slot = "data", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  vlnplot.lineweight = 0.5)
dittoPlot(BM1_query, "nCount_RNA", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  vlnplot.lineweight = 0.5)
dittoPlot(BM1_query, "nFeature_RNA", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  vlnplot.lineweight = 0.5)
dittoPlot(BM1_query, "percent.mt", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  vlnplot.lineweight = 0.5)

dittoRidgePlot(BM1_query, "CD34", assay = "RNA", slot = "data", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), ridgeplot.lineweight = 0.5, max = 2)
dittoRidgePlot(BM1_query, "nCount_RNA", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  ridgeplot.lineweight = 0.5, max = 75000)
dittoRidgePlot(BM1_query, "nFeature_RNA", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  ridgeplot.lineweight = 0.5)
dittoRidgePlot(BM1_query, "percent.mt", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  ridgeplot.lineweight = 0.5)

# 8.4.3- dittoBarPlot & dittoFreqPlot

# A couple of very handy visualizations missing from some other major single-cell visualization toolsets, 
# these functions quantify and display frequencies of clusters or cell types (or other discrete data) per sample (or other discrete groupings). 
# Such visualizations are quite useful for QC-ing clustering for batch effects and generally assessing cell type fluctuations.
# For both, data can be represented as percentages or counts, and this is controlled by the scale input.
# Visualisation of  frequency of cell types in bar plots

dittoBarPlot(BM1_query, var = "predicted.celltype", group.by = "orig.ident")
dittoBarPlot(BM1_query, var = "predicted.celltype", group.by = "orig.ident", scale = "count")

# dittoFreqPlot separates each cell type into its own facet, and thus puts more emphasis on individual cells. 
# An additional sample.by input controls splitting of cells within group.by-groups into individual samples.

# 8.4.4-  dittoHeatmap

# This function is essentially a wrapper for generating heatmaps with pheatmap, but with the same automatic, user-friendly, data extraction, (subsetting,) and metadata integration common to other dittoSeq functions.
# For large, many cell, single-cell datasets, it can be necessary to turn off clustering by cells in generating the heatmap because the process is very memory intensive. 
# As an alternative, dittoHeatmap offers the ability to order columns in functional ways using the order.by input. 
# This input will default to the first annotation provided to annot.by for single cell datasets, but can also be controlled separately.

# Pick Genes
genes <- c("TFRC", "ITGA4", "GYPA", "CD34", "MS4A3", "CD247", "CCL5", "CD69",
           "GZMB", "PF4", "SDC1", "CSF3R", "CD14", "FCGR3A","MME", "ARG1",
           "PECAM1", "ICAM1", "ANGPT1", "SELE", "KITLG", "SPP1")

# Annotating and ordering cells by some meaningful feature(s):
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
display.brewer.all(n=10, exact.n=FALSE)
# col = colorRampPalette(rev(brewer.pal(n = 9, name = "Spectral")))(3)
# col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), brewer.pal(n = 5, name = "YlOrRd"), space = "RGB")
# colorRamp2(c(-0.5, 0, 0.5), c("blue", "yellow", "red"), space = "RGB"),
# col = colorRamp2(c(-1,-0.5, 0, 0.5, 1), magma(5, direction=-1))
# col = colorRamp2(c(-1, -0.5, 0, 0.5, 1), viridis(256) 

dittoHeatmap(BM1_query, genes, assay = "RNA", slot = "scale.data", annot.by = c("predicted.celltype.l2"), 
             heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50), scale = "column", heatmap.colors.max.scaled = inferno(10), cluster_rows = FALSE, scaled.to.max = FALSE)

dittoHeatmap(BM1_query, genes, assay = "RNA", slot = "scale.data", annot.by = c("predicted.celltype.l2"), 
             heatmap.colors = inferno(3), scale = "column", heatmap.colors.max.scaled = inferno(10), cluster_rows = FALSE, scaled.to.max = FALSE)

BM1_RNA_data_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
# Plot the top10 marker genes per cluster
dittoHeatmap(BM1_query, top10$gene, assay = "RNA", slot = "scale.data", annot.by = c("predicted.celltype.l2"), 
             heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50), scale = "column", 
             heatmap.colors.max.scaled = inferno(10), cluster_rows = FALSE, scaled.to.max = FALSE, cluster_cols = FALSE,
             fontsize_row = 2)

# 8.4.5-  Multi-Plotters

# These create either multiple plots or create plots that summarize data for multiple variables all in one plot. 
# They make it easier to create summaries for many genes or many cell types without the need for writing loops.

# 8.4.5.1- dittoDotPlot

# A very succinct representation that is useful for showing differences between groups. 
# The plot uses differently colored and sized dots to summarizes both expression level (color) and percent of cells/samples with non-zero expression (size) for multiple genes (or values of metadata) within different groups of cells/samples.
# By default, expression values for all groups are centered and scaled to ensure a similar range of values for all vars displayed and to emphasize differences between groups.
dittoDotPlot(BM1_query, vars = genes, 
             assay = "RNA", slot = "data",
             group.by = "predicted.celltype.l2", scale = FALSE) + theme_light()
# String which sets whether the values shown with color (default: mean non-zero expression) should be centered and scaled.
# scale = FALSE - plots average expression
# scale = TRUE (default) - plots relative expression


# 9- Predict and remove doublets with DoubletFinder - should be done after UMAP and before clustering and run this per sample before merging or integration!

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
sweep.res.list_BM1 <- paramSweep_v3(BM1_query, PCs = 1:10, sct = FALSE)
sweep.stats_BM1 <- summarizeSweep(sweep.res.list_BM1, GT = FALSE)
bcmvn_BM1 <- find.pK(sweep.stats_BM1)
barplot(bcmvn_BM1$BCmetric, names.arg = bcmvn_BM1$pK, las=2) # pK = 0.005

## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
sweep.res.list_BM1 <- paramSweep_v3(BM1_query, PCs = 1:10, sct = FALSE)
gt.calls <- BM1_query@meta.data[rownames(sweep.res.list_BM1[[1]]), "GT"]   ## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
sweep.stats_BM1 <- summarizeSweep(sweep.res.list_BM1, GT = TRUE, GT.calls = gt.calls)
bcmvn_BM1 <- find.pK(sweep.stats_BM1)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(BM1_query@meta.data$predicted.celltype.l2)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.004*nrow(BM1_query@meta.data))  ## Assuming 4% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
BM1_query <- doubletFinder_v3(BM1_query, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi) #Run 1
BM1_query <- doubletFinder_v3(BM1_query, PCs = 1:10, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_25") #Run 2

# Plot results - name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(BM1_query@meta.data)[grepl("DF.classification", colnames(BM1_query@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(BM1_query, group.by = "orig.ident") + NoAxes(),
                   DimPlot(BM1_query, group.by = DF.name) + NoAxes())
VlnPlot(BM1_query, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

# Remove all predicted doublets from our data.
BM1_query = BM1_query[, BM1_query@meta.data[, DF.name] == "Singlet"]
dim(BM1_query)

# 10- Detecting and removing empty droplets with DoubletUtils - run this per sample before merging or integration!

#R eference: https://bioconductor.org/packages/devel/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html#detecting-empty-droplets
# https://rockefelleruniversity.github.io/scRNA-seq/exercises/answers/exercise2_answers.html

# Empty droplets often contain RNA from the ambient solution, resulting in non-zero counts after debarcoding. 
# The emptyDrops function is designed to distinguish between empty droplets and cells. It does so by testing each barcode’s expression profile for significant deviation from the ambient profile. 
# Given a matrix my.counts containing UMI counts for all barcodes, we call:

# To use DropletUtils convert Seurat object to SCE object
BM1.sce <- as.SingleCellExperiment(BM1_query, assay = "RNA")

# Draw a knee plot and identify inflection point and knee point.
bcrank <- barcodeRanks(counts(BM1.sce))
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
e.out <- emptyDrops(counts(BM1.sce),lower=limit, test.ambient=TRUE) # it did not identified empty droplets or ambient RNA

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
BM1.sce2 <- BM1.sce[,which(e.out$FDR <= 0.01)]

################################################################################################################################################################################

# 11- CITE-seq analysis - Normalisation using dsb package

# References: https://github.com/niaid/dsb
# https://cran.r-project.org/web/packages/dsb/vignettes/end_to_end_workflow.html
# Normalizing and denoising protein expression data from droplet-based single cell proﬁling
# https://github.com/niaid/dsb_manuscript/#instructions

library(dsb)

DefaultAssay(BM1_query) <- "ADT"

# 11.1- read raw data using the Seurat function "Read10X" 
raw = Seurat::Read10X("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/scRNAseq/2586_Test_run_PBMC_BM_lonza_exp1/Files/Sample_BM-1/outs/raw_feature_bc_matrix")
cells = Seurat::Read10X("/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/scRNAseq/2586_Test_run_PBMC_BM_lonza_exp1/Files/Sample_BM-1/outs/filtered_feature_bc_matrix/")

# define cell-containing barcodes and separate cells and empty drops
stained_cells = colnames(cells$`Gene Expression`)
background = setdiff(colnames(raw$`Gene Expression`), stained_cells)

# split the data into separate matrices for RNA and ADT
prot = raw$`Antibody Capture`
rna = raw$`Gene Expression`

# Now calculate some standard meta data for cells that we will use for quality control using standard approaches used in scRNAseq analysis pipelines.

# create metadata of droplet QC stats used in standard scRNAseq processing
mtgene = grep(pattern = "^MT-", rownames(rna), value = TRUE) # used below

md = data.frame(
  rna.size = log10(Matrix::colSums(rna)), 
  prot.size = log10(Matrix::colSums(prot)), 
  n.gene = Matrix::colSums(rna > 0), 
  mt.prop = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
)
# add indicator for barcodes Cell Ranger called as cells
md$drop.class = ifelse(rownames(md) %in% stained_cells, 'cell', 'background')

# remove barcodes with no evidence of capture in the experiment
md = md[md$rna.size > 0 & md$prot.size > 0, ]

background_drops = rownames(
  md[ md$prot.size > 1.5 & 
        md$prot.size < 3 & 
        md$rna.size < 2.5, ]
) 
background.adt.mtx = as.matrix(prot[ , background_drops])

# 11.2- calculate statistical thresholds for droplet filtering.
cellmd = md[md$drop.class == 'cell', ]

# filter drops with + / - 3 median absolute deviations from the median library size
rna.mult = (3*mad(cellmd$rna.size))
prot.mult = (3*mad(cellmd$prot.size))
rna.lower = median(cellmd$rna.size) - rna.mult
rna.upper = median(cellmd$rna.size) + rna.mult
prot.lower = median(cellmd$prot.size) - prot.mult
prot.upper = median(cellmd$prot.size) + prot.mult

# filter rows based on droplet qualty control metrics
qc_cells = rownames(
  cellmd[cellmd$prot.size > prot.lower & 
           cellmd$prot.size < prot.upper & 
           cellmd$rna.size > rna.lower & 
           cellmd$rna.size < rna.upper & 
           cellmd$mt.prop < 0.14, ]
)

# Check: are the number of cells passing QC in line with the expected recovery from the experiment?
length(qc_cells)

cell.adt.raw = as.matrix(prot[ , qc_cells]) 
cell.rna.raw = rna[ ,qc_cells]
cellmd = cellmd[qc_cells, ]

# 11.3- Optional step; remove proteins without staining

# Proteins without raw data in stained cells (below, the maximum UMI count for one protein was 4, an order of magnitude below even the isotype controls) can be removed from both matrices prior to normalization. In many cases, removing proteins is not necessary, but we recommend checking your raw data.

# flter 
pm = sort(apply(cell.adt.raw, 1, max))
pm2 = apply(background.adt.mtx, 1, max)
head(pm)

# 11.4- Normalize protein data with the DSBNormalizeProtein Function

# We are now ready to use dsb to normalize and denoise the ADTs. 
# We normalize the raw ADT matrix before creating a Seurat object since we also use the empty droplets.

# The method is carried out in a single step with a call to the DSBNormalizeProtein() function.
# cells_citeseq_mtx - the raw ADT UMI count matrix containing cells
# empty_drop_citeseq_mtx - the raw ADT UMI count matrix containing background droplets
# denoise.counts - we set to TRUE (the recommended default) to define and remove technical cell to cell variations. 
# isotype.control.name.vec - a vector of the isotype control antibodies from the rownames of the raw ADT matrix defined from rownames(cell.adt.raw) (see vignettes for data without isotype controls). 
# For data without isotype controls, see the vignette section Using dsb with data lacking isotype controls.

# The function returns a matrix of normalized protein values which can be integrated with any single cell analysis software. 
# Advanced users may want to examine internal stats used in dsb, in that case use return.stats = TRUE. 
# If the range of values after normalization is large, this is often due to a few low or high magnitude outliers. 
# The simplest way to address these outliers is to clip the maximum values by setting quantile.clipping = TRUE. 
# Finally a different pseudocount can be used with define.pseudocount = TRUE and pseudocount.use. 
# Please see the vignettes on CRAN to learn more about different parameters that can be used in dsb

# define isotype controls 
isotype.controls = c("IgG1_TotalSeqB", "IgG2a_TotalSeqB")

# normalize and denoise with dsb with 
cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = cell.adt.raw, 
  empty_drop_matrix = background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

# 11.5- Integrating dsb with Seurat - should run dsb before creating Seurat object. 

# 11.5.1- Create a Seurat object. 

# Make sure to add the dsb normalized matrix cell.adt.dsb to the data slot, not the counts slot.
# Seurat workflow 
# library(Seurat)
# integrating with Seurat
# stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.adt.raw))))
# stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.rna.raw))))

# create Seurat object note: min.cells is a gene filter, not a cell filter
# s = Seurat::CreateSeuratObject(counts = cell.rna.raw, meta.data = cellmd, assay = "RNA", min.cells = 20)

# add dsb normalized matrix "cell.adt.dsb" to the "CITE" data (not counts!) slot
# s[["CITE"]] = Seurat::CreateAssayObject(data = cells.dsb.norm)

# 11.6- Integrating dsb with Seurat object after was created

# add dsb normalized matrix "cell.adt.dsb" to the "CITE" data (not counts!) slot 
# define isotype controls 
isotype.controls = c("IgG1-TotalSeqB", "IgG2a-TotalSeqB")

# renaming rows of the backgroung matrix to mtch the rownames of the ADT counts matrix of Seurat object

new_names <- c("CD3-TotalSeqB","CD11b-TotalSeqB","CD20-TotalSeqB","CD31-TotalSeqB","CD34-TotalSeqB","CD44-TotalSeqB",
               "CD71-TotalSeqB", "CD49d-TotalSeqB", "IgG1-TotalSeqB", "IgG2a-TotalSeqB")
rownames(background.adt.mtx) <- new_names
background.adt.mtx <- as.matrix(background.adt.mtx)

# normalize and denoise with dsb with 
cells.dsb.norm = DSBNormalizeProtein(
  cell_protein_matrix = BM1_query@assays[["ADT"]]@counts, 
  empty_drop_matrix = background.adt.mtx, 
  denoise.counts = TRUE, 
  use.isotype.control = TRUE, 
  isotype.control.name.vec = isotype.controls
)

library("Matrix")
cells.dsb.norm <- as(cells.dsb.norm, "dgCMatrix")
cells.dsb.norm
BM1_query[["dsb"]] = Seurat::CreateAssayObject(data = cells.dsb.norm)

# 12- CITE-seq analysis - Normalisation using CLR from Seurat
# The resulting values can be interpreted as either a natural log ratio of the count for a given protein relative to the other proteins in the cell (CLR “across proteins”, as implemented in the original report of CITEseq 2 ) 
# or relative to other cells (CLR “across cells”, a modiﬁcation used in later work by the authors 19 , which renders CLR less dependent on the composition of the antibody panel). 
# The CLR transformation helps to better separate cell populations, but it does not directly estimate and correct for speciﬁc sources of technical noise including the apparent background noise.

# Normalize ADT data
DefaultAssay(BM1_query) <- "ADT"
BM1_query <- NormalizeData(BM1_query, normalization.method = "CLR", margin = 2) # If performing CLR normalization, normalize across features (1) or cells (2)
BM1_query@assays[["ADT"]]@data <- as(BM1_query@assays[["ADT"]]@data, "dgCMatrix")

p1 <- FeaturePlot(BM1_query, "CD20-TotalSeqB", slot = "data", cols = c("lightgrey", "darkgreen")) + ggtitle("CD20 protein")
p1


# 13- Visualisations of the CITE-seq normalisation methods - dsb vs CLR - and with dittoSeq

DefaultAssay(BM1_query) <- "dsb"
DefaultAssay(BM1_query) <- "ADT"

# Pick proteins
proteins <- c("CD3-TotalSeqB","CD11b-TotalSeqB","CD20-TotalSeqB","CD31-TotalSeqB","CD34-TotalSeqB","CD44-TotalSeqB",
              "CD71-TotalSeqB", "CD49d-TotalSeqB", "IgG1-TotalSeqB", "IgG2a-TotalSeqB")

# dittoDimPlot
# Visualisation of specific proteins in UMAP
multi_dittoDimPlot(BM1_query, vars = proteins, reduction.use = "proj.umap", size = 0.25,
             do.label = TRUE, assay = "ADT",
             slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression based on dsb")

# dittoScatterPlot
# Visualisation of specific combination of proteins in Scatter plot, colored by cell types
dittoScatterPlot(
  object = BM1_query,
  x.var = "CD3-TotalSeqB", y.var = "CD11b-TotalSeqB",
  color.var = "predicted.celltype.l2")

dittoScatterPlot(
  object = BM1_query,
  x.var = "CD71-TotalSeqB", y.var = "CD49d-TotalSeqB",
  color.var = "predicted.celltype.l2")

# dittoDotPlot
dittoDotPlot(BM1_query, vars = proteins, 
             assay = "ADT", slot = "data",min.color = "grey90",
             max.color = "#C51B7D",
             group.by = "predicted.celltype.l2", scale = FALSE) + theme_classic(base_size = 12)

# Violin and RidgePlot with Seurat
proteins1 <- c("CD3-TotalSeqB","CD11b-TotalSeqB","CD20-TotalSeqB","CD31-TotalSeqB")
proteins2 <- c("CD34-TotalSeqB","CD44-TotalSeqB","CD71-TotalSeqB", "CD49d-TotalSeqB")
proteins3 <- c("IgG1-TotalSeqB", "IgG2a-TotalSeqB")

VlnPlot(BM1_query, features = proteins1, assay = "ADT", slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()
VlnPlot(BM1_query, features = proteins2, assay = "ADT", slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()
VlnPlot(BM1_query, features = proteins3, assay = "ADT", slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()

RidgePlot(BM1_query, features = proteins1, assay = "ADT", slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
RidgePlot(BM1_query, features = proteins2, assay = "ADT", slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
RidgePlot(BM1_query, features = proteins3, assay = "ADT", slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 

# 13.1- Visualisations of the imputed CITE-seq data from the mapped Seurat ref dataset

DefaultAssay(BM1_query) <- "predicted_ADT"
rownames(BM1_query@assays[["predicted_ADT"]])
#[1] "CD11a"       "CD11c"       "CD123"       "CD127-IL7Ra" "CD14"        "CD16"        "CD161"       "CD19"       
#[9] "CD197-CCR7"  "CD25"        "CD27"        "CD278-ICOS"  "CD28"        "CD3"         "CD34"        "CD38"       
#[17] "CD4"         "CD45RA"      "CD45RO"      "CD56"        "CD57"        "CD69"        "CD79b"       "CD8a"       
#[25] "HLA.DR" 

# Pick proteins
proteins <- c("CD34", "CD11a","CD11c","CD123","CD14","CD16","CD161", "CD56", "CD57", "CD3", "CD4", "CD8a", "CD127-IL7Ra", "CD197-CCR7", "CD25","CD27", "CD278-ICOS", "CD28", "CD69", "CD19","CD79b","HLA.DR","CD45RA","CD45RO", "CD38")
proteins1 <- c("CD34", "CD11a","CD11c","CD123","CD14","CD16","CD161", "CD56", "CD57")
proteins2 <- c("CD3", "CD4", "CD8a", "CD127-IL7Ra", "CD197-CCR7", "CD25","CD27", "CD278-ICOS", "CD28", "CD69")
proteins3 <- c("CD19","CD79b","HLA.DR","CD45RA","CD45RO", "CD38")

# dittoDimPlot
# Visualisation of specific proteins in UMAP
multi_dittoDimPlot(BM1_query, vars = proteins1, reduction.use = "proj.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

multi_dittoDimPlot(BM1_query, vars = proteins2, reduction.use = "proj.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

multi_dittoDimPlot(BM1_query, vars = proteins3, reduction.use = "proj.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

# dittoScatterPlot
# Visualisation of specific combination of proteins in Scatter plot, colored by cell types
dittoScatterPlot(
  object = BM1_query,
  x.var = "CD3", y.var = "CD19",
  color.var = "predicted.celltype.l2")

# dittoDotPlot
dittoDotPlot(BM1_query, vars = proteins1, 
             slot = "data",min.color = "grey90",
             max.color = "#C51B7D",
             group.by = "predicted.celltype.l2", scale = FALSE) + theme_classic(base_size = 12)

dittoDotPlot(BM1_query, vars = proteins2, 
             slot = "data",min.color = "grey90",
             max.color = "#C51B7D",
             group.by = "predicted.celltype.l2", scale = FALSE) + theme_classic(base_size = 12)

dittoDotPlot(BM1_query, vars = proteins3, 
             slot = "data",min.color = "grey90",
             max.color = "#C51B7D",
             group.by = "predicted.celltype.l2", scale = FALSE) + theme_classic(base_size = 12)

# Violin and RidgePlot with Seurat
proteins1 <- c("CD34", "CD11a","CD11c","CD123")
proteins2 <- c("CD14","CD16", "CD161", "CD56", "CD57")
proteins3 <- c("CD3", "CD4", "CD8a", "CD127-IL7Ra", "CD197-CCR7")
proteins4 <- c("CD25", "CD27", "CD278-ICOS", "CD28", "CD69")
proteins5 <- c("CD19","CD79b","HLA.DR","CD45RA","CD45RO", "CD38")

VlnPlot(BM1_query, features = proteins1, slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()
VlnPlot(BM1_query, features = proteins2, slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()
VlnPlot(BM1_query, features = proteins3, slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()
VlnPlot(BM1_query, features = proteins4, slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()
VlnPlot(BM1_query, features = proteins5, slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()

RidgePlot(BM1_query, features = proteins1, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
RidgePlot(BM1_query, features = proteins2, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
RidgePlot(BM1_query, features = proteins3, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
RidgePlot(BM1_query, features = proteins4, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
RidgePlot(BM1_query, features = proteins5, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 

# 14- Clustering cells based on dsb normalized protein using Seurat

# Cluster cells based on dsb normalized protein levels. 
# Similar to workflow used in our paper Kotliarov et al. 2020 we don’t cluster based on principal components from ADT, instead directly using the normalized values.
DefaultAssay(BM1_query) <- "dsb"
# define proteins to use in clustering (non-isptype controls)
prots = rownames(BM1_query@assays[["dsb"]])[1:8]

# cluster and run umap 
BM1_query = FindNeighbors(object = BM1_query, dims = NULL, assay = 'dsb', 
                          features = prots, k.param = 30, 
                          verbose = FALSE)

# direct graph clustering 
BM1_query = FindClusters(object = BM1_query, resolution = 1, 
                         algorithm = 3,
                         graph.name = 'dsb_snn', 
                         verbose = FALSE)
# umap (optional)
BM1_query = RunUMAP(object = BM1_query, assay = "dsb", features = prots,
                     seed.use = 1990, min.dist = 0.2, n.neighbors = 30, reduction.key = "umap.dsb",
                     verbose = FALSE)

# make results dataframe 
d = cbind(BM1_query@meta.data, 
          as.data.frame(t(BM1_query@assays$dsb@data)),
          BM1_query@reductions$umap@cell.embeddings)

# To see if we recovered the expected cell populations, it is often more informative and interpretable to first look at the summarized (median or mean) protein expression in each cluster.
# we recommend doing this before trying to look at cells in 2-d visualization plots like umap.

# 14.1- dsb derived cluster interpretation

# dsb values are interpretable as the number of standard deviations of each protein from the expected noise (if using the default settings) with additional correction for cell to cell technical variations. 
# One can use this to set a threshold across all proteins for positivity, e.g. expression below 3 or 4 can be deemed unexpressed. If normalized with the same parameters from dsb, the same threshold can be applied across multiple datasets.
library(magrittr)
# calculate the median protein expression separately for each cluster 
adt_plot = d %>% 
  dplyr::group_by(dsb_snn_res.1) %>% 
  dplyr::summarize_at(.vars = prots, .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("dsb_snn_res.1") 
# plot a heatmap of the average dsb normalized values for each cluster
pheatmap::pheatmap(t(adt_plot), 
                   color = viridis::viridis(25, option = "B"), 
                  fontsize_row = 8, border_color = NA)

dittoHeatmap(BM1_query, var = "ident", assay = "dsb", slot = "data", annot.by = c("dsb_snn_res.1"), 
             heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50), scale = "column", 
             heatmap.colors.max.scaled = inferno(10), cluster_rows = FALSE, scaled.to.max = FALSE, cluster_cols = FALSE,
             fontsize_row = 10)

dittoDotPlot(BM1_query, vars = proteins, 
             assay = "dsb", slot = "data",min.color = "grey90",
             max.color = "#C51B7D",
             group.by = "dsb_snn_res.1", scale = FALSE) + theme_classic(base_size = 12)

# Visualisation of clusters in UMAP using different projections
dittoDimPlot(BM1_query, var = "dsb_snn_res.1",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from dsb on UMAP")

dittoDimPlot(BM1_query, var = "dsb_snn_res.1",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from dsb on UMAP")

# Frequency and total number of dsb clusters
dittoBarPlot(BM1_query, var = "dsb_snn_res.1", group.by = "orig.ident")
dittoBarPlot(BM1_query, var = "dsb_snn_res.1", group.by = "orig.ident", scale = "count")

# 14.2- Assigning cell type identity to dsb clusters

# Set idents from a value in object metadata - in this case the dsb clusters

colnames(BM1_query[[]])
Idents(BM1_query) <- 'dsb_snn_res.1'
levels(BM1_query)

# Rename idents of dsb clusters 
dsb.cluster.ids <- c("B cells", "T cells", "T cells", "T cells", "CD31+CD11b- cells", "B cells",
                     "CD71HighCD49d+CD44+ cells", "CD71LowCD49d+CD44+ cells", "CD31+CD11b+ cells", "CD34+ cells")

names(dsb.cluster.ids) <- levels(BM1_query@meta.data[["dsb_snn_res.1"]])
BM1_query <- RenameIdents(BM1_query, dsb.cluster.ids)

# Adding object metadata with cluster names in order - to be used in the bar plots
BM1_query$dsb_cluster_order <- factor(BM1_query@active.ident,
                                      levels = c("CD71HighCD49d+CD44+ cells", "CD71LowCD49d+CD44+ cells", "B cells", "T cells", "CD31+CD11b+ cells", "CD31+CD11b- cells", "CD34+ cells"))

# Visualise renamed/merged clusters
dittoDimPlot(BM1_query, var = "ident",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from dsb on UMAP")

dittoDimPlot(BM1_query, var = "ident",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from dsb on UMAP")

# Plot a heatmap of the average dsb normalized values for each cluster
e = cbind(BM1_query@active.ident, 
          as.data.frame(t(BM1_query@assays$dsb@data)))

adt_plot = e %>% 
  dplyr::group_by(BM1_query@active.ident) %>% 
  dplyr::summarize_at(.vars = prots, .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("BM1_query@active.ident") 

pheatmap::pheatmap(t(adt_plot), 
                   color = viridis::viridis(25, option = "B"), 
                   fontsize_row = 8, border_color = NA)


# Frequency and total number of dsb clusters
dittoBarPlot(BM1_query, var = "dsb_cluster_order", group.by = "orig.ident", retain.factor.levels = TRUE)
dittoBarPlot(BM1_query, var = "dsb_cluster_order", group.by = "orig.ident", retain.factor.levels = TRUE, scale = "count")

# 15- Clustering cells based on predicted ADT expression

# Cluster cells based on predicted ADT protein levels. 
DefaultAssay(BM1_query) <- "predicted_ADT"
# define proteins to use in clustering (non-isptype controls)
prots = rownames(BM1_query@assays[["predicted_ADT"]])

# cluster and run umap 
BM1_query = FindNeighbors(object = BM1_query, dims = NULL, assay = 'predicted_ADT', 
                          features = prots, k.param = 30, 
                          verbose = FALSE)

# direct graph clustering 
BM1_query = FindClusters(object = BM1_query, resolution = 1, 
                         algorithm = 3,
                         graph.name = 'predicted_ADT_snn', 
                         verbose = FALSE)
# umap (optional)
BM1_query = RunUMAP(object = BM1_query, assay = "predicted_ADT", features = prots,
                    seed.use = 1990, min.dist = 0.2, n.neighbors = 30, reduction.key = "predADT",
                    verbose = FALSE)

# make results dataframe 
d = cbind(BM1_query@meta.data, 
          as.data.frame(t(BM1_query@assays$predicted_ADT@data)))
          #BM1_query@reductions$umap@cell.embeddings)

# 15.1- Cluster interpretation

# calculate the median protein expression separately for each cluster 
adt_plot = d %>% 
  dplyr::group_by(predicted_ADT_snn_res.1) %>% 
  dplyr::summarize_at(.vars = prots, .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("predicted_ADT_snn_res.1") 
# plot a heatmap of the average dsb normalized values for each cluster
pheatmap::pheatmap(t(adt_plot), 
                   color = viridis::viridis(25, option = "B"), 
                   fontsize_row = 8, border_color = NA)

# Visualisation of clusters in UMAP using different projections
dittoDimPlot(BM1_query, var = "predicted_ADT_snn_res.1",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from predicted_ADT on UMAP")

dittoDimPlot(BM1_query, var = "predicted_ADT_snn_res.1",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from predicted_ADT on UMAP")

dittoBarPlot(BM1_query, var = "predicted_ADT_snn_res.1", group.by = "orig.ident")
dittoBarPlot(BM1_query, var = "predicted_ADT_snn_res.1", group.by = "orig.ident", scale = "count")

dittoDotPlot(BM1_query, vars = proteins, 
             assay = "predicted_ADT", slot = "data",min.color = "grey90",
             max.color = "#C51B7D",
             group.by = "predicted_ADT_snn_res.1", scale = FALSE) + theme_classic(base_size = 12)


# 16- Assigning major cell types to cell identities

# Set idents from a value in object metadata - in this case the dsb clusters

colnames(BM1_query[[]])
Idents(BM1_query) <- 'predicted.celltype.l2'
levels(BM1_query)

# Rename idents of dsb clusters 
major.cluster.ids <- c("B cell", "HSPC", "DC", "CD4 T cell", "CD8 T cell", "Erythroid", "CD4 T cell",
                     "B cell", "Myeloid", "B cell", "CD8 T cell", "B cell", "HSPC", "Erythroid",
                     "B cell", "HSPC", "HSCs", "DC", "DC", "HSPC", "CD4 T cell",
                     "CD8 T cell", "HSPC", "Natural Killer", "HSPC", "Myeloid")

names(major.cluster.ids) <- levels(BM1_query)
BM1_query <- RenameIdents(BM1_query, major.cluster.ids)

# Adding object metadata with cluster names in order - to be used in the bar plots
BM1_query$major_celltype.l2 <- factor(BM1_query@active.ident)

# Visualise renamed/merged clusters
dittoDimPlot(BM1_query, var = "major_celltype.l2",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Major cell types on UMAP")

# Frequency and total number of dsb clusters
dittoBarPlot(BM1_query, var = "major_celltype.l2", group.by = "orig.ident", retain.factor.levels = TRUE)
dittoBarPlot(BM1_query, var = "major_celltype.l2", group.by = "orig.ident", retain.factor.levels = TRUE, scale = "count")

# 17- Weighted Nearest Neighbor multimodal clustering using dsb normalized values

# The dsb normalized values can be used with the Weighted Nearest Neighbor multimodal clustering method. WNN is an excellent way to find fine-grained cell states

# 17.1- Method 1: Seurat WNN default with PCA on dsb normalized protein
# set up dsb values to use in WNN analysis (do not normalize with CLR, use dsb normalized values)
DefaultAssay(BM1_query) = "dsb"
prots = rownames(BM1_query@assays[["dsb"]])[1:8]
VariableFeatures(BM1_query) = prots
BM1_query = BM1_query %>% 
  ScaleData() %>% 
  RunPCA(reduction.name = 'apca', verbose = FALSE)

# run WNN 
BM1_query = FindMultiModalNeighbors(
  BM1_query, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:7), 
  modality.weight.name = "RNA.weight", 
  verbose = FALSE
)

# cluster 
BM1_query <- FindClusters(BM1_query, graph.name = "wsnn", 
                  algorithm = 3, 
                  resolution = 1.5, 
                  verbose = FALSE, 
                  random.seed = 1990)

# Visualise renamed/merged clusters
dittoDimPlot(BM1_query, var = "wsnn_res.1.5",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

# 17.2- Method 2– Seurat WNN with dsb normalized protein directly without PCA - for small ADT panels

# Below we show a version of WNN where we directly use normalized protein values without PCA compression. 
# We have found this procedure works well for smaller (e.g. less than 30) protein panels, 
# whereas datasets with many cells generated with recently available pre-titrated panels consisting of more than 100 or 200 proteins may benefit more from dimensionality reduction with PCA.

# set up dsb values to use in WNN analysis 
DefaultAssay(BM1_query) = "dsb"
prots = rownames(BM1_query@assays[["dsb"]])[1:8]
VariableFeatures(BM1_query) = prots

BM1_query = RunPCA(BM1_query, reduction.name = 'pdsb', features = VariableFeatures(BM1_query), verbose = FALSE)

# make matrix of normalized protein values to add as dr embeddings
pseudo = t(GetAssayData(BM1_query, slot = 'data', assay = 'dsb'))[,1:8]
colnames(pseudo) = paste('pseudo', 1:8, sep = "_")
pseudo <- as.matrix(pseudo)
BM1_query@reductions$pdsb@cell.embeddings = pseudo

# run WNN directly using dsb normalized values. 
BM1_query = FindMultiModalNeighbors(
  object = BM1_query,
  reduction.list = list("pca", "pdsb"),
  weighted.nn.name = "dsb_wnn", 
  knn.graph.name = "dsb_knn",
  modality.weight.name = "dsb_weight",
  snn.graph.name = "dsb_snn",
  dims.list = list(1:30, 1:8), 
  verbose = FALSE
)

# cluster based on the join RNA and protein graph
BM1_query = FindClusters(BM1_query, graph.name = "dsb_knn", 
                 algorithm = 3, 
                 resolution = 1.5,
                 random.seed = 1990, 
                 verbose = FALSE)

# Visualise renamed/merged clusters
dittoDimPlot(BM1_query, var = "dsb_knn_res.1.5",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

# Exploring the expression of canonical marker genes
Idents(BM1_query) <- 'dsb_knn_res.1.5'
VlnPlot(BM1_query, features = c("TFRC", "ITGA4", "GYPA"), assay = "RNA", slot = "data", sort = TRUE, ncol=1) + NoLegend()
VlnPlot(BM1_query, features = c("TFRC", "ITGA4", "GYPA"), assay = "SCT", slot = "data", sort = TRUE, ncol=1) + NoLegend()
VlnPlot(BM1_query, features = c("CD71-TotalSeqB", "CD44-TotalSeqB", "CD34-TotalSeqB"), assay = "dsb", slot = "data", sort = TRUE, ncol=1, pt.size=0.5) + NoLegend()

# 17.3- Defining marker genes for each cluster and visualizing aggregated per cluster expression

# create multimodal heatmap 
vf = BM1_query@assays[["SCT"]]@var.features

# find marker genes for the joint clusters 
Idents(BM1_query) = "dsb_knn_res.1.5"

DefaultAssay(BM1_query)  = "RNA"

rnade = FindAllMarkers(BM1_query, features = vf, only.pos = TRUE, verbose = FALSE)

gene_plot = rnade %>% 
  dplyr::filter(avg_log2FC > 1 ) %>%  
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(10) %$% gene %>% unique 

cite_data = GetAssayData(BM1_query,slot = 'data',assay = 'dsb') %>% t()
rna_subset = GetAssayData(BM1_query,assay = 'RNA',slot = 'data')[gene_plot, ] %>%
  as.data.frame() %>% 
  t() %>% 
  as.matrix()

# combine into dataframe 
d = cbind(BM1_query@meta.data, cite_data, rna_subset) 

# define proteins to use in analysis (non-isptype controls)
prots = rownames(BM1_query@assays[["dsb"]])[1:8]

# calculate the median protein expression per cluster
dat_plot = d %>% 
  dplyr::group_by(dsb_knn_res.1.5) %>% 
  dplyr::summarize_at(.vars = c(prots, gene_plot), .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("dsb_knn_res.1.5") 

# Next we can make a joint heatmap of protein and mRNA expression. 
# It is helpful to visualize protein levels on the normalized dsb scale which was also used to cluster cells with relative mRNA expression as a z score.
# make a combined plot 
suppressMessages(library(ComplexHeatmap)); ht_opt$message = FALSE

# protein heatmap 
prot_col = circlize::colorRamp2(breaks = seq(-1,25, by = 1), 
                                colors = viridis::viridis(n = 27, option = "B"))



# plot a heatmap of the average dsb normalized values for each cluster
p1 = pheatmap::pheatmap(t(dat_plot)[prots, ], name = "protein", 
                   color = viridis::viridis(25, option = "B"), 
                   fontsize_row = 8, border_color = NA)
# mRNA heatmap 
p2 = pheatmap::pheatmap(t(dat_plot)[gene_plot, ], name = "mRNA", 
                        color = viridis::viridis(25, option = "B"), 
                        fontsize_row = 5, border_color = NA)

# combine heatmaps 
library(ComplexHeatmap)

prot_col = circlize::colorRamp2(breaks = seq(-1,25, by = 1), 
                                colors = viridis::viridis(n = 27, option = "B"))
p3 = Heatmap(t(dat_plot)[prots, ], 
             name = "protein", 
             col = prot_col, 
             use_raster = T,
             row_names_gp = gpar(color = "black", fontsize = 5))


# mRNA heatmap 
mrna = t(dat_plot)[gene_plot, ]
rna_col = circlize::colorRamp2(breaks = c(-2,-1,0,1,2), 
                               colors = colorspace::diverge_hsv(n = 5))
p4 = Heatmap(t(dat_plot)[gene_plot, ], 
             name = "mRNA", 
             col = rna_col,
             use_raster = T, 
             column_names_gp = gpar(color = "black", fontsize = 7), 
             row_names_gp = gpar(color = "black", fontsize = 5))

ht_list = p3 %v% p4
saveRDS(ht_list, "ht_list.RDS")
ht_list <- readRDS("ht_list.RDS")
draw(ht_list)

# 18- Weighted Nearest Neighbor multimodal clustering using dsb normalized values - Seurat vignette

# Reference: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-cite-seq-rna-adt-1
# Similar to the method 1 above
DefaultAssay(BM1_query) = "dsb"

prots = rownames(BM1_query@assays[["dsb"]])[1:8]
VariableFeatures(BM1_query) = prots

BM1_query = BM1_query %>% 
  ScaleData() %>% 
  RunPCA(reduction.name = 'apca', verbose = FALSE)

# For each cell, we calculate its closest neighbors in the dataset based on a weighted combination of RNA and protein similarities. 
# The cell-specific modality weights and multimodal neighbors are calculated in a single function, which takes ~2 minutes to run on this dataset. 
# We specify the dimensionality of each modality (similar to specifying the number of PCs to include in scRNA-seq clustering), but you can vary these settings to see that small changes have minimal effect on the overall results.

# Identify multimodal neighbors. These will be stored in the neighbors slot, 
# and can be accessed using BM1_query[['weighted.nn']]
# The WNN graph can be accessed at BM1_query[["wknn"]], 
# and the SNN graph used for clustering at BM1_query[["wsnn"]]

# Cell-specific modality weights can be accessed at BM1_query$RNA.weight
BM1_query = FindMultiModalNeighbors(
  BM1_query, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:7), 
  modality.weight.name = "RNA.weight", 
  verbose = FALSE
)

# Warning message:
# In FindMultiModalNeighbors(BM1_query, reduction.list = list("pca",  :
# The number of provided modality.weight.name is not equal to the number of modalities. 
# SCT.weight dsb.weight are used to store the modality weights

# We can now use these results for downstream analysis, such as visualization and clustering. 
# For example, we can create a UMAP visualization of the data based on a weighted combination of RNA and protein data.
# We can also perform graph-based clustering and visualize these results on the UMAP, alongside a set of cell annotations.
BM1_query <- RunUMAP(BM1_query, nn.name = "weighted.nn", 
                     reduction.name = "wnn.umap", 
                     reduction.key = "wnnUMAP_")

BM1_query <- FindClusters(BM1_query, graph.name = "wsnn", 
                          algorithm = 3, 
                          resolution = 1.5, 
                          verbose = FALSE, 
                          random.seed = 1990)

# Visualise WNN clusters in the HCA reference dataset UMAP
dittoDimPlot(BM1_query, var = "wsnn_res.1.5",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

# Visualise WNN clusters in the WNN (ADT+RNA) UMAP projection
dittoDimPlot(BM1_query, var = "wsnn_res.1.5",
             reduction.use = "wnn.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())


# Visualise predicted clusters from the HCA reference dataset in the WNN (ADT+RNA) UMAP projection
dittoDimPlot(BM1_query, var = "predicted.celltype.l2",
             reduction.use = "wnn.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

p1 <- DimPlot(BM1_query, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(BM1_query, reduction = 'wnn.umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2

# We can also compute UMAP visualization based on only the RNA and protein data and compare. 
# We find that the RNA analysis is more informative than the ADT analysis in identifying progenitor states (the ADT panel contains markers for differentiated cells),
# while the converse is true of T cell states (where the ADT analysis outperforms RNA).
BM1_query <- RunUMAP(BM1_query, reduction = 'pca', dims = 1:30, assay = 'RNA', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
BM1_query <- RunUMAP(BM1_query, reduction = 'apca', dims = 1:7, assay = 'dsb', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

dittoDimPlot(BM1_query, var = "predicted.celltype.l2",
             reduction.use = "rna.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

dittoDimPlot(BM1_query, var = "predicted.celltype.l2",
             reduction.use = "adt.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

p3 <- DimPlot(BM1_query, reduction = 'rna.umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(BM1_query, reduction = 'adt.umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p3 + p4

# Visualisation of proteins in WNN UMAP
DefaultAssay(BM1_query) <- "dsb"
# Pick proteins
proteins <- c("CD3-TotalSeqB","CD11b-TotalSeqB","CD20-TotalSeqB","CD31-TotalSeqB","CD34-TotalSeqB","CD44-TotalSeqB",
              "CD71-TotalSeqB", "CD49d-TotalSeqB")
multi_dittoDimPlot(BM1_query, vars = proteins, reduction.use = "wnn.umap", size = 0.25,
                   do.label = TRUE, assay = "dsb",
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN UMAP")

# WNN RNA UMAP
multi_dittoDimPlot(BM1_query, vars = proteins, reduction.use = "rna.umap", size = 0.25,
                   do.label = TRUE, assay = "dsb",
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN RNA UMAP")

# WNN ADT UMAP
multi_dittoDimPlot(BM1_query, vars = proteins, reduction.use = "adt.umap", size = 0.25,
                   do.label = TRUE, assay = "dsb",
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN ADT UMAP")

# We can visualize the expression of canonical marker genes and proteins on the multimodal UMAP, which can assist in verifying the provided annotations:

# Finally, we can visualize the modality weights that were learned for each cell. 
# Each of the populations with the highest RNA weights represent progenitor cells, while the populations with the highest protein weights represent T cells. 
# This is in line with our biological expectations, as the antibody panel does not contain markers that can distinguish between different progenitor populations.
ditto_VlnPlot(BM1_query, features = "SCT.weight", group.by = 'predicted.celltype.l2', sort = TRUE, pt.size = 0.1) + NoLegend()
VlnPlot(BM1_query, features = "dsb.weight", group.by = 'predicted.celltype.l2', sort = TRUE, pt.size = 0.1) + NoLegend()

# 19- Weighted Nearest Neighbor multimodal clustering using predicted ADT values
DefaultAssay(BM1_query) <- 'predicted_ADT'

BM1_query = BM1_query %>% 
  ScaleData() %>% 
  RunPCA(reduction.name = 'apca.predicted.adt', verbose = FALSE)

# Cell-specific modality weights can be accessed at BM1_query$RNA.weight
BM1_query = FindMultiModalNeighbors(
  BM1_query, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:7), knn.graph.name = "predicted.ADT.wknn",
  snn.graph.name = "predicted.ADT.wsnn",
  weighted.nn.name = "predicted.ADT.weighted.nn",
  modality.weight.name = "predicted.ADT.weight", 
  verbose = FALSE
)

BM1_query <- RunUMAP(BM1_query, nn.name = "predicted.ADT.weighted.nn", 
                     reduction.name = "predicted.ADT.wnn.umap", 
                     reduction.key = "predicted.ADT.wnnUMAP_")

BM1_query <- FindClusters(BM1_query, graph.name = "predicted.ADT.wsnn", 
                          algorithm = 3, 
                          resolution = 1.5, 
                          verbose = FALSE, 
                          random.seed = 1990)

# Visualise WNN clusters in the HCA reference dataset UMAP
dittoDimPlot(BM1_query, var = "predicted.ADT.wsnn_res.1.5",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

# Visualise WNN clusters in the WNN (ADT+RNA) UMAP projection
dittoDimPlot(BM1_query, var = "predicted.ADT.wsnn_res.1.5",
             reduction.use = "predicted.ADT.wnn.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

# Visualise predicted clusters from the HCA reference dataset in the WNN (ADT+RNA) UMAP projection
dittoDimPlot(BM1_query, var = "predicted.celltype.l2",
             reduction.use = "predicted.ADT.wnn.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

# We can also compute UMAP visualization based on only the RNA and protein data and compare. 
# We find that the RNA analysis is more informative than the ADT analysis in identifying progenitor states (the ADT panel contains markers for differentiated cells),
# while the converse is true of T cell states (where the ADT analysis outperforms RNA).
BM1_query <- RunUMAP(BM1_query, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                     reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
BM1_query <- RunUMAP(BM1_query, reduction = 'apca', dims = 1:7, assay = 'predicted_ADT', 
                     reduction.name = 'predicted.adt.umap', reduction.key = 'predictedADTUMAP_')

dittoDimPlot(BM1_query, var = "predicted.celltype.l2",
             reduction.use = "rna.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

dittoDimPlot(BM1_query, var = "predicted.celltype.l2",
             reduction.use = "predicted.adt.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

p3 <- DimPlot(BM1_query, reduction = 'rna.umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(BM1_query, reduction = 'predicted.adt.umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p3 + p4

# Visualisation of proteins in WNN UMAP
# Pick proteins

# Pick proteins
proteins <- c("CD34", "CD11a","CD11c","CD123","CD14","CD16","CD161", "CD56", "CD57", "CD3", "CD4", "CD8a", "CD127-IL7Ra", "CD197-CCR7", "CD25","CD27", "CD278-ICOS", "CD28", "CD69", "CD19","CD79b","HLA.DR","CD45RA","CD45RO", "CD38")
proteins1 <- c("CD34", "CD11a","CD11c","CD123","CD14","CD16","CD161", "CD56", "CD57")
proteins2 <- c("CD3", "CD4", "CD8a", "CD127-IL7Ra", "CD197-CCR7", "CD25","CD27", "CD278-ICOS", "CD28", "CD69")
proteins3 <- c("CD19","CD79b","HLA.DR","CD45RA","CD45RO", "CD38")

# WNN UMAP
multi_dittoDimPlot(BM1_query, vars = proteins1, reduction.use = "predicted.ADT.wnn.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

multi_dittoDimPlot(BM1_query, vars = proteins2, reduction.use = "predicted.ADT.wnn.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

multi_dittoDimPlot(BM1_query, vars = proteins3, reduction.use = "predicted.ADT.wnn.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

# WNN RNA UMAP
multi_dittoDimPlot(BM1_query, vars = proteins1, reduction.use = "rna.umap", size = 0.25,
                   do.label = TRUE,
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN RNA UMAP")

multi_dittoDimPlot(BM1_query, vars = proteins2, reduction.use = "rna.umap", size = 0.25,
                   do.label = TRUE,
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN RNA UMAP")

multi_dittoDimPlot(BM1_query, vars = proteins3, reduction.use = "rna.umap", size = 0.25,
                   do.label = TRUE,
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN RNA UMAP")

# WNN ADT UMAP
multi_dittoDimPlot(BM1_query, vars = proteins1, reduction.use = "adt.umap", size = 0.25,
                   do.label = TRUE,
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN RNA UMAP")

multi_dittoDimPlot(BM1_query, vars = proteins2, reduction.use = "adt.umap", size = 0.25,
                   do.label = TRUE,
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN RNA UMAP")

multi_dittoDimPlot(BM1_query, vars = proteins3, reduction.use = "adt.umap", size = 0.25,
                   do.label = TRUE,
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN RNA UMAP")




# 20- Subset Seurat object to cluster major cell types and re-add the info to the original seurat object

# reference: https://github.com/satijalab/seurat/issues/1748
# Subset BM1_query object based on the clusters and sub-cluster
BM1_eryth <- subset(BM1_query, idents = c("Late Eryth", "Early Eryth"))
BM1_eryth <- FindNeighbors(BM1_eryth, dims = 1:10, k.param = 5)
BM1_eryth <- FindClusters(BM1_eryth) 

# Generate a new column called sub_cluster in the metadata
BM1_query$sub_cluster <- as.character(Idents(BM1_query))

# Change the information of cells containing sub-cluster information
WhichCells(BM1_query, c("Late Eryth", "Early Eryth"))
rownames(BM1_eryth@meta.data) #to get cells names from BM1_eryth

BM1_query$sub_cluster[Cells(BM1_eryth)] <- paste("BM1_eryth",Idents(BM1_eryth))
DimPlot(BM1_query, group.by = "sub_cluster")





