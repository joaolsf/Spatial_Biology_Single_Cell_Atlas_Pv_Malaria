#start by reading in the data. The Read10X() function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix. 
#The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).
#We next use the count matrix to create a Seurat object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset. 
#For a technical discussion of the Seurat object structure, check out our GitHub Wiki. For example, the count matrix is stored in pbmc[["RNA"]]@counts.
setwd("~/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/Seurat_analysis")

library(usethis) 
usethis::edit_r_environ()
rm(large_df, large_list, large_vector, temp_variables)

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

# 1- Loading data from 10X published experiments
# 1- Loading data from 10X multi-modal experiments

BM9.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_A/")
BM9 <- CreateSeuratObject(counts = BM9.data, project = "HD_BM9", min.cells = 3, min.features = 0)

BM10.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_B/")
BM10 <- CreateSeuratObject(counts = BM10.data, project = "HD_BM19", min.cells = 3, min.features = 0)

BM11.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_C1/")
BM11 <- CreateSeuratObject(counts = BM11.data, project = "HD_BM11", min.cells = 3, min.features = 0)

BM12.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_C2/")
BM12 <- CreateSeuratObject(counts = BM12.data, project = "HD_BM12", min.cells = 3, min.features = 0)

BM13.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_Ck/")
BM13 <- CreateSeuratObject(counts = BM13.data, project = "HD_BM13", min.cells = 3, min.features = 0)

BM14.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_E/")
BM14 <- CreateSeuratObject(counts = BM14.data, project = "HD_BM14", min.cells = 3, min.features = 0)

BM15.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_F/")
BM15 <- CreateSeuratObject(counts = BM15.data, project = "HD_BM15", min.cells = 3, min.features = 0)

BM16.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_G/")
BM16 <- CreateSeuratObject(counts = BM16.data, project = "HD_BM16", min.cells = 3, min.features = 0)

BM17.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_H/")
BM17 <- CreateSeuratObject(counts = BM17.data, project = "HD_BM17", min.cells = 3, min.features = 0)

BM18.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_J/")
BM18 <- CreateSeuratObject(counts = BM18.data, project = "HD_BM18", min.cells = 3, min.features = 0)

BM19.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_K/")
BM19 <- CreateSeuratObject(counts = BM19.data, project = "HD_BM19", min.cells = 3, min.features = 0)

BM20.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_L/")
BM20 <- CreateSeuratObject(counts = BM20.data, project = "HD_BM20", min.cells = 3, min.features = 0)

BM21.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_M/")
BM21 <- CreateSeuratObject(counts = BM21.data, project = "HD_BM21", min.cells = 3, min.features = 0)

BM22.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_N/")
BM22 <- CreateSeuratObject(counts = BM22.data, project = "HD_BM22", min.cells = 3, min.features = 0)

BM23.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_O/")
BM23 <- CreateSeuratObject(counts = BM23.data, project = "HD_BM23", min.cells = 3, min.features = 0)

BM24.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_P/")
BM24 <- CreateSeuratObject(counts = BM24.data, project = "HD_BM24", min.cells = 3, min.features = 0)

BM25.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_Q/")
BM25 <- CreateSeuratObject(counts = BM25.data, project = "HD_BM25", min.cells = 3, min.features = 0)

BM26.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_R/")
BM26 <- CreateSeuratObject(counts = BM26.data, project = "HD_BM26", min.cells = 3, min.features = 0)

BM27.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_S1/")
BM27 <- CreateSeuratObject(counts = BM27.data, project = "HD_BM27", min.cells = 3, min.features = 0)

BM28.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_S2/")
BM28 <- CreateSeuratObject(counts = BM28.data, project = "HD_BM28", min.cells = 3, min.features = 0)

BM29.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_Sk1/")
BM29 <- CreateSeuratObject(counts = BM29.data, project = "HD_BM29", min.cells = 3, min.features = 0)

BM30.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_Sk2/")
BM30 <- CreateSeuratObject(counts = BM30.data, project = "HD_BM30", min.cells = 3, min.features = 0)

BM31.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_T/")
BM31 <- CreateSeuratObject(counts = BM31.data, project = "HD_BM31", min.cells = 3, min.features = 0)

BM32.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_U/")
BM32 <- CreateSeuratObject(counts = BM32.data, project = "HD_BM32", min.cells = 3, min.features = 0)

BM33.data <- Read10X(data.dir = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/JCI_BM_samples/BM_W/")
BM33 <- CreateSeuratObject(counts = BM33.data, project = "HD_BM33", min.cells = 3, min.features = 0)


# 1.3- Metadata addition ----

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats - calculate and add the percentage of MT genes

BM34 <- readRDS("./RDS_files/Seurat_Multimodal_Reference.rds")

BM9[["percent.mt"]] <- PercentageFeatureSet(BM9, pattern = "^MT-")
BM10[["percent.mt"]] <- PercentageFeatureSet(BM10, pattern = "^MT-")
BM11[["percent.mt"]] <- PercentageFeatureSet(BM11, pattern = "^MT-")
BM12[["percent.mt"]] <- PercentageFeatureSet(BM12, pattern = "^MT-")
BM13[["percent.mt"]] <- PercentageFeatureSet(BM13, pattern = "^MT-")
BM14[["percent.mt"]] <- PercentageFeatureSet(BM14, pattern = "^MT-")
BM15[["percent.mt"]] <- PercentageFeatureSet(BM15, pattern = "^MT-")
BM16[["percent.mt"]] <- PercentageFeatureSet(BM16, pattern = "^MT-")
BM17[["percent.mt"]] <- PercentageFeatureSet(BM17, pattern = "^MT-")
BM18[["percent.mt"]] <- PercentageFeatureSet(BM18, pattern = "^MT-")
BM20[["percent.mt"]] <- PercentageFeatureSet(BM20, pattern = "^MT-")
BM21[["percent.mt"]] <- PercentageFeatureSet(BM21, pattern = "^MT-")
BM22[["percent.mt"]] <- PercentageFeatureSet(BM22, pattern = "^MT-")
BM23[["percent.mt"]] <- PercentageFeatureSet(BM23, pattern = "^MT-")
BM24[["percent.mt"]] <- PercentageFeatureSet(BM24, pattern = "^MT-")
BM25[["percent.mt"]] <- PercentageFeatureSet(BM25, pattern = "^MT-")
BM26[["percent.mt"]] <- PercentageFeatureSet(BM26, pattern = "^MT-")
BM27[["percent.mt"]] <- PercentageFeatureSet(BM27, pattern = "^MT-")
BM28[["percent.mt"]] <- PercentageFeatureSet(BM28, pattern = "^MT-")
BM29[["percent.mt"]] <- PercentageFeatureSet(BM29, pattern = "^MT-")
BM30[["percent.mt"]] <- PercentageFeatureSet(BM30, pattern = "^MT-")
BM31[["percent.mt"]] <- PercentageFeatureSet(BM31, pattern = "^MT-")
BM32[["percent.mt"]] <- PercentageFeatureSet(BM32, pattern = "^MT-")
BM33[["percent.mt"]] <- PercentageFeatureSet(BM33, pattern = "^MT-")
BM34[["percent.mt"]] <- PercentageFeatureSet(BM34, pattern = "^MT-")


# Calculation of the percentage and absolute numbers of human transcripts/UMIs
BM9[["UMI"]] <- BM9$nCount_RNA
BM10[["UMI"]] <- BM10$nCount_RNA
BM11[["UMI"]] <- BM11$nCount_RNA
BM12[["UMI"]] <- BM12$nCount_RNA
BM13[["UMI"]] <- BM13$nCount_RNA
BM14[["UMI"]] <- BM14$nCount_RNA
BM15[["UMI"]] <- BM15$nCount_RNA
BM16[["UMI"]] <- BM16$nCount_RNA
BM17[["UMI"]] <- BM17$nCount_RNA
BM18[["UMI"]] <- BM18$nCount_RNA
BM19[["UMI"]] <- BM19$nCount_RNA
BM20[["UMI"]] <- BM20$nCount_RNA
BM21[["UMI"]] <- BM21$nCount_RNA
BM22[["UMI"]] <- BM22$nCount_RNA
BM23[["UMI"]] <- BM23$nCount_RNA
BM24[["UMI"]] <- BM24$nCount_RNA
BM25[["UMI"]] <- BM25$nCount_RNA
BM26[["UMI"]] <- BM26$nCount_RNA
BM27[["UMI"]] <- BM27$nCount_RNA
BM28[["UMI"]] <- BM28$nCount_RNA
BM29[["UMI"]] <- BM29$nCount_RNA
BM30[["UMI"]] <- BM30$nCount_RNA
BM31[["UMI"]] <- BM31$nCount_RNA
BM32[["UMI"]] <- BM32$nCount_RNA
BM33[["UMI"]] <- BM33$nCount_RNA
BM34[["UMI"]] <- BM34$nCount_RNA

# Calculation of the number of  genes detected per cell
BM19[["Genes"]] <- BM9$nFeature_RNA 
BM10[["Genes"]] <- BM10$nFeature_RNA 
BM11[["Genes"]] <- BM11$nFeature_RNA 
BM12[["Genes"]] <- BM12$nFeature_RNA 
BM13[["Genes"]] <- BM13$nFeature_RNA 
BM14[["Genes"]] <- BM14$nFeature_RNA 
BM15[["Genes"]] <- BM15$nFeature_RNA 
BM16[["Genes"]] <- BM16$nFeature_RNA 
BM17[["Genes"]] <- BM17$nFeature_RNA 
BM18[["Genes"]] <- BM18$nFeature_RNA 
BM19[["Genes"]] <- BM19$nFeature_RNA 
BM20[["Genes"]] <- BM20$nFeature_RNA 
BM21[["Genes"]] <- BM21$nFeature_RNA 
BM22[["Genes"]] <- BM22$nFeature_RNA 
BM23[["Genes"]] <- BM23$nFeature_RNA 
BM24[["Genes"]] <- BM24$nFeature_RNA 
BM25[["Genes"]] <- BM25$nFeature_RNA 
BM26[["Genes"]] <- BM26$nFeature_RNA 
BM27[["Genes"]] <- BM27$nFeature_RNA 
BM28[["Genes"]] <- BM28$nFeature_RNA 
BM29[["Genes"]] <- BM29$nFeature_RNA 
BM30[["Genes"]] <- BM30$nFeature_RNA 
BM31[["Genes"]] <- BM31$nFeature_RNA 
BM32[["Genes"]] <- BM32$nFeature_RNA 
BM33[["Genes"]] <- BM33$nFeature_RNA 
BM34[["Genes"]] <- BM34$nFeature_RNA 

# add metadata - sample_id and experiment_date
BM9$sample_id = "HD_BM9"
BM10$sample_id = "HD_BM10"
BM11$sample_id = "HD_BM11"
BM12$sample_id = "HD_BM12"
BM13$sample_id = "HD_BM13"
BM14$sample_id = "HD_BM14"
BM15$sample_id = "HD_BM15"
BM16$sample_id = "HD_BM16"
BM17$sample_id = "HD_BM17"
BM18$sample_id = "HD_BM18"
BM19$sample_id = "HD_BM19"
BM20$sample_id = "HD_BM20"
BM21$sample_id = "HD_BM21"
BM22$sample_id = "HD_BM22"
BM23$sample_id = "HD_BM23"
BM24$sample_id = "HD_BM24"
BM25$sample_id = "HD_BM25"
BM26$sample_id = "HD_BM26"
BM27$sample_id = "HD_BM27"
BM28$sample_id = "HD_BM28"
BM29$sample_id = "HD_BM29"
BM30$sample_id = "HD_BM30"
BM31$sample_id = "HD_BM31"
BM32$sample_id = "HD_BM32"
BM33$sample_id = "HD_BM33"
BM34$sample_id = "HD_BM34"

BM9$experiment_date = "12_18"
BM10$experiment_date = "12_18"
BM11$experiment_date = "12_18"
BM12$experiment_date = "12_18"
BM13$experiment_date = "12_18"
BM14$experiment_date = "12_18"
BM15$experiment_date = "12_18"
BM16$experiment_date = "12_18"
BM17$experiment_date = "12_18"
BM18$experiment_date = "12_18"
BM19$experiment_date = "12_18"
BM20$experiment_date = "12_18"
BM21$experiment_date = "12_18"
BM22$experiment_date = "12_18"
BM23$experiment_date = "12_18"
BM24$experiment_date = "12_18"
BM25$experiment_date = "12_18"
BM26$experiment_date = "12_18"
BM27$experiment_date = "12_18"
BM28$experiment_date = "12_18"
BM29$experiment_date = "12_18"
BM30$experiment_date = "12_18"
BM31$experiment_date = "12_18"
BM32$experiment_date = "12_18"
BM33$experiment_date = "12_18"
BM34$experiment_date = "06_19"

# 1.4- Note that all operations below are performed on the RNA assay Set and verify that the
# default assay is RNA
DefaultAssay(BM1) <- "RNA"

# Data information ----
# Different ways to look at the data. Doesn't change it
head(BM1)
dim(BM1) # Gives you the numbers of genes (first number) and cells (second number) in your seurat object

saveRDS(BM9, "BM9.rds")
saveRDS(BM10, "BM10.rds")
saveRDS(BM11, "BM11.rds")
saveRDS(BM12, "BM12.rds")
saveRDS(BM13, "BM13.rds")
saveRDS(BM14, "BM14.rds")
saveRDS(BM15, "BM15.rds")
saveRDS(BM16, "BM16.rds")
saveRDS(BM17, "BM17.rds")
saveRDS(BM18, "BM18.rds")
saveRDS(BM19, "BM19.rds")
saveRDS(BM20, "BM20.rds")
saveRDS(BM21, "BM21.rds")
saveRDS(BM22, "BM22.rds")
saveRDS(BM23, "BM23.rds")
saveRDS(BM24, "BM24.rds")
saveRDS(BM25, "BM25.rds")
saveRDS(BM26, "BM26.rds")
saveRDS(BM27, "BM27.rds")
saveRDS(BM28, "BM28.rds")
saveRDS(BM29, "BM29.rds")
saveRDS(BM30, "BM30.rds")
saveRDS(BM31, "BM31.rds")
saveRDS(BM32, "BM32.rds")
saveRDS(BM33, "BM33.rds")
saveRDS(BM34, "BM34.rds")

BM1 <- readRDS("./RDS_files/BM1_dsb_normalised.rds")
BM2 <- readRDS("./RDS_files/BM2_dsb_normalised.rds")
BM3 <- readRDS("./RDS_files/BM3_dsb_normalised.rds")
BM4 <- readRDS("./RDS_files/BM4_dsb_normalised.rds")
BM5 <- readRDS("./RDS_files/BM5_dsb_normalised.rds")
BM6 <- readRDS("./RDS_files/BM6_dsb_normalised.rds")
BM7 <- readRDS("./RDS_files/BM7_dsb_normalised.rds")
BM8 <- readRDS("./RDS_files/BM8_dsb_normalised.rds")

BM9 <- readRDS("./RDS_files/BM9.rds")
BM10 <- readRDS("./RDS_files/BM10.rds")
BM11 <- readRDS("./RDS_files/BM11.rds")
BM12 <- readRDS("./RDS_files/BM12.rds")
BM13 <- readRDS("./RDS_files/BM13.rds")
BM14 <- readRDS("./RDS_files/BM14.rds")
BM15 <- readRDS("./RDS_files/BM15.rds")
BM16 <- readRDS("./RDS_files/BM16.rds")
BM17 <- readRDS("./RDS_files/BM17.rds")
BM18 <- readRDS("./RDS_files/BM18.rds")
BM19 <- readRDS("./RDS_files/BM19.rds")
BM20 <- readRDS("./RDS_files/BM20.rds")
BM21 <- readRDS("./RDS_files/BM21.rds")
BM22 <- readRDS("./RDS_files/BM22.rds")
BM23 <- readRDS("./RDS_files/BM23.rds")
BM24 <- readRDS("./RDS_files/BM24.rds")
BM25 <- readRDS("./RDS_files/BM25.rds")
BM26 <- readRDS("./RDS_files/BM26.rds")
BM27 <- readRDS("./RDS_files/BM27.rds")
BM28 <- readRDS("./RDS_files/BM28.rds")
BM29 <- readRDS("./RDS_files/BM29.rds")
BM30 <- readRDS("./RDS_files/BM30.rds")
BM31 <- readRDS("./RDS_files/BM31.rds")
BM32 <- readRDS("./RDS_files/BM32.rds")
BM33 <- readRDS("./RDS_files/BM33.rds")


# 3- Create one merged Seurat object before pre-processing and QC as we will use the same parameters across all samples

# Reference 1: https://satijalab.org/seurat/articles/merge_vignette.html
# Reference 2: https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets

# Add dataset labels as cell.ids just in case you have overlapping barcodes between the datasets. 
# Merge datasets into one single seurat object
BM_merged_data <- merge(BM1, c(BM2, BM3, BM4, BM5, BM6, BM7, BM8, BM9, BM10, BM11, BM12, BM13, BM14, BM15, BM16,
                        BM17, BM18, BM19, BM20, BM21, BM22, BM23, BM24, BM25, BM26, BM27, BM28, BM29, BM30, BM31, BM32, BM33), 
                        add.cell.ids = c("HD_BM1.1", "HD_BM2.1", "HD_BM3.1", "HD_BM4.1", "HD_BM1.2", "HD_BM2.2", "HD_BM3.2", "HD_BM4.2","HD_BM9",
                                         "HD_BM10", "HD_BM11", "HD_BM12", "HD_BM13", "HD_BM14", "HD_BM15", "HD_BM16", "HD_BM17", "HD_BM18", "HD_BM19", "HD_BM20", "HD_BM21", "HD_BM22", "HD_BM23","HD_BM24", "HD_BM25", "HD_BM26", "HD_BM27", "HD_BM28", "HD_BM29", "HD_BM30", "HD_BM31", "HD_BM32", "HD_BM33"), project = "HD_BM")

saveRDS(BM_merged_data, "BM_merged_data.rds")
BM_merged_data <- readRDS("./RDS_files/Merged_datasets/BM_merged_data.rds")

# 3- Standard pre-processing workflow

# The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. 
# These represent the selection and filtration of cells based on QC metrics, 
# data normalization and scaling, and the detection of highly variable features.

# 3.1- Calculate QC

# Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. 
# A few QC metrics commonly used by the community include:
# The number of unique genes detected in each cell.
# Low-quality cells or empty droplets will often have very few genes
# Cell doublets or multiplets may exhibit an aberrantly high gene count

# Calculate and add the percentage of ribossomal genes
BM_merged_data <- PercentageFeatureSet(BM_merged_data, "^RP[SL]", col.name = "percent_ribo")

# Calculate proportion hemoglobin genes
BM_merged_data <- PercentageFeatureSet(BM_merged_data, "^HB[^(P)]", col.name = "percent_hb")

# Visualize QC metrics as a violin plot before filtering
feats <- c("Genes", "UMI", "percent.mt", "percent_ribo", "percent_hb")
VlnPlot(BM_merged_data, features = feats, pt.size = 0.01, ncol = 3, raster=FALSE) +
  NoLegend()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(BM_merged_data, feature1 = "UMI", feature2 = "percent.mt", raster=FALSE)
plot2 <- FeatureScatter(BM_merged_data, feature1 = "UMI", feature2 = "Genes", raster=FALSE)
plot1 + plot2

# 3.2- Filtering

# 3.2.1- Detection-based filtering and removal of cells with high % of mitochondrial genes (optional: high % of ribossomal genes)

# A standard approach is to filter cells with low amount of reads as well as genes that are present in at least a certain amount of cells. 
# Here we will only consider cells with at least 200 detected genes. 
# Please note that those values are highly dependent on the library preparation method used.
BM_merged_subset <- subset(BM_merged_data, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
head(BM_merged_subset)
dim(BM_merged_subset) # Gives you the numbers of genes (first number) and cells (second number) in your seurat object
head(BM_subset) 
unique(sapply(X = strsplit(colnames(BM_subset), split = "_"), FUN = "[", 1))
table(BM_subset$orig.ident)

# 3.2.2- Filter genes

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
# BM_merged_subset <- BM_merged_subset[!grepl("MALAT1", rownames(BM_merged_subset)), ]

# Filter Mitocondrial
# BM_subset2 <- BM_subset2[!grepl("^MT-", rownames(BM_subset)), ]
# Filter Ribossomal gene (optional if that is a problem on your data) data.filt
# <- data.filt[ ! grepl('^RP[SL]', rownames(data.filt)), ]
# Filter Hemoglobin gene (optional if that is a problem on your data)
#data.filt <- data.filt[!grepl("^HB[^(P)]", rownames(data.filt)), ]
dim(BM_merged_subset)

# Plot filtered QC metrics

feats <- c("Genes", "UMI", "percent.mt", "percent_ribo", "percent_hb")
VlnPlot(BM_merged_subset, features = feats, pt.size = 0.1, ncol = 3, raster=FALSE) + NoLegend()
plot1 <- FeatureScatter(BM_merged_subset, feature1 = "UMI", feature2 = "percent.mt", raster=FALSE)
plot2 <- FeatureScatter(BM_merged_subset, feature1 = "UMI", feature2 = "Genes", raster=FALSE)
plot1 + plot2

saveRDS(BM_merged_subset, "BM_merged_subset.rds")

# 4- Normalize the data

# After removing unwanted cells from the dataset, the next step is to normalize the data. 
# By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
# Normalized values are stored in BM1_subset[["RNA"]]@data.

BM_merged_subset <- NormalizeData(BM_merged_subset, normalization.method = "LogNormalize", scale.factor = 10000)

# 4.1- Identification of highly variable features (feature selection)

# We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others).
# We and others have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.
# By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.

BM_merged_subset <- FindVariableFeatures(BM_merged_subset, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(BM_merged_subset), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(BM_merged_subset)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# 4.2- Scaling the data

# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:
# Shifts the expression of each gene, so that the mean expression across cells is 0
# Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
# The results of this are stored in pbmc[["RNA"]]@scale.data
(future.globals.maxSize = 4000 * 1024^2) 
(future.globals.maxSize = 4000 * 1024^8)

#all.genes <- rownames(BM_merged_subset)
BM_merged_subset <- ScaleData(BM_merged_subset) #scale/center only the variable features.

# 8- SCTransform

# Different algorithm to normalise and scale the data, the results are stored in a different assay, so our previous normalisation is not overwritten
# Biological heterogeneity in single-cell RNA-seq data is often confounded by technical factors including sequencing depth.
# This procedure omits the need for heuristic steps including pseudocount addition or log-transformation and improves common downstream analytical tasks 
# such as variable gene selection, dimensional reduction, and differential expression.
# Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
# Transformed data will be available in the SCT assay, which is set as the default after running sctransform
# During normalization, we can also remove confounding sources of variation, for example, mitochondrial mapping percentage
# The latest version of sctransform also supports using glmGamPoi package which substantially improves the speed of the learning procedure. It can be invoked by specifying method="glmGamPoi".

library(sctransform)
BM_merged_subset <- SCTransform(BM_merged_subset, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

# 8.1- Perform linear dimensional reduction

# Calculate PCs using variable features determined by SCTransform (3000 by default)
BM_merged_subset<- RunPCA(BM_merged_subset, assay = "SCT", npcs = 50)

# Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction(), DimPlot(), and DimHeatmap()
VizDimLoadings(BM_merged_subset, dims = 1:2, reduction = "pca")
DimPlot(BM_merged_subset, reduction = "pca")

# In particular DimHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. 
# Both cells and features are ordered according to their PCA scores. 
# Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. 
# Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated feature sets.
DimHeatmap(BM_merged_subset, dims = 1, cells = 1000, balanced = TRUE)
DimHeatmap(BM_merged_subset, dims = 1:15, cells = 500, balanced = TRUE)

# 8.3- Run non-linear dimensional reduction (UMAP/tSNE) before integration

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
BM_merged_subset <- RunUMAP(BM_merged_subset, dims = 1:30)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual samples
DimPlot(BM_merged_subset, reduction = "umap", label = FALSE) #try other ploting packages and explore the function

dittoDimPlot(BM_merged_subset, var = "sample_id",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

# individual 10X runs
dittoDimPlot(BM_merged_subset, var = "experiment_date",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("10X runs on UMAP")

saveRDS(BM_merged_subset, "./RDS_files/BM_merged_subset.rds")

BM_merged_subset <- readRDS("./RDS_files/Merged_datasets/BM_merged_subset.rds")



# 9- Harmony integration
BM_merged_integrated <- RunHarmony(BM_merged_subset, 
                            group.by.vars = c("sample_id", "experiment_date"), 
                            reduction = "pca", assay.use = "SCT", reduction.save = "harmony", kmeans_init_nstart=20, kmeans_init_iter_max=5000, epsilon.cluster = -Inf,
                            epsilon.harmony = -Inf)

BM_merged_integrated <- RunUMAP(BM_merged_integrated, reduction = "harmony", assay = "SCT", dims = 1:50)

# individual samples
DimPlot(BM_merged_integrated, reduction = "umap", label = FALSE) #try other ploting packages and explore the function

dittoDimPlot(BM_merged_integrated, var = "sample_id",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Samples on UMAP")

# individual 10X runs
dittoDimPlot(BM_merged_integrated, var = "experiment_date",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("10X runs on UMAP")

saveRDS(BM_merged_integrated, "BM_merged_integrated.rds")








