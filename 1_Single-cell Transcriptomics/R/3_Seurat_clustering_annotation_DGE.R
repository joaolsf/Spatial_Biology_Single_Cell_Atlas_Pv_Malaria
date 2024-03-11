
#start by reading in the data. The Read10X() function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix. 
#The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).
#We next use the count matrix to create a Seurat object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset. 
#For a technical discussion of the Seurat object structure, check out our GitHub Wiki. For example, the count matrix is stored in pbmc[["RNA"]]@counts.
setwd("~/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/Seurat_analysis")

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
library(EnhancedVolcano)

BM_integrated <- readRDS("./RDS_files/Lonza_dataset/BM_integrated.rds")

saveRDS(BM_integrated, "./RDS_files/BM_integrated.rds")

saveRDS(BM_query, "./RDS_files/Lonza_dataset/BM_query.rds")

BM_query <- readRDS("./RDS_files/Lonza_dataset/BM_query.rds")


# 7- Clustering analysis

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

saveRDS(BM_integrated, "./RDS_files/Merged_datasets/BM_merged_integrated.rds")

# 7.1- Finding differentially expressed features (cluster biomarkers) - marker genes using the SCT slot

# find markers that define clusters via differential expression. 
# By default, it identifies positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. 
# FindAllMarkers() automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells.

# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, 
# and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. 
# You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are unlikely to be highly discriminatory. 
# As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. 
# While there is generally going to be a loss in power, the speed increases can be significant and the most highly differentially expressed features will likely still rise to the top.

# find markers for every cluster compared to all remaining cells, report only the positive ones
DefaultAssay(BM_integrated) <- "SCT" 
BM_integrated.markers <- FindAllMarkers(BM_integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
clipr::write_clip(BM_integrated.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste


# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
BM_integrated.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap(BM_integrated, features = top10$gene, size = 2) + NoLegend() + theme(axis.text.y = element_text(size = 3))
# changing the default color
# DoHeatmap(BM_integrated, features = top10$gene, size = 2) + NoLegend() + theme(axis.text.y = element_text(size = 3)) +
   # scale_fill_viridis()

# Heatmap visualization - DittoHeatmap
dittoHeatmap(BM_integrated, genes = top10$gene,
             assay = "SCT", order.by = c("seurat_clusters"), cluster_rows = FALSE,
             cluster_cols = FALSE, scale = "none",
             heatmap.colors = viridis(100), 
             annot.by = c("seurat_clusters"), fontsize_row = 2)

# Quantifications
# Code to quantify the numbers of cells in each cluster and the proportion of cells 
table(BM_integrated$seurat_clusters)
prop.table(table(BM_integrated$seurat_clusters))

# median number of RNA molecules per cluster
tapply(BM_integrated$nCount_RNA, BM_integrated$seurat_clusters,  median)

saveRDS(BM_integrated.markers, "BM_integrated.markers_SCT.rds")

# 7.2- Marker gene Identification using RNA slot 

# Calculates the genes upregulated in each cluster. It's important that you switch to the RNA assay for that   

DefaultAssay(BM_integrated) <- "RNA" 
BM_integrated.marker2 <- FindAllMarkers(object = BM_integrated, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
clipr::write_clip(BM_integrated.marker2) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# Generate a list of marker genes that you would like to visualise
Marker <- c("TFRC", "GYPA")  

# Different plotting options. Feel free to play around
DotPlot(BM_integrated, features = Marker)
VlnPlot(BM_integrated, features = Marker)
FeaturePlot(BM_integrated, Marker, cols = viridis(100, direction = -1)) # cols changes the colors of the plot. I really like viridis for those heatmaps


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

# 8.2.1- Download the Azimuth reference and extract the archive
# Load the reference - it has to be a list object!
# Change the file path based on where the reference is located on your system.
reference <- LoadReference(path = "/Users/joaoluizsfilho/Library/CloudStorage/Dropbox/Work_Files/Matthias_Lab/Projects/Vivax_human_project/Results/Prospective Study 1/BM_lonza/scRNAseq/Seurat_analysis/Reference_datasets/HCA") 

# 8.2.2- Load the query object for mapping
# Change the file path based on where the query file is located on your system.
# query <- LoadFileInput(path = "BM1.rds") # it is already loaded
query <- BM_integrated

query <- ConvertGeneNames(
  object = query,
  reference.names = rownames(x = reference$map),
  homolog.table = 'https://seurat.nygenome.org/azimuth/references/homologs.rds'
)

# 8.2.3- Calculate nCount_RNA and nFeature_RNA if the query does not contain them already
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

# 8.2.4- Calculate percent mitochondrial genes if the query contains genes
# matching the regular expression "^MT-"
if (any(grepl(pattern = '^MT-', x = rownames(x = query)))) {
  query <- PercentageFeatureSet(
    object = query,
    pattern = '^MT-',
    col.name = 'percent.mt',
    assay = "RNA"
  )
}

# 8.2.5- Filter cells based on the thresholds for nCount_RNA and nFeature_RNA
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

# 8.2.6- Remove filtered cells from the query - DO NOT SKIP THIS STEP
query <- query[, cells.use]

# 8.2.7- Preprocess with SCTransform
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

# 8.2.8- Find anchors between query and reference
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

# 8.2.9- Transfer cell type labels and impute protein expression
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

# Plot the samples in the UMAP projection from the ref dataset
p1 <- DimPlot(query, reduction = "proj.umap", group.by = "orig.ident")
p1

# DimPlot of the reference
DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE) + NoLegend()

p2 <- dittoDimPlot(reference$plot, var = "celltype.l2",
                   reduction.use = "refUMAP", size = 0.25,
                   do.label = TRUE, labels.size = 3, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")
p2

# DimPlot of the query, colored by predicted cell type
DimPlot(object = query, reduction = "proj.umap", group.by = predicted.id, label = TRUE) + NoLegend()

p3 <- dittoDimPlot(query, var = "predicted.celltype.l2",
                   reduction.use = "proj.umap", size = 0.75,
                   do.label = TRUE, labels.size = 3, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")
p3

# Plot cell numbers and frequencies by sample
p4 <- dittoBarPlot(BM_query, var = "predicted.celltype.l2", group.by = "orig.ident",  scale = c("percent")) 
#+ theme_classic(base_size = 8) + theme(axis.text.x = element_text(size = 5))
p4

p4 <- dittoBarPlot(query, var = "predicted.celltype.l2", group.by = "orig.ident",  scale = c("count")) 
#+ theme_classic(base_size = 8) + theme(axis.text.x = element_text(size = 5))
p4

# Plot the score for the predicted cell type of the query
FeaturePlot(object = query, features = paste0(predicted.id, ".score"), reduction = "proj.umap")
VlnPlot(object = query, features = paste0(predicted.id, ".score"), group.by = predicted.id) + NoLegend()

# Plot the mapping score
FeaturePlot(object = query, features = "mapping.score", reduction = "proj.umap")
VlnPlot(object = query, features = "mapping.score", group.by = predicted.id) + NoLegend()

# cell type prediction scores
DefaultAssay(query) <- 'prediction.score.celltype.l2'
p5 <- FeaturePlot(query, features = c("CD14 Mono", "CD4 Naive", "CD4 Memory", "Memory B"), reduction = 'proj.umap', ncol = 2, 
                  cols = c("lightgrey", "darkred"))
p5

p6 <- FeaturePlot(query, features = c("MAIT", "NK", "HSC", "GMP"), reduction = 'proj.umap', ncol =2, 
                  cols = c("lightgrey", "darkred"))
p6

p7 <- FeaturePlot(query, features = c("EMP", "Early Eryth", "Late Eryth", "Plasma"), reduction = 'proj.umap', ncol =2, 
                  cols = c("lightgrey", "darkred"))
p7

# Plot an RNA feature
FeaturePlot(object = query, features = "TFRC", reduction = "proj.umap")
VlnPlot(object = query, features = "TFRC", group.by = predicted.id, sort = TRUE) + NoLegend()

FeaturePlot(object = query, features = "GYPA", reduction = "proj.umap")
VlnPlot(object = query, features = "GYPA", group.by = predicted.id, sort = TRUE) + NoLegend()

# Visualisation of dsb normalized proteins 
DefaultAssay(BM_query) <- "dsb"

# Pick proteins
proteins <- c("CD34", "CD31", "CD3", "CD20", "CD11b", "CD44", "CD71", "CD49d")

# Visualisation of specific proteins in UMAP with multidittoDimPlot
multi_dittoDimPlot(BM_query, vars = proteins, reduction.use = "proj.umap", size = 0.25,
                   do.label = TRUE,
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

BM_query <- query

saveRDS(BM_query, "./RDS_files/Merged_datasets/BM_merged_query")

BM_query <- readRDS("./RDS_files/Merged_datasets/BM_merged_query")

################################################################################################################################################################################

# 8.3- Map to any published dataset 

# 8.3.1- Multimodal reference mapping - Option 1
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
#InstallData("bmcite")
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

# Query dataset mapping/projection onto the reference dataset

# We then find anchors between each donor query dataset and the multimodal reference. 
# This command is optimized to minimize mapping time, by passing in a pre-computed set of reference neighbors, and turning off anchor filtration.
DefaultAssay(BM_query) <- "RNA" # make sure to change the default assat to match the reference default assay

# Command below if mapping several samples in a loop - when they are not integrated but only merged in the Seurat object - SKIP THIS
#anchors <- list()
#for (i in 1:length(hcabm40k.batches)) {
 # anchors[[i]] <- FindTransferAnchors(
  #  reference = bm,
   # query = hcabm40k.batches[[i]],
  #  k.filter = NA,
  #  reference.reduction = "spca", 
  #  reference.neighbors = "spca.annoy.neighbors", 
  #  dims = 1:50
#  )
# }
# Adapt command to run each sample individually or integrated Seurat object
anchors <- FindTransferAnchors(
    reference = bm,
    query = BM_query,
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )

# We then individually map each of the datasets.
# Command below if mapping several samples in a loop
# for (i in 1:length(hcabm40k.batches)) {
  # hcabm40k.batches[[i]] <- MapQuery(
   # anchorset = anchors[[i]], 
  #  query = hcabm40k.batches[[i]],
  #  reference = bm, 
  #  refdata = list(
  #    celltype = "celltype.l2", 
  #    predicted_ADT = "ADT"),
#    reference.reduction = "spca",
#    reduction.model = "wnn.umap"
#  )
# }

# Adapt command to run each sample individually or integrated Seurat object
BM_query <- MapQuery(
    anchorset = anchors, 
    query = BM_query,
    reference = bm, 
    refdata = list(
      celltype = "celltype.l2", 
      predicted_ADT = "ADT"),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )

rownames(BM_query@assays[["predicted_ADT"]]) #full list of ADT features in the imputed ADT SLOT
 
# VISUALISATIONS 

# Project annotated cell types on the reference dataset UMAP projection
p1 <- DimPlot(BM_query, reduction = 'ref.umap', group.by = 'predicted.celltype', label = TRUE) + NoLegend()
p1

p2 <- dittoDimPlot(bm, var = "celltype.l2",
                   reduction.use = "wnn.umap", size = 0.25,
                   do.label = TRUE, labels.size = 3, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")
p2

p3 <- dittoDimPlot(BM_query, var = "predicted.celltype",
                   reduction.use = "ref.umap", size = 0.75,
                   do.label = TRUE, labels.size = 3, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")
p3

# If mapping several samples at once: We can also merge all the objects into one dataset. Note that they have all been integrated into a common space, defined by the reference. 
# We can then visualize the results together.
# Merge the batches 
#hcabm40k <- merge(hcabm40k.batches[[1]], hcabm40k.batches[2:length(hcabm40k.batches)], merge.dr = "ref.umap")
#DimPlot(hcabm40k, reduction = "ref.umap", group.by =  "predicted.celltype", label = TRUE, repel = TRUE, label.size = 3) + NoLegend()

# cell type prediction scores
# Plot the score for the predicted cell type of the query
FeaturePlot(object = BM_query, features = "predicted.celltype.score", reduction = "ref.umap")
VlnPlot(object = BM_query, features = "predicted.celltype.score", group.by = "predicted.celltype", pt.size=0.2) + NoLegend()

# cell type prediction scores
DefaultAssay(BM_query) <- 'prediction.score.celltype'
p5 <- FeaturePlot(BM_query, features = c("CD14 Mono", "CD4 Naive", "CD4 Memory", "Memory B"), reduction = 'ref.umap', ncol = 2, 
                  cols = c("lightgrey", "darkred"))
p5

p6 <- FeaturePlot(BM_query, features = c("MAIT", "NK", "Prog-RBC", "HSC"), reduction = 'ref.umap', ncol =2, 
                  cols = c("lightgrey", "darkred"))
p6

# Plot cell numbers and frequencies by sample
p7 <- dittoBarPlot(BM_query, var = "predicted.celltype", group.by = "orig.ident",  scale = c("percent"))
p7

p8 <- dittoBarPlot(BM_query, var = "predicted.celltype", group.by = "orig.ident",  scale = c("count"))
p8

# imputed protein levels from predicted_ADT slot
DefaultAssay(BM_query) <- 'predicted_ADT'

# Pick proteins
proteins <- c("CD34", "CD11a","CD11c","CD123","CD14","CD16","CD161", "CD56", "CD57", "CD3", "CD4", "CD8a", "CD127-IL7Ra", "CD197-CCR7", "CD25","CD27", "CD278-ICOS", "CD28", "CD69", "CD19","CD79b","HLA.DR","CD45RA","CD45RO", "CD38")
proteins1 <- c("CD34", "CD11a","CD11c","CD123","CD14","CD16","CD161", "CD56", "CD57")
proteins2 <- c("CD3", "CD4", "CD8a", "CD127-IL7Ra", "CD197-CCR7", "CD25","CD27", "CD278-ICOS", "CD28", "CD69")
proteins3 <- c("CD19","CD79b","HLA.DR","CD45RA","CD45RO", "CD38")

# Visualisation of specific proteins in UMAP with multidittoDimPlot
multi_dittoDimPlot(BM_query, vars = proteins1, reduction.use = "ref.umap", size = 0.25,
                   do.label = TRUE,
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

multi_dittoDimPlot(BM_query, vars = proteins2, reduction.use = "ref.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

multi_dittoDimPlot(BM_query, vars = proteins3, reduction.use = "ref.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

# Visualisation with Seurat
p6 <- FeaturePlot(BM_query, features = c("CD45RA","CD45RO", "CD14", "CD16"), reduction = 'proj.umap',
                  min.cutoff = 'q10', max.cutoff = 'q99', cols = c("lightgrey", "darkgreen") ,
                  ncol = 2)
p6

p7 <- FeaturePlot(BM_query, features = c("CD19", "CD3", "CD8a", "CD127-IL7Ra"), reduction = 'proj.umap',
                  min.cutoff = 'q10', max.cutoff = 'q99', cols = c("lightgrey", "darkgreen") ,
                  ncol = 2)
p7

p8 <- FeaturePlot(BM_query, features = c("CD56", "CD57", "CD161", "CD11c"), reduction = 'proj.umap',
                  min.cutoff = 'q10', max.cutoff = 'q99', cols = c("lightgrey", "darkgreen") ,
                  ncol = 2)
p8

# Compare this mapping with the HCA ref mapping, use the Azimuth ref reduction to generate the UMAPs

p1 <- DimPlot(BM_query, reduction = 'proj.umap', group.by = 'predicted.celltype', label = TRUE) + NoLegend()
p1

p2 <- dittoDimPlot(BM_query, var = "predicted.celltype",
                   reduction.use = "proj.umap", size = 0.75,
                   do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")
p2

# Computing a new UMAP visualization - use harmony here?

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
BM_query$id <- 'query'
refquery <- merge(bm, BM_query)
refquery[["spca"]] <- merge(bm[["spca"]], BM_query[["ref.spca"]])
refquery <- RunUMAP(refquery, reduction = 'spca', dims = 1:50)
DimPlot(refquery, group.by = 'id', shuffle = TRUE)
DimPlot(refquery, group.by = 'predicted.celltype', label = TRUE) + NoLegend()
DimPlot(refquery, group.by = 'predicted.celltype.l2', label = TRUE) + NoLegend()

table(BM_query$predicted.celltype)
table(BM_query$predicted.celltype.l2)

# 9- Assigning major cell types to cell identities

# Set idents from a value in object metadata - in this case the dsb clusters

colnames(BM_query[[]])
Idents(BM_query) <- 'predicted.celltype.l2'
levels(BM_query)

# Rename idents of dsb clusters 
major.cluster.ids <- c("B cell", "Hematopoietic Stem Cell", "DC", "CD4 T cell", "CD4 T cell", "Erythroid",
                       "CD8 T cell", "B cell",  "Myeloid", "HSPC", "CD8 T cell", "B cell",
                       "Erythroid", "HSPC", "B cell", "HSPC", "HSPC", "DC",
                       "DC", "HSPC", "CD4 T cell", "CD8 T cell", "Natural Killer", "DC",
                       "HSPC", "B cell", "DC", "DC", "Invariant T cell", "Myeloid",
                       "Myeloid", "Stromal cell", "B cell", "Innate Lymphocyte", "CD8 T cell", "CD8 T cell",
                       "Natural Killer", "Proliferating T cell", "Platelets", "Natural Killer")

names(major.cluster.ids) <- levels(BM_query)
BM_query <- RenameIdents(BM_query, major.cluster.ids)

# Adding object metadata with cluster names in order - to be used in the bar plots
BM_query$major_celltype.l2 <- factor(BM_query@active.ident)

# Visualise renamed/merged clusters
dittoDimPlot(BM_query, var = "major_celltype.l2",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Major cell types on UMAP")

# Frequency and total number of dsb clusters
dittoBarPlot(BM_query, var = "major_celltype.l2", group.by = "orig.ident", retain.factor.levels = TRUE)
dittoBarPlot(BM_query, var = "major_celltype.l2", group.by = "orig.ident", retain.factor.levels = TRUE, scale = "count")

################################################################################################################################################################################

# 10- Differential expression testing

# 10.1- Perform default differential expression tests - comparisons between clusters within the integrated/mapped Seurat object

# The bulk of Seurat’s differential expression features can be accessed through the FindMarkers() function. 
# As a default, Seurat performs differential expression based on the non-parametric Wilcoxon rank sum test. 
# This replaces the previous default test (‘bimod’). 
# To test for differential expression between two specific groups of cells, specify the ident.1 and ident.2 parameters. 

# which assay should be used for DGE comparisons? 
# Reference: Current best practices in single-cell RNA-seq analysis: a tutorial" described to use the raw pre-normalised measured values.

# list options for groups to perform differential expression on
DefaultAssay(BM_query) <- "SCT"
Idents(BM_query) <- 'predicted.celltype.l2'
levels(BM_query)

# 10.1.1- Prefilter features or cells to increase the speed of DE testing

# Reference: https://satijalab.org/seurat/articles/de_vignette.html
# To increase the speed of marker discovery, particularly for large datasets, Seurat allows for pre-filtering of features or cells. 
# For example, features that are very infrequently detected in either group of cells, or features that are expressed at similar average levels, are unlikely to be differentially expressed. 
# Example use cases of the min.pct, logfc.threshold, min.diff.pct, and max.cells.per.ident parameters are demonstrated below.

# Use min.pct to pre-filter features that are detected at <X% frequency;
# Use logfc.threshold to pre-filter features that have less than a X-fold change between the average expression;
# Use min.diff.pct to pre-filter features whose detection percentages across the two groups are similar;
# Use max.cells.per.ident to subsample each group to a maximum of 200 cells. Can be very useful for large clusters, or computationally-intensive DE tests;
# Increasing min.pct, logfc.threshold, and min.diff.pct, will increase the speed of DE testing, but could also miss features that are prefiltered.

# 10.1.2- Using the RNA normalised slot

BM_RNA_data_markers1 <- FindAllMarkers(BM_query, assay = "RNA", slot = "data", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
clipr::write_clip(BM_RNA_data_markers1) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste
saveRDS(BM_RNA_data_markers1, "BM_RNA_data_markers.rds")

# top 5 upregulated genes per cluster
BM1_RNA_data_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>%
  print(n=100)

# 10.1.3- Using the SCT normalised slot

BM_SCT_data_markers <- FindAllMarkers(BM_query, assay = "SCT", slot = "data", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
clipr::write_clip(BM_SCT_data_markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste
saveRDS(BM_SCT_data_markers, "BM_SCT_data_markers.rds")
BM_SCT_data_markers <- readRDS("./Plots/12_DGE_after_annotation/BM_SCT_data_markers.rds")

BM_SCT_data_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>% 
  print(n=100)

BM_SCT_data_markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

DoHeatmap(BM_query, features = top10$gene, size = 2) + NoLegend() + theme(axis.text.y = element_text(size = 3))
# Subset Seurat object top plot heatmap of specific clusters
DoHeatmap(subset(BM_query, downsample = 500, idents = c("Late Eryth", "Early Eryth")), features = top10$gene, size = 2, assay = "RNA", combine = TRUE) + theme(axis.text.y = element_text(size = 3))

# Heatmap visualization - DittoHeatmap
dittoHeatmap(BM_query, genes = top10$gene,
             assay = "SCT", order.by = c("predicted.celltype.l2"), cluster_rows = FALSE,
             cluster_cols = FALSE, scale = "none",
             heatmap.colors = inferno(100), 
             annot.by = c("predicted.celltype.l2"), fontsize_row = 2)


#  Find differentially expressed features between Late Eryth and all other cells, only search for positive markers
Late_Eryth_markers <- FindMarkers(BM_query, assay = "RNA", slot = "data", ident.1 = "Late Eryth", only.pos = TRUE, logfc.threshold = 0.1)
print((Late_Eryth_markers))
clipr::write_clip(Late_Eryth_markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# Find all markers distinguishing late vs early erythrocytes
Late_markers <- FindMarkers(BM_query, assay = "RNA", slot = "data", ident.1 = "Late Eryth", ident.2 = "Early Eryth", min.pct = 0.25, only.pos = FALSE, logfc.threshold = 0.1)
print(Late_markers)
clipr::write_clip(Late_markers)

# Find all markers distinguishing early erythrocytes vs EMP
Early_markers <- FindMarkers(BM_query, assay = "RNA", slot = "data", ident.1 = "Early Eryth", ident.2 = "EMP", min.pct = 0.25, only.pos = FALSE, logfc.threshold = 0.1)
print(Early_markers)
clipr::write_clip(Early_markers)

# Find all markers distinguishing CD8 Effector 1 erythrocytes vs CD8 Effector 2
CD8_effector_markers <- FindMarkers(BM_query, assay = "RNA", slot = "data", ident.1 = "CD8 Effector_2", ident.2 = "CD8 Effector_1", min.pct = 0.25, only.pos = FALSE, logfc.threshold = 0.1)
print(CD8_effector_markers)
clipr::write_clip(CD8_effector_markers)


DGE1 <- read.csv("./Plots/12_DGE_after_annotation/DGE_CD8Effector2_CD8Effector1.csv")
options(ggrepel.max.overlaps = 50000)
EnhancedVolcano(toptable = DGE1,
                x = "avg_log2FC",
                y = "p_val_adj",
                lab = DGE1$Features,
                pCutoff = 0.001,
                FCcutoff = 0.5,
                pointSize = 1,
                labSize = 3.0,
                labCol = 'black',
                boxedLabels = FALSE,
                title = "DGE CD8 Effector_2 vs CD8 Effector_1",
                legendLabels = c(
                  'Not significant',
                  'Fold change (but do not pass padj cutoff)',
                  'Pass padj cutoff',
                  'Pass both padj & fold change'),
                legendPosition = 'right',
                legendLabSize = 6.0,
                legendIconSize = 2.0,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                typeConnectors = "open",
                colConnectors = 'black',
) + theme_light()


# Highlight specific genes
celltype1 <- c('TFRC','GYPA', 'ITGA4', 'CD34')
celltype2 <- c('ITGA4')

EnhancedVolcano(toptable = DGE1,
                x = "avg_log2FC",
                y = "p_val_adj",
                selectLab = c(celltype1, celltype2),
                lab = DGE1$Features,
                pCutoff = 0.001,
                FCcutoff = 0.5,
                pointSize = 1,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                shape = 42,
                boxedLabels = FALSE,
                title = "DGE late eryth vs early eryth",
                legendLabels = c(
                  'Not significant',
                  'Fold change (but do not pass padj cutoff)',
                  'Pass padj cutoff',
                  'Pass both padj & fold change'),
                legendPosition = 'right',
                legendLabSize = 6.0,
                legendIconSize = 2.0,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                typeConnectors = "open",
                colConnectors = 'black',
                # encircle
                encircle = celltype1,
                encircleCol = 'skyblue',
                encircleSize = 2.5,
                encircleFill = 'pink',
                encircleAlpha = 1/2,
                # shade
                shade = celltype2,
                shadeAlpha = 1/2,
                shadeFill = 'skyblue',
                shadeSize = 1,
                shadeBins = 5
) + theme_bw()

# Exploring the expression of canonical marker genes
Idents(BM1_query) <- 'predicted.celltype.l2'
VlnPlot(BM1_query, features = c("TFRC", "ITGA4", "GYPA"), assay = "RNA", slot = "data", sort = TRUE, ncol=1) + NoLegend()
VlnPlot(BM1_query, features = c("TFRC", "ITGA4", "GYPA"), assay = "SCT", slot = "data", sort = TRUE, ncol=1) + NoLegend()

# 10.2- Subset to compare and plot specific clusters within the same condition

# Reference: https://satijalab.org/seurat/articles/de_vignette.html
# Reference: https://satijalab.org/seurat/articles/visualization_vignette.html

BM_query_eryth <- subset(BM_query, idents = c("Late Eryth", "Early Eryth"))
Late_markers <- FindAllMarkers(BM_query_eryth, assay = "RNA", slot = "data", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.1)

Late_markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC) -> top50

DoHeatmap(subset(BM_query_eryth, downsample = 500), assay = "RNA", features = top50$gene, size = 4) + NoLegend() + theme(axis.text.y = element_text(size = 6))

dittoHeatmap(subset(BM_query_eryth, downsample = 500), genes = top50$gene,
             assay = "RNA", slot = "scale.data",order.by = c("predicted.celltype.l2"), cluster_rows = FALSE,
             cluster_cols = FALSE, scale = "none",
             heatmap.colors = inferno(50), 
             annot.by = c("predicted.celltype.l2"), fontsize_row = 6)

# 10.3- DGE between same clusters different samples/conditions

# Reference: https://satijalab.org/seurat/articles/sctransform_v2_vignette.html#identify-differential-expressed-genes-across-conditions-1
# Using the normalized datasets with known celltype annotation, we can ask what genes change in different conditions for cells of the same type. 
# First, we create a column in the meta.data slot to hold both the cell type and stimulation information and switch the current ident to that column.
BM_query$celltype.sample <- paste(BM_query$predicted.celltype.l2, BM_query$sample_id,
                                  sep = "_")
Idents(BM_query) <- "celltype.sample"
levels(BM_query)

# Prior to performing differential expression, we first run PrepSCTFindMarkers, which ensures that the fixed value is set properly. 
# Then we use FindMarkers(assay="SCT") to find differentially expressed genes. 
BM_query <- PrepSCTFindMarkers(BM_query)

CD8_markers <- FindMarkers(BM_query, assay = "SCT", ident.1 = "CD8 Effector_2_HD_BM1.1", ident.2 = "CD8 Effector_2_HD_BM4.1", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.1,
                                     verbose = FALSE)
head(CD8_markers, n = 20)
clipr::write_clip(CD8_markers)
DGE1 <- read.csv("./Plots/12_DGE_after_annotation/DGE_CD8_BM1_BM4.csv")
options(ggrepel.max.overlaps = 50000)
EnhancedVolcano(toptable = DGE1,
                x = "avg_log2FC",
                y = "p_val_adj",
                lab = DGE1$Features,
                pCutoff = 0.001,
                FCcutoff = 0.5,
                pointSize = 1,
                labSize = 3.0,
                labCol = 'black',
                boxedLabels = FALSE,
                title = "DGE CD8 Effector_2 vs CD8 Effector_1",
                legendLabels = c(
                  'Not significant',
                  'Fold change (but do not pass padj cutoff)',
                  'Pass padj cutoff',
                  'Pass both padj & fold change'),
                legendPosition = 'right',
                legendLabSize = 6.0,
                legendIconSize = 2.0,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                typeConnectors = "open",
                colConnectors = 'black',
) + theme_light()


# If running on a subset of the original object after running PrepSCTFindMarkers(), 
# FindMarkers() should be invoked with recorrect_umi = FALSE to use the existing corrected counts:
BM_query_CD8 <- subset(BM_query, idents = c("CD8 Effector_2_HD_BM1.1", "CD8 Effector_2_HD_BM4.1"))
CD8_markers <- FindAllMarkers(BM_query_CD8, assay = "SCT", ident.1 = "CD8 Effector_2_HD_BM1.1",
                                            ident.2 = "CD8 Effector_2_HD_BM4.1", verbose = FALSE, recorrect_umi = FALSE)

# 10.4- Identify conserved cell type markers

# To identify canonical cell type marker genes that are conserved across conditions, we provide the FindConservedMarkers() function. 
# This function performs differential gene expression testing for each dataset/group and combines the p-values using meta-analysis methods from the MetaDE R package. 
# For example, we can identify genes that are conserved markers irrespective of stimulation condition in NK cells. 
# Note that the PrepSCTFindMarkers command does not to be rerun here.
Idents(BM_query) <- "predicted.celltype.l2"
LE.markers <- FindConservedMarkers(BM_query, assay = "SCT", ident.1 = "Late Eryth", grouping.var = "sample_id",
                                   verbose = FALSE)
head(LE.markers)

################################################################################################################################################################################

# 11- Visualisations in Seurat

# Reference: https://satijalab.org/seurat/articles/visualization_vignette.html
# Five visualizations of marker feature expression

# 11.1- # Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
features <- c("TFRC", "ITGA4", "GYPA", "CD34")
features2 <- c("GATA1", "GATA2", "SP1", "CEBPB", "PAX5", "EBF1")
features3 <- c("MS4A3", "CD247", "CCL5", "CD69")
features4 <- c("GZMB", "PF4", "SDC1", "CSF3R")
features5 <- c("CXCL12", "CXCR4","KITLG", "SPP1", "JAG1")
features6 <- c("PECAM1", "ICAM1", "ANGPT1", "SELE")
features7 <- c("ALAS1", "ALAS2", "NFE2", "TAL1", "ARID3A")

p1 <- RidgePlot(BM_query, features = features, assay = "SCT", slot = "data", ncol = 2, sort = TRUE, combine = TRUE, group.by = "predicted.celltype.l2") 
#+ theme(plot.title = element_text(size = rel(0.5)), axis.text = element_text(size = rel(0.5)))
p1

p1.5 <- RidgePlot(BM_query, features = c("TFRC", "ITGA4", "GYPA"), assay = "SCT", slot = "data", idents = c("Late Eryth", "Early Eryth"), ncol = 2, sort = TRUE, combine = TRUE) 
p1.5

p2 <- RidgePlot(BM_query, features = features2, assay = "SCT", slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
p2

p3 <- RidgePlot(BM_query, features = features3, assay = "SCT", slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
p3

p4 <- RidgePlot(BM_query, features = features4, assay = "SCT", slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
p4

p5 <- RidgePlot(BM_query, features = features5, assay = "SCT", slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
p5

p6 <- RidgePlot(BM_query, features = features6, assay = "SCT", slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
p6

# 11.2- Violin plot - Visualize single cell expression distributions in each cluster
p7 <- VlnPlot(BM_query, features = features, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE)
p7
p8 <- VlnPlot(BM_query, features = features2, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE)
p8
p9 <- VlnPlot(BM_query, features = features3, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE)
p9
p10 <- VlnPlot(BM_query, features = features4, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE)
p10
p11 <- VlnPlot(BM_query, features = features5, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE)
p11
p12 <- VlnPlot(BM_query, features = features6, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE)
p12
p13 <- VlnPlot(BM_query, features = features7, assay = "RNA", slot = "data", ncol = 2, sort = TRUE, combine = TRUE)
p13

p14 <- VlnPlot(BM_query, features = c("UMI"), sort = TRUE, combine = TRUE)
p14
p15 <- VlnPlot(BM_query, features = c("Genes"), sort = TRUE, combine = TRUE)
p15
p16 <- VlnPlot(BM_query, features = c("percent.mt"), sort = TRUE, combine = TRUE)
p16
p17 <- VlnPlot(BM_query, features = c("percent_ribo"), sort = TRUE, combine = TRUE)
p17
p18 <- VlnPlot(BM_query, features = c("percent_hb"), sort = TRUE, combine = TRUE)
p18

# 11.3- Feature plot - visualize feature expression in low-dimensional space
DefaultAssay(BM_query) <- 'RNA'
FeaturePlot(BM_query, features = features, ncol = 2)
FeaturePlot(BM_query, features = features2, ncol = 2)
FeaturePlot(BM_query, features = features3, ncol = 2)
FeaturePlot(BM_query, features = features4, ncol = 2)
FeaturePlot(BM_query, features = features5, ncol = 2)
FeaturePlot(BM_query, features = features6, ncol = 2)

# 11.4- Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
colnames(BM_query[[]])
Idents(BM_query) <- 'predicted.celltype.l2'
levels(BM_query)
levels(BM_query) <- c("Stromal", "HSC", "LMPP", "EMP", "Prog Mk", 
                      "Early Eryth", "Late Eryth", "BaEoMa",  "GMP", "pre-mDC", "cDC2", "ASDC", "pDC", "CD14 Mono", "CD16 Mono", 
                      "Macrophage", "CLP", "pro B", "pre B", "transitional B", "Naive B", "Memory B", "Plasma", "CD4 Naive", 
                      "CD4 Effector", "CD4 Memory", "CD8 Effector_1", "CD8 Effector_2", "CD8 Memory","NK", "MAIT")


BM_query$predicted.celltype.l2.ordered <- factor(BM_query@active.ident, 
                            levels=c("Stromal", 
                                     "HSC",
                                     "LMPP", 
                                     "EMP", 
                                     "Prog Mk", 
                                     "Early Eryth", 
                                     "Late Eryth", "BaEoMa", 
                                     "GMP", 
                                     "pre-mDC", "cDC2", "ASDC", "pDC", "CD14 Mono", "CD16 Mono", "Macrophage", 
                                     "CLP", "pro B", "pre B", "transitional B", "Naive B", "Memory B", "Plasma",
                                     "CD4 Naive", "CD4 Effector", "CD4 Memory", "CD8 Effector_1", "CD8 Effector_2", "CD8 Memory","NK", "MAIT"))

Idents(BM_query) <- 'predicted.celltype.l2.ordered'
levels(BM_query)

features7 <- c("CXCL12", "KITLG", "SPP1", "JAG1", "CD34", "GATA1", "GATA2", "TFRC", "ITGA4", "GYPA", "PF4", "SP1", "CEBPB", "CSF3R", "PAX5", "EBF1", "SDC1", "CD247", "CCL5", "GZMB", "MS4A3", "PECAM1", "ICAM1", "ANGPT1", "SELE")
DotPlot(BM_query, features = features7, assay = "RNA", cols = c("lightgrey", "blue")) + RotatedAxis()
DotPlot(BM_query, features = features) + RotatedAxis()
DotPlot(BM_query, features = features2) + RotatedAxis()
DotPlot(BM_query, features = features3) + RotatedAxis()
DotPlot(BM_query, features = features4) + RotatedAxis()

# 11.5- Single cell heatmap of feature expression
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

################################################################################################################################################################################

# 12- Visualisations in DittoSeq

# Reference: https://bioconductor.org/packages/devel/bioc/vignettes/dittoSeq/inst/doc/dittoSeq.html
# All plotting functions use these colors, stored in dittoColors(), by default.

# 12.1- Functions:

# DimPlot/ (I)FeaturePlot / UMAPPlot / etc.	dittoDimPlot / multi_dittoDimPlot
# VlnPlot / RidgePlot	dittoPlot / multi_dittoPlot
# DotPlot	dittoDotPlot
# FeatureScatter / GenePlot	dittoScatterPlot
# DoHeatmap	dittoHeatmap*
# [No Seurat Equivalent]	dittoBarPlot / dittoFreqPlot
# [No Seurat Equivalent]	dittoDimHex / dittoScatterHex
# [No Seurat Equivalent]	dittoPlotVarsAcrossGroups
# SpatialDimPlot, SpatialFeaturePlot, etc.	dittoSpatial (coming soon!)

# 12.2- Input:

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

# 12.3- Helper Functions

# dittoSeq’s helper functions make it easy to determine the metadata, gene, and dimensionality reduction options for plotting.
getMetas(BM_query)
# Query for the presence of a metadata slot
isMeta("nCount_RNA", BM_query)
# Retrieve metadata values:
meta("predicted.celltype.l2", BM_query)[1:10]
# Retrieve unique values of a metadata
metaLevels("predicted.celltype.l2", BM_query)
metaLevels("predicted.celltype", BM_query)
metaLevels("seurat_clusters", BM_query)
# Retrieve all gene names
DefaultAssay(BM_query) <- "SCT"
getGenes(BM_query)[1:10]
# Query for the presence of a gene(s)
isGene("CD3E", BM_query)
isGene(c("CD3E","ENO1","INS","non-gene"), BM_query, return.values = TRUE)
# Retrieve gene expression values:
gene("CD3E", BM_query)[1:10]
# Retrieve all dimensionality reductions
getReductions(BM_query)

# 12.4- Visualisations

#There are many different types of dittoSeq visualizations. 
# Each has intuitive defaults which allow creation of immediately usable plots. 
# Each also has many additional tweaks available through discrete inputs that can help ensure you can create precisely-tuned, deliberately-labeled, publication-quality plots out-of-the-box.

# 12.4.1- dittoDimPlot & dittoScatterPlot
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
dittoDimPlot(BM_query, "nCount_RNA", reduction.use = "proj.umap", size = 0.5,
             do.label = TRUE, assay = "RNA",
             slot = "data", do.contour = FALSE,
             contour.color = "lightblue", # Optional, black by default
             contour.linetype = "dashed", max = 50000) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("UMI counts")

dittoDimPlot(BM_query, "nFeature_RNA", reduction.use = "proj.umap", size = 0.5,
             do.label = TRUE, assay = "RNA",
             slot = "data", do.contour = FALSE,
             contour.color = "lightblue", # Optional, black by default
             contour.linetype = "dashed") + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Number of Genes")

dittoDimPlot(BM_query, "percent.mt", reduction.use = "proj.umap", size = 0.5,
             do.label = TRUE, assay = "RNA",
             slot = "data", do.contour = FALSE,
             contour.color = "lightblue", # Optional, black by default
             contour.linetype = "dashed") + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Percentage of mitochondrial genes")

dittoScatterPlot(
  object = BM_query,
  x.var = "nCount_RNA", y.var = "nFeature_RNA",
  color.var = "percent.mt")

# Visualisation of prediction scores in UMAP
dittoDimPlot(BM_query, "predicted.celltype.score")
dittoDimPlot(BM_query, "predicted.celltype.l2.score")

# Visualisation of cell types in UMAP
dittoDimPlot(BM_query, var = "predicted.celltype",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")

dittoDimPlot(BM_query, var = "predicted.celltype.l2",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities on UMAP")

# 12.4.2- dittoPlot (and dittoRidgePlot + dittoBoxPlot wrappers)

# These display continuous cells/samples’ data on a y-axis (or x-axis for ridgeplots) grouped on the x-axis by sample, age, condition, or any discrete grouping metadata. 
# Data can be represented with violin plots, box plots, individual points for each cell/sample, and/or ridge plots. 
# The plots input controls which data representations are used. 
# The group.by input controls how the data are grouped in the x-axis. And the color.by input controls the colors that fill in violin, box, and ridge plots.

# dittoPlot() is the main function, but dittoRidgePlot() and dittoBoxPlot() are wrappers which essentially just adjust the default for the plots input from c(“jitter”, “vlnplot”) to c(“ridgeplot”) or c(“boxplot”,“jitter”), respectively.

dittoPlot(BM_query, "CD34", assay = "RNA", slot = "data", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  vlnplot.lineweight = 0.5)
dittoPlot(BM_query, "nCount_RNA", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  vlnplot.lineweight = 0.5)
dittoPlot(BM_query, "nFeature_RNA", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  vlnplot.lineweight = 0.5)
dittoPlot(BM_query, "percent.mt", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  vlnplot.lineweight = 0.5)

dittoRidgePlot(BM_query, "CD34", assay = "RNA", slot = "data", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), ridgeplot.lineweight = 0.5, max = 2)
dittoRidgePlot(BM_query, "nCount_RNA", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  ridgeplot.lineweight = 0.5, max = 75000)
dittoRidgePlot(BM_query, "nFeature_RNA", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  ridgeplot.lineweight = 0.5)
dittoRidgePlot(BM_query, "percent.mt", group.by = "predicted.celltype.l2", 
          theme = theme_classic(), jitter.size = 0.5,  ridgeplot.lineweight = 0.5)

# 12.4.3- dittoBarPlot & dittoFreqPlot

# A couple of very handy visualizations missing from some other major single-cell visualization toolsets, 
# these functions quantify and display frequencies of clusters or cell types (or other discrete data) per sample (or other discrete groupings). 
# Such visualizations are quite useful for QC-ing clustering for batch effects and generally assessing cell type fluctuations.
# For both, data can be represented as percentages or counts, and this is controlled by the scale input.
# Visualisation of  frequency of cell types in bar plots

dittoBarPlot(BM_query, var = "predicted.celltype", group.by = "orig.ident")
dittoBarPlot(BM_query, var = "predicted.celltype", group.by = "orig.ident", scale = "count")

# dittoFreqPlot separates each cell type into its own facet, and thus puts more emphasis on individual cells. 
# An additional sample.by input controls splitting of cells within group.by-groups into individual samples.

# 12.4.4-  dittoHeatmap

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

dittoHeatmap(BM_query, genes, assay = "RNA", slot = "scale.data", annot.by = c("predicted.celltype.l2"), 
             heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50), scale = "column", heatmap.colors.max.scaled = inferno(10), cluster_rows = FALSE, scaled.to.max = FALSE)

dittoHeatmap(BM_query, genes, assay = "RNA", slot = "scale.data", annot.by = c("predicted.celltype.l2"), 
             heatmap.colors = inferno(3), scale = "column", heatmap.colors.max.scaled = inferno(10), cluster_rows = FALSE, scaled.to.max = FALSE)

BM_RNA_data_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

# Plot the top10 marker genes per cluster
dittoHeatmap(BM_query, top10$gene, assay = "RNA", slot = "scale.data", annot.by = c("predicted.celltype.l2"), 
             heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50), scale = "column", 
             heatmap.colors.max.scaled = inferno(10), cluster_rows = FALSE, scaled.to.max = FALSE, cluster_cols = FALSE,
             fontsize_row = 2)

# 12.4.5-  Multi-Plotters

# These create either multiple plots or create plots that summarize data for multiple variables all in one plot. 
# They make it easier to create summaries for many genes or many cell types without the need for writing loops.

# 12.4.5.1- dittoDotPlot

# A very succinct representation that is useful for showing differences between groups. 
# The plot uses differently colored and sized dots to summarizes both expression level (color) and percent of cells/samples with non-zero expression (size) for multiple genes (or values of metadata) within different groups of cells/samples.
# By default, expression values for all groups are centered and scaled to ensure a similar range of values for all vars displayed and to emphasize differences between groups.
dittoDotPlot(BM_query, vars = genes, 
             assay = "RNA", slot = "data",
             group.by = "predicted.celltype.l2", scale = FALSE) + theme_light()
# String which sets whether the values shown with color (default: mean non-zero expression) should be centered and scaled.
# scale = FALSE - plots average expression
# scale = TRUE (default) - plots relative expression


################################################################################################################################################################################

# 13- CITE-seq analysis

# 13.1- Visualisations of the imputed CITE-seq data from the mapped Seurat ref dataset

DefaultAssay(BM_query) <- "predicted_ADT"
rownames(BM_query@assays[["predicted_ADT"]])
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
multi_dittoDimPlot(BM_query, vars = proteins1, reduction.use = "proj.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

multi_dittoDimPlot(BM_query, vars = proteins2, reduction.use = "proj.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

multi_dittoDimPlot(BM_query, vars = proteins3, reduction.use = "proj.umap", size = 0.25,
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
dittoDotPlot(BM_query, vars = proteins1, 
             slot = "data",min.color = "grey90",
             max.color = "#C51B7D",
             group.by = "predicted.celltype.l2", scale = FALSE) + theme_classic(base_size = 12)

dittoDotPlot(BM_query, vars = proteins2, 
             slot = "data",min.color = "grey90",
             max.color = "#C51B7D",
             group.by = "predicted.celltype.l2", scale = FALSE) + theme_classic(base_size = 12)

dittoDotPlot(BM_query, vars = proteins3, 
             slot = "data",min.color = "grey90",
             max.color = "#C51B7D",
             group.by = "predicted.celltype.l2", scale = FALSE) + theme_classic(base_size = 12)

# Violin and RidgePlot with Seurat
proteins1 <- c("CD34", "CD11a","CD11c","CD123")
proteins2 <- c("CD14","CD16", "CD161", "CD56", "CD57")
proteins3 <- c("CD3", "CD4", "CD8a", "CD127-IL7Ra", "CD197-CCR7")
proteins4 <- c("CD25", "CD27", "CD278-ICOS", "CD28", "CD69")
proteins5 <- c("CD19","CD79b","HLA.DR","CD45RA","CD45RO", "CD38")

VlnPlot(BM_query, features = proteins1, slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()
VlnPlot(BM_query, features = proteins2, slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()
VlnPlot(BM_query, features = proteins3, slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()
VlnPlot(BM_query, features = proteins4, slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()
VlnPlot(BM_query, features = proteins5, slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()

RidgePlot(BM_query, features = proteins1, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
RidgePlot(BM_query, features = proteins2, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
RidgePlot(BM_query, features = proteins3, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
RidgePlot(BM_query, features = proteins4, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
RidgePlot(BM_query, features = proteins5, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 

# 14- Visualisation of cells based on dsb normalized protein using Seurat

DefaultAssay(BM_query) <- "dsb"

# Pick proteins
proteins <- c("CD34", "CD31", "CD3", "CD20", "CD11b", "CD44", "CD71", "CD49d")

# Visualisation of specific proteins in UMAP with multidittoDimPlot
multi_dittoDimPlot(BM_query, vars = proteins, reduction.use = "proj.umap", size = 0.25,
                   do.label = TRUE,
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

multi_dittoDimPlot(BM_query, vars = proteins, reduction.use = "ref.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

multi_dittoDimPlot(BM_query, vars = proteins, reduction.use = "ref.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

# dittoScatterPlot
# Visualisation of specific combination of proteins in Scatter plot, colored by cell types
dittoScatterPlot(
  object = BM_query,
  x.var = "CD3", y.var = "CD20",
  color.var = "predicted.celltype.l2")

# dittoDotPlot
dittoDotPlot(BM_query, vars = proteins, 
             slot = "data",min.color = "grey90",
             max.color = "#C51B7D",
             group.by = "predicted.celltype.l2.ordered", scale = FALSE) + theme_classic(base_size = 12)

# Violin and RidgePlot with Seurat
proteins1 <- c("CD34", "CD31", "CD3", "CD20")
proteins2 <- c("CD11b", "CD44", "CD71", "CD49d")

VlnPlot(BM_query, features = proteins1, slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()
VlnPlot(BM_query, features = proteins2, slot = "data", ncol = 2, pt.size = 0.25, sort = TRUE, combine = TRUE) + NoLegend()

RidgePlot(BM_query, features = proteins1, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
RidgePlot(BM_query, features = proteins2, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 


################################################################################################################################################################################

# 15- Clustering cells based on dsb normalized protein using Seurat

# Cluster cells based on dsb normalized protein levels. 
# Similar to workflow used in our paper Kotliarov et al. 2020 we don’t cluster based on principal components from ADT, instead directly using the normalized values.
DefaultAssay(BM_query) <- "dsb"
# define proteins to use in clustering (non-isptype controls)
prots = rownames(BM_query@assays[["dsb"]])[1:8]
proteins <- c("CD3","CD11b","CD20","CD31","CD34","CD44",
              "CD71", "CD49d")

# cluster and run umap 
BM_query = FindNeighbors(object = BM_query, dims = NULL, assay = 'dsb', 
                         features = prots, k.param = 30, 
                         verbose = FALSE)

# direct graph clustering 
BM_query = FindClusters(object = BM_query, resolution = 1, 
                        algorithm = 3,
                        graph.name = 'dsb_snn', 
                        verbose = FALSE)
# umap (optional)
BM_query = RunUMAP(object = BM_query, assay = "dsb", features = prots,
                   seed.use = 1990, min.dist = 0.2, n.neighbors = 30, reduction.key = "umap.dsb",
                   verbose = FALSE)

# make results dataframe 
d = cbind(BM_query@meta.data, 
          as.data.frame(t(BM_query@assays$dsb@data)),
          BM_query@reductions$umap@cell.embeddings)

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

dittoHeatmap(BM_query, var = "ident", assay = "dsb", slot = "data", annot.by = c("dsb_snn_res.1"), 
             heatmap.colors = colorRampPalette(c("blue", "white", "red"))(50), scale = "column", 
             heatmap.colors.max.scaled = inferno(10), cluster_rows = FALSE, scaled.to.max = FALSE, cluster_cols = FALSE,
             fontsize_row = 10)

dittoDotPlot(BM_query, vars = proteins, 
             assay = "dsb", slot = "data",min.color = "grey90",
             max.color = "#C51B7D",
             group.by = "dsb_snn_res.1", scale = FALSE) + theme_classic(base_size = 12)

# Visualisation of clusters in UMAP using different projections
dittoDimPlot(BM_query, var = "dsb_snn_res.1",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from dsb on UMAP")

dittoDimPlot(BM_query, var = "dsb_snn_res.1",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from dsb on UMAP")

# Frequency and total number of dsb clusters
dittoBarPlot(BM_query, var = "dsb_snn_res.1", group.by = "orig.ident")
dittoBarPlot(BM_query, var = "dsb_snn_res.1", group.by = "orig.ident", scale = "count")

# 14.2- Assigning cell type identity to dsb clusters

# Set idents from a value in object metadata - in this case the dsb clusters

colnames(BM_query[[]])
Idents(BM_query) <- 'dsb_snn_res.1'
levels(BM_query)

# Rename idents of dsb clusters 
dsb.cluster.ids <- c("CD34+ cells", "T cells", "T cells" ,"CD34+ cells", "B cells", "CD31+CD11b+ cells", 
                     "CD71IntCD49dIntCD44Int cells", "CD71LowCD49dLowCD44Low cells", "CD31+CD11b- cells", 
                     "CD31+CD11b+ cells", "T cells", "CD34+CD31+ cells", "CD71HighCD49d+CD44+ cells", "T cells",
                     "B cells", "CD31+CD11b+ cells", "T cells")

names(dsb.cluster.ids) <- levels(BM_query@meta.data[["dsb_snn_res.1"]])
BM_query <- RenameIdents(BM_query, dsb.cluster.ids)

# Adding object metadata with cluster names in order - to be used in the bar plots
BM_query$dsb_cluster_order <- factor(BM_query@active.ident,
                                     levels = c("CD71HighCD49d+CD44+ cells", "CD71IntCD49dIntCD44Int cells", "CD71LowCD49d+CD44+ cells", "B cells", "T cells",
                                                "CD31+CD11b+ cells", "CD31+CD11b- cells", "CD34+CD31+ cells", "CD34+ cells"))

# Visualise renamed/merged clusters
dittoDimPlot(BM_query, var = "dsb_cluster_order",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 6, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from dsb on UMAP")

dittoDimPlot(BM_query, var = "ident",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 6, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from dsb on UMAP")

# Plot a heatmap of the average dsb normalized values for each cluster
e = cbind(BM_query@active.ident, 
          as.data.frame(t(BM_query@assays$dsb@data)))

adt_plot = e %>% 
  dplyr::group_by(BM_query@active.ident) %>% 
  dplyr::summarize_at(.vars = prots, .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("BM_query@active.ident") 

pheatmap::pheatmap((adt_plot), 
                   color = viridis::viridis(25, option = "B"), 
                   fontsize_row = 14, border_color = NA, fontsize_col = 14, angle_col = "0")


# Frequency and total number of dsb clusters
dittoBarPlot(BM_query, var = "dsb_cluster_order", group.by = "orig.ident", retain.factor.levels = TRUE)
dittoBarPlot(BM_query, var = "dsb_cluster_order", group.by = "orig.ident", retain.factor.levels = TRUE, scale = "count")

# 15- Clustering cells based on predicted ADT expression

# Cluster cells based on predicted ADT protein levels. 
DefaultAssay(BM_query) <- "predicted_ADT"

# define proteins to use in clustering (non-isptype controls)
prots = rownames(BM_query@assays[["predicted_ADT"]])

proteins <- c("CD34", "CD11a","CD11c","CD123","CD14","CD16","CD161", "CD56", "CD57", "CD3", "CD4", "CD8a", "CD127-IL7Ra", "CD197-CCR7", "CD25","CD27", "CD278-ICOS", "CD28", "CD69", "CD19","CD79b","HLA.DR","CD45RA","CD45RO", "CD38")
proteins1 <- c("CD34", "CD11a","CD11c","CD123","CD14","CD16","CD161", "CD56", "CD57")
proteins2 <- c("CD3", "CD4", "CD8a", "CD127-IL7Ra", "CD197-CCR7", "CD25","CD27", "CD278-ICOS", "CD28", "CD69")
proteins3 <- c("CD19","CD79b","HLA.DR","CD45RA","CD45RO", "CD38")

# cluster and run umap 
BM_query = FindNeighbors(object = BM_query, dims = NULL, assay = 'predicted_ADT', 
                         features = prots, k.param = 30, 
                         verbose = FALSE)

# direct graph clustering 
BM_query = FindClusters(object = BM_query, resolution = 1, 
                        algorithm = 3,
                        graph.name = 'predicted_ADT_snn', 
                        verbose = FALSE)
# umap (optional)
BM_query = RunUMAP(object = BM_query, assay = "predicted_ADT", features = prots,
                   seed.use = 1990, min.dist = 0.2, n.neighbors = 30, reduction.key = "predADT",
                   verbose = FALSE)

colnames(BM_query[[]])
Idents(BM_query) <- 'predicted.celltype.l2.score'
levels(BM_query)

# make results dataframe 
d = cbind(BM_query@meta.data, 
          as.data.frame(t(BM_query@assays$predicted_ADT@data)))
#BM1_query@reductions$umap@cell.embeddings)

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
dittoDimPlot(BM_query, var = "predicted_ADT_snn_res.1",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from predicted_ADT on UMAP")

dittoDimPlot(BM_query, var = "predicted_ADT_snn_res.1",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from predicted_ADT on UMAP")

dittoBarPlot(BM_query, var = "predicted_ADT_snn_res.1", group.by = "orig.ident")
dittoBarPlot(BM_query, var = "predicted_ADT_snn_res.1", group.by = "orig.ident", scale = "count")

dittoDotPlot(BM_query, vars = proteins, 
             assay = "predicted_ADT", slot = "data",min.color = "grey90",
             max.color = "#C51B7D",
             group.by = "predicted_ADT_snn_res.1", scale = FALSE) + theme_classic(base_size = 12)

# 15.1- Cluster interpretation

colnames(BM_query[[]])
Idents(BM_query) <- 'predicted_ADT_snn_res.1'
levels(BM_query)

# Rename idents of dsb clusters 
predicted.ADT.cluster.ids <- c("CD4 Naive T cell", "CD4 Naive T cell","CD45RAHighCD38- B cell", "CD4 Naive T cell", "CD8 Naive T cell",
      "CD14High Monocyte", "CD27Low CD4 T cell", "CD27Low CD4 T cell", "Other", "CD45RAIntCD38Int B cell", "Other", "CD8 Naive T cell", "CD11aHighCD123Low cell",
      "CD4 Memory T cell", "CD34+ cell", "Other", "CD45RALowCD38Int B cell", "CD4 Memory T cell", "CD8 Effector T cell", "CD38High Plasma cell", 
      "Other", "CD45RAIntCD38Int B cell", "CD45RAIntCD38Int B cell", "CD45RAIntCD38Int B cell", "Other", "CD14Low Monocyte", "CD45RALowCD38Int B cell",
      "CD123HighCD4+ DC", "CD34+ cell", "Other", "CD11aHighCD11cHigh cell", "CD34+ cell", "Other", "CD45RAHighCD16HighCD11aHigh cell")

names(predicted.ADT.cluster.ids) <- levels(BM_query@meta.data[["predicted_ADT_snn_res.1"]])
BM_query <- RenameIdents(BM_query, predicted.ADT.cluster.ids)

# Adding object metadata with cluster names in order - to be used in the bar plots
BM_query$predicted.ADT.cluster.ids <- factor(BM_query@active.ident)

# Visualise renamed/merged clusters
dittoDimPlot(BM_query, var = "predicted.ADT.cluster.ids",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 6, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from predicted ADT on UMAP")

dittoDimPlot(BM_query, var = "predicted.ADT.cluster.ids",
             reduction.use = "umap", size = 0.75,
             do.label = TRUE, labels.size = 6, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("Cell identities from predicted ADT on UMAP")

# Plot a heatmap of the average dsb normalized values for each cluster
e = cbind(BM_query@active.ident, 
          as.data.frame(t(BM_query@assays$predicted_ADT@data)))

adt_plot = e %>% 
  dplyr::group_by(BM_query@active.ident) %>% 
  dplyr::summarize_at(.vars = prots, .funs = median) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("BM_query@active.ident") 

pheatmap::pheatmap((adt_plot), 
                   color = viridis::viridis(25, option = "B"), 
                   fontsize_row = 12, border_color = NA, fontsize_col = 12, angle_col = "45")


# Frequency and total number of dsb clusters
dittoBarPlot(BM_query, var = "predicted.ADT.cluster.ids", group.by = "orig.ident", retain.factor.levels = TRUE)
dittoBarPlot(BM_query, var = "predicted.ADT.cluster.ids", group.by = "orig.ident", retain.factor.levels = TRUE, scale = "count")

################################################################################################################################################################################

# 16- Weighted Nearest Neighbor multimodal clustering using dsb normalized values - Seurat vignette

# The dsb normalized values can be used with the Weighted Nearest Neighbor multimodal clustering method. WNN is an excellent way to find fine-grained cell states
# Reference: https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-cite-seq-rna-adt-1

DefaultAssay(BM_query) = "dsb"

prots = rownames(BM_query@assays[["dsb"]])[1:8]
VariableFeatures(BM_query) = prots

BM_query = BM_query %>% 
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
BM_query = FindMultiModalNeighbors(
  BM_query, reduction.list = list("pca", "apca"), 
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
BM_query <- RunUMAP(BM_query, nn.name = "weighted.nn", 
                     reduction.name = "wnn.umap", 
                     reduction.key = "wnnUMAP_")

BM_query <- FindClusters(BM_query, graph.name = "wsnn", 
                          algorithm = 3, 
                          resolution = 1.5, 
                          verbose = FALSE, 
                          random.seed = 1990)

# Visualise WNN clusters in the HCA reference dataset UMAP
dittoDimPlot(BM_query, var = "wsnn_res.1.5",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

# Visualise WNN clusters in the WNN (ADT+RNA) UMAP projection
dittoDimPlot(BM_query, var = "wsnn_res.1.5",
             reduction.use = "wnn.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())


# Visualise predicted clusters from the HCA reference dataset in the WNN (ADT+RNA) UMAP projection
dittoDimPlot(BM_query, var = "predicted.celltype.l2",
             reduction.use = "wnn.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

p1 <- DimPlot(BM1_query, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p2 <- DimPlot(BM1_query, reduction = 'wnn.umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p1 + p2

# We can also compute UMAP visualization based on only the RNA and protein data and compare. 
# We find that the RNA analysis is more informative than the ADT analysis in identifying progenitor states (the ADT panel contains markers for differentiated cells),
# while the converse is true of T cell states (where the ADT analysis outperforms RNA).
BM_query <- RunUMAP(BM_query, reduction = 'pca', dims = 1:30, assay = 'RNA', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
BM_query <- RunUMAP(BM_query, reduction = 'apca', dims = 1:7, assay = 'dsb', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

dittoDimPlot(BM_query, var = "predicted.celltype.l2",
             reduction.use = "rna.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

dittoDimPlot(BM_query, var = "predicted.celltype.l2",
             reduction.use = "adt.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

p3 <- DimPlot(BM_query, reduction = 'rna.umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p4 <- DimPlot(BM_query, reduction = 'adt.umap', group.by = 'predicted.celltype.l2', label = TRUE, repel = TRUE, label.size = 2.5) + NoLegend()
p3 + p4

# Visualisation of proteins in WNN UMAP
DefaultAssay(BM_query) <- "dsb"
# Pick proteins
proteins <- c("CD3","CD11b","CD20","CD31","CD34","CD44",
              "CD71", "CD49d")

multi_dittoDimPlot(BM_query, vars = proteins, reduction.use = "wnn.umap", size = 0.25,
                   do.label = TRUE, assay = "dsb",
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN UMAP")

# WNN RNA UMAP
multi_dittoDimPlot(BM_query, vars = proteins, reduction.use = "rna.umap", size = 0.25,
                   do.label = TRUE, assay = "dsb",
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN RNA UMAP")

# WNN ADT UMAP
multi_dittoDimPlot(BM_query, vars = proteins, reduction.use = "adt.umap", size = 0.25,
                   do.label = TRUE, assay = "dsb",
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN ADT UMAP")

# We can visualize the expression of canonical marker genes and proteins on the multimodal UMAP, which can assist in verifying the provided annotations:

# Finally, we can visualize the modality weights that were learned for each cell. 
# Each of the populations with the highest RNA weights represent progenitor cells, while the populations with the highest protein weights represent T cells. 
# This is in line with our biological expectations, as the antibody panel does not contain markers that can distinguish between different progenitor populations.
ditto_VlnPlot(BM_query, features = "SCT.weight", group.by = 'predicted.celltype.l2', sort = TRUE, pt.size = 0.1) + NoLegend()
VlnPlot(BM_query, features = "SCT.weight", group.by = 'predicted.celltype.l2', sort = TRUE, pt.size = 0.1) + NoLegend()

################################################################################################################################################################################

# 17- Weighted Nearest Neighbor multimodal clustering using predicted ADT values
DefaultAssay(BM_query) <- 'predicted_ADT'

prots = rownames(BM_query@assays[["predicted_ADT"]])
VariableFeatures(BM_query) = prots

BM_query = BM_query %>% 
  ScaleData() %>% 
  RunPCA(reduction.name = 'apca.predicted.adt', verbose = FALSE)

# Cell-specific modality weights can be accessed at BM1_query$RNA.weight
BM_query = FindMultiModalNeighbors(
  BM_query, reduction.list = list("pca", "apca.predicted.adt"), 
  dims.list = list(1:30, 1:24), knn.graph.name = "predicted.ADT.wknn",
  snn.graph.name = "predicted.ADT.wsnn",
  weighted.nn.name = "predicted.ADT.weighted.nn",
  modality.weight.name = "predicted.ADT.weight", 
  verbose = FALSE
)

BM_query <- RunUMAP(BM_query, nn.name = "predicted.ADT.weighted.nn", 
                     reduction.name = "predicted.ADT.wnn.umap", 
                     reduction.key = "predicted.ADT.wnnUMAP_")

BM_query <- FindClusters(BM_query, graph.name = "predicted.ADT.wsnn", 
                          algorithm = 3, 
                          resolution = 1.5, 
                          verbose = FALSE, 
                          random.seed = 1990)

# Visualise WNN clusters in the HCA reference dataset UMAP
dittoDimPlot(BM_query, var = "predicted.ADT.wsnn_res.1.5",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

# Visualise WNN clusters in the WNN (ADT+RNA) UMAP projection
dittoDimPlot(BM_query, var = "predicted.ADT.wsnn_res.1.5",
             reduction.use = "predicted.ADT.wnn.umap", size = 0.75,
             do.label = TRUE, labels.size = 4, labels.highlight = TRUE) +
  theme(legend.title = element_blank())

# Visualise predicted clusters from the HCA reference dataset in the WNN (ADT+RNA) UMAP projection
dittoDimPlot(BM_query, var = "predicted.celltype.l2",
             reduction.use = "predicted.ADT.wnn.umap", size = 0.75,
             do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
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
multi_dittoDimPlot(BM_query, vars = proteins1, reduction.use = "predicted.ADT.wnn.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

multi_dittoDimPlot(BM_query, vars = proteins2, reduction.use = "predicted.ADT.wnn.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

multi_dittoDimPlot(BM_query, vars = proteins3, reduction.use = "predicted.ADT.wnn.umap", size = 0.25,
                   do.label = TRUE, 
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Predicted protein expression")

# WNN RNA UMAP
multi_dittoDimPlot(BM_query, vars = proteins1, reduction.use = "rna.umap", size = 0.25,
                   do.label = TRUE,
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN RNA UMAP")

multi_dittoDimPlot(BM_query, vars = proteins2, reduction.use = "rna.umap", size = 0.25,
                   do.label = TRUE,
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN RNA UMAP")

multi_dittoDimPlot(BM_query, vars = proteins3, reduction.use = "rna.umap", size = 0.25,
                   do.label = TRUE,
                   slot = "data", do.contour = FALSE) + # Optional, solid by default
  theme(legend.title = element_blank()) +
  ggtitle("Protein expression on WNN RNA UMAP")

# WNN ADT UMAP
multi_dittoDimPlot(BM_query, vars = proteins1, reduction.use = "adt.umap", size = 0.25,
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

################################################################################################################################################################################

# 18- Automated cell annotation using SingleR

# Reference: https://github.com/dviraran/SingleR
# Reference: https://github.com/LTLA/SingleR/blob/master/README.md
# Reference: Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage.
# Reference: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html
# Reference: http://bioconductor.org/books/devel/SingleRBook/
# Reference: https://nbisweden.github.io/excelerate-scRNAseq/session-celltypeid/celltypeid.html

# Example of the use of SingleR() with the Human Primary Cell Atlas dataset (Mabbott et al. 2013) as the reference.
# Loading reference data with Ensembl annotations.
# Each celldex dataset actually has three sets of labels that primarily differ in their resolution. 
# For the purposes of this demonstration, we will use the “fine” labels in the label.fine metadata field, which represents the highest resolution of annotation available for this dataset.
library(celldex)
ref.data <- HumanPrimaryCellAtlasData(ensembl=FALSE)
saveRDS(ref.data, "HumanPrimaryCellAtlasData.rds")

hemato <- NovershternHematopoieticData(ensembl=FALSE)
saveRDS(hemato, "./Reference_datasets/NovershternHematopoieticData/NovershternHematopoieticData.rds")
hemato <- readRDS("./Reference_datasets/NovershternHematopoieticData/NovershternHematopoieticData.rds")

blueprint <- BlueprintEncodeData(ensembl=FALSE)
saveRDS(hemato, "BlueprintEncodeData.rds")

# Human dtabases in the library celldex:
# HumanPrimaryCellAtlasData - ownload and cache the normalized expression values of the data stored in the Human Primary Cell Atlas.
# BlueprintEncodeData - Download and cache the normalized expression values of 259 RNA-seq samples of pure stroma and immune cells as generated and supplied by Blueprint and ENCODE.
# DatabaseImmuneCellExpressionData - Download and cache the normalized expression values of 1561 bulk RNA-seq samples of sorted cell populations from the Database of Immune Cell Expression (DICE).
# NovershternHematopoieticData - ownload and cache the normalized expression values of 211 bulk human microarray samples of sorted hematopoietic cell populations that can be found in GSE24759.
# Samples were additionally annotated to 38 fine cell types ("label.fine").
# MonacoImmuneData - Download and cache the normalized expression values of 114 bulk RNA-seq samples of sorted immune cell populations


# 20.1- Performing predictions. Convert Seurat object to SingleCellExperiment object
library(SingleCellExperiment)
library(scater)
library(SingleR)

DefaultAssay(BM_query) <- "RNA" 
Idents(BM_query) <- 'predicted.celltype.l2'
levels(BM_query)
BM_sce <- as.SingleCellExperiment(BM_query, assay = "RNA")
saveRDS(BM_sce, "BM_query_sce.rds")

# use the sce created for the analysis with milo
#BM_sce <- readRDS("./RDS_files/Merged_datasets/BM_merged_sce.rds")

# We perform annotation by calling SingleR() on our test dataset and the reference (ImmGen) dataset, 
# leaving the default of de.method="classic" to use the original marker detection scheme. 
# This applies the algorithm described in Section 1.2, returning a DataFrame where each row contains prediction results for a single cell in the sce object. 
# Labels are provided in the labels column.
# See 'Choices of assay data' for 'assay.type.test=' explanation.
pred <- SingleR(test = BM_sce, ref = hemato, 
                labels = hemato$label.fine, assay.type.test=1)
colnames(pred)
saveRDS(pred, "./RDS_files/Lonza_dataset/SingleR_pred.rds")
pred <- readRDS("./RDS_files/Lonza_dataset/SingleR_pred.rds")

# 20.2- Quality control
# Upon examining the distribution of assigned labels, we see that many of them are related to cells. 
head(sort(table(pred$labels), decreasing=TRUE))

# SingleR results labels can be easily added back to the metadata of these objects as well:
BM_query[["SingleR.labels"]] <- pred$labels

# Visualise Single R labels in the UMAP from one of the BM reference
dittoDimPlot(BM_query, var = "SingleR.labels",
             reduction.use = "proj.umap", size = 0.75,
             do.label = TRUE, labels.size = 2, labels.highlight = TRUE) +
  theme(legend.title = element_blank()) +
  ggtitle("SingleR cell types on UMAP")

# Frequency and total number of clusters
dittoBarPlot(BM_query, var = "SingleR.labels", group.by = "sample_id", retain.factor.levels = TRUE)
dittoBarPlot(BM_query, var = "SingleR.labels", group.by = "sample_id", retain.factor.levels = TRUE, scale = "count")

# Or if `method="cluster"` was used:
# seurat.obj[["SingleR.cluster.labels"]] <- singler.results$labels[match(seurat.obj[[]][["my.input.clusters"]], rownames(singler.results))]

# Plot normalised protein expression per cluster
DefaultAssay(BM_query) <- "dsb"
Idents(BM_query) <- "SingleR.labels"

# Violin and RidgePlot with Seurat
proteins1 <- c("CD3", "CD11b","CD20","CD31")
proteins2 <- c("CD34","CD44", "CD71", "CD49d")

VlnPlot(BM_query, features = proteins1, layer = "data", ncol = 2, pt.size = 0.1, sort = TRUE, combine = TRUE) + NoLegend()
VlnPlot(BM_query, features = proteins2, slot = "data", ncol = 2, pt.size = 0.1, sort = TRUE, combine = TRUE) + NoLegend()

RidgePlot(BM_query, features = proteins1, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 
RidgePlot(BM_query, features = proteins2, slot = "data", ncol = 2, sort = TRUE, combine = TRUE) 

DefaultAssay(BM_query) <- "RNA" 
Idents(BM_query) <- "SingleR.labels"
SingleR.markers <- FindAllMarkers(BM_query, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.25) #min.pct = 0.25,
clipr::write_clip(SingleR.markers) # Copies the table to the clipboard. Go to excel and press cmd+V (mac) or ctrl+V (windows) to paste

# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
SingleR.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

# changing the default color
DoHeatmap((subset(BM_query, downsample = 1000)), features = top10$gene, size = 3) + NoLegend() + theme(axis.text.y = element_text(size = 8))

dittoHeatmap(BM_query, top10$gene, 
             cluster_rows = FALSE, 
             annot.by = "SingleR.labels", 
             fontsize_row = 2)

avgexp <- AverageExpression(subset(BM_query, features = top10$gene))
avgexp <- as.data.frame(avgexp$RNA)
#write.csv(avgexp, "Average_expression_Top_expressed_genes_CIPR_immune_clusters_RNA.csv")
#avgexp$gene <- rownames(avgexp)
#avgexp <- avgexp[,1:13]

#avgexp2 <- scale(avgexp)
avgexp3 <- t(avgexp)
avgexp4 <- scale(avgexp3)
avgexp5 <- t(avgexp4)

library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
heat <- Heatmap(avgexp5, cluster_columns = FALSE, cluster_rows = FALSE,
                col = colorRamp2(c(1, 0.6, 0, -0.6, -1), brewer.pal(n = 5, name = "RdYlBu"), space = "RGB"),
                heatmap_legend_param = list(color_bar = "continuous"),
                show_row_dend = TRUE, show_column_dend = TRUE,
                row_names_gp = gpar(fontsize = 2),
                column_names_gp = gpar(fontsize = 6))

print(heat)

# Plot TFRC, GYPA, gene counts, UMI counts and nuclear transcripts in SingleR and CIPR erythroid clusters

DefaultAssay(BM_query) <- "RNA" 
Idents(BM_query) <- "SingleR.labels"
levels(BM_query)

# median number of RNA molecules per cluster
tapply(BM_query$nCount_RNA, BM_query$SingleR.labels,  median)

BM_erythro <- subset(BM_query, idents = c("Megakaryocyte/erythroid progenitors","Erythroid_CD34+ CD71+ GlyA-", "Erythroid_CD34- CD71+ GlyA-", "Erythroid_CD34- CD71+ GlyA+",
                                          "Erythroid_CD34- CD71lo GlyA+", "Erythroid_CD34- CD71- GlyA+"))

Idents(BM_erythro) <- "SingleR.labels"
levels(BM_erythro)
my_levels <- c("Megakaryocyte/erythroid progenitors", "Erythroid_CD34+ CD71+ GlyA-", "Erythroid_CD34- CD71+ GlyA-", "Erythroid_CD34- CD71+ GlyA+", "Erythroid_CD34- CD71lo GlyA+", "Erythroid_CD34- CD71- GlyA+")
# Relevel object@ident
BM_erythro@active.ident <- factor(x = BM_erythro@active.ident, levels = my_levels)
BM_erythro$SingleR.labels_ordered <- factor(BM_erythro@active.ident)

# Visualize number of RNA molecules per cluster as a violin plot 
feats <- c("Genes", "nCount_RNA", "percent_hb", "TFRC", "GYPA", "ITGA4")
feats2 <- c("POLR1A", "POLR2A", "POLR3A")
VlnPlot(BM_erythro, features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
VlnPlot(BM_erythro, features = feats2, pt.size = 0.1, ncol = 3) + NoLegend()

DefaultAssay(BM_erythro) <- "dsb"
Idents(BM_erythro) <- "SingleR.labels_ordered"

# Violin and RidgePlot with Seurat
proteins1 <- c("CD3", "CD11b","CD20","CD31")
proteins2 <- c("CD34","CD44", "CD71", "CD49d")

VlnPlot(BM_erythro, features = proteins2, slot = "data", ncol = 2, pt.size = 0.1) + NoLegend()

################################################################################################################################################################################


# 21- Convert Seurat object to anndata for scvelo and CellRank analysis in Python

library(sceasy)
DefaultAssay(BM_query) <- "RNA" 
sceasy::convertFormat(BM_query, from="seurat", to="anndata",
                      outFile='BM_query.h5ad')

################################################################################################################################################################################


# 22- Subset Seurat object to cluster major cell types and re-add the info to the original seurat object

# reference: https://github.com/satijalab/seurat/issues/1748

DefaultAssay(BM_query) <- "RNA" 
Idents(BM_query) <- "major_celltype.l2"
levels(BM_query)

BM_eryth <- subset(BM_query, idents = c("Hematopoietic Stem Cell", "HSPC", "Erythroid"))
BM_myelo <- subset(BM_query, idents = c("Hematopoietic Stem Cell", "HSPC", "Myeloid", "DC"))
BM_bcell <- subset(BM_query, idents = c("Hematopoietic Stem Cell", "HSPC", "B cell"))
BM_lympho <- subset(BM_query, idents = c("CD4 T cell", "CD8 T cell", "Invariant T cell", "Natural Killer"))
BM_stem <- subset(BM_query, idents = c("Hematopoietic Stem Cell", "HSPC"))

saveRDS(BM_eryth, "./RDS_files/Lonza_dataset/BM_eryth.rds")
saveRDS(BM_myelo, "./RDS_files/Lonza_dataset/BM_myelo.rds")
saveRDS(BM_bcell, "./RDS_files/Lonza_dataset/BM_bcell.rds")
saveRDS(BM_lympho, "./RDS_files/Lonza_dataset/BM_lympho.rds")
saveRDS(BM_stem, "./RDS_files/Lonza_dataset/BM_stem.rds")

sceasy::convertFormat(BM_eryth, from="seurat", to="anndata",
                      outFile='BM_eryth.h5ad')
sceasy::convertFormat(BM_myelo, from="seurat", to="anndata",
                      outFile='BM_myelo.h5ad')
sceasy::convertFormat(BM_bcell, from="seurat", to="anndata",
                      outFile='BM_bcell.h5ad')
sceasy::convertFormat(BM_lympho, from="seurat", to="anndata",
                      outFile='BM_lympho.h5ad')
sceasy::convertFormat(BM_stem, from="seurat", to="anndata",
                      outFile='BM_stem.h5ad')

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





