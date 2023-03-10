# Script: "sc-ATACSeq Analysis"
# Author: Keshav Prasad Gubbi

# SECTION1: Initialization
suppressPackageStartupMessages(library("ArchR"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("patchwork"))
set.seed(1234)

# Set the number of threads
addArchRThreads(threads = 8)

here::i_am("scATACseqAnalysis.R")

# Input data

day1 <- "/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day12Sample/2_Processed_data/V16_d12_230105/fragments.tsv.gz"
day5 <- "/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day5Sample/2_Processed_data/AJ_v12/fragments.tsv.gz"
day10 <- "/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day10Sample/compressed_tars/2_Processed_data/AJ_V12_d10_L129_UI/fragments.tsv.gz"
day25 <- "/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day25Sample/fragments.tsv.gz"
day4dpi <- "/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/4dpiSample/fragments.tsv.gz"

# Samples Input Vector
samples <- c(day1, day5, day10, day25, day4dpi)
# Set understandable names
names(samples) <- c("d01", "d05", "d10", "d25", "Inf-d4")
# inputFiles <- c('scATAC' = samples)
inputFiles <- c(samples)
print(names(inputFiles))

# Add respective Genome
addArchRGenome("mm10")

## Creating Arrow Files: #Each Arrow file stores all of the data associated with an individual sample (i.e. metadata, accessible fragments, and data matrices).
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  # filterFrags = 2500,
  minTSS = 5, # 4 was the default value
  # Dont set this too high because you can always increase later
  minFrags = 2500,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  subThreading = TRUE,
  verbose = TRUE,
  cleanTmp = TRUE,
  force = TRUE, # will recreate arrow files anew each time
  logFile = createLogFile("All5")
)

##############################################################################################
#                             ## SECTION 2: Quality Control
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            ## SECTION2.1 : Doublet Score                 ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## Inferring scATAC-seq Doublets with ArchR : A doublet refers to a single droplet that received a single barcoded bead and more than one nucleus. This causes the reads from more than one cell to appear as a single cell that is effectively the average of the two cells. These are removed computationally!
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 30, # Refers to how many cells near a "pseudo-doublet" to count
  knnMethod = "LSI", # Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force = TRUE,
  verbose = TRUE
)

##############################################################################################
##                     SECTION3 : Setting up ArchR project
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          # 3.1 ## Creating ArchR Project                   ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

ATACSeq_project_All5 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "/home/keshavprasad/Documents/scATACseq/All5/",
  copyArrows = FALSE # This is recommended so that you maintain an unaltered copy for later usage.
)
paste0("Memory Size = ", round(object.size(ATACSeq_project_All5) / 10^6, digits = 3), " MB")
getAvailableMatrices(ATACSeq_project_All5)

# Inspect the newly created Project.
print(ATACSeq_project_All5)
# Save the newly created
saveArchRProject(
  ArchRProj = ATACSeq_project_All5,
  outputDirectory = getOutputDirectory(ATACSeq_project_All5),
  overwrite = TRUE,
  load = TRUE,
  logFile = createLogFile("saveArchRProject_All5"),
)
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          # 3.2 Manipulating An ArchRProject                ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Cell names
# head(ATACSeq_project_All5$cellNames)
# # Sample names
# head(ATACSeq_project_All5$Sample)
# # TSS Enrichment
# head(ATACSeq_project_All5$TSSEnrichment)
# # One can access the TSS Enrichment Scores for each cell:
# quantile(ATACSeq_project_All5$TSSEnrichment)
# # doublet Enrichment
# head(ATACSeq_project_All5$DoubletEnrichment)
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          # 3.3 Obtaining columns from cellColData         ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
DoubletEnrichment_df <- getCellColData(ATACSeq_project_All5, select = "DoubletEnrichment")
head(DoubletEnrichment_df, 30)

## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## #3.3 Plotting QC metrics - log10(Unique Fragments) vs TSS enrichment score      ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
QCMetric_df <- getCellColData(ATACSeq_project_All5, select = c("log10(nFrags)", "TSSEnrichment"))
(scatterplot_UniqueFragsVsEnrichment <- ggPoint(
  x = QCMetric_df[, 1],
  y = QCMetric_df[, 2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(QCMetric_df[, 1], probs = 0.99)),
  ylim = c(0, quantile(QCMetric_df[, 2], probs = 0.99))
) +
  geom_hline(yintercept = 5, lty = "dashed") +
  geom_vline(xintercept = 3.4, lty = "dashed"))

plotPDF(scatterplot_UniqueFragsVsEnrichment,
  name = "TSS-vs-Frags_PreFiltering.pdf",
  ArchRProj = ATACSeq_project_All5,
  addDOC = TRUE
)
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## #3.4 Plotting QC metrics - PreFiltering                                     ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Make a violin plot for each sample based on a column from cellcolldata using plotGroups().
create_GroupPlot <- function(projectname, cellColData_colname){
  Group_Plot <- plotGroups(
    ArchRProj = projectname,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = cellColData_colname,
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
  )
  return(Group_Plot)
}

### violin plot for each sample for TSSEnrichment.
(Group_plot_TSS_violin_PreFiltering <- create_GroupPlot(ATACSeq_project_All5, "TSSEnrichment") + ggtitle("TSSEnrichment"))
### Plots (per sample) for log10 (unique nuclear fragments) #### log10 (unique nuclear fragments)
(Group_plot_nFrags_violin_PreFiltering <- create_GroupPlot(ATACSeq_project_All5, "log10(nFrags)") + ggtitle("nFrags"))
### Plots (per sample) for BlacklistRatio
(Group_plot_BLR_violin_PreFiltering <- create_GroupPlot(ATACSeq_project_All5, "BlacklistRatio") + ggtitle("BlacklistRatio"))
### Plots (per sample) for NucleosomeRatio
(Group_plot_NR_violin_PreFiltering <- create_GroupPlot(ATACSeq_project_All5, "NucleosomeRatio") + ggtitle("NucleosomeRatio"))
### Plots (per sample) for DoubletScore
(Group_plot_DS_violin_PreFiltering <- create_GroupPlot(ATACSeq_project_All5, "DoubletScore") + ggtitle("DoubletScore"))
### Plots (per sample) for DoubletEnrichment
(Group_plot_DE_violin_PreFiltering <- create_GroupPlot(ATACSeq_project_All5, "DoubletEnrichment") + ggtitle("DoubletEnrichment"))

plotPDF(Group_plot_TSS_violin_PreFiltering, Group_plot_nFrags_violin_PreFiltering, Group_plot_BLR_violin_PreFiltering, Group_plot_NR_violin_PreFiltering, Group_plot_DS_violin_PreFiltering, Group_plot_DE_violin_PreFiltering, 
        name = "QC-Metrics_PreFiltering.pdf", ArchRProj = ATACSeq_project_All5, addDOC = TRUE, width = 12, height = 12)

## Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
(FragSizePlot <- plotFragmentSizes(ArchRProj = ATACSeq_project_All5))
(TSSEnrichmentPlot <- plotTSSEnrichment(ArchRProj = ATACSeq_project_All5))
plotPDF(FragSizePlot, TSSEnrichmentPlot, name = "QC-Sample-FragSizes-TSSProfile_PreFiltering.pdf", ArchRProj = ATACSeq_project_All5, addDOC = TRUE, width = 10, height = 10)

# Filter out low Quality Cells from ArchR Project, based on the violin plots of QC Metrics
ATACSeq_project_All5 <- ATACSeq_project_All5[ATACSeq_project_All5$NucleosomeRatio < 2.5 &
                                             ATACSeq_project_All5$TSSEnrichment > 6 &
                                             ATACSeq_project_All5$BlacklistRatio < 0.05 &
                                             ATACSeq_project_All5$DoubletScore < 50 &
                                             ATACSeq_project_All5$DoubletEnrichment < 3.5
                                            ]

# QC Plots PostFiltering
### violin plot for each sample for TSSEnrichment.
(Group_plot_TSS_violin_PostFiltering <- create_GroupPlot(ATACSeq_project_All5, "TSSEnrichment") + ggtitle("TSSEnrichment"))
### Plots (per sample) for log10 (unique nuclear fragments) #### log10 (unique nuclear fragments)
(Group_plot_nFrags_violin_PostFiltering <- create_GroupPlot(ATACSeq_project_All5, "log10(nFrags)") + ggtitle("nFrags"))
### Plots (per sample) for BlacklistRatio
(Group_plot_BLR_violin_PostFiltering <- create_GroupPlot(ATACSeq_project_All5, "BlacklistRatio") + ggtitle("BlacklistRatio"))
### Plots (per sample) for NucleosomeRatio
(Group_plot_NR_violin_PostFiltering <- create_GroupPlot(ATACSeq_project_All5, "NucleosomeRatio") + ggtitle("NucleosomeRatio"))
### Plots (per sample) for DoubletScore
(Group_plot_DS_violin_PostFiltering <- create_GroupPlot(ATACSeq_project_All5, "DoubletScore") + ggtitle("DoubletScore"))
### Plots (per sample) for DoubletEnrichment
(Group_plot_DE_violin_PostFiltering <- create_GroupPlot(ATACSeq_project_All5, "DoubletEnrichment") + ggtitle("DoubletEnrichment"))

plotPDF(Group_plot_TSS_violin_PostFiltering, Group_plot_nFrags_violin_PostFiltering, Group_plot_BLR_violin_PostFiltering, Group_plot_NR_violin_PostFiltering, Group_plot_DS_violin_PostFiltering, Group_plot_DE_violin_PostFiltering, name = "QC-Metrics_PostFiltering.pdf", ArchRProj = ATACSeq_project_All5, addDOC = TRUE, width = 12, height = 12)

## Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
(FragSizePlot_PostFiltering <- plotFragmentSizes(ArchRProj = ATACSeq_project_All5))
(TSSEnrichmentPlot_PostFiltering <- plotTSSEnrichment(ArchRProj = ATACSeq_project_All5))
plotPDF(FragSizePlot_PostFiltering, TSSEnrichmentPlot_PostFiltering, name = "QC-Sample-FragSizes-TSSProfile_PostFiltering.pdf", ArchRProj = ATACSeq_project_All5, addDOC = TRUE, width = 10, height = 10)
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## # 3.6 Filtering Doublets from an ArchRProject                       ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# We filter putative doublets based on the previously determined doublet scores using the filterDoublets() function. This doesn’t physically remove data from the Arrow files but rather tells the ArchRProject to ignore these cells for downstream analysis.
# The higher the filterRatio, the greater the number of cells potentially removed as doublets.

## Finding doublets
ATACSeq_project_All5 <- filterDoublets(ArchRProj = ATACSeq_project_All5)

# Save the newly created
saveArchRProject(
  ArchRProj = ATACSeq_project_All5,
  outputDirectory = getOutputDirectory(ATACSeq_project_All5),
  overwrite = TRUE,
  load = TRUE,
  logFile = createLogFile("saveArchRProject_All5_PostQC"),
)

# Load the saved ArchR project
# readRDS(file.path("/media/keshavprasad/HornefLab_Data3/scRNA_AnnotationData_Johannes",
# "scrna_with_day25.Rds"))

# loadArchRProject(path = "./", force = TRUE, showLogo = TRUE)
##############################################################################################
##                     SECTION4 : Dimensional Reduction, INTEGRATION, CLUSTERING & EMBEDDING
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

# It is important to note that LSI is not deterministic.
# This means that even if you run LSI in exactly the same way with exactly the same parameters,
# you will not get exactly the same results. Of course, they will be highly similar,
# but not identical. So make sure to save your ArchRProject or the relevant LSI information
# once you’ve settled on an ideal dimensionality reduction. For the purposes of this project,
# we will create a reducedDims object called “IterativeLSI_all5”.

# In the context of scATAC-seq and scRNA-seq data ArchR performs LSI following these steps:
#   
#   scATAC-seq: documents=samples, words=regions/peaks. scRNA-seq: documents=samples, words=genes.
# Calculate word frequency by depth normalization per single cell.
# Normalize word frequency by the inverse document frequency which weights features by how often they occur
# Results in a word frequency-inverse document frequency (TF-IDF) matrix, which reflects how important a word (aka region/peak) is to a document (aka sample).
# Perform singular value decomposition (SVD) on the TF-IDF matrix. 

ATACSeq_project_All5 <- addIterativeLSI(# 4.1 ## Reducing Dims via Iterative LSI.# The most common parameters to tweak are iterations, varFeatures, and resolution.
    ArchRProj = ATACSeq_project_All5,
    useMatrix = "TileMatrix",
    name = "IterativeLSI_all5",
    # iterations = 2,
    saveIterations = FALSE,
    clusterParams = list(
      resolution = c(0.2),
      sampleCells = 10000,
      n.start = 10
    ),
    varFeatures = 15000,
    dimsToUse = 1:30,
    force = TRUE
) 

ATACSeq_project_All5 <- addHarmony(# 4.4 ## Batch Effect Correction with HARMONY                   ##
    ArchRProj = ATACSeq_project_All5,
    reducedDims = "IterativeLSI_all5",
    name = "Harmony_all5",
    groupBy = "Sample",
    verbose = TRUE,
    force = TRUE
) 

ATACSeq_project_All5 <- addUMAP(# 5.1 ## Visualising clusters via UMAP
    ArchRProj = ATACSeq_project_All5,
    reducedDims = "Harmony_all5", # Need to use the Harmony object for clustering and UMAP for proper integration
    name = "UMAP_all5",
    nNeighbors = 30,
    minDist = 0.4,
    metric = "cosine",
    force = TRUE,
    verbose = TRUE,
    saveModel = TRUE
)

# Create a Split view UMAP so that the UMAP for each sample can be visualized separately.
create_splitUMAP <- function(ATACSeq_project_All5, sampleObject) {
  split_UMAP <- plotEmbedding(
    ArchRProj = ATACSeq_project_All5,
    colorBy = "cellColData",
    name = "Sample",
    embedding = "UMAP_all5",
    highlightCells = getCellNames(ArchRProj = ATACSeq_project_All5)[which(ATACSeq_project_All5@cellColData$Sample == sampleObject)],
    keepAxis = TRUE,
    baseSize = 10,
    plotAs = "points"
)
  return(split_UMAP)
}

# Get the the Samples from the ArchR project
samplelist <- ATACSeq_project_All5@cellColData$Sample
# Call the function to create UMAP for each sample
(d1_UMAP = create_splitUMAP(ATACSeq_project_All5, samplelist@values[1]))
(d5_UMAP = create_splitUMAP(ATACSeq_project_All5, samplelist@values[5]))
(d10_UMAP = create_splitUMAP(ATACSeq_project_All5, samplelist@values[4]))
(d25_UMAP = create_splitUMAP(ATACSeq_project_All5, samplelist@values[3]))
(d4pi_UMAP = create_splitUMAP(ATACSeq_project_All5, samplelist@values[2]))

# (ggAlignPlots(p1_UMAP, p2_UMAP, type = "h"))
plotPDF(d1_UMAP, d5_UMAP, d10_UMAP, d25_UMAP, d4pi_UMAP,
  name = "Plot-Splitview UMAP-Sample-Clusters.pdf",
  ArchRProj = ATACSeq_project_All5,
  addDOC = TRUE,
  width = 10,
  height = 10
)

# Create UMAP for all samples together
p1_UMAP <- plotEmbedding(
  ArchRProj = ATACSeq_project_All5,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_all5",
  keepAxis = TRUE
)
plotPDF(p1_UMAP, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = ATACSeq_project_All5, addDOC = TRUE, width = 10, height = 10)

# Adding Clusters 
ATACSeq_project_All5 <- addClusters(# 5.1 ## Creating clusters based on Iterative LSI object      ##
    input = ATACSeq_project_All5,
    reducedDims = "Harmony_all5", # Harmony_all5
    method = "Seurat",
    name = "Clusters_all5",
    resolution = 0.5,
    force = TRUE
)
ATACSeq_project_All5 <- addUMAP(# 5.1 ## Visualising clusters via UMAP
    ArchRProj = ATACSeq_project_All5,
    reducedDims = "Harmony_all5", # Need to used the HArmony object for cluistering and UMAP for proper integration
    name = "ClustersUMAP_all5",
    nNeighbors = 30,
    minDist = 0.4,
    metric = "cosine",
    force = TRUE,
    verbose = TRUE,
    saveModel = TRUE
)

(p2_UMAP <- plotEmbedding(ArchRProj = ATACSeq_project_All5, colorBy = "cellColData", name = "Clusters_all5", embedding = "ClustersUMAP_all5", keepAxis = TRUE))
plotPDF(p2_UMAP, name = "Plot-UMAP-Clusters.pdf", ArchRProj = ATACSeq_project_All5, addDOC = TRUE, width = 10, height = 10)

# 
# # Accessing the clusters
head(ATACSeq_project_All5$Clusters_all5)
# tabulate the number of cells present in each cluster:
table(ATACSeq_project_All5$Clusters_all5)

# To better understand which samples reside in which clusters, we can create a
# cluster confusion matrix across each sample using the confusionMatrix() function.
cM <- confusionMatrix(
  paste0(ATACSeq_project_All5$Clusters_all5),
  paste0(ATACSeq_project_All5$Sample)
)
print(cM)

# plotting the Confusion Matrix
cM <- cM / Matrix::rowSums(cM)
(p <- pheatmap::pheatmap(
  mat = as.matrix(cM),
  color = paletteContinuous("whiteBlue"),
  border_color = "black"
))


## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##            # 5.3 ## Dimensionality Reduction after Harmony                         ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

# can assess the effects of Harmony by visualizing the embedding using UMAP or t-SNE and
# comparing this to the embeddings visualized in the previous sections for iterative LSI.

# ATACSeq_project_All5 <- addUMAP(ArchRProj = ATACSeq_project_All5,
#                                 reducedDims = "Harmony_all5",
#                                 name = "UMAPHarmony_all5",
#                                 nNeighbors = 30,
#                                 minDist = 0.5,
#                                 metric = "cosine"
# )
#
# p3_UMAP <- plotEmbedding(ArchRProj = ATACSeq_project_All5,
#                          colorBy = "cellColData",
#                          name = "Sample",
#                          embedding = "UMAPHarmony_all5")
# p4_UMAP <- plotEmbedding(ArchRProj = ATACSeq_project_All5,
#                          colorBy = "cellColData",
#                          name = "Clusters",
#                          embedding = "UMAPHarmony_all5")
# (ggAlignPlots(p3_UMAP, p4_UMAP, type = "h"))
#
# plotPDF(p1_UMAP,p2_UMAP, p3_UMAP,p4_UMAP,
#         name = "Plot-UMAP2Harmony-Sample-Clusters.pdf",
#         ArchRProj = ATACSeq_project_All5,
#         addDOC = TRUE,
#         width = 10,
#         height = 10)

##############################################################################################
##                     SECTION7 : Gene Scores and Marker Genes
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##            # 6.3 ## Identifying Marker Genes                               ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#
markersGS <- getMarkerFeatures(
  ArchRProj = ATACSeq_project_All5,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters_all5",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

## Markers List
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1.0")
head(markerList@listData, 19)
# cluster1_df <- as.data.frame(markerList$C1)
# cluster2_df <- as.data.frame(markerList$C2)
# cluster3_df <- as.data.frame(markerList$C3)
# cluster4_df <- as.data.frame(markerList$C4)
# cluster5_df <- as.data.frame(markerList$C5)
# cluster6_df <- as.data.frame(markerList$C6)
# cluster7_df <- as.data.frame(markerList$C7)
# cluster8_df <- as.data.frame(markerList$C8)
# cluster9_df <- as.data.frame(markerList$C9)
# cluster10_df <- as.data.frame(markerList$C10)

# SECTION3. Marker Genes
##  visualize all of the marker features simultaneously
markerGenesList <- c("Plag2g2a", "Defa-rs1", "Mmp7", "Lyz1", "Spdef", "Tcf7l2", "Ephb3", "Sis", "Ada", "Lct")

(heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.0",
  limits = c(-3, 3),
  returnMatrix = FALSE,
  plotLog2FC = TRUE,
  labelMarkers = markerGenesList,
  transpose = TRUE,
  labelRows = TRUE, clusterCols = TRUE, nPrint = 10
))

plotPDF(heatmapGS,
  name = "GeneScores-Marker-Heatmap", width = 10, height = 10,
  ArchRProj = ATACSeq_project_All5, addDOC = TRUE
)



## –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##          SECTION 7: Annotating Cell types with a Reference Dataset       ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# ArchR includes a function to align a reference scRNA-seq dataset, and impute cell type annotations based on the reference annotation (addGeneIntegrationMatrix). 
# As a reference, we will use a pre-processed scRNA-seq dataset for human PBMCs.

##           7.1 Unconstrained Integration       ##
# Read-in the reference
annotate_reference <- readRDS(file.path(
  "/media/keshavprasad/HornefLab_Data3/scRNA_AnnotationData_Johannes",
  "scrna_with_day25.Rds"
))

# add gene integration matrix
ATACSeq_project_All5 <- addGeneIntegrationMatrix(
  ArchRProj = ATACSeq_project_All5,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI_all5",
  seRNA = annotate_reference,
  addToArrow = FALSE,
  groupRNA = "int_0.3_broad_tuft",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

# Plot UMAP with predicted cell types
(annotated_UMAP_scRNAseqRef <- plotEmbedding(ATACSeq_project_All5,
  name = "predictedGroup_Un",
  embedding = "UMAP_all5",
  size = 1.5,
  labelAsFactors = F,
  labelMeans = F
))

plotPDF(annotated_UMAP_scRNAseqRef,
  name = "annotated_UMAP_scRNAseqRef_unconstrained",
  width = 10, height = 10, ArchRProj = ATACSeq_project_All5, addDOC = TRUE
)

##           7.2 Constrained Integration       ##
## Now that we have our preliminary unconstrained integration, we will identify general cell types to profide a framework to further refine the integration results. 
# First, we will identify which cell types from the scRNA-seq data are most abundant in each of our scATAC-seq clusters and , we would ideally constrain the integration to associate similar cell types togther. 

cM_cells <- as.matrix(confusionMatrix(ATACSeq_project_All5$Clusters_all5, ATACSeq_project_All5$predictedGroup_Un))
preClust <- colnames(cM_cells)[apply(cM_cells, 1 , which.max)]
cbind(preClust, rownames(cM_cells)) #Assignments

# First, lets look at the cell type labels from our scRNA-seq data that were used in our unconstrained integration:
unique(unique(ATACSeq_project_All5$predictedGroup_Un))
#From scRNA
cTNK <- paste0(paste0(1:12), collapse="|")
cTNK

#Assign scATAC to these categories
clustTNK <- rownames(cM)[grep(cTNK, preClust)]
clustTNK
##                            Calling Peaks                         ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
