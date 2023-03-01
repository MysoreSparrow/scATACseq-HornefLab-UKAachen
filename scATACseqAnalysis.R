#Script: "sc-ATACSeq Analysis"
#Author: Keshav Prasad Gubbi

# SECTION1: Initialization
suppressPackageStartupMessages(library("ArchR"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("patchwork"))
set.seed(1234)

# Set the number of threads
addArchRThreads(threads = 8)

here::i_am("scATACseqAnalysis.R")

# Input data 

day1 <- "/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day12Sample/2_Processed_data/V16_d12_230105/fragments.tsv.gz"
day5 <- '/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day5Sample/2_Processed_data/AJ_v12/fragments.tsv.gz'
day10 <- '/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day10Sample/compressed_tars/2_Processed_data/AJ_V12_d10_L129_UI/fragments.tsv.gz'
day25 <- '/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day25Sample/fragments.tsv.gz'
day4dpi <- '/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/4dpiSample/fragments.tsv.gz'

inputFiles <- c('scATAC' = c(day1, day5, day10, day25, day4dpi))
print(names(inputFiles))

# Add respective Genome
addArchRGenome("mm10")

## Creating Arrow Files: #Each Arrow file stores all of the data associated with an individual sample (i.e. metadata, accessible fragments, and data matrices).
ArrowFiles <- createArrowFiles(inputFiles = inputFiles,
                               sampleNames = names(inputFiles),
                               # filterFrags = 2500, 
                               minTSS = 5, #4 was the default value
                               #Dont set this too high because you can always increase later
                               minFrags = 2500,
                               addTileMat = TRUE,
                               addGeneScoreMat = TRUE,
                               subThreading = TRUE,
                               verbose = TRUE,
                               cleanTmp = TRUE,
                               force = FALSE, # will recreate arrow files anew each time
                               logFile = createLogFile("All5")
                               )

##############################################################################################
#                             ## SECTION 2: Quality Control
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            ## SECTION2.1 : Doublet Score                 ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##Inferring scATAC-seq Doublets with ArchR : A doublet refers to a single droplet that received a single barcoded bead and more than one nucleus. This causes the reads from more than one cell to appear as a single cell that is effectively the average of the two cells. These are removed computationally!
doubScores <- addDoubletScores(input = ArrowFiles,
                               k = 30, #Refers to how many cells near a "pseudo-doublet" to count
                               knnMethod = "LSI", #Refers to the embedding to use for nearest neighbor search.
                               LSIMethod = 1,
                               force = FALSE,
                               verbose = TRUE
                               )

##############################################################################################
##                     SECTION3 : Setting up ArchR project
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          # 3.1 ## Creating ArchR Project                   ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

ATACSeq_project_All5 <- ArchRProject(ArrowFiles = ArrowFiles,
                                    outputDirectory = "/home/keshavprasad/Documents/scATACseq/All5/",
                                    copyArrows = FALSE #This is recommended so that you maintain an unaltered copy for later usage.
)
paste0("Memory Size = ", round(object.size(ATACSeq_project_All5) / 10^6, digits = 3), " MB")
getAvailableMatrices(ATACSeq_project_All5)

# Inspect the newly created Project.
print(ATACSeq_project_All5)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          # 3.2 Manipulating An ArchRProject                ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Cell names
head(ATACSeq_project_All5$cellNames)
# Sample names
head(ATACSeq_project_All5$Sample)
# TSS Enrichment
head(ATACSeq_project_All5$TSSEnrichment)
# One can access the TSS Enrichment Scores for each cell:
quantile(ATACSeq_project_All5$TSSEnrichment)
# doublet Enrichment 
head(ATACSeq_project_All5$DoubletEnrichment)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          # 3.3 Obtaining columns from cellColData         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
DoubletEnrichment_df <- getCellColData(ATACSeq_project_All5, select = "DoubletEnrichment")
head(DoubletEnrichment_df, 30)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## #3.3 Plotting QC metrics - log10(Unique Fragments) vs TSS enrichment score      ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
QCMetric_df <- getCellColData(ATACSeq_project_All5, select = c("log10(nFrags)", "TSSEnrichment"))
(scatterplot_UniqueFragsVsEnrichment <- ggPoint(x = QCMetric_df[,1], 
                                                y = QCMetric_df[,2],
                                                colorDensity = TRUE,
                                                continuousSet = "sambaNight",
                                                xlabel = "Log10 Unique Fragments",
                                                ylabel = "TSS Enrichment",
                                                xlim = c(log10(500), quantile(QCMetric_df[,1], probs = 0.99)),
                                                ylim = c(0, quantile(QCMetric_df[,2], probs = 0.99))) + 
    geom_hline(yintercept = 5, lty = "dashed") +
    geom_vline(xintercept = 3, lty = "dashed"))

plotPDF(scatterplot_UniqueFragsVsEnrichment, 
        name = "TSS-vs-Frags.pdf", 
        ArchRProj = ATACSeq_project_All5, 
        addDOC = TRUE)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## #3.4 Plotting Sample Statistics from an ArchRProject                       ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Make a violin plot for each sample for the TSS enrichment scores.

### violin plot for each sample for TSSEnrichment.
(Group_plot_TSS_violin <- plotGroups(ArchRProj = ATACSeq_project_All5,
                                     groupBy = "Sample",
                                     colorBy = "cellColData",
                                     name = "TSSEnrichment",
                                     plotAs = "violin",
                                     alpha = 0.4,
                                     addBoxPlot = TRUE) + ggtitle("TSSEnrichment"))

### Plots (per sample) for log10 (unique nuclear fragments)
#### log10 (unique nuclear fragments)
(Group_plot_nFrags_violin <- plotGroups(ArchRProj = ATACSeq_project_All5,
                                        groupBy = "Sample",
                                        colorBy = "cellColData",
                                        name = "log10(nFrags)",
                                        plotAs = "violin",
                                        alpha = 0.4,
                                        addBoxPlot = TRUE) + ggtitle("nFrags"))
### violin plot for each sample for TSSEnrichment.
(Group_plot_TSS_violin <- plotGroups(ArchRProj = ATACSeq_project_All5,
                                     groupBy = "Sample",
                                     colorBy = "cellColData",
                                     name = "TSSEnrichment",
                                     plotAs = "violin",
                                     alpha = 0.4,
                                     addBoxPlot = TRUE) + ggtitle("TSSEnrichment"))

(Group_plot_BLR_violin <- plotGroups(ArchRProj = ATACSeq_project_All5,
                                     groupBy = "Sample",
                                     colorBy = "cellColData",
                                     name = "BlacklistRatio",
                                     plotAs = "violin",
                                     alpha = 0.4,
                                     addBoxPlot = TRUE) + ggtitle("BlacklistRatio"))

(Group_plot_NR_violin <- plotGroups(ArchRProj = ATACSeq_project_All5,
                                    groupBy = "Sample",
                                    colorBy = "cellColData",
                                    name = "NucleosomeRatio",
                                    plotAs = "violin",
                                    alpha = 0.4,
                                    addBoxPlot = TRUE) + ggtitle("NucleosomeRatio"))

(Group_plot_DS_violin <- plotGroups(ArchRProj = ATACSeq_project_All5,
                                    groupBy = "Sample",
                                    colorBy = "cellColData",
                                    name = "DoubletScore",
                                    plotAs = "violin",
                                    alpha = 0.4,
                                    addBoxPlot = TRUE) + ggtitle("DoubletScore"))

plotPDF(Group_plot_nFrags_violin, 
        Group_plot_TSS_violin, 
        Group_plot_BLR_violin ,
        Group_plot_NR_violin,
        Group_plot_DS_violin, 
        name = "QC-Sample-Statistics.pdf", 
        ArchRProj = ATACSeq_project_All5, 
        addDOC = TRUE, width = 12, height = 12)

## Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
(FragSizePlot = plotFragmentSizes(ArchRProj = ATACSeq_project_All5))
(TSSEnrichmentPlot = plotTSSEnrichment(ArchRProj = ATACSeq_project_All5))
plotPDF(FragSizePlot, 
        TSSEnrichmentPlot, 
        name = "QC-Sample-FragSizes-TSSProfile.pdf", 
        ArchRProj = ATACSeq_project_All5, 
        addDOC = TRUE, width = 10, height = 10)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## # 3.6 Filtering Doublets from an ArchRProject                       ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# We filter putative doublets based on the previously determined doublet scores using the filterDoublets() function. This doesn’t physically remove data from the Arrow files but rather tells the ArchRProject to ignore these cells for downstream analysis. 
# The higher the filterRatio, the greater the number of cells potentially removed as doublets.

## Finding doublets
ATACSeq_project_All5 <- filterDoublets(ArchRProj = ATACSeq_project_All5)

##############################################################################################
##                     SECTION4 : Dimensional Reduction
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          # 4.1 ## Reducing Dims via Iterative LSI                   ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# The most common parameters to tweak are iterations, varFeatures, and resolution. 

ATACSeq_project_All5 <- addIterativeLSI(ArchRProj = ATACSeq_project_All5,
                                        useMatrix = "TileMatrix", 
                                        name = "IterativeLSI_all5", 
                                        iterations = 6, 
                                        clusterParams = list( #See Seurat::FindClusters
                                          resolution = c(0.8), 
                                          sampleCells = 10000, 
                                          n.start = 10), 
                                        varFeatures = 25000, 
                                        dimsToUse = 1:30,
                                        force = TRUE
                                      ) 

# It is important to note that LSI is not deterministic. 
# This means that even if you run LSI in exactly the same way with exactly the same parameters, 
# you will not get exactly the same results. Of course, they will be highly similar, 
# but not identical. So make sure to save your ArchRProject or the relevant LSI information 
# once you’ve settled on an ideal dimensionality reduction. For the purposes of this project, 
# we will create a reducedDims object called “IterativeLSI_all5”.
  
#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##           # 4.4 ## Batch Effect Correction with HARMONY                   ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

ATACSeq_project_All5 <- addHarmony(ArchRProj = ATACSeq_project_All5, 
                        reducedDims = "IterativeLSI_all5", 
                        name = "Harmony_all5",
                        groupBy = "Sample",
                        force = TRUE
)

##############################################################################################
##                     SECTION5 : Clustering
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##            # 5.1 ## Creating clusters based on Iterative LSI object      ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Uses graph clustering approach from Seurat, as used for scRNAseq
ATACSeq_project_All5 <- addClusters(input = ATACSeq_project_All5,
                                    reducedDims = "IterativeLSI_all5",
                                    method = "Seurat",
                                    name = "Clusters",
                                    resolution = 0.8
)

# Accessing the clusters
head(ATACSeq_project_All5$Clusters)
# tabulate the number of cells present in each cluster:
table(ATACSeq_project_All5$Clusters)

# To better understand which samples reside in which clusters, we can create a 
# cluster confusion matrix across each sample using the confusionMatrix() function.
cM <- confusionMatrix(paste0(ATACSeq_project_All5$Clusters), paste0(ATACSeq_project_All5$Sample))
print(cM)

# plotting the Confusion Matrix
cM <- cM / Matrix::rowSums(cM)
(p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
))


# Adjusting the clustering method and parameters via scran
ATACSeq_project_All5_scran <- addClusters(
  input = ATACSeq_project_All5,
  reducedDims = "IterativeLSI_all5",
  method = "scran",
  name = "ScranClusters",
  k = 15
)

##############################################################################################
##                     SECTION5 : Single Cell Embeddings
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##            # 5.1 ## Visualising clusters via UMAP                          ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

# In ArchR, embeddings, such as Uniform Manifold Approximation and Projection (UMAP) or 
# t-distributed stochastic neighbor embedding (t-SNE), are used to visualize single cells in 
# reduced dimension space. These embeddings each have distinct advantages and disadvantages. 
# We call these “embeddings” because they are strictly used to visualize the clusters and are 
# not used to identify clusters which is done in an LSI sub-space as mentioned in previous chapters.

ATACSeq_project_All5 <- addUMAP(ArchRProj = ATACSeq_project_All5, 
                                reducedDims = "IterativeLSI_all5", 
                                name = "UMAP_all5", 
                                nNeighbors = 30, 
                                minDist = 0.5, 
                                metric = "cosine"
)

p1_UMAP <- plotEmbedding(ArchRProj = ATACSeq_project_All5, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP_all5")
p2_UMAP <- plotEmbedding(ArchRProj = ATACSeq_project_All5, 
                    colorBy = "cellColData",
                    name = "Clusters", 
                    embedding = "UMAP_all5")
(ggAlignPlots(p1_UMAP, p2_UMAP, type = "h"))
plotPDF(p1_UMAP,p2_UMAP, 
        name = "Plot-UMAP-Sample-Clusters.pdf", 
        ArchRProj = ATACSeq_project_All5, 
        addDOC = TRUE, 
        width = 10, 
        height = 10)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##            # 5.2 ## Visualising clusters via t-SNE                         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

ATACSeq_project_All5 <- addTSNE(ArchRProj = ATACSeq_project_All5, 
                               reducedDims = "IterativeLSI_all5", 
                               name = "TSNE_all5",
                               perplexity = 30
)

p1_tSNE <- plotEmbedding(ArchRProj = ATACSeq_project_All5, 
                         colorBy = "cellColData", 
                         name = "Sample", 
                         embedding = "TSNE_all5")
p2_tSNE <- plotEmbedding(ArchRProj = ATACSeq_project_All5, 
                         colorBy = "cellColData",
                         name = "Clusters", 
                         embedding = "TSNE_all5")
(ggAlignPlots(p1_tSNE, p2_tSNE, type = "h"))
plotPDF(p1_tSNE,p2_tSNE, 
        name = "Plot-tSNE-Sample-Clusters.pdf", 
        ArchRProj = ATACSeq_project_All5, 
        addDOC = TRUE, 
        width = 10, 
        height = 10)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##            # 5.3 ## Dimensionality Reduction after Harmony                         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

# can assess the effects of Harmony by visualizing the embedding using UMAP or t-SNE and 
# comparing this to the embeddings visualized in the previous sections for iterative LSI.

ATACSeq_project_All5 <- addUMAP(ArchRProj = ATACSeq_project_All5, 
                                reducedDims = "Harmony_all5", 
                                name = "UMAPHarmony_all5", 
                                nNeighbors = 30, 
                                minDist = 0.5, 
                                metric = "cosine"
)

p3_UMAP <- plotEmbedding(ArchRProj = ATACSeq_project_All5, 
                         colorBy = "cellColData", 
                         name = "Sample", 
                         embedding = "UMAPHarmony_all5")
p4_UMAP <- plotEmbedding(ArchRProj = ATACSeq_project_All5, 
                         colorBy = "cellColData",
                         name = "Clusters", 
                         embedding = "UMAPHarmony_all5")
(ggAlignPlots(p3_UMAP, p4_UMAP, type = "h"))

plotPDF(p1_UMAP,p2_UMAP, p3_UMAP,p4_UMAP, 
        name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", 
        ArchRProj = ATACSeq_project_All5, 
        addDOC = TRUE, 
        width = 10, 
        height = 10)

##############################################################################################
##                     SECTION7 : Gene Scores and Marker Genes
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##            # 6.3 ## Identifying Marker Genes                               ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

markersGS <- getMarkerFeatures(
  ArchRProj = ATACSeq_project_All5, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
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
markerGenesList  <- c("Plag2g2a", "Defa-rs1", "Mmp7", "Lyz1", "Spdef", "Tcf7l2", "Ephb3", "Sis", "Ada", "Lct" )

(heatmapGS <- plotMarkerHeatmap(seMarker = markersGS, 
                                cutOff = "FDR <= 0.01 & Log2FC >= 1.0",
                                limits = c(-3, 3), 
                                returnMatrix = FALSE,
                                plotLog2FC = TRUE, 
                                labelMarkers = markerGenesList, 
                                transpose = FALSE,
                                labelRows = TRUE, clusterCols = TRUE, nPrint = 10))

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 10, height = 10,
        ArchRProj = ATACSeq_project_All5, addDOC = TRUE)



##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##          SECTION 7: Annotating Cell types with a Reference Dataset       ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# ArchR includes a function to align a reference scRNA-seq dataset, and impute cell type annotations based on the reference annotation (addGeneIntegrationMatrix). As a reference, we will use a pre-processed scRNA-seq dataset for human PBMCs.

# Read-in the reference
annotate_reference <- readRDS(file.path("/media/keshavprasad/HornefLab_Data3/scRNA_AnnotationData_Johannes", 
                                        "scrna_with_day25.Rds"))

# add gene integration matrix
ATACSeq_project_All5 <- addGeneIntegrationMatrix(ArchRProj   = ATACSeq_project_All5,
                                            useMatrix   = "GeneScoreMatrix",
                                            matrixName  = "GeneIntegrationMatrix",
                                            reducedDims = "IterativeLSI_all5",
                                            seRNA       = annotate_reference,
                                            addToArrow  = FALSE,
                                            groupRNA    = "int_0.3_broad_tuft",
                                            nameCell    = "predictedCell_Un",
                                            nameGroup   = "predictedGroup_Un",
                                            nameScore   = "predictedScore_Un"
)

# Plot UMAP with predicted cell types
(annotated_UMAP_scRNAseqRef <- plotEmbedding(ATACSeq_project_All5, 
              name = "predictedGroup_Un",
              embedding = "UMAP_all5",
              size = 1.5,
              labelAsFactors = F,
              labelMeans = F))
  
plotPDF(annotated_UMAP_scRNAse, name = "annotated_UMAP_scRNAseqRef_unconstrained",
        width = 10, height = 10, ArchRProj = ATACSeq_project_All5, addDOC = TRUE)
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            Calling Peaks                         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##