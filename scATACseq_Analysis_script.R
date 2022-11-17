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
suppressPackageStartupMessages(library("dplyr"))
set.seed(1234)

# addArchRThreads(threads = 6)
here::i_am("scATACseq_Analysis_script.R")

# Data Input
day5 <- "D:/ATACSeq/Data/Sample1/2_Processed_data_S1/AJ_v12/fragments.tsv.gz"

day10 <- "D:/ATACSeq/Data/Sample2/220928_Juergens_Hornef_microbiology_scATACseq/compressed_tars/_2_Processed_data_S2/2_Processed_data/AJ_V12_d10_L129_UI/fragments.tsv.gz"

inputFiles <- c('scATAC' = c(day5, day10))

# Add respective Genome
addArchRGenome("mm10")

## Creating Arrow Files: #Each Arrow file stores all of the data associated with an individual sample (i.e. metadata, accessible fragments, and data matrices).
ArrowFiles <- createArrowFiles(inputFiles = inputFiles,
                                sampleNames = names(inputFiles),
                                minTSS = 6, #4 was the default value
                               #Dont set this too high because you can always increase later
                                minFrags = 1000,
                                addTileMat = TRUE,
                                addGeneScoreMat = TRUE,
                                subThreading = TRUE,
                                verbose = TRUE,
                                cleanTmp = TRUE,
                                #force = TRUE, # will recreate arrow files anew each time
                                logFile = createLogFile("createArrows_Nov")
)

## Creating ArchR Project
ATACSeq_project <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "D:/ATACSeq/",
  copyArrows = TRUE #This is recommended so that you maintain an unaltered copy for later usage.
)
getAvailableMatrices(ATACSeq_project)


##############################################################################################
# SECTION 2: Quality Control

## SECTION2.1 : Doublet Removal
##Inferring scATAC-seq Doublets with ArchR : A doublet refers to a single droplet that received a single barcoded bead and more than one nucleus. This causes the reads from more than one cell to appear as a single cell that is effectively the average of the two cells. These are removed computationally!
doubScores <- addDoubletScores(input = ArrowFiles,
                               k = 30, #Refers to how many cells near a "pseudo-doublet" to count
                               knnMethod = "LSI", #Refers to the embedding to use for nearest neighbor search.
                               LSIMethod = 1,
                               force = TRUE
)
# We filter putative doublets based on the previously determined doublet scores using the filterDoublets() function. This doesn’t physically remove data from the Arrow files but rather tells the ArchRProject to ignore these cells for downstream analysis. The higher the filterRatio, the greater the number of cells potentially removed as doublets.
# # Finding doublets
ATACSeq_project <- filterDoublets(ArchRProj = ATACSeq_project)

### Doublet score
Doublet_score_df <- as.data.frame(ATACSeq_project$DoubletScore)
quantile(ATACSeq_project$DoubletScore)
### TSS Enrichment Scores for each cell:
quantile(ATACSeq_project$TSSEnrichment)

# Filter out Low Quality Cells
ATACSeq_project <- ATACSeq_project[ATACSeq_project$TSSEnrichment > 6 &
                                   ATACSeq_project$nFrags > 2000 &
                                   ATACSeq_project$NucleosomeRatio < 2]

### Plotting QC metrics - log10(Unique Fragments) vs TSS enrichment score
df <- getCellColData(ATACSeq_project, select = c("log10(nFrags)", "TSSEnrichment"))

(scatterplot_FragsVsEnrichment <- ggPoint(x = df[,1], y = df[,2],
                                          colorDensity = TRUE,
                                          continuousSet = "sambaNight",
                                          xlabel = "Log10 Unique Fragments",
                                          ylabel = "TSS Enrichment",
                                          xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
                                          ylim = c(0, quantile(df[,2], probs = 0.99))
                                      ) +
                                    geom_hline(yintercept = 4, lty = "dashed") +
                                    geom_vline(xintercept = 3, lty = "dashed"))

## Plotting Sample Statistics from an ArchRProject

### Plots (per sample) for log10 (unique nuclear fragments)
#### log10 (unique nuclear fragments)
(Group_plot_nFrags_violin <- plotGroups(ArchRProj = ATACSeq_project,
                                       groupBy = "Sample",
                                       colorBy = "cellColData",
                                       name = "log10(nFrags)",
                                       plotAs = "violin",
                                       alpha = 0.4,
                                       addBoxPlot = TRUE
)) + ggtitle("nFrags")
### violin plot for each sample for TSSEnrichment.
(Group_plot_TSS_violin <- plotGroups(ArchRProj = ATACSeq_project,
                                        groupBy = "Sample",
                                        colorBy = "cellColData",
                                        name = "TSSEnrichment",
                                        plotAs = "violin",
                                        alpha = 0.4,
                                        addBoxPlot = TRUE
)) + ggtitle("TSSEnrichment")

(Group_plot_BLR_violin <- plotGroups(ArchRProj = ATACSeq_project,
                                     groupBy = "Sample",
                                     colorBy = "cellColData",
                                     name = "BlacklistRatio",
                                     plotAs = "violin",
                                     alpha = 0.4,
                                     addBoxPlot = TRUE
)) + ggtitle("BlacklistRatio")

(Group_plot_NR_violin <- plotGroups(ArchRProj = ATACSeq_project,
                                     groupBy = "Sample",
                                     colorBy = "cellColData",
                                     name = "NucleosomeRatio",
                                     plotAs = "violin",
                                     alpha = 0.4,
                                     addBoxPlot = TRUE
)) + ggtitle("NucleosomeRatio")

(Group_plot_DS_violin <- plotGroups(ArchRProj = ATACSeq_project,
                                    groupBy = "Sample",
                                    colorBy = "cellColData",
                                    name = "DoubletScore",
                                    plotAs = "violin",
                                    alpha = 0.4,
                                    addBoxPlot = TRUE
)) + ggtitle("DoubletScore")

Group_plot_nFrags_violin + Group_plot_TSS_violin + Group_plot_BLR_violin +
  Group_plot_NR_violin + Group_plot_DS_violin + patchwork::plot_layout(nrow = 1)

## Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
(FragSizePlot <- plotFragmentSizes(ArchRProj = ATACSeq_project))
(TSSEnrichmentPlot <- plotTSSEnrichment(ArchRProj = ATACSeq_project))
FragSizePlot + TSSEnrichmentPlot + patchwork::plot_layout(nrow = 1)

#Saving both these plots and ArchR Project
plotPDF(FragSizePlot, TSSEnrichmentPlot,name = "QC-Sample-FragSizes-TSSProfile.pdf",
        ArchRProj = ATACSeq_project, addDOC = FALSE, width = 8, height = 8)
saveArchRProject(ArchRProj = ATACSeq_project, outputDirectory = "D:/ATACSeq/", load = FALSE)


## SECTION3 : Dimensionality Reduction and Clustering
ATACSeq_project <- addIterativeLSI(
    ArchRProj = ATACSeq_project,
    useMatrix = "TileMatrix",
    name = "ATACSeq_LSI",
    iterations = 2,
    clusterParams = list( #See Seurat::FindClusters
        resolution = 0.2,
        sampleCells = 10000,
        n.start = 10
    ),
    varFeatures = 10000,
    dimsToUse = 1:30,
    force = TRUE
)


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            UMAP on the LSI results                         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
archrproj <- addUMAP(ATACSeq_project,
                     reducedDims = "ATACSeq_LSI",
                     name = "UMAP_ATAC",
                     minDist = 0.8,
                     force = TRUE
                     )

# Clustering
ATACSeq_project <- addClusters(
    input = ATACSeq_project,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force = TRUE
)


# ### Clustering using scran
# ATACSeq_project <- addClusters(
#     input = ATACSeq_project,
#     reducedDims = "IterativeLSI",
#     method = "scran",
#     name = "ScranClusters",
#     k = 15,
#     force = TRUE
# )

## Clusters:
### number of cells present in each cluster:
table(ATACSeq_project$Clusters)

### Cluster Confusion Matrix
confusion_Matrix <- confusionMatrix(paste0(ATACSeq_project$Clusters),
                                    paste0(ATACSeq_project$Sample))

confusion_Matrix <- confusion_Matrix / Matrix::rowSums(confusion_Matrix)
(pheatmap::pheatmap(
   mat = as.matrix(confusion_Matrix),
   color = paletteContinuous("whiteBlue"),
   border_color = "black"))

p1 <- (plotEmbedding(ArchRProj = ATACSeq_project,
                    colorBy = "cellColData",
                    name = "Sample",
                    embedding = "UMAP"))
p2 <- (plotEmbedding(ArchRProj = ATACSeq_project,
              colorBy = "cellColData",
              name = "Clusters",
              embedding = "UMAP"))

ggAlignPlots(p1, p2, type = "h")
(plotPDF(p1,p2,
         name = "Plot-UMAP-Sample-Clusters.pdf",
         ArchRProj = ATACSeq_project,
         addDOC = TRUE,
         width = 9,
         height = 9))

## t-Stocastic Neighbor Embedding (t-SNE)
# ATACSeq_project <- addTSNE(
#     ArchRProj = ATACSeq_project,
#     reducedDims = "IterativeLSI",
#     name = "TSNE",
#     perplexity = 30
# )
#
# ## Plotting tSNE
# q1 <- plotEmbedding(ArchRProj = ATACSeq_project, colorBy = "cellColData",
#                     name = "Sample", embedding = "TSNE")
# q2 <- plotEmbedding(ArchRProj = ATACSeq_project, colorBy = "cellColData",
#                     name = "Clusters", embedding = "TSNE")
# ggAlignPlots(q1, q2, type = "h")
# plotPDF(q1,q2, name = "Plot-TSNE-Sample-Clusters.pdf", ArchRProj = ATACSeq_project,
#         addDOC = FALSE, width = 10, height = 10)

# Gene Scores and Marker Genes with ArchR
## Identifying Marker Genes
markersGS <- getMarkerFeatures(
    ArchRProj = ATACSeq_project,
    useMatrix = "GeneScoreMatrix",
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
## Markers List
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.0")
cluster1_df <- as.data.frame(markerList$C1)
cluster2_df <- as.data.frame(markerList$C2)
cluster3_df <- as.data.frame(markerList$C3)
cluster4_df <- as.data.frame(markerList$C4)
cluster5_df <- as.data.frame(markerList$C5)
cluster6_df <- as.data.frame(markerList$C6)
cluster7_df <- as.data.frame(markerList$C7)
cluster8_df <- as.data.frame(markerList$C8)
cluster9_df <- as.data.frame(markerList$C9)
cluster10_df <- as.data.frame(markerList$C10)
# cluster11_df <- as.data.frame(markerList$C11)
# cluster12_df <- as.data.frame(markerList$C12)
# cluster13_df <- as.data.frame(markerList$C13)

# 3. Marker Genes
##  visualize all of the marker features simultaneously
markerGenes  <- c("Plag2g2a", "Defa-rs1", "Mmp7",
                  "Lyz1", "Spdef", "Tcf7l2", "Ephb3", "Sis", "Ada", "Lct" )
heatmapGS <- markerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.0",
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

