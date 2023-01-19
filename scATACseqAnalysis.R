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

# Set the number of threads
addArchRThreads(threads = 8)

here::i_am("scATACseqAnalysis.R")

# Input data 
day12 <- "/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day12Sample/2_Processed_data/V16_d12_230105/fragments.tsv.gz"
inputFiles <- c('scATAC' = day12)

# Add respective Genome
addArchRGenome("mm10")

## Creating Arrow Files: #Each Arrow file stores all of the data associated with an individual sample (i.e. metadata, accessible fragments, and data matrices).
ArrowFiles <- createArrowFiles(inputFiles = inputFiles,
                               sampleNames = names(inputFiles),
                               minTSS = 4, #4 was the default value
                               #Dont set this too high because you can always increase later
                               minFrags = 1000,
                               addTileMat = TRUE,
                               addGeneScoreMat = TRUE,
                               subThreading = TRUE,
                               verbose = TRUE,
                               cleanTmp = TRUE,
                               #force = TRUE, # will recreate arrow files anew each time
                               logFile = createLogFile("d12"))

##############################################################################################
#                             ## SECTION 2: Quality Control
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            ## SECTION2.1 : Doublet Score                 ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##Inferring scATAC-seq Doublets with ArchR : A doublet refers to a single droplet that received a single barcoded bead and more than one nucleus. This causes the reads from more than one cell to appear as a single cell that is effectively the average of the two cells. These are removed computationally!
doubScores <- addDoubletScores(input = ArrowFiles,
                               k = 30, #Refers to how many cells near a "pseudo-doublet" to count
                               knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
                               LSIMethod = 1,
                               force = FALSE)

##############################################################################################
##                     SECTION3 : Setting up ArchR project
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          # 3.1 ## Creating ArchR Project                   ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

ATACSeq_project_d12 <- ArchRProject(ArrowFiles = ArrowFiles,
                                    outputDirectory = "/home/keshavprasad/Documents/scATACseq/d12/",
                                    copyArrows = FALSE #This is recommended so that you maintain an unaltered copy for later usage.
)
paste0("Memory Size = ", round(object.size(ATACSeq_project_d12) / 10^6, digits = 3), " MB")
getAvailableMatrices(ATACSeq_project_d12)

# Inspect the newly created Project.
print(ATACSeq_project_d12)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          # 3.2 Manipulating An ArchRProject                ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Cell names
head(ATACSeq_project_d12$cellNames)
# Sample names
head(ATACSeq_project_d12$Sample)
# TSS Enrichment
head(ATACSeq_project_d12$TSSEnrichment)
# One can access the TSS Enrichment Scores for each cell:
quantile(ATACSeq_project_d12$TSSEnrichment)
# doublet Enrichment 
head(ATACSeq_project_d12$DoubletEnrichment)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                          # 3.3 Obtaining columns from cellColData         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
DoubletEnrichment_df <- getCellColData(ATACSeq_project_d12, select = "DoubletEnrichment")
head(DoubletEnrichment_df, 30)


##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## #3.3 Plotting QC metrics - log10(Unique Fragments) vs TSS enrichment score      ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
QCMetric_df <- getCellColData(ATACSeq_project_d12, select = c("log10(nFrags)", "TSSEnrichment"))
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
        ArchRProj = ATACSeq_project_d12, 
        addDOC = TRUE)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## #3.4 Plotting Sample Statistics from an ArchRProject                       ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# Make a violin plot for each sample for the TSS enrichment scores.

### violin plot for each sample for TSSEnrichment.
Group_plot_TSS_violin <- plotGroups(ArchRProj = ATACSeq_project_d12,
                                     groupBy = "Sample",
                                     colorBy = "cellColData",
                                     name = "TSSEnrichment",
                                     plotAs = "violin",
                                     alpha = 0.4,
                                     addBoxPlot = TRUE) + ggtitle("TSSEnrichment")

### Plots (per sample) for log10 (unique nuclear fragments)
#### log10 (unique nuclear fragments)
Group_plot_nFrags_violin <- plotGroups(ArchRProj = ATACSeq_project_d12,
                                        groupBy = "Sample",
                                        colorBy = "cellColData",
                                        name = "log10(nFrags)",
                                        plotAs = "violin",
                                        alpha = 0.4,
                                        addBoxPlot = TRUE) + ggtitle("nFrags")
### violin plot for each sample for TSSEnrichment.
Group_plot_TSS_violin <- plotGroups(ArchRProj = ATACSeq_project_d12,
                                     groupBy = "Sample",
                                     colorBy = "cellColData",
                                     name = "TSSEnrichment",
                                     plotAs = "violin",
                                     alpha = 0.4,
                                     addBoxPlot = TRUE) + ggtitle("TSSEnrichment")

Group_plot_BLR_violin <- plotGroups(ArchRProj = ATACSeq_project_d12,
                                     groupBy = "Sample",
                                     colorBy = "cellColData",
                                     name = "BlacklistRatio",
                                     plotAs = "violin",
                                     alpha = 0.4,
                                     addBoxPlot = TRUE) + ggtitle("BlacklistRatio")

Group_plot_NR_violin <- plotGroups(ArchRProj = ATACSeq_project_d12,
                                    groupBy = "Sample",
                                    colorBy = "cellColData",
                                    name = "NucleosomeRatio",
                                    plotAs = "violin",
                                    alpha = 0.4,
                                    addBoxPlot = TRUE) + ggtitle("NucleosomeRatio")

Group_plot_DS_violin <- plotGroups(ArchRProj = ATACSeq_project_d12,
                                    groupBy = "Sample",
                                    colorBy = "cellColData",
                                    name = "DoubletScore",
                                    plotAs = "violin",
                                    alpha = 0.4,
                                    addBoxPlot = TRUE) + ggtitle("DoubletScore")

plotPDF(Group_plot_nFrags_violin, 
        Group_plot_TSS_violin, 
        Group_plot_BLR_violin ,
        Group_plot_NR_violin,
        Group_plot_DS_violin, 
        name = "QC-Sample-Statistics.pdf", 
        ArchRProj = ATACSeq_project_d12, 
        addDOC = TRUE, width = 12, height = 12)

## Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
(FragSizePlot = plotFragmentSizes(ArchRProj = ATACSeq_project_d12))
(TSSEnrichmentPlot = plotTSSEnrichment(ArchRProj = ATACSeq_project_d12))
plotPDF(FragSizePlot, 
        TSSEnrichmentPlot, 
        name = "QC-Sample-FragSizes-TSSProfile.pdf", 
        ArchRProj = ATACSeq_project_d12, 
        addDOC = TRUE, width = 10, height = 10)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## # 3.6 Filtering Doublets from an ArchRProject                       ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# We filter putative doublets based on the previously determined doublet scores using the filterDoublets() function. This doesn’t physically remove data from the Arrow files but rather tells the ArchRProject to ignore these cells for downstream analysis. 
# The higher the filterRatio, the greater the number of cells potentially removed as doublets.

## Finding doublets
ATACSeq_project_d12 <- filterDoublets(ArchRProj = ATACSeq_project_d12)



#  Saving both these plots and ArchR Project
# plotPDF(FragSizePlot, TSSEnrichmentPlot,
#         name = "QC-Sample-FragSizes-TSSProfile.pdf",
#         ArchRProj = ATACSeq_project_d12,
#         addDOC = TRUE, width = 8, height = 8)

# saveArchRProject(ArchRProj = ATACSeq_project_d12, 
#                  outputDirectory = "/home/keshavprasad/Documents/scATACseq/d12/", 
#                  load = FALSE)




# 
# ### Doublet score
# Doublet_score_df <- as.data.frame(ATACSeq_project_d12$DoubletScore)
# quantile(ATACSeq_project_d12$DoubletScore)
# ### TSS Enrichment Scores for each cell:
# quantile(ATACSeq_project_d12$TSSEnrichment)
# 
# # Filter out Low Quality Cells
# ATACSeq_project_d12 <- ATACSeq_project_d12[ATACSeq_project_d12$TSSEnrichment > 5 &
#                                            ATACSeq_project_d12$nFrags > 2000 &
#                                            ATACSeq_project_d12$NucleosomeRatio < 2]







ATACSeq_project <- addIterativeLSI(
  ArchRProj = ATACSeq_project,
  useMatrix = "TileMatrix",
  name = "ATACSeq_LSI",
  iterations = 2,
  clusterParams = list(resolution = 0.2, sampleCells = 10000, n.start = 10),#See Seurat::FindClusters
  varFeatures = 10000,
  dimsToUse = 1:30,
  #force = TRUE
)
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            UMAP on the LSI results                         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
ATACSeq_project <- addUMAP(ATACSeq_project,
                           reducedDims = "ATACSeq_LSI",
                           name = "UMAP_ATAC",
                           minDist = 0.8,
                           #force = TRUE
)

# Clustering
ATACSeq_project <- addClusters(
  input = ATACSeq_project,
  reducedDims = "ATACSeq_LSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  #force = TRUE
)

plotEmbedding(ATACSeq_project,
              name = "Clusters",
              embedding = "UMAP_ATAC",
              size = 1.5,
              labelAsFactors = F,
              labelMeans=F)

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
                     embedding = "UMAP_ATAC"))
p2 <- (plotEmbedding(ArchRProj = ATACSeq_project,
                     colorBy = "cellColData",
                     name = "Clusters",
                     embedding = "UMAP_ATAC"))

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
markersGS <- getMarkerFeatures(ArchRProj = ATACSeq_project, useMatrix = "GeneScoreMatrix",
                               groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon")
## Markers List
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.0")
cluster1_df <- as.data.frame(markerList$C1)
write.csv(cluster1_df, file.path(here(), "/cluster1_df.csv"))
cluster2_df <- as.data.frame(markerList$C2)
write.csv(cluster2_df, file.path(here(), "/cluster2_df.csv"))
cluster3_df <- as.data.frame(markerList$C3)
write.csv(cluster3_df, file.path(here(), "/cluster3_df.csv"))
cluster4_df <- as.data.frame(markerList$C4)
write.csv(cluster4_df, file.path(here(), "/cluster4_df.csv"))
cluster5_df <- as.data.frame(markerList$C5)
write.csv(cluster5_df, file.path(here(), "/cluster5_df.csv"))
cluster6_df <- as.data.frame(markerList$C6)
write.csv(cluster6_df, file.path(here(), "/cluster6_df.csv"))
cluster7_df <- as.data.frame(markerList$C7)
write.csv(cluster7_df, file.path(here(), "/cluster7_df.csv"))
cluster8_df <- as.data.frame(markerList$C8)
write.csv(cluster8_df, file.path(here(), "/cluster8_df.csv"))
cluster9_df <- as.data.frame(markerList$C9)
write.csv(cluster9_df, file.path(here(), "/cluster9_df.csv"))
cluster10_df <- as.data.frame(markerList$C10)
write.csv(cluster10_df, file.path(here(), "/cluster10_df.csv"))
##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##          SECTION3. Marker Genes                                           ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##  visualize all of the marker features simultaneously
markerGenesList  <- c("Plag2g2a", "Defa-rs1", "Mmp7", "Lyz1", "Spdef", "Tcf7l2", "Ephb3", "Sis", "Ada", "Lct" )
markerGenesList2 <- c("Mmp7", "Lyz1", "Spdef", "Tcf7l2", "Ephb3", "Sis", "Ada", "Lct" )
(heatmapGS <- plotMarkerHeatmap(seMarker = markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.0",
                                #limits = c(-3, 3), #returnMatrix = TRUE
                                plotLog2FC = TRUE, labelMarkers = markerGenesList2, transpose = FALSE,
                                labelRows = TRUE, clusterCols = TRUE, nPrint = 10))
# ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 8,
        ArchRProj = ATACSeq_project, addDOC = FALSE)

# The gene activities can be used to visualize the expression of marker genes on the scATAC-seq clusters.
# features <- getFeatures(ArchRProj = ATACSeq_project, useMatrix = "GeneScoreMatrix")

# Visualizing Marker Genes on an Embedding :
(MarkerGeneEmbeddingPlot <- plotEmbedding(ArchRProj = ATACSeq_project, colorBy = "GeneScoreMatrix",
                                          name = markerGenesList2, embedding = "UMAP_ATAC",
                                          quantCut = c(0.01, 0.95), imputeWeights = NULL))
# Plot all marker genes via cow plot
MarkerGeneEmbedding_CowPlot <- lapply(MarkerGeneEmbeddingPlot, function(x){ x + guides(color = FALSE,
                                                                                       fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
})
#do.call(cowplot::plot_grid, c(list(ncol = 3),MarkerGeneEmbedding_CowPlot))
patchwork::wrap_plots(MarkerGeneEmbedding_CowPlot)
plotPDF(plotList = MarkerGeneEmbedding_CowPlot, name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf",
        ArchRProj = ATACSeq_project, addDOC = FALSE, width = 8, height = 8)

##–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##          SECTION 4: Annotating Cell types with a Reference Dataset       ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
# ArchR includes a function to align a reference scRNA-seq dataset, and impute cell type annotations based on the reference annotation (addGeneIntegrationMatrix). As a reference, we will use a pre-processed scRNA-seq dataset for human PBMCs.

# Read-in the reference
annotate_reference <- readRDS(file.path("D:/scRNA_AnnotationData_Johannes", "scrna_with_day25.Rds"))

# add gene integration matrix
ATACSeq_project <- addGeneIntegrationMatrix(ArchRProj   = ATACSeq_project,
                                            useMatrix   = "GeneScoreMatrix",
                                            matrixName  = "GeneIntegrationMatrix",
                                            reducedDims = "ATACSeq_LSI",
                                            seRNA       = annotate_reference,
                                            addToArrow  = FALSE,
                                            groupRNA    = "int_0.3_broad_tuft",
                                            nameCell    = "predictedCell_Un",
                                            nameGroup   = "predictedGroup_Un",
                                            nameScore   = "predictedScore_Un"
)

# Plot UMAP with predicted cell types
plotEmbedding(ATACSeq_project, name = "predictedGroup_Un",
              embedding = "UMAP_ATAC",
              size = 1.5,
              labelAsFactors = F,
              labelMeans = F)

##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            Calling Peaks                         ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##



##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
##                            ## SECTION : Per Cell Quality Control        ##
##––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
## Plotting Sample Statistics from an ArchRProject



Group_plot_nFrags_violin + Group_plot_TSS_violin + Group_plot_BLR_violin +
  Group_plot_NR_violin +
  # Group_plot_DS_violin +
  patchwork::plot_layout(nrow = 2)