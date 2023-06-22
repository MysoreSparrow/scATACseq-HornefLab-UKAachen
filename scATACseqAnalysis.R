#' ---
#' title: "scATAC-seq Analysis"
#' author: "KeshavPrasadGubbi"
#' output:
#'    html_document:
#'      toc: yes
#'      toc_depth: 5
#'      toc_float: yes
#'      fig_width: 8
#'      fig_height: 8
#'      fig_caption: yes
#'      number_sections: yes
#' ---

#' SECTION1: Initialization

Package_List <- c("ArchR", "ComplexHeatmap", "writexl", "here", "patchwork", "tidyverse", 
                  "pheatmap", "ArchRtoSignac", "Signac", "Seurat", "biovizBase", "stringr")
not_installed <- Package_List[!(Package_List %in% installed.packages()[, "Package"])] 
# Extract not installed packages
if (length(not_installed)) install.packages(not_installed) # Install the uninstalled packages
invisible(lapply(Package_List, suppressPackageStartupMessages(library), character.only = TRUE))

# Set the number of threads
# addArchRThreads(threads = 8)

# File Path Declarations
here::i_am(path = "scATACseqAnalysis.R")
paste0(here())

# set seed
set.seed(1989)

# Function to save generic plots
saveplot <- function(plot, plotname) {
  # Function to save the plots
  extension <- ".jpeg"
  ggsave(
    filename = file.path(getOutputDirectory(ATACSeq_project_All5), 
                         paste(plotname, extension), sep = ""),
    plot = plot,
    dpi = 300,
    width = 8,
    height = 8,
    units = "in"
  )
  # dev.off()
  while (!is.null(dev.list())) {
    dev.off()
  }
}


#' Load the ArchR Project
ATACSeq_project_All5 <- loadArchRProject(path = "F:/scATACseq/All5/",
                                         force = FALSE, 
                                         showLogo = FALSE)
# Aesthetics
my_gg_theme <- theme(axis.title = element_text(size = 15),
                    axis.text = element_text(size = 15),
                    legend.text = element_text(size = 10)
)
print(getOutputDirectory(ATACSeq_project_All5))


#' Dimensional Reduction, INTEGRATION, CLUSTERING & EMBEDDING
# 4.1 ## Reducing Dims via Iterative LSI.# The most common parameters to tweak are iterations, varFeatures, and resolution

ATACSeq_project_All5 <- addIterativeLSI(
  ArchRProj = ATACSeq_project_All5,
  useMatrix = "TileMatrix",
  name = "IterativeLSI_all5",
  iterations = 3,
  saveIterations = TRUE,
  depthCol = "nFrags",
  # corCutOff = 0.75,
  clusterParams = list(resolution = c(0.2, 0.1, 0.05) , sampleCells = 10000, n.start = 10),
  varFeatures = 15000,
  dimsToUse = 1:30,
  seed = 1989,
  force = TRUE, 
  verbose = TRUE
)

#############################
# Compute CorCutOff Values
CorCutOffValues <- ATACSeq_project_All5@reducedDims@listData[["IterativeLSI_all5"]]@listData[["corToDepth"]][["none"]] %>% as.data.frame() %>% tibble::rownames_to_column(., "LSI_Dim") %>% rename("Cor_Values" = ".")
print(head(CorCutOffValues, 10))

################################

##' Using LSI object as Input to correct for Batch Effects via Harmony
ATACSeq_project_All5 <- addHarmony( # 4.4 ## Batch Effect Correction with HARMONY   ##
  ArchRProj = ATACSeq_project_All5,
  reducedDims = "IterativeLSI_all5",
  name = "Harmony_all5",
  groupBy = "Sample",
  corCutOff = 0.75,
  # dimsToUse = 2:30,
  verbose = TRUE,
  force = TRUE
) 

##' Create a UMAP
ATACSeq_project_All5 <- addUMAP( # 5.1 ## Visualising clusters via UMAP
  ArchRProj = ATACSeq_project_All5,
  reducedDims = "Harmony_all5", # Need to use the Harmony object for clustering and UMAP for proper integration
  name = "UMAP_all5",
  nNeighbors = 30,
  minDist = 0.4,
  # dimsToUse = 2:30,
  metric = "cosine",
  force = TRUE,
  verbose = TRUE,
  saveModel = TRUE
)

##' Create UMAP Embedding for all samples together
(p1_UMAP <- plotEmbedding(
  ArchRProj = ATACSeq_project_All5,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_all5",
  keepAxis = TRUE
) + my_gg_theme)

plotPDF(p1_UMAP, name = "Plot-UMAP-Samplewise.pdf", ArchRProj = ATACSeq_project_All5,
        width = 10, height = 10)

saveplot(plot = p1_UMAP, plotname = "V2_UMAP-Samplewise")

#' Create SplitView UMAP
my_gg_theme <- theme(axis.title = element_text(size = 15),
                     axis.text = element_text(size = 15),
                     legend.text = element_text(size = 10))
# Create a Split view UMAP so that the UMAP for each sample can be visualized separately.
create_splitUMAP <- function(ATACSeq_project_All5, sampleObject) {
  split_UMAP <- plotEmbedding(
    ArchRProj = ATACSeq_project_All5,
    colorBy = "cellColData",
    name = "Sample",
    embedding = "UMAP_all5",
    highlightCells = getCellNames(ArchRProj = ATACSeq_project_All5)[which(ATACSeq_project_All5@cellColData$Sample == sampleObject)],
    keepAxis = TRUE,
    baseSize = 15,
    plotAs = "points",
    size = 0.05
  )
  return(split_UMAP)
}

# Get the the Samples from the ArchR project
samplelist <- ATACSeq_project_All5@cellColData$Sample
# Call the function to create UMAP for each sample
(d1_UMAP <- create_splitUMAP(ATACSeq_project_All5, samplelist@values[1]) + my_gg_theme)
(d5_UMAP <- create_splitUMAP(ATACSeq_project_All5, samplelist@values[5]) + my_gg_theme)
(d10_UMAP <- create_splitUMAP(ATACSeq_project_All5, samplelist@values[4]) + my_gg_theme)
(d25_UMAP <- create_splitUMAP(ATACSeq_project_All5, samplelist@values[3]) + my_gg_theme)
(d4pi_UMAP <- create_splitUMAP(ATACSeq_project_All5, samplelist@values[2]) + my_gg_theme)

(ggAlignPlots(d1_UMAP, d5_UMAP, d10_UMAP, d25_UMAP, d4pi_UMAP, type = "h"))
plotPDF(d1_UMAP, d5_UMAP, d10_UMAP, d25_UMAP, d4pi_UMAP,
        name = "Plot-Splitview UMAP-Sample-Clusters.pdf",
        ArchRProj = ATACSeq_project_All5, 
        width = 10,
        height = 10
) 
saveplot(plot = d1_UMAP, plotname = "V2_d1_UMAP_SampleWise")
saveplot(plot = d5_UMAP, plotname = "V2_d5_UMAP_SampleWise")
saveplot(plot = d10_UMAP, plotname = "V2_d10_UMAP_SampleWise")
saveplot(plot = d25_UMAP, plotname = "V2_d25_UMAP_SampleWise")
saveplot(plot = d4pi_UMAP, plotname = "V2_d4pi_UMAP_SampleWise")

#' Adding Clusters

# Adding Clusters
ATACSeq_project_All5 <- addClusters( # 5.1 ## Creating clusters based on Iterative LSI object      ##
  input = ATACSeq_project_All5,
  reducedDims = "Harmony_all5", # Harmony_all5
  method = "Seurat",
  name = "Clusters_all5",
  resolution = 0.5,
  force = TRUE
)
ATACSeq_project_All5 <- addUMAP( # 5.1 ## Visualising clusters via UMAP
  ArchRProj = ATACSeq_project_All5,
  reducedDims = "Harmony_all5", # Need to used the HArmony object for clustering and UMAP to ensure proper integration
  name = "ClustersUMAP_all5",
  nNeighbors = 30,
  minDist = 0.4,
  metric = "cosine",
  force = TRUE,
  verbose = TRUE,
  saveModel = TRUE
)

(p2_UMAP <- plotEmbedding(ArchRProj = ATACSeq_project_All5, 
                          colorBy = "cellColData", 
                          name = "Clusters_all5", 
                          embedding = "ClustersUMAP_all5", 
                          plotAs = "points",
                          size = 0.05,
                          keepAxis = TRUE) + my_gg_theme)
plotPDF(p2_UMAP, name = "Plot-UMAP-Clusters.pdf", ArchRProj = ATACSeq_project_All5,  width = 10, height = 10)
saveplot(plot = p2_UMAP, plotname = "V2_UMAP-Clusters_CorCutOff0.75")

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
  mat = as.matrix(cM), color = paletteContinuous("whiteBlue"),
  border_color = "black"
) + my_gg_theme)



## –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
#'          SECTION 7: Annotating Cell types with a Reference Dataset       ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

##           7.1 Unconstrained Integration       ##
# Read-in the reference
annotate_reference <- readRDS(file.path("F:/scRNA_AnnotationData_Johannes", "scrna_with_day25.Rds"))

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
                                             baseSize = 15,
                                             plotAs = "points",
                                             size = 0.05,
                                             labelAsFactors = F,
                                             labelMeans = F))

plotPDF(annotated_UMAP_scRNAseqRef, name = "annotated_UMAP_scRNAseqRef_unconstrained", width = 10, height = 10, ArchRProj = ATACSeq_project_All5, addDOC = TRUE)
saveplot(plot = annotated_UMAP_scRNAseqRef, plotname = "V2_annotated_UMAP_scRNAseqRef")

#' Create Split Annotated UMAP for each cellType
# SampleWise
create_SampleWise_splitAnnotatedUMAP <- function(ATACSeq_project_All5, CellType) {
  split_AnnotatedUMAP <- plotEmbedding(
    ArchRProj = ATACSeq_project_All5,
    colorBy = "cellColData",
    name = "Sample", # ATACSeq_project_All5@cellColData@listData$predictedGroup_Un
    embedding = "UMAP_all5",
    highlightCells = getCellNames(ArchRProj = ATACSeq_project_All5)[which(ATACSeq_project_All5@cellColData@listData$predictedGroup_Un == CellType)],
    keepAxis = TRUE,
    baseSize = 15,
    plotAs = "points",
    size = 0.05
  )
  return(split_AnnotatedUMAP)
}

# Call the function to create UMAP for each sample
(Enterocyte_UMAP <- create_SampleWise_splitAnnotatedUMAP(ATACSeq_project_All5, "Enterocyte") + my_gg_theme + ggtitle("Enterocyte_UMAP"))
(GP_UMAP <- create_SampleWise_splitAnnotatedUMAP(ATACSeq_project_All5, "Goblet+Paneth") + my_gg_theme + ggtitle("Goblet+Paneth_UMAP"))
(EEC_UMAP <- create_SampleWise_splitAnnotatedUMAP(ATACSeq_project_All5, "EEC") + my_gg_theme + ggtitle("EEC_UMAP"))
(Stem_UMAP <- create_SampleWise_splitAnnotatedUMAP(ATACSeq_project_All5, "Stem") + my_gg_theme + ggtitle("Stem_UMAP"))
(Tuft_UMAP <- create_SampleWise_splitAnnotatedUMAP(ATACSeq_project_All5, "Tuft") + my_gg_theme + ggtitle("Tuft_UMAP"))

plotPDF(Enterocyte_UMAP, GP_UMAP, EEC_UMAP, Stem_UMAP, Tuft_UMAP,
        name = "Plot- SampleWise_Splitview-CellTypeAnnotatedUMAP.pdf",
        ArchRProj = ATACSeq_project_All5, 
        width = 10,
        height = 10
) 
saveplot(plot = Enterocyte_UMAP, plotname = "V2_Enterocyte_UMAP")
saveplot(plot = GP_UMAP, plotname = "V2_GP_UMAP")
saveplot(plot = EEC_UMAP, plotname = "V2_EEC_UMAP")
saveplot(plot = Stem_UMAP, plotname = "V2_Stem_UMAP")
saveplot(plot = Tuft_UMAP, plotname = "V2_Tuft_UMAP")

# CellTypeWise
# create_CellTypeWise_splitAnnotatedUMAP <- function(ATACSeq_project_All5, CellType) {
#   ct = ATACSeq_project_All5@cellColData@listData$predictedGroup_Un == CellType
#   ct_split_AnnotatedUMAP <- plotEmbedding(
#     ArchRProj = ATACSeq_project_All5,
#     colorBy = "cellColData",
#     name = "ct", # ATACSeq_project_All5@cellColData@listData$predictedGroup_Un
#     embedding = "UMAP_all5",
#     highlightCells = getCellNames(ArchRProj = ATACSeq_project_All5)[which(ATACSeq_project_All5@cellColData@listData$predictedGroup_Un == CellType)],
#     keepAxis = TRUE,
#     baseSize = 10,
#     plotAs = "points"
#   )
#   return(ct_split_AnnotatedUMAP)
# }
# 
# # Call the function to create UMAP for each sample
# (Enterocyte_UMAP1 <- create_CellTypeWise_splitAnnotatedUMAP(ATACSeq_project_All5, "Enterocyte") + my_gg_theme)
# (GP_UMAP1 <- create_CellTypeWise_splitAnnotatedUMAP(ATACSeq_project_All5, "Goblet+Paneth") + my_gg_theme)
# (EEC_UMAP1 <- create_CellTypeWise_splitAnnotatedUMAP(ATACSeq_project_All5, "EEC") + my_gg_theme)
# (Stem_UMAP1 <- create_CellTypeWise_splitAnnotatedUMAP(ATACSeq_project_All5, "Stem") + my_gg_theme)
# (Tuft_UMAP1 <- create_CellTypeWise_splitAnnotatedUMAP(ATACSeq_project_All5, "Tuft") + my_gg_theme)
# plotPDF(Enterocyte_UMAP1, GP_UMAP1, EEC_UMAP1, Stem_UMAP1, Tuft_UMAP1,
#         name = "Plot- SampleWise_Splitview-CellTypeAnnotatedUMAP.pdf",
#         ArchRProj = ATACSeq_project_All5, 
#         width = 10,
#         height = 10
# ) 
##############################################################################################
#'                     SECTION7 : Gene Scores and Marker Genes
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

##'            Identifying Marker Genes                               ##
markersGS <- getMarkerFeatures(
  ArchRProj = ATACSeq_project_All5,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters_all5",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

#' Markers List
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0")
# Create an empty list to store dataframes
df_list <- list()

# Loop over the list and convert each element to a dataframe
for (i in 1:length(markerList)) {
  cluster_name <- names(markerList)[i]
  cluster_df <- as.data.frame(markerList[[i]])
  write.csv(cluster_df, "cluster_df.csv")
  
  # Add the dataframe to the list with a unique name
  df_list[[paste0(cluster_name, "_sheet")]] <- cluster_df
}

# Write the dataframes to different sheets of an Excel file
write_xlsx(df_list, path = "Clusters_ATACseq_AllGenes.xlsx")

#' SECTION3. Marker Genes
# Human_and_mouse_cell_markers_df <- as.data.frame(read.delim("D:/scATACseq/Human_and_mouse_cell_markers-Markers.tsv"))
# 
# Haber_marker_genes <- Human_and_mouse_cell_markers_df %>% 
#   filter(Tissue.of.Origin == 'Intestine') %>% select(Mouse.Gene) %>% pull(Mouse.Gene)
##  visualize all of the marker features simultaneously
# Provided by Johannes
markerGenesList <- c("Plag2g2a", "Defa-rs1", "Mmp7", "Lyz1", "Spdef", "Tcf7l2", "Ephb3", "Sis", "Ada", "Lct") #%>% as.data.frame()
# colnames(markerGenesList) <- c('Mouse.Gene')
# 
# Haber_marker_genes <- full_join(Haber_marker_genes, markerGenesList, by = join_by(Mouse.Gene)) %>% distinct()
# Haber_marker_genes <-  c(Haber_marker_genes, markerGenesList)
markerGenesList <- c("Plag2g2a", "Defa-rs1", "Mmp7", "Lyz1", "Spdef", "Tcf7l2", "Ephb3", "Sis", "Ada", "Lct")
(heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.05 & Log2FC >= 0.0",
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

mg <- c("Mmp7", "Lyz1", "Spdef", "Tcf7l2", "Ephb3", "Sis", "Ada", "Lct")
##'            Visualizing Marker Genes                               ##
(markerGeneEmbedding_Object <- plotEmbedding(ArchRProj = ATACSeq_project_All5, 
                                             colorBy = "GeneScoreMatrix", 
                                             name = mg, 
                                             embedding = "UMAP_all5", 
                                             quantCut = c(0.01, 0.95), 
                                             imputeWeights = NULL))

markerGenes_UMAP_Plot <- lapply(markerGeneEmbedding_Object, function(x){
  x + guides(colour = guide_colorbar(), size = guide_legend(),
             shape = guide_legend(), fill = FALSE) +
    theme_ArchR(baseSize = 6.5) 
  # theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  # theme(
  #   axis.text.x = element_blank(), 
  #   axis.ticks.x = element_blank(), 
  #   axis.text.y = element_blank(), 
  #   axis.ticks.y = element_blank()
  # )
})
do.call(cowplot::plot_grid, c(list(ncol = 2),markerGenes_UMAP_Plot))

plotPDF(plotList = markerGenes_UMAP_Plot, 
        name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = ATACSeq_project_All5, 
        addDOC = TRUE, width = 8, height = 8)


##'            Marker Genes Imputation with MAGIC                              ##
ATACSeq_project_All5 <- addImputeWeights(ATACSeq_project_All5,
                                         reducedDims = "Harmony_all5"
)

(ImputedmarkerGeneEmbedding_Object <- plotEmbedding(ArchRProj = ATACSeq_project_All5, 
                                                    colorBy = "GeneScoreMatrix", 
                                                    name = mg, 
                                                    embedding = "UMAP_all5", 
                                                    quantCut = c(0.01, 0.95), 
                                                    imputeWeights = getImputeWeights(ATACSeq_project_All5)
))

ImputedmarkerGenes_UMAP_Plot <- lapply(ImputedmarkerGeneEmbedding_Object, function(x){
  x + 
    # guides(colour = guide_colorbar(), size = guide_legend(), shape = guide_legend(), fill = FALSE) +
    theme_ArchR(baseSize = 6.5) 
  # theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  # theme(
  #   axis.text.x = element_blank(), 
  #   axis.ticks.x = element_blank(), 
  #   axis.text.y = element_blank(), 
  #   axis.ticks.y = element_blank()
  # )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),ImputedmarkerGenes_UMAP_Plot))

plotPDF(plotList = ImputedmarkerGenes_UMAP_Plot, 
        name = "Plot-UMAP-Marker-Genes-With-Imputation.pdf", 
        ArchRProj = ATACSeq_project_All5, 
        addDOC = TRUE, width = 8, height = 8)

library(ggpubr)

# (listplot <- ggarrange(plotlist = ImputedmarkerGenes_UMAP_Plot, ncol = 2, nrow = 3,
#                       common.legend = FALSE, align = "h"))
#   
# plotPDF(plotList = listplot, 
#         name = "Plot-UMAP-Marker-Genes-With-Imputation.pdf", 
#         ArchRProj = ATACSeq_project_All5, 
#         addDOC = TRUE, width = 8, height = 8)

## All genes UMAP with Imputed objects. 

AllGenesList <- as.character(cluster_df$name)
ATACSeq_project_All5 <- addImputeWeights(ATACSeq_project_All5,
                                         reducedDims = "Harmony_all5"
)

(ImputedmarkerGeneEmbedding_Object <- plotEmbedding(ArchRProj = ATACSeq_project_All5, 
                                                    colorBy = "GeneScoreMatrix", 
                                                    name = AllGenesList, #mg, #
                                                    embedding = "UMAP_all5", 
                                                    quantCut = c(0.01, 0.95), 
                                                    imputeWeights = getImputeWeights(ATACSeq_project_All5)
))

(ImputedmarkerGenes_UMAP_Plot <- lapply(ImputedmarkerGeneEmbedding_Object, function(x){
  x + 
    theme_ArchR(baseSize = 6, legendPosition = "right") +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    
  }))
# do.call(cowplot::plot_grid, c(list(ncol = 3),ImputedmarkerGenes_UMAP_Plot))

(listplot <- ggarrange(plotlist = ImputedmarkerGenes_UMAP_Plot, ncol = 2, nrow = 3,
                       common.legend = FALSE, align = "hv"))

plotPDF(plotList = listplot,
         name = "Plot-UMAP-Marker-Genes-With-Imputation.pdf",
         ArchRProj = ATACSeq_project_All5,
         addDOC = TRUE, width = 16, height = 16)

##'            Track Plotting                              ##

Trackplotter <- plotBrowserTrack(ArchRProj = ATACSeq_project_All5, groupBy = "Clusters_all5", 
                                 geneSymbol = mg, upstream = 50000, downstream = 50000)

plotPDF(plotList = Trackplotter, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = ATACSeq_project_All5, 
        addDOC = FALSE, width = 5, height = 5)

# saveplot(plot = Trackplotter, plotname = "Plot-Tracks-Marker-Genes")


