
# Script to Perform LSI, Batch Correction, Clustering, scRNAseq Integration, Gene Scores and Marker Genes Heatmaps, Imputation, Peak Calling

Package_List <- c("ArchR", "ComplexHeatmap", "writexl", "here", "patchwork", "tidyverse", "magick",
                  "pheatmap","Signac", "Seurat", "stringr", "BSgenome.Mmusculus.UCSC.mm10", 
                  "ggpubr", "chromVARmotifs")
not_installed <- Package_List[!(Package_List %in% installed.packages()[, "Package"])] 
# Extract not installed packages
if (length(not_installed)) install.packages(not_installed) # Install the uninstalled packages
invisible(lapply(Package_List, suppressMessages(library), character.only = TRUE))

# Set the number of threads
addArchRThreads(threads = 15)

# File Path Declarations
here::i_am(path = "scATACseqAnalysis.R")
paste0(here())

# set seed
set.seed(1)

# Set Data Paths
PreProcessed_Data_Path = file.path("/media", "horneflablinux", "HornefLab_Data3", "scATACseq", "All5/")
scRNAseq_Johannes_Reference_Path = file.path("/media", "horneflablinux", "HornefLab_Data3",   
                                             "scRNA_AnnotationData_Johannes", "scrna_with_day25.Rds")

# Load the PreProcessed ATACseq Project
ATACSeq_project_All5 <- loadArchRProject(path = PreProcessed_Data_Path, force = TRUE, showLogo = TRUE)

# Set the Aesthetics
my_gg_theme <- theme(axis.title = element_text(size = 15),
                     axis.text = element_text(size = 15),
                     legend.text = element_text(size = 10))
print(getOutputDirectory(ATACSeq_project_All5))

# Function to save generic plots in jpeg format in jpeg folder.
saveplot <- function(plot, plotname) {
  extension <- ".jpeg"
  filename <- paste(plotname, extension, sep = "")
  ggsave(
    filename = file.path(getOutputDirectory(ATACSeq_project_All5), "jpeg", filename),
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

#' Dimensional Reduction, INTEGRATION, CLUSTERING & EMBEDDING
# 4.1 ## Reducing Dims via Iterative LSI. #The most common parameters to tweak are iterations, varFeatures, and resolution

ATACSeq_project_All5 <- ATACSeq_project_All5 %>%
  addIterativeLSI(
    useMatrix = "TileMatrix",
    name = "IterativeLSI_all5",
    iterations = 5,# no. of iterations must be 1 less than length of list of resolution numbers for clusterParams below. Else it will fail!
    depthCol = "nFrags",
    corCutOff = 0.75, # A numeric cutoff for the correlation of each dimension to the sequencing depth. If the dimension has a correlation to sequencing depth that is greater than the corCutOff, it will be excluded from analysis.
    firstSelection = "Top", # for scATACseq data
    clusterParams = list(resolution = c(04, 03, 02, 01, 0.5, 0.25, 0.125, 0.05), 
                         # Bigger the value of resolutions, the clusters will be more pulled apart. Smaller the value, more closer they are!
                         sampleCells = 10000, 
                         n.start = 1),
    varFeatures = 15000,
    dimsToUse = 1:30,
    # seed = 1,
    scaleDims = FALSE, # A boolean that indicates whether to z-score the reduced dimensions for each cell. This is useful for minimizing the contribution of strong biases (dominating early PCs) and lowly abundant populations. However, this may lead to stronger sample-specific biases since it is over-weighting latent PCs. If set to NULL this will scale the dimensions based on the value of scaleDims when the reducedDims were originally created during dimensionality reduction.
    filterBias = TRUE,# bias clusters often refer to clusters of cells that are formed due to technical artifacts or biases rather than biological variation. These biases could be introduced during the experimental process, such as during cell capture, library preparation, or sequencing. They could also arise from the computational processing of the data.
    saveIterations = TRUE,
    force = TRUE, 
    verbose = TRUE
  ) %>%
  addHarmony( ## Using LSI object as Input to correct for Batch Effects via Harmony with default values
    reducedDims = "IterativeLSI_all5",
    name = "Harmony_all5",
    groupBy = "Sample",
    corCutOff = 0.75,
    dimsToUse = 1:30,
    verbose = TRUE,
    force = TRUE
  ) %>%
  addUMAP(# 5.1 ## Visualising clusters via UMAP
    reducedDims = "Harmony_all5", # Need to use the Harmony object for clustering and UMAP for proper integration
    name = "UMAP_all5",
    nNeighbors = 30,
    minDist = 0.4,
    dimsToUse = 1:30,
    corCutOff = 0.75,
    metric = "cosine",
    force = TRUE,
    verbose = TRUE,
    saveModel = TRUE,
    # seed = 1
  )

## Create UMAP Embedding for all samples together
(SampleWise_UMAP <- plotEmbedding(
  ArchRProj = ATACSeq_project_All5,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_all5",
  plotAs = "points",
  size = 0.05,
  keepAxis = TRUE
) + my_gg_theme)

########## Compute CorCutOff Values####################
CorCutOffValues <- ATACSeq_project_All5@reducedDims@listData[["IterativeLSI_all5"]]@listData[["corToDepth"]][["none"]] %>% as.data.frame() %>% tibble::rownames_to_column(., "LSI_Dim") %>% rename("Cor_Values" = ".")
print(head(CorCutOffValues, 10))
###############################################

plotPDF(SampleWise_UMAP, 
        name = "UMAP-Samplewise_All5.pdf", 
        ArchRProj = ATACSeq_project_All5,
        width = 10, height = 10, addDOC = FALSE)
saveplot(plot = SampleWise_UMAP, plotname = "UMAP-Samplewise_All5")

## Create SplitView UMAP
create_splitUMAP <- function(ATACSeq_project_All5, sampleObject) {# Create a Split view UMAP so that the UMAP for each sample can be visualized separately.
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
d1_UMAP <- create_splitUMAP(ATACSeq_project_All5, samplelist@values[1]) + my_gg_theme
d5_UMAP <- create_splitUMAP(ATACSeq_project_All5, samplelist@values[5]) + my_gg_theme
d10_UMAP <- create_splitUMAP(ATACSeq_project_All5, samplelist@values[4]) + my_gg_theme
d25_UMAP <- create_splitUMAP(ATACSeq_project_All5, samplelist@values[3]) + my_gg_theme
d4pi_UMAP <- create_splitUMAP(ATACSeq_project_All5, samplelist@values[2]) + my_gg_theme

plotPDF(d1_UMAP, d5_UMAP, d10_UMAP, d25_UMAP, d4pi_UMAP,
        name = "SplitviewUMAP_SampleWise_All5.pdf",
        ArchRProj = ATACSeq_project_All5, 
        width = 10,
        height = 10, addDOC = FALSE
)
saveplot(plot = d1_UMAP, plotname = "d1_UMAP")
saveplot(plot = d5_UMAP, plotname = "d5_UMAP")
saveplot(plot = d10_UMAP, plotname = "d10_UMAP")
saveplot(plot = d25_UMAP, plotname = "d25_UMAP")
saveplot(plot = d4pi_UMAP, plotname = "d4pi_UMAP")

#' Adding Clusters
ATACSeq_project_All5 <- ATACSeq_project_All5 %>% 
  addClusters( # 5.1 ## Creating clusters based on Iterative LSI object      ##
    reducedDims = "Harmony_all5", # Harmony_all5
    method = "Seurat",
    name = "Clusters_all5",
    resolution = c(0.2, 0.5, 0.8),
    # seed = 1,
    corCutOff = 0.75,
    force = TRUE
  ) %>% 
  addUMAP( # 5.1 ## Visualising clusters via UMAP
    # ArchRProj = ATACSeq_project_All5,
    reducedDims = "Harmony_all5", # Need to used the Harmony object for clustering and UMAP to ensure proper integration
    name = "ClustersUMAP_all5",
    nNeighbors = 30,
    minDist = 0.4,
    metric = "cosine",
    corCutOff = 0.75,
    # seed = 1,
    dimsToUse = 1:30,
    force = TRUE,
    verbose = TRUE,
    saveModel = TRUE
  )

ClusterWise_UMAP <- plotEmbedding(ArchRProj = ATACSeq_project_All5, 
                                  colorBy = "cellColData", 
                                  name = "Clusters_all5", 
                                  embedding = "ClustersUMAP_all5", 
                                  plotAs = "points",
                                  size = 0.05,
                                  keepAxis = TRUE) + my_gg_theme
plotPDF(ClusterWise_UMAP, name = "ClustersUMAP_All5.pdf", ArchRProj = ATACSeq_project_All5,  
        width = 10, height = 10, addDOC = FALSE)
saveplot(plot = ClusterWise_UMAP, plotname = "ClustersUMAP_All5")

# # Accessing the clusters
head(ATACSeq_project_All5$Clusters_all5)
# tabulate the number of cells present in each cluster:
table(ATACSeq_project_All5$Clusters_all5)

# To better understand which samples reside in which clusters, we can create a
# cluster confusion matrix across each sample using the confusionMatrix() function.
cM <- confusionMatrix(paste0(ATACSeq_project_All5$Clusters_all5),
                      paste0(ATACSeq_project_All5$Sample))
print(cM)

# plotting the Confusion Matrix
cM <- cM / Matrix::rowSums(cM)
cM_heatmap <- pheatmap::pheatmap(mat = as.matrix(cM), color = paletteContinuous("whiteBlue"), 
                                 border_color = "black") + my_gg_theme

## –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##
#'          SECTION 7: Annotating Cell types with a Reference Dataset       ##
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

##           7.1 Unconstrained Integration       ##
# Read-in the reference
scRNA_AnnotationData_Johannes <- readRDS(scRNAseq_Johannes_Reference_Path)

# add gene integration matrix
ATACSeq_project_All5 <- ATACSeq_project_All5 %>% 
  addGeneIntegrationMatrix(
    # ArchRProj = ATACSeq_project_All5,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI_all5",
    seRNA = scRNA_AnnotationData_Johannes,
    addToArrow = TRUE,
    groupRNA = "int_0.3_broad_tuft",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    dimsToUse = 1:30,
    corCutOff = 0.75,
    plotUMAP = TRUE,
    force = TRUE
  )

# Plot UMAP with predicted cell types
annotated_UMAP_scRNAseqRef <- plotEmbedding(ATACSeq_project_All5,
                                            name = "predictedGroup_Un",
                                            embedding = "UMAP_all5",
                                            baseSize = 15,
                                            plotAs = "points",
                                            size = 0.05,
                                            labelAsFactors = F,
                                            labelMeans = F)

plotPDF(annotated_UMAP_scRNAseqRef, 
        name = "CellTypeAnnotatedUMAP_scRNAseqRef_All5",
        width = 10, height = 10, 
        ArchRProj = ATACSeq_project_All5, 
        addDOC = FALSE)
saveplot(plot = annotated_UMAP_scRNAseqRef, plotname = "Annotated_UMAP_scRNAseqRef")

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

# Call the function to create UMAP for each Major CellType
Enterocyte_UMAP <- create_SampleWise_splitAnnotatedUMAP(ATACSeq_project_All5, "Enterocyte") + my_gg_theme + ggtitle("Enterocyte_UMAP")
GP_UMAP <- create_SampleWise_splitAnnotatedUMAP(ATACSeq_project_All5, "Goblet+Paneth") + my_gg_theme + ggtitle("Goblet+Paneth_UMAP")
EEC_UMAP <- create_SampleWise_splitAnnotatedUMAP(ATACSeq_project_All5, "EEC") + my_gg_theme + ggtitle("EEC_UMAP")
Stem_UMAP <- create_SampleWise_splitAnnotatedUMAP(ATACSeq_project_All5, "Stem") + my_gg_theme + ggtitle("Stem_UMAP")
Tuft_UMAP <- create_SampleWise_splitAnnotatedUMAP(ATACSeq_project_All5, "Tuft") + my_gg_theme + ggtitle("Tuft_UMAP")

plotPDF(Enterocyte_UMAP, GP_UMAP, EEC_UMAP, Stem_UMAP, Tuft_UMAP,
        name = "SampleWise_Splitview-CellTypeAnnotatedUMAP.pdf",
        ArchRProj = ATACSeq_project_All5, 
        width = 10,
        height = 10,
        addDOC = FALSE
)
saveplot(plot = Enterocyte_UMAP, plotname = "Enterocyte_UMAP")
saveplot(plot = GP_UMAP, plotname = "GP_UMAP")
saveplot(plot = EEC_UMAP, plotname = "EEC_UMAP")
saveplot(plot = Stem_UMAP, plotname = "Stem_UMAP")
saveplot(plot = Tuft_UMAP, plotname = "Tuft_UMAP")

## Constrained Integration

# Now that we have our preliminary unconstrained integration, we will identify general cell types to profide a framework to further refine the integration results.

# we would ideally constrain the integration to associate similar cell types togther. First, we will identify which cell types from the scRNA-seq data are most abundant in each of our scATAC-seq clusters. 
cM <- as.matrix(confusionMatrix(ATACSeq_project_All5$Clusters_all5, ATACSeq_project_All5$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments. This list shows which scRNA-seq cell type is most abundant in each of the 12 scATAC-seq clusters.

# lets look at the cell type labels from our scRNA-seq data that were used in our unconstrained integration:
unique(unique(ATACSeq_project_All5$predictedGroup_Un))


# with the results of our scATAC-seq and scRNA-seq integration, we can re-run the integration with 
# addToArrow = TRUE to add the linked gene expression data to each of the Arrow files. 

# The other key parameters for this function provide column names in cellColData where certain information will be stored: nameCell will store the cell ID from the matched scRNA-seq cell, nameGroup will store the group ID from the scRNA-seq cell, and nameScore will store the cross-platform integration score.
# ATACSeq_project_All5 <- ATACSeq_project_All5 %>%
#   addGeneIntegrationMatrix(
#     # ArchRProj = ATACSeq_project_All5,
#     useMatrix = "GeneScoreMatrix",
#     matrixName = "GeneIntegrationMatrix",
#     reducedDims = "IterativeLSI_all5",
#     seRNA = scRNA_AnnotationData_Johannes,
#     addToArrow = TRUE,
#     groupRNA = "int_0.3_broad_tuft",
#     nameCell = "predictedCell_Un",
#     nameGroup = "predictedGroup_Un",
#     nameScore = "predictedScore_Un",
#     dimsToUse = 1:30,
#     corCutOff = 0.75,
#     plotUMAP = TRUE,
#     force = TRUE
#   )

## Labeling scATAC-seq clusters with scRNA-seq information
cM <- confusionMatrix(ATACSeq_project_All5$Clusters_all5, ATACSeq_project_All5$predictedGroup_Un)
labelOld <- rownames(cM)
print(labelOld)
# for each of our scATAC-seq clusters, we identify the cell type from predictedGroup which best defines that cluster
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
print(labelNew)

# Next we need to reclassify these new cluster labels to make a simpler categorization system. For each scRNA-seq cluster
remapClust <- c(
  "Enterocyte" = "Aline/Enterocyte",
  "Tuft" = "Kaiyi/Tuft",
  "Stem" = "Matthias/Stem",
  "Goblet+Paneth" = "Fabian/Goblet+Paneth",
  "EEC" = "Annalena/EEC"
)
remapClust <- remapClust[names(remapClust) %in% labelNew]
labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
labelNew2
# use the mapLabels() function again to create new cluster labels in cellColData.
ATACSeq_project_All5$Clusters2 <- mapLabels(ATACSeq_project_All5$Clusters_all5, 
                                            newLabels = labelNew2, 
                                            oldLabels = labelOld)
Renamed_UMAP <- plotEmbedding(ATACSeq_project_All5, 
                              colorBy = "cellColData", 
                              embedding = "ClustersUMAP_all5",
                              name = "Clusters2")
plotPDF(Renamed_UMAP, name = "Renamed_UMAP-Clusters.pdf", ArchRProj = ATACSeq_project_All5, 
        addDOC = FALSE, width = 8, height = 8)
saveplot(plot = Renamed_UMAP, plotname = "Renamed_UMAP")
##############################################################################################
#'                     SECTION7 : Gene Scores and Marker Genes
## ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––##

#' ##'            Identifying Marker Genes                               ##
# Identifying Marker Genes    
markersGS <- getMarkerFeatures(ArchRProj = ATACSeq_project_All5,
                               useMatrix = "GeneScoreMatrix",
                               groupBy = "Clusters_all5",
                               bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon"
)

# Markers List
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
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

## Marker Genes

# Human_and_mouse_cell_markers_df <- as.data.frame(read.delim("D:/scATACseq/Human_and_mouse_cell_markers-Markers.tsv"))
# 
# Haber_marker_genes <- Human_and_mouse_cell_markers_df %>% 
#   filter(Tissue.of.Origin == 'Intestine') %>% select(Mouse.Gene) %>% pull(Mouse.Gene)
##  visualize all of the marker features simultaneously
##  
# Provided by Johannes
markerGenesList <- c("Plag2g2a", "Defa-rs1", "Mmp7", "Lyz1", "Spdef", "Tcf7l2", "Ephb3", "Sis", "Ada", "Lct") #%>% as.data.frame()
# colnames(markerGenesList) <- c('Mouse.Gene')
# 
# Haber_marker_genes <- full_join(Haber_marker_genes, markerGenesList, by = join_by(Mouse.Gene)) %>% distinct()
# Haber_marker_genes <-  c(Haber_marker_genes, markerGenesList)
# markerGenesList <- c("Plag2g2a", "Defa-rs1", "Mmp7", "Lyz1", "Spdef", "Tcf7l2", "Ephb3", "Sis", "Ada", "Lct")
# 
# ### Create Heatmap for Marker Genes
heatmap_MarkerGene <- plotMarkerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.05 & Log2FC >= 1.0",
  limits = c(-3, 3),
  returnMatrix = FALSE,
  plotLog2FC = TRUE,
  labelMarkers = markerGenesList,
  transpose = TRUE,
  labelRows = TRUE, clusterCols = TRUE, nPrint = 10)

plotPDF(heatmap_MarkerGene,
        name = "GeneScores-Marker-Heatmap",
        width = 10, height = 10,
        ArchRProj = ATACSeq_project_All5, addDOC = FALSE)

mg <- c("Plag2g2a", "Defa-rs1", "Mmp7", "Lyz1", "Spdef", "Tcf7l2", "Ephb3", "Sis", "Ada", "Lct")

# #            Visualizing Marker Genes 
markerGeneEmbedding_Object <- plotEmbedding(ArchRProj = ATACSeq_project_All5,
                                            colorBy = "GeneScoreMatrix",
                                            name = mg,
                                            embedding = "UMAP_all5",
                                            quantCut = c(0.01, 0.95),
                                            imputeWeights = NULL)

markerGenes_UMAP_Plot <- lapply(markerGeneEmbedding_Object, function(x){
  x + guides(colour = guide_colorbar(), size = guide_legend(),
             shape = guide_legend(), fill = FALSE) +
    theme_ArchR(baseSize = 6.5)+ 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
})

plotPDF(plotList = markerGenes_UMAP_Plot,
        name = "UMAP-Marker-Genes-Without-Imputation.pdf",
        ArchRProj = ATACSeq_project_All5,
        addDOC = FALSE, width = 8, height = 8)

# ### Marker Genes with Imputation with MAGIC
ATACSeq_project_All5 <- addImputeWeights(ATACSeq_project_All5,
                                         reducedDims = "Harmony_all5")

ImputedmarkerGeneEmbedding_Object <- plotEmbedding(ArchRProj = ATACSeq_project_All5,
                                                   colorBy = "GeneScoreMatrix",
                                                   name = mg,
                                                   embedding = "UMAP_all5",
                                                   quantCut = c(0.01, 0.95),
                                                   imputeWeights = getImputeWeights(ATACSeq_project_All5)
)

ImputedmarkerGenes_UMAP_Plot <- lapply(ImputedmarkerGeneEmbedding_Object, function(x){
  x + guides(colour = guide_colorbar(), fill = FALSE) +
    theme_ArchR(baseSize = 6.5)+ 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
})

ImputedmarkerGenes_UMAP_listplot_withImputation <- ggarrange(plotlist = ImputedmarkerGenes_UMAP_Plot,
                                                             ncol = 2, nrow = 3, common.legend = FALSE, align = "h")

plotPDF(plotList = ImputedmarkerGenes_UMAP_listplot_withImputation,
        name = "Marker-Genes-UMAP-With-Imputation.pdf",
        ArchRProj = ATACSeq_project_All5,
        addDOC = FALSE, width = 8, height = 8)
## All genes UMAP with Imputed objects.

AllGenesList <- as.character(cluster_df$name)
ATACSeq_project_All5 <- addImputeWeights(ATACSeq_project_All5, reducedDims = "Harmony_all5")

ImputedmarkerGeneEmbedding_Object <- plotEmbedding(ArchRProj = ATACSeq_project_All5,
                                                   colorBy = "GeneScoreMatrix", name = AllGenesList, 
                                                   embedding = "UMAP_all5", 
                                                   quantCut = c(0.01, 0.95), 
                                                   imputeWeights = getImputeWeights(ATACSeq_project_All5))

ImputedmarkerGenes_UMAP_Plot <- lapply(ImputedmarkerGeneEmbedding_Object, function(x){
  x + theme_ArchR(baseSize = 6, legendPosition = "right") })
ImputedmarkerGenes_UMAP_Plot_listplot <- ggarrange(plotlist = ImputedmarkerGenes_UMAP_Plot,
                                                   ncol = 3, nrow = 3,
                                                   common.legend = FALSE,
                                                   align = "hv")

plotPDF(plotList = ImputedmarkerGenes_UMAP_Plot_listplot,
        name = "ListPlot-Marker-Genes-UMAP-With-Imputation.pdf",
        ArchRProj = ATACSeq_project_All5,
        addDOC = FALSE, width = 10, height = 10)

# Track Plotting
Trackplotter = plotBrowserTrack(ArchRProj = ATACSeq_project_All5,
                                groupBy = "Clusters_all5",
                                geneSymbol = mg,
                                upstream = 5000, downstream = 5000)
plotPDF(plotList = Trackplotter,
        name = "TracksPlot-Marker-Genes.pdf",  
        ArchRProj = ATACSeq_project_All5,
        addDOC = FALSE, width = 8, height = 8)

# PseudoBulk Replicates in ArchR

## Defining pseudo-bulk replicates
# First we have to define pseudo-bulk replicates to call peaks on them, ArchR merges cells within each designated cell group. The key parameter here is groupBy which defines the groups for which pseudo-bulk replicates should be made.
ATACSeq_project_All5 <- addGroupCoverages(ArchRProj = ATACSeq_project_All5, 
                                          groupBy = "Clusters_all5", force = TRUE)

# Calling Peaks

# One of the Most important bits of scATACseq Analysis. Because per-cell scATAC-seq data is essentially binary (accessible or not accessible), we cannot call peaks on an individual cell basis. For this reason, we defined groups of cells, typically clusters. Moreover, we created pseudo-bulk replicates to allow us to assess the reproducibility of our peak calls.

## Calling Peaks with MACS2 via Peak Matrix
pathToMacs2 <- findMacs2()#
ATACSeq_project_All5 <- addReproduciblePeakSet(ArchRProj = ATACSeq_project_All5,
                                               groupBy = "Clusters_all5",
                                               peaksPerCell = 500, # To avoid bias from pseudo-bulk replicates with very few cells, provide a cutoff for the upper limit of the number of peaks called per cell via the peaksPerCell parameter.
                                               pathToMacs2 = pathToMacs2)
PeakSet <- getPeakSet(ATACSeq_project_All5)
ATACSeq_project_All5 <- addPeakMatrix(ATACSeq_project_All5)

## Calling Peaks with ArchR Native Peak Caller via TileMatrix.
# ATACSeq_project_All5 <- addReproduciblePeakSet(
#   ArchRProj = ATACSeq_project_All5, 
#   groupBy = "Clusters_all5",
#   peakMethod = "Tiles",
#   method = "p"
# )
# getPeakSet(ATACSeq_project_All5)

# # save Project after Peak Calling
saveArchRProject(
  ArchRProj = ATACSeq_project_All5,
  outputDirectory = "/home/horneflablinux/postpeakcalling/",
  overwrite = TRUE,
  load = TRUE,
  logFile = createLogFile("saveArchRProject_All5_postpeakcalling")
)
