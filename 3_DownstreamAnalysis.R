#' SECTION1: Initialization

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

# Set Data Paths /home/horneflablinux//
PostPeakCalling_Path = file.path("/home", "horneflablinux", "postpeakcalling/")
scRNAseq_Johannes_Reference_Path = file.path("/media", "horneflablinux", "HornefLab_Data3",   
                                             "scRNA_AnnotationData_Johannes", "scrna_with_day25.Rds")

# Load the PreProcessed ATACseq Project
ATACSeq_project_All5 <- loadArchRProject(path = PostPeakCalling_Path, force = TRUE, showLogo = TRUE)

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

mg <- c("Plag2g2a", "Defa-rs1", "Mmp7", "Lyz1", "Spdef", "Tcf7l2", "Ephb3", "Sis", "Ada", "Lct")

#  Identifying Marker Peaks with ArchR

# Marker features are unique to a specific cell grouping. These can be very useful in understanding cluster- or cell type-specific biology.
# we are interested to know which peaks are unique to an individual cluster or a small group of clusters. We can do this in an unsupervised fashion in ArchR using the addMarkerFeatures() function in combination with useMatrix = "PeakMatrix". 
markersPeaks <- getMarkerFeatures(ArchRProj = ATACSeq_project_All5, 
                                  useMatrix = "PeakMatrix", 
                                  groupBy = "Clusters_all5", 
                                  testMethod = "wilcoxon",
                                  bias = c("TSSEnrichment", "log10(nFrags)")# to account for differences in data quality amongst the cell groups by setting the bias parameter to account for TSS enrichment and the number of unique fragments per cell.
)
print(markersPeaks)
## retrieve particular slices of this SummarizedExperiment that we are interested in. The default behavior of this function is to return a list of DataFrame objects, one for each cell group.
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)
print(markerList)

# Plotting Marker Peaks 
heatmap_MarkerPeaks <- plotMarkerHeatmap(seMarker = markersPeaks,
                                         cutOff = "FDR <= 0.1 & Log2FC >= 0.5", 
                                         plotLog2FC = FALSE,
                                         limits = c(-2, 2),
                                         nLabel = 15,
                                         nPrint = 15,
                                         transpose = TRUE)

# ComplexHeatmap::draw(heatmap_MarkerPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmap_MarkerPeaks, name = "MarkerPeaks-Heatmap", width = 8, height = 8,
        ArchRProj = ATACSeq_project_All5, addDOC = FALSE)

## Volcano Plots for Marker Peaks
MarkerPeak_VolcanoPlot <- markerPlot(seMarker = markersPeaks, name = "C1", 
                                     cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
plotPDF(MarkerPeak_VolcanoPlot, name = "MarkerPeaks-VolcanoPlot", width = 8, height = 8, 
        ArchRProj = ATACSeq_project_All5, addDOC = FALSE)
saveplot(plot = MarkerPeak_VolcanoPlot, plotname = "MarkerPeak_VolcanoPlot")

## Marker Peaks in Browser Tracks
# visualise peak regions overlayed on our browser tracks by passing the relevant peak regions to the features parameterin the plotBrowserTrack() function. This will add an additional BED-style track of marker peak regions to the bottom of our ArchR track plot.

MarkerPeak_TrackPlot <- plotBrowserTrack(ArchRProj = ATACSeq_project_All5, 
                                         groupBy = "Clusters2", 
                                         geneSymbol = mg,
                                         features =  getMarkers(markersPeaks, 
                                                                cutOff = "FDR <= 0.1 & Log2FC >= 1", 
                                                                returnGR = TRUE)["C1"],
                                         upstream = 5000, 
                                         downstream = 5000)
plotPDF(MarkerPeak_TrackPlot, name = "MarkerPeak_TrackPlot-Tracks-With-Features", 
        width = 8, height = 8, ArchRProj = ATACSeq_project_All5, addDOC = FALSE)

# Pairwise Testing Between Groups
markerTest <- getMarkerFeatures(ArchRProj = ATACSeq_project_All5, useMatrix = "PeakMatrix",
                                groupBy = "Clusters2", testMethod = "wilcoxon", 
                                bias = c("TSSEnrichment", "log10(nFrags)"), 
                                useGroups = "Matthias/Stem", bgdGroups = "Aline/Enterocyte")
markerTest_volcanoPlot <- markerPlot(seMarker = markerTest, name = "Matthias/Stem", 
                                     cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
plotPDF(markerTest_volcanoPlot, name = "Matthias_Stem-vs-Aline_Enterocyte-Markers-VolcanoPlot", 
        width = 8, height = 8, ArchRProj = ATACSeq_project_All5, addDOC = FALSE)
saveplot(plot = markerTest_volcanoPlot, plotname = "PairWisemarkerTest_volcanoPlot")

# Motif and Feature Enrichment

# to predict what transcription factors may be mediating the binding events that create those accessible chromatin sites. This can be helpful in assessing marker peaks or differential peaks to understand if these groups of peaks are enriched for binding sites of specific transcription factors. 

# look for motifs that are enriched in peaks that are up or down in various cell types. To do this, we must first add these motif annotations to our ArchRProject. This effectively creates a binary matrix where the presence of a motif in each peak is indicated numerically.
ATACSeq_project_All5 <- addMotifAnnotations(ArchRProj = ATACSeq_project_All5, motifSet = "cisbp", name = "Motif")

## Motif Enrichment in Differential Peaks
motifsUp <- peakAnnoEnrichment(seMarker = markerTest, # use the differential testing SummarizedExperiment object markerTest which was generated in the previous chapter to define the set of significantly differential peaks that we are interested in testing for motif enrichment.
                               ArchRProj = ATACSeq_project_All5, 
                               peakAnnotation = "Motif", 
                               cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
print(motifsUp) # is a SummarizedExperiment object containing multiple assays that store the results of enrichment testing with the hypergeometric test.

# Prepare this data for plotting with ggplot we can create a simplified data.frame object containing the motif names, the corrected p-values, and the significance rank.
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)

# plot the rank-sorted TF motifs and color them by the significance of their enrichment. 
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

plotPDF(ggUp, name = "ggUp-Motifs-Enriched", width = 8, height = 8, ArchRProj = ATACSeq_project_All5,
        addDOC = FALSE)
saveplot(plot = ggUp, plotname = "ggUp-Motifs-Enriched")

# We can perform the same analyses for the peaks that are more accessible in the “Progenitor” cells by using peaks with Log2FC <= -0.5.

motifsDo <- peakAnnoEnrichment(seMarker = markerTest, ArchRProj = ATACSeq_project_All5,
                               peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC <= -0.5")
print(motifsDo)

df_DO <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df_DO <- df[order(df_DO$mlog10Padj, decreasing = TRUE),]
df_DO$rank <- seq_len(nrow(df_DO))
head(df_DO)

ggDo <- ggplot(df_DO, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

plotPDF(ggDo, name = "ggDo-Motifs-Enriched", width = 8, height = 8, ArchRProj = ATACSeq_project_All5,
        addDOC = FALSE)
saveplot(plot = ggDo, plotname = "ggDo-Motifs-Enriched")

## Motif Enrichment in Marker Peaks
enrichMotifs_MarkerPeaks <- peakAnnoEnrichment(seMarker = markersPeaks, ArchRProj = ATACSeq_project_All5,
                                               peakAnnotation = "Motif", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
print(enrichMotifs_MarkerPeaks)

# Plot this enrich Motifs directly as a Heatmap
heatmap_enrichMotifs_MarkerPeaks <- plotEnrichHeatmap(enrichMotifs_MarkerPeaks, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmap_enrichMotifs_MarkerPeaks, heatmap_legend_side = "bot", 
                     annotation_legend_side = "bot")

plotPDF(heatmap_enrichMotifs_MarkerPeaks, name = "EnrichMotifs_MarkerPeaks-Heatmap", 
        width = 8, height = 8, ArchRProj = ATACSeq_project_All5, addDOC = FALSE)

# ChromVar Deviations Enrichment
# TF motif enrichments can help us predict which regulatory factors are most active in our cell type of interest. These enrichments, however, are not calculated on a per-cell basis and they do not take into account the insertion sequence bias of the Tn5 transposase.

# ChromVAR is designed for predicting enrichment of TF activity on a per-cell basis from sparse chromatin accessibility data. The two primary outputs of chromVAR are: "Deviations" and "Zscore/DeviationScore".

# check if the motifs are added to teh respective ArchR Project
if("Motif" %ni% names(ATACSeq_project_All5@peakAnnotation)){
  ATACSeq_project_All5 <- addMotifAnnotations(ArchRProj = ATACSeq_project_All5, motifSet = "cisbp", name = "Motif")}

# add a set of background peaks which are used in computing deviations. Background peaks are chosen using the chromVAR::getBackgroundPeaks() function which samples peaks based on similarity in GC-content and number of fragments across all samples using the Mahalanobis distance.
ATACSeq_project_All5 <- addBgdPeaks(ATACSeq_project_All5)

# now ready to compute per-cell deviations accross all of our motif annotations using the addDeviationsMatrix() function. 
ATACSeq_project_All5 <- addDeviationsMatrix(ArchRProj = ATACSeq_project_All5, 
                                            peakAnnotation = "Motif", force = TRUE)

# To access these deviations, we use the getVarDeviations() function. If we want this function to return a ggplot object, we set plot = TRUE. 
plotVarDev <- getVarDeviations(ATACSeq_project_All5, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 8, height = 8, 
        ArchRProj = ATACSeq_project_All5, addDOC = FALSE)
saveplot(plot = plotVarDev, "Variable-Motif-Deviation-Scores")

# Extract Subset of Motives
motifs <- mg 
markerMotifs <- getFeatures(ATACSeq_project_All5,
                            select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
print(markerMotifs)

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
# markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs

###' plot the distribution of chromVAR deviation scores for each cluster.
ChromVarDeviations_plot <- plotGroups(ArchRProj = ATACSeq_project_All5, 
                                      groupBy = "Clusters2", colorBy = "MotifMatrix", 
                                      name = markerMotifs, 
                                      imputeWeights = getImputeWeights(ATACSeq_project_All5))

# use cowplot to plot the distributions of all of these motifs in a single plot.
ChromVarDeviations_Allplots <- lapply(seq_along(ChromVarDeviations_plot), function(x){
  if(x != 1){
    ChromVarDeviations_plot[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    ChromVarDeviations_plot[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})

plotPDF(ChromVarDeviations_plot, 
        name = "ChromVar-Groups-Deviations-w-Imputation", width = 8, height = 8, 
        ArchRProj = ATACSeq_project_All5, addDOC = FALSE)

# Instead of looking at the distributions of these z-scores, we can overlay the z-scores on our UMAP embedding as we’ve done previously for gene scores.

CV_on_GeneScore <- plotEmbedding(ArchRProj = ATACSeq_project_All5, colorBy = "MotifMatrix", 
                                 name = sort(markerMotifs), embedding = "UMAP_all5",
                                 imputeWeights = getImputeWeights(ATACSeq_project_All5))
## Motif UMAPs

Motif_UMAPs_GeneScore <- lapply(CV_on_GeneScore, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})

# To see how these TF deviation z-scores compare to the inferred gene expression via gene scores of the corresponding TF genes, we can overlay the gene scores for each of these TFs on the UMAP embedding.
markerRNA_TFDeviation <- getFeatures(ATACSeq_project_All5, select = paste(motifs, collapse="|"), 
                                     useMatrix = "GeneScoreMatrix")
#markerRNA_TFDeviation <- markerRNA_TFDeviation#[markerRNA_TFDeviation %ni% c("Spdef","Ephb3", "Ada", "Lct", "Sis", "Lyz1")]
markerRNA_TFDeviation_UMAP <- plotEmbedding(ArchRProj = ATACSeq_project_All5, colorBy = "GeneScoreMatrix",  
                                            name = sort(markerRNA_TFDeviation), embedding = "UMAP_all5", 
                                            imputeWeights = getImputeWeights(ATACSeq_project_All5))

markerRNA_TFDeviation_UMAP_All <- lapply(markerRNA_TFDeviation_UMAP, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

markerRNA_TFDeviation_UMAP_All_listplot <- ggarrange(plotlist = markerRNA_TFDeviation_UMAP_All, 
                                                     ncol = 2, nrow = 2, common.legend = TRUE, 
                                                     align = "hv")

plotPDF(markerRNA_TFDeviation_UMAP_All_listplot, 
        name = "ChromVar-Groups-Deviations-w-Imputation", width = 8, height = 8, 
        ArchRProj = ATACSeq_project_All5, addDOC = FALSE)

###################13.1 ArchR incomplete#################################

########MOTIF FOOTPRINTING##################

# obtain the positions of the relevant motifs.
motifpositions <- getPositions(ATACSeq_project_All5) # This creates a GRangesList object where each TF motif is represented by a separate GRanges object.
print(motifpositions)
# subset this GRangesList to a few TF motifs that we are interested in. To accurately profile TF footprints, a large number of reads are required. Therefore, cells are grouped to create pseudo-bulk ATAC-seq profiles that can be then used for TF footprinting. These pseudo-bulk profiles are stored as group coverage files which we originally created in a previous chapter to perform peak calling.

ATACSeq_project_All5 <- addGroupCoverages(ArchRProj = ATACSeq_project_All5, groupBy = "Clusters2")

# now compute footprints for the subset of marker motifs that we previously selected using the getFootprints() function
seFoot <- getFootprints(ArchRProj = ATACSeq_project_All5, positions = motifpositions, 
                        groupBy = "Clusters2")

# Once we have retrieved these footprints, we can plot them using the plotFootprints() function. This function can simultaneously normalize the footprints in various ways. 

# Normalization of Footprints for Tn5 Bias
# One major challenge with TF footprinting using ATAC-seq data is the insertion sequence bias of the Tn5 transposase which can lead to misclassification of TF footprints. To account for Tn5 insertion bias, ArchR identifies the k-mer (user-defined length, default length 6) sequences surrounding each Tn5 insertion site.


# Strategy1: Subtracting the Tn5 Bias
TF_footprintg_S1 <- plotFootprints(
  seFoot = seFoot,
  ArchRProj = ATACSeq_project_All5, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

# Strategy2: Dividing by the Tn5 Bias
TF_footprintg_S2 <-plotFootprints(
  seFoot = seFoot,
  ArchRProj = ATACSeq_project_All5, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

# Strategy3: Footprinting Without Normalization for Tn5 Bias
# TF_footprintg_S3 <-plotFootprints(
#   seFoot = seFoot,
#   ArchRProj = ATACSeq_project_All5, 
#   normMethod = "None",
#   plotName = "Footprints-No-Normalization",
#   addDOC = FALSE,
#   smoothWindow = 5
# )

# Feature Footprinting

# Create a TSS insertion profile . It is just a specialized sub-case of footprinting. create TSS insertion profiles without normalization for Tn5 bias. The main difference from our previous analyses is that we specify flank = 2000 to extend these footprints 2000 bp on either side of each TSS.

seTSS <- getFootprints(ArchRProj = ATACSeq_project_All5, 
                       positions = GRangesList(TSS = getTSS(ATACSeq_project_All5)), 
                       groupBy = "Clusters2", flank = 2000)
# then plot the TSS insertion profiles for each cell group using plotFootprints().

FeatureFootprint_TSSInsertion <- plotFootprints(seFoot = seTSS, ArchRProj = ATACSeq_project_All5, 
                                                normMethod = "None", 
                                                plotName = "FeatureFootprint_TSSInsertion_TSS-No-Normalization",
                                                addDOC = FALSE, flank = 2000, flankNorm = 100)

# Integrative Analysis with ArchR
# 1. ONLY ATACseq data to identify Peak Co-accesibility to predict regulatory interactions
# 2. analyses that integrate scRNA-seq data such as prediction of enhancer activity through peak-to-gene linkage analysis. 

# Creating Low-Overlapping Aggregates of Cells

# Integrative Analysis #1. ONLY ATACseq data to identify Peak Co-accesibility to predict regulatory interactions

# Co-accessibility with ArchR
ATACSeq_project_All5 <- addCoAccessibility(ArchRProj = ATACSeq_project_All5, reducedDims = "IterativeLSI_all5")
cA <- getCoAccessibility(ArchRProj = ATACSeq_project_All5, corCutOff = 0.75, resolution = 1000, returnLoops = TRUE)
print(cA)

# Plotting browser tracks of Co-accessibility
markerGenes  <- mg
CoAccessibility_BrowserTrack <- plotBrowserTrack(ArchRProj = ATACSeq_project_All5, groupBy = "Clusters2", 
                                                 geneSymbol = markerGenes, 
                                                 upstream = 5000,
                                                 downstream = 5000,
                                                 loops = getCoAccessibility(ATACSeq_project_All5))
plotPDF(plotList = CoAccessibility_BrowserTrack, 
        name = "Marker-Genes-with-CoAccessibility.pdf", 
        ArchRProj = ATACSeq_project_All5, 
        addDOC = FALSE, width = 8, height = 8)


# Peak2GeneLinkage with ArchR
ATACSeq_project_All5 <- addPeak2GeneLinks(ArchRProj = ATACSeq_project_All5, 
                                          corCutOff = 0.75,
                                          predictionCutoff = 0.4, # A numeric describing the cutoff for RNA integration to use when picking cells for groupings.
                                          reducedDims = "IterativeLSI_all5")
# retrieve these peak-to-gene links in a similar fashion to how we retrieved co-accessibility interactions by using the getPeak2GeneLinks() function. 
Peak2GeneLinkage <- getPeak2GeneLinks(ArchRProj = ATACSeq_project_All5,
                                      corCutOff = 0.75,
                                      resolution = 1000,
                                      returnLoops = TRUE)
print(Peak2GeneLinkage)

# Plotting browser tracks with peak-to-gene links
Peak2GeneLinkage_BrowserTrack <- plotBrowserTrack(ArchRProj = ATACSeq_project_All5, groupBy = "Clusters2", 
                                                  geneSymbol = markerGenes, 
                                                  upstream = 5000,
                                                  downstream = 5000,
                                                  loops = getPeak2GeneLinks(ATACSeq_project_All5))
plotPDF(plotList = Peak2GeneLinkage_BrowserTrack, 
        name = "Marker-Genes-with-Peak2GeneLinkage_BrowserTrack.pdf", 
        ArchRProj = ATACSeq_project_All5, 
        addDOC = FALSE, width = 8, height = 8)

# Heatmap of peak2Gene Links
Heatmap_peak2Gene <- plotPeak2GeneHeatmap(ArchRProj = ATACSeq_project_All5, groupBy = "Clusters2")

################################## Identification of Positive TF-Regulators #################################
# Step 1. Identify Deviant TF Motifs
seGroupMotif <- getGroupSE(ArchRProj = ATACSeq_project_All5, useMatrix = "MotifMatrix", groupBy = "Clusters2")
# this SummarizedExperiment object comes from the MotifMatrix is has two seqnames - “deviations” and “z” - corresponding to the raw deviations and deviation z-scores from chromVAR.

# subset this SummarizedExperiment to just the deviation z-scores.
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
# identify the maximum delta in z-score between all clusters. This will be helpful in stratifying motifs based on the degree of variation observed across clusters.
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){rowMaxs(assay(seZ) - assay(seZ)[,x])}) %>% 
  Reduce("cbind", .) %>% rowMaxs

# Step 2. Identify Correlated TF Motifs and TF Gene Score/Expression

# To identify TFs whose motif accessibility is correlated with with their own gene activity (either by gene score or gene expression), we use the correlateMatrices() function and provide the two matrices that we are interested in, in this case the GeneScoreMatrix and the MotifMatrix. As mentioned previously, these correlations are determined across many low-overlapping cell aggregates determined in the lower dimension space specified in the reducedDims parameter.
corGSM_MM <- correlateMatrices(ArchRProj = ATACSeq_project_All5, useMatrix1 = "GeneScoreMatrix", 
                               useMatrix2 = "MotifMatrix", reducedDims = "IterativeLSI_all5")
# perform the same analysis using the GeneIntegrationMatrix instead of the GeneScoreMatrix.
corGSM_MM <- correlateMatrices(ArchRProj = ATACSeq_project_All5, useMatrix1 = "GeneIntegrationMatrix", 
                               useMatrix2 = "MotifMatrix", reducedDims = "IterativeLSI_all5")
# Step 3. Add Maximum Delta Deviation to the Correlation Data Frame
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
# Step 4. Identify Positive TF Regulators
# we consider positive regulators as those TFs whose correlation between motif and gene score (or gene expression) is greater than 0.5 with an adjusted p-value less than 0.01 and a maximum inter-cluster difference in deviation z-score that is in the top quartile. We apply these selection criteria and do a little text juggling to isolate the TF names.
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

# Step5. Plotting
# Having identified these positive TF regulators from gene scores and motif deviation z-scores, we can highlight them in a dot plot.
PositiveTFRegulator_basedonGSM <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )

# Same Analysis based on Gene Integration Matrix
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])

PositiveTFRegulator_basedonGIM <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )

