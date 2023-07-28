# Initialization

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

# Input data
day1 <- "/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day12Sample/2_Processed_data/V16_d12_230105/fragments.tsv.gz"
day5 <- "/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day5Sample/2_Processed_data/AJ_v12/fragments.tsv.gz"
day10 <- "/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day10Sample/compressed_tars/2_Processed_data/AJ_V12_d10_L129_UI/fragments.tsv.gz"
day25 <- "/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/Day25Sample/fragments.tsv.gz"
day4dpi <- "/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/4dpiSample/fragments.tsv.gz"

# Samples Input Vector
# samples <- c(day1, day5, day10, day25, day4dpi)
# # Set understandable names
# names(samples) <- c("d01", "d05", "d10", "d25", "Inf-d4")
# # inputFiles <- c('scATAC' = samples)
# inputFiles <- c(samples)
# print(names(inputFiles))

# Add respective Genome
addArchRGenome("mm10")

# Set Different ( higher) thresholds for Day1 sample separately.
names(day1) <- c("d01")
inputFile1 <- c('scATAC' = day1)
ArrowFiles1 <- createArrowFiles(
  inputFiles = inputFile1,
  sampleNames = names(day1),
  # filterFrags = 2500,
  minTSS = 8, # 4 was the default value. Dont set this too high;can always increase later
  minFrags = 5500,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  subThreading = TRUE,
  cleanTmp = TRUE,
  force = TRUE, # will recreate arrow files anew each time
  logFile = createLogFile("onlyDay1")
)

# Samples for all other time points
samples <- c(day5, day10, day25, day4dpi)
# # Set understandable names
names(samples) <- c("d05", "d10", "d25", "Inf-d4")
inputFiles <- c('scATAC' = samples)
print(names(inputFiles))

# Create Arrow Files

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(samples),
  # filterFrags = 2500,
  minTSS = 5, # 4 was the default value. Dont set this too high;can always increase later
  minFrags = 2500,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  subThreading = TRUE,
  cleanTmp = TRUE,
  force = TRUE, # will recreate arrow files anew each time
  logFile = createLogFile("All5")
)

# Quality Control
# Compute Doublet Score
doubScores <- addDoubletScores(
  input = c(ArrowFiles1, ArrowFiles),
  k = 30, # Refers to how many cells near a "pseudo-doublet" to count
  knnMethod = "LSI", # Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force = FALSE,
  verbose = TRUE
)


# Setup ArchR project
ATACSeq_project_All5 <- ArchRProject(
  ArrowFiles = c(ArrowFiles1, ArrowFiles),
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

# Plotting QC metrics
# UniqueFrags Vs Enrichment
QCMetric_df <- getCellColData(ATACSeq_project_All5, select = c("log10(nFrags)", "TSSEnrichment"))
(scatterplot_UniqueFragsVsEnrichment_PreFiltering <- ggPoint(
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
    geom_vline(xintercept = 3.4, lty = "dashed") +
    theme(axis.title = element_text(size = 15), axis.text = element_text(size = 15)))
plotPDF(scatterplot_UniqueFragsVsEnrichment_PreFiltering, 
        name = "TSS-vs-Frags_PreFiltering.pdf", ArchRProj = ATACSeq_project_All5)


create_GroupPlot <- function(projectname, cellColData_colname) {
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

my_gg_theme <- theme(
  axis.title = element_text(size = 15),
  axis.text = element_text(size = 15),
  legend.text = element_text(size = 10)
)


# violin plot for each sample for TSSEnrichment.
Group_plot_TSS_violin_PreFiltering <- create_GroupPlot(ATACSeq_project_All5, "TSSEnrichment") + ggtitle("TSSEnrichment") + my_gg_theme

### Plots (per sample) for log10 (unique nuclear fragments)
Group_plot_nFrags_violin_PreFiltering <- create_GroupPlot(ATACSeq_project_All5, "log10(nFrags)") + ggtitle("nFrags") + my_gg_theme

### Plots (per sample) for BlacklistRatio
Group_plot_BLR_violin_PreFiltering <- create_GroupPlot(ATACSeq_project_All5, "BlacklistRatio") + ggtitle("BlacklistRatio") + my_gg_theme

### Plots (per sample) for NucleosomeRatio
Group_plot_NR_violin_PreFiltering <- create_GroupPlot(ATACSeq_project_All5, "NucleosomeRatio") + ggtitle("NucleosomeRatio") + my_gg_theme 

### Plots (per sample) for DoubletScore
Group_plot_DS_violin_PreFiltering <- create_GroupPlot(ATACSeq_project_All5, "DoubletScore") + ggtitle("DoubletScore") + my_gg_theme

### Plots (per sample) for DoubletEnrichment

Group_plot_DE_violin_PreFiltering <- create_GroupPlot(ATACSeq_project_All5, "DoubletEnrichment") + ggtitle("DoubletEnrichment") + my_gg_theme

plotPDF(Group_plot_TSS_violin_PreFiltering, Group_plot_nFrags_violin_PreFiltering,
        Group_plot_BLR_violin_PreFiltering, Group_plot_NR_violin_PreFiltering,
        Group_plot_DS_violin_PreFiltering, Group_plot_DE_violin_PreFiltering, 
        name = "QC-Metrics_PreFiltering.pdf", 
        ArchRProj = ATACSeq_project_All5, 
        width = 8, height = 8)

FragSizePlot <- plotFragmentSizes(ArchRProj = ATACSeq_project_All5) + my_gg_theme
TSSEnrichmentPlot <- plotTSSEnrichment(ArchRProj = ATACSeq_project_All5) + my_gg_theme
plotPDF(FragSizePlot, TSSEnrichmentPlot, name = "QC-Sample-FragSizes-TSSProfile_PreFiltering.pdf", ArchRProj = ATACSeq_project_All5, width = 8, height = 8)

# Filtering
# Filter out low Quality Cells from ArchR Project, based on the violin plots of QC Metrics.

ATACSeq_project_All5 <- ATACSeq_project_All5[ATACSeq_project_All5$NucleosomeRatio < 2.5 &
                                               ATACSeq_project_All5$TSSEnrichment > 6 &
                                               ATACSeq_project_All5$BlacklistRatio < 0.05 &
                                               ATACSeq_project_All5$DoubletScore < 50 &
                                               ATACSeq_project_All5$DoubletEnrichment < 3.5]

# QC Plots PostFiltering
# violin plot for each sample for TSSEnrichment.
(Group_plot_TSS_violin_PostFiltering <- create_GroupPlot(ATACSeq_project_All5, "TSSEnrichment") + ggtitle("TSSEnrichment")) + my_gg_theme
### Plots (per sample) for log10 (unique nuclear fragments) #### log10 (unique nuclear fragments)
Group_plot_nFrags_violin_PostFiltering <- create_GroupPlot(ATACSeq_project_All5, "log10(nFrags)") + ggtitle("nFrags") + my_gg_theme
# Plots (per sample) for BlacklistRatio
Group_plot_BLR_violin_PostFiltering <- create_GroupPlot(ATACSeq_project_All5, "BlacklistRatio") + ggtitle("BlacklistRatio") + my_gg_theme
# Plots (per sample) for NucleosomeRatio
Group_plot_NR_violin_PostFiltering <- create_GroupPlot(ATACSeq_project_All5, "NucleosomeRatio") + ggtitle("NucleosomeRatio") + my_gg_theme
# Plots (per sample) for DoubletScore
Group_plot_DS_violin_PostFiltering <- create_GroupPlot(ATACSeq_project_All5, "DoubletScore") + ggtitle("DoubletScore") + my_gg_theme
# Plots (per sample) for DoubletEnrichment
(Group_plot_DE_violin_PostFiltering <- create_GroupPlot(ATACSeq_project_All5, "DoubletEnrichment") + ggtitle("DoubletEnrichment")) + my_gg_theme
plotPDF(Group_plot_TSS_violin_PostFiltering, Group_plot_nFrags_violin_PostFiltering, Group_plot_BLR_violin_PostFiltering, Group_plot_NR_violin_PostFiltering, Group_plot_DS_violin_PostFiltering, Group_plot_DE_violin_PostFiltering, name = "QC-Metrics_PostFiltering.pdf", ArchRProj = ATACSeq_project_All5, width = 8, height = 8)

# Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
FragSizePlot_PostFiltering <- plotFragmentSizes(ArchRProj = ATACSeq_project_All5) + my_gg_theme

TSSEnrichmentPlot_PostFiltering <- plotTSSEnrichment(ArchRProj = ATACSeq_project_All5) + my_gg_theme
plotPDF(FragSizePlot_PostFiltering, TSSEnrichmentPlot_PostFiltering, name = "QC-Sample-FragSizes-TSSProfile_PostFiltering.pdf", ArchRProj = ATACSeq_project_All5, width = 8, height = 8)

# Save the ArchR project with the Preprocessed and filtered Project.
saveArchRProject(
  ArchRProj = ATACSeq_project_All5,
  outputDirectory = getOutputDirectory(ATACSeq_project_All5),
  overwrite = TRUE,
  load = TRUE,
  logFile = createLogFile("saveArchRProject_All5"),
)
