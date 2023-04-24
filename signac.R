#' SECTION1: Initialization

Package_List <- c("ComplexHeatmap", "writexl", "here", "patchwork", "tidyverse", "hdf5r", 
                  "EnsDb.Mmusculus.v79", "GenomeInfoDb", "pheatmap", "ArchRtoSignac", "Signac",
                  "Seurat", "biovizBase", "stringr")
not_installed <- Package_List[!(Package_List %in% installed.packages()[, "Package"])] # Extract not installed packages
if (length(not_installed)) install.packages(not_installed) # Install the uninstalled packages
invisible(lapply(Package_List, suppressPackageStartupMessages(library), character.only = TRUE))

set.seed(1234)

# File Path Declarations
here::i_am(path = "signac.R")
paste0(here())

data_folder <- '/media/keshavprasad/HornefLab_Data3/ATACSeq/Data/'

# now for Day01
# d01_counts <- Read10X_h5(filename = file.path(data_folder,
#                                           "Day12Sample/2_Processed_data/V16_d12_230105",
#                                           "filtered_peak_bc_matrix.h5"))
# d01_meta <- read.csv(file = file.path(data_folder,
#                                   "Day12Sample/2_Processed_data/V16_d12_230105",
#                                   'singlecell.csv'), header = TRUE, row.names = 1)
# d01_chrom_assay <- CreateChromatinAssay(counts = d01_counts,
#                                     sep = c(":", "-"),
#                                     genome = 'mm10',
#                                     fragments = file.path(data_folder, "Day12Sample/2_Processed_data/V16_d12_230105", 'fragments.tsv.gz'), min.cells = 10, min.features = 200)
# d01_data <- CreateSeuratObject(counts = d01_chrom_assay, assay = "peaks", meta.data = d01_meta)
# 
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79) # EnsDb.Mmusculus.v79 is a specific version of the Ensembl Database for the mouse genome assembly GRCm38 (mm10).
# seqlevelsStyle(annotations) <- 'UCSC'
# 
# 
# Annotation(d01_data) <- annotations
# 
# d01_data <- NucleosomeSignal(object = d01_data)
# d01_data <- TSSEnrichment(object = d01_data, fast = FALSE)
# d01_data$blacklist_ratio <- d01_data$blacklist_region_fragments / d01_data$peak_region_fragments
# d01_data$pct_reads_in_peaks <- d01_data$peak_region_fragments / d01_data$passed_filters * 100
# 
# 
# (VlnPlot(
#   object = d01_data,
#   features = c('peak_region_fragments', 'pct_reads_in_peaks',
#                'blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'),
#   pt.size = 0.1,
#   ncol = 5
# ))



# low_prf <- quantile(d01_data[["peak_region_fragments"]]$peak_region_fragments, probs = 0.02)
# hig_prf <- quantile(d01_data[["peak_region_fragments"]]$peak_region_fragments, probs = 0.98)
# low_prp <- quantile(d01_data[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.02)
# high_blr <- quantile(d01_data[["blacklist_ratio"]]$blacklist_ratio, probs = 0.98)
# hig_ns <- quantile(d01_data[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98)
# low_ts <- quantile(d01_data[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)
# 
# print(low_prf)
# print(hig_prf)
# print(low_prp)
# print(high_blr)
# print(hig_ns)
# print(low_ts)
# 
# d01_data <- subset(
#   x = d01_data,
#   subset = peak_region_fragments > low_prf &
#     peak_region_fragments < hig_prf &
#     pct_reads_in_peaks > low_prp &
#     blacklist_ratio < high_blr &
#     nucleosome_signal < hig_ns &
#     TSS.enrichment > low_ts
# )

# Normalization, dimension reduction

# d01_data <- RunTFIDF(d01_data)
# 
# d01_data <- FindTopFeatures(d01_data, min.cutoff = 'q0')
# 
# d01_data <- RunSVD(d01_data)
# 
# DepthCor(d01_data)

###############################################
# Read in the files :
import_atac <- function(count_path, meta_path, fragment_path){
  counts <- Read10X_h5(filename = count_path)
  meta <- read.csv(file = meta_path, header = TRUE, row.names = 1)
  chrom_assay <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), genome = 'mm10', fragments = fragment_path, min.cells = 10, min.features = 200)

  data <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = meta)

  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79) # EnsDb.Mmusculus.v79 is a specific version of the Ensembl Database for the mouse genome assembly GRCm38 (mm10).
  seqlevelsStyle(annotations) <- 'UCSC'
  Annotation(data) <- annotations

  data <- NucleosomeSignal(object = data) #fragment ratio 147-294: <147  ---  mononucleosome:nucleosome-free

  data <- TSSEnrichment(object = data, fast = FALSE)

  data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments

  data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100

  # low_prf <- quantile(data[["peak_region_fragments"]]$peak_region_fragments, probs = 0.02)
  # hig_prf <- quantile(data[["peak_region_fragments"]]$peak_region_fragments, probs = 0.98)
  # low_prp <- quantile(data[["pct_reads_in_peaks"]]$pct_reads_in_peaks, probs = 0.02)
  # 
  # high_blr <- quantile(data[["blacklist_ratio"]]$blacklist_ratio, probs = 0.98)
  # 
  # hig_ns <- quantile(data[["nucleosome_signal"]]$nucleosome_signal, probs = 0.98)
  # 
  # low_ts <- quantile(data[["TSS.enrichment"]]$TSS.enrichment, probs = 0.02)
  # 
  # data <- subset(x = data,
  #                subset = peak_region_fragments > low_prf & peak_region_fragments < hig_prf & pct_reads_in_peaks > low_prp & blacklist_ratio < high_blr &
  #                  nucleosome_signal < hig_ns & TSS.enrichment > low_ts)
  # data <- RunTFIDF(data)
  # data <- FindTopFeatures(data, min.cutoff = 'q0')
  # data <- RunSVD(data)
  return(data)
}

(d01 <- import_atac(count_path = file.path(data_folder, "Day12Sample/2_Processed_data/V16_d12_230105", "filtered_peak_bc_matrix.h5"),
                    meta_path =  file.path(data_folder, "Day12Sample/2_Processed_data/V16_d12_230105", "singlecell.csv"),
                    fragment_path = file.path(data_folder, "Day12Sample/2_Processed_data/V16_d12_230105", "fragments.tsv.gz")
))
(d05 <- import_atac(count_path = file.path(data_folder, "Day5Sample/2_Processed_data/AJ_v12", "filtered_peak_bc_matrix.h5"),
                    meta_path =  file.path(data_folder, "Day5Sample/2_Processed_data/AJ_v12", "singlecell.csv"),
                    fragment_path = file.path(data_folder, "Day5Sample/2_Processed_data/AJ_v12", "fragments.tsv.gz")
))
(d10 <- import_atac(count_path = file.path(data_folder, "Day10Sample/compressed_tars/2_Processed_data/AJ_V12_d10_L129_UI", "filtered_peak_bc_matrix.h5"),
                    meta_path =  file.path(data_folder, "Day10Sample/compressed_tars/2_Processed_data/AJ_V12_d10_L129_UI", "singlecell.csv"),
                    fragment_path = file.path(data_folder, "Day10Sample/compressed_tars/2_Processed_data/AJ_V12_d10_L129_UI", "fragments.tsv.gz")
))
(d25 <- import_atac(count_path = file.path(data_folder, "Day25Sample/221123_Juergens_Hornef_microbiology_scATACseq/compressed_tars/_2_Processed_data/2_Processed_data/d25_221123", "filtered_peak_bc_matrix.h5"),
                    meta_path =  file.path(data_folder, "Day25Sample/221123_Juergens_Hornef_microbiology_scATACseq/compressed_tars/_2_Processed_data/2_Processed_data/d25_221123", "singlecell.csv"),
                    fragment_path = file.path(data_folder, "Day25Sample/221123_Juergens_Hornef_microbiology_scATACseq/compressed_tars/_2_Processed_data/2_Processed_data/d25_221123", "fragments.tsv.gz")
))
(Inf_dpi4 <- import_atac(count_path = file.path(data_folder, "4dpiSample/4dpi/genomics.rwth-aachen.de/data/221215_Juergens_Hornef_microbiology_scATACseq/2_Processed_data/V15_4dpi_221215", "filtered_peak_bc_matrix.h5"),
                         meta_path =  file.path(data_folder, "4dpiSample/4dpi/genomics.rwth-aachen.de/data/221215_Juergens_Hornef_microbiology_scATACseq/2_Processed_data/V15_4dpi_221215", "singlecell.csv"),
                         fragment_path = file.path(data_folder, "4dpiSample/4dpi/genomics.rwth-aachen.de/data/221215_Juergens_Hornef_microbiology_scATACseq/2_Processed_data/V15_4dpi_221215", "fragments.tsv.gz")
))

(d01_violinPlot <- VlnPlot(object = d01,
                           features = c('peak_region_fragments', 'pct_reads_in_peaks','blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'), pt.size = 0.1, ncol = 5))
(d05_violinPlot <- VlnPlot(object = d05,
                           features = c('peak_region_fragments', 'pct_reads_in_peaks','blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'), pt.size = 0.1, ncol = 5))
(d10_violinPlot <- VlnPlot(object = d10,
                           features = c('peak_region_fragments', 'pct_reads_in_peaks','blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'), pt.size = 0.1, ncol = 5))
(d25_violinPlot <- VlnPlot(object = d25,
                           features = c('peak_region_fragments', 'pct_reads_in_peaks','blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'), pt.size = 0.1, ncol = 5))
(Inf_dpi4_violinPlot <- VlnPlot(object = Inf_dpi4,
                           features = c('peak_region_fragments', 'pct_reads_in_peaks','blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'), pt.size = 0.1, ncol = 5))

# Merge Datasets
atac_data <- Reduce(function(x, y) merge(x, y), list(d01, d05, d10, d25, Inf_dpi4))
# Save an object to a file
saveRDS(atac_data, file = "merged_signacobject_all5_unfiltered.rds")


atac_data <- RunTFIDF(atac_data)
atac_data <- FindTopFeatures(atac_data, min.cutoff = 'q0')
atac_data <- RunSVD(atac_data)
DepthCorPlot_all5 <- DepthCor(atac_data)

# Determine Quality metrics for the merged dataset
(atac_data_all5_plot <- VlnPlot(object = atac_data, features = c('peak_region_fragments', 'pct_reads_in_peaks','blacklist_ratio', 'nucleosome_signal', 'TSS.enrichment'), pt.size = 0.1, ncol = 5))

# Subsetting the data based on Violin plot insights
atac_data <- subset(x = atac_data, subset = peak_region_fragments > 3000 & peak_region_fragments < 100000 & pct_reads_in_peaks > 30 & blacklist_ratio < 0.05 & nucleosome_signal < 3 & TSS.enrichment > 3)


# NonLinear Dimensional Reduction and Clustering
atac_data <- RunUMAP(object = atac_data, reduction = 'lsi', dims = 2:30)
(DimPlot_all5 <- DimPlot(object = atac_data, label = TRUE) + NoLegend())
atac_data <- FindNeighbors(object = atac_data, reduction = 'lsi', dims = 2:30)
atac_data <- FindClusters(object = atac_data, verbose = FALSE, algorithm = 3)
(DimPlot_all5 <- DimPlot(object = atac_data, label = TRUE) + NoLegend())
atac_data <- RunUMAP(object = atac_data, reduction = 'lsi', dims = 2:30)
