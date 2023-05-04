# R script to perform Dimensional Reduction via peakVI model from scvi tools
# Author: Keshav Prasad Gubbi

#' SECTION1: Initialization

Package_List <- c("ComplexHeatmap", "writexl", "here", "patchwork", "tidyverse", "hdf5r", 'sceasy', "reticulate", "EnsDb.Mmusculus.v79", "GenomeInfoDb", "pheatmap", "ArchRtoSignac", "Signac", "Seurat", "biovizBase", "stringr", "anndata", "SeuratDisk")
not_installed <- Package_List[!(Package_List %in% installed.packages()[, "Package"])] # Extract not installed packages
if (length(not_installed)) install.packages(not_installed) # Install the uninstalled packages
invisible(lapply(Package_List, suppressPackageStartupMessages(library), character.only = TRUE))

set.seed(1234)
# Install the reticulate package if you haven't already
# if (!requireNamespace("reticulate", quietly = TRUE)) {
#   install.packages("reticulate")
# }

# Load the reticulate package
# File Path Declarations
here::i_am(path = "peak6.R")
paste0(here())

use_python("/home/keshavprasad/miniconda3/bin/python", required = TRUE)
# /home/keshavprasad/miniconda3/bin/python
# /home/keshavprasad/miniconda3/envs/scvi-env/bin/python
sc <- import("scanpy") 
# Need to fix this. this is where the problem lies!
scvi <- import("scvi") 

 # Load the merged object 
atac_data <- readRDS(file = 'merged_signacobject_all5_unfiltered.rds')
  
# Convert the seurat object into a anndata object
adata <- convertFormat(atac_data, from = "seurat", to = "anndata", 
              main_layer = "counts", 
              assay = "peaks", drop_single_values = FALSE, 
              outFile = "converted_signacobject_all5_unfiltered.h5ad"
              )
print(adata) # Note generally in Python, dataset conventions are obs x var

# Save the 

# write_h5ad(anndata = adata, filename = 'anndata_object_from_signac.h5ad', compression = "gzip")

# Run the standard PeakVI workflow
# scvi$models$PEAKVI$setup_anndata(adata)
# pvi <- scvi$model$PEAKVI(adata)
# pvi$train()

