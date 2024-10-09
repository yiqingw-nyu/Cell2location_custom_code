# This script converts the single cell Seurat object to AnnData object 
# that can be used as input to Cell2location.
# Author: Yiqing Wang
# Date: 2024-6-16

library(Seurat)
library(SeuratDisk)

# Read in the Seurat object
test_data <- "path/to/data/"
ref_folder <- file.path(test_data, "reference")
results_folder <- file.path(test_data, "test_results")
sc <- readRDS(file.path(ref_folder, "reference_file.rds"))

# Save as H5 Seurat
SaveH5Seurat(sc, filename = file.path(results_folder, "reference_file.h5Seurat"))

# Convert H5 Seurat to AnnData (h5ad file)
# It is important to specify RNA assay, as it is not the default assay (the default assay is "integrated").
# Cell2location needs untransformed, unnormalized data.
Convert(file.path(results_folder, "reference_file.h5Seurat"), 
        dest = "h5ad", assay = "RNA")

# The resulting h5ad file was moved to the "RNA_assay" folder.

# In the resulting AnnData, "X" attribute is filled with "data" of the RNA assay;
# "raw" attribute is filled with "counts" of the RNA assay.
# See https://mojaveazure.github.io/seurat-disk/reference/Convert.html for details.
