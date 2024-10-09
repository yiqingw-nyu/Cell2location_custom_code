# This script conducts estimation of cell type signatures. 
# It is run on HPC, after preprocessing is done on the single cell dataset.

# Author: Yiqing Wang
# Date: 2024-6-18

# INPUT: preprocessed, filtered single cell data (an AnnData object)
# OUTPUT: 
# 1) NB regression model 
# 2) AnnData object with estimated cell type signatures 
# 3) CSV file with estimated cell type signatures 
# 4) ELBO loss history of model training

import scanpy as sc  # single cell library
import os
import matplotlib.pyplot as plt

# Tell theano (Cell2location dependency) to use GPU; this line is to be run before importing cell2location.
# This line is taken from Cell2location common error page. 
os.environ["THEANO_FLAGS"] = "device=cuda,floatX=float32,force_device=True"

import cell2location

### 1. Loading the filtered single cell data

dir = "path/to/data/"
os.chdir(dir)
results_folder = "./test_results/"
ref_run_name = f"{results_folder}/run_name"

adata_ref = sc.read_h5ad(f"{ref_run_name}/adata_ref_filtered.h5ad")

### 2. Estimating reference cell type signatures via NB regression

# Prepare anndata for the regression model
# All fields (in quotation marks) are columns in adata.obs, which is cell metadata.
cell2location.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    # 10X reaction / sample / batch to correct for
    batch_key="ident",
    # cell type annotation, covariate used for constructing signatures
    labels_key="label",
    # other categorical covariates to correct for
    categorical_covariate_keys=["external_donor_name"],
)

# Create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# View Anndata setup as a sanity check
mod.view_anndata_setup()

# Train the model
# max_epochs, batch_size, and lr can be customized to optimize the model.
# Two parameter settings were attempted, and the second was adopted for further analysis.
# mod.train(max_epochs=250, train_size=1, batch_size=1024, lr=0.01)
mod.train(max_epochs=500, train_size=1, batch_size=1024, lr=0.001)

### 3. Model QC and exporting posterior distributions of model parameters to Anndata

# Plot ELBO loss history, excluding first 100 epochs
mod.plot_history(iter_start=100)
plt.savefig(f"{ref_run_name}/loss_history_100.png")

# Export the estimated cell type signatures (summary of the posterior distribution).
# num_samples is the number of samples to draw from the posterior distribution.
# batch size is data batch size (if you run out of memory, reduce this number).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={"num_samples": 1000, "batch_size": 2048}
)

# Save the model
mod.save(f"{ref_run_name}", overwrite=True)

# Save the anndata object with results
adata_file = f"{ref_run_name}/sc_trained.h5ad"
adata_ref.write(adata_file)

### 4. Exporting estimated expression of each cell type to CSV

# Export estimated expression in each cluster (cell type) for spatial mapping
# Factor names are cell type names.
# "means_per_cluster_mu_fg_{cell type}" is the estimated expression of each gene of that cell type (mean of the posterior)
if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
    inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
        [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
    ].copy()
else:
    inf_aver = adata_ref.var[
        [f"means_per_cluster_mu_fg_{i}" for i in adata_ref.uns["mod"]["factor_names"]]
    ].copy()

inf_aver.columns = adata_ref.uns["mod"][
    "factor_names"
]  # rename columns to cell type names
inf_aver.iloc[0:5, 0:5]  # display the first 5 rows and columns

# Save the estimated expressions to a CSV
inf_aver.to_csv(f"{ref_run_name}/inf_aver.csv", index=True)