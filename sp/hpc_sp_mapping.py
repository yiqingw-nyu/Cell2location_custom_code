# This script conducts spatial mapping and model QC.
# This script is adapted to run on HPC.

# Author: Yiqing Wang
# Date: 2024-6-14

# INPUT:
# 1) CSV of cell type signatures from the single cell pipeline
# 2) spatial data (AnnData) that had MT genes removed

# OUTPUT:
# 1) trained spatial model
# 2) spatial AnnData with posterior distributions of model parameters, such as cell abundance
# 3) ELBO loss history plots with some epochs skipped
# 4) reconstruction accuracy plot

import os
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

import cell2location

print("import complete")

### 1. Load cell type signature and Visium spatial data

dir = "path/to/data"
os.chdir(dir)
sample = "D1"
results_folder = "./test_results/"
ref_run_name = f"{results_folder}/run_name"
run_name = f"{results_folder}/{sample}_run_name"

# Load cell type signature inferred by the NB regression model
inf_aver = pd.read_csv(f"{ref_run_name}/inf_aver.csv", index_col=0)

print("single-cell cell type average expression data loaded")
print(f"number of genes in cell type signatures: {inf_aver.shape[0]}")

# Load Visium spatial data that had MT gene removed
adata_vis = sc.read_h5ad(f"{run_name}/adata_vis_MTremoved.h5ad")

print("spatial data loaded")
print(f"number of genes in Visium: {adata_vis.shape[1]}")

### 2. Convert gene row names from gene ids to gene symbols in the spatial data
# This step is optional. The purpose is so that the gene identifiers are consistent between spatial and single cell data.

# Since the symbols have duplicates, we need to drop duplicates (keep first) before converting.
duplicates = adata_vis.var.duplicated(subset="SYMBOL", keep="first")
print(f"Number of duplicated entries: {sum(duplicates)}")
print("\nDuplicated entries:")
print(adata_vis.var["SYMBOL"][duplicates])
adata_vis = adata_vis[:, ~duplicates].copy()  # drop duplicates

# Convert adata_vis index from ensembl id to gene symbol
adata_vis.var["gene_id"] = adata_vis.var.index
adata_vis.var = adata_vis.var.set_index("SYMBOL", drop=True)
adata_vis.var_names = adata_vis.var.index

### 3. Subset the datasets to shared genes

# Find shared genes from both spatial AnnData and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
print(f"number of shared genes: {len(intersect)}")

# subset datasets to shared genes
adata_vis = adata_vis[
    :, intersect
].copy()  # AnnData objects are designed to be subsetted like NumPy arrays
inf_aver = inf_aver.loc[
    intersect, :
].copy()  # .loc is used to select rows and columns by labels

print("datasets subsetted")

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(
    adata=adata_vis, batch_key="sample"
)  # label_key needs not be set as all spots come from the same sample

### 4. Map single cell signatures to spatial data

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis,
    cell_state_df=inf_aver,

    # The expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    # This number was based on manual counting.

    N_cells_per_location=8.6,

    # Hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    # 200 is for high regularization for low variability, 20 is vice versa.

    detection_alpha=200,
)

mod.view_anndata_setup() # model sanity check

# Training parameters epochs and learning rate can be optimized.
mod.train(
    max_epochs=30000,

    # train using full data (batch_size=None)
    # If None, no minibath is used, and the whole dataset is used at once.
    batch_size=None,

    # use all data points in training because
    # we need to estimate cell abundance at all locations
    train_size=1,
)

print("model trained")

### 5. Plot ELBO loss history during training

# include all epochs
plt.figure()
mod.plot_history()
plt.savefig(f"{run_name}/loss_history_all.png")
plt.close()

# skip first 100 epochs
plt.figure()
mod.plot_history(iter_start=100)
plt.savefig(f"{run_name}/loss_history_100.png")
plt.close()

# skip first 1000 epochs
plt.figure()
mod.plot_history(iter_start=1000)
plt.savefig(f"{run_name}/loss_history_1000.png")
plt.close()

print("loss history plots done")

### 6. Export parameter posterior distributions to AnnData and save trained files

# Export poseterior, parameters from tutorial
adata_vis = mod.export_posterior(
    adata_vis,
    sample_kwargs={"num_samples": 1000, "batch_size": mod.adata.n_obs},
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis) # how the model can be loaded later

# Save anndata object with results
adata_file = f"{run_name}/sp_mapped.h5ad"
adata_vis.write(adata_file)

print("model and Anndata object saved")

## Note on where the results are saved: ##
#
# The mean, std, 5% quantile, and 95% quantile of the posterior distribution of cell abundance (Ws,f) are saved in adata_vis.obsm.
# For example, adata_vis.obsm["q05_cell_abundance_w_sf"] is the 5% quantile of the posterior distribution of cell abundance.
# 
# The mean, sts, 5% quantile, and 95% quantile of the posterior distribution of all model parameters are saved in adata_vis.uns["mod"].
# For example, adata_vis_new.uns["mod"]["post_sample_means"]["detection_mean_y_e"] is 
# the mean of the posterior distribution of the latent average detection efficiency of each batch (refer to supplementary methods). 
##########################################

### 7. Reconstruction accuracy for model QC

# Plot the QC plot
plt.figure()
mod.plot_QC()
plt.savefig(f"{run_name}/reconstruction_accuracy.png")
plt.close()

print("QC plot done")
