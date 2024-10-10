# Custom Codes for Cell2location

This repository contains custom scripts to utilize [Cell2location](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html) for performing cell type deconvolution on Visium spatial transcriptomics (ST) data. 
These scripts are inspired by, and complement, the official Cell2location tutorials.

## Repository Structure

### `sc` Folder

This folder includes scripts for estimating cell type signatures from reference scRNA-seq data. The order of script execution and the expected outputs are outlined below:
![Cell2location Flowchart-Single Cell drawio](https://github.com/user-attachments/assets/2fcc6a99-04e4-449d-92bd-0c3d7af5a434)

### `sp` Folder

This folder contains scripts for spatial mapping of cell type abundances using reference signatures and the ST data. The order of script execution and the expected outputs are illustrated below:
![Cell2location Flowchart-Spatial drawio](https://github.com/user-attachments/assets/349ea6cc-177d-4e2d-bb33-891193c8d5a0)

## References

Kleshchevnikov, V., Shmatko, A., Dann, E. et al. Cell2location maps fine-grained cell types in spatial transcriptomics. Nat Biotechnol 40, 661–671 (2022). https://doi.org/10.1038/s41587-021-01139-4
