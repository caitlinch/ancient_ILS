# ancient_ILS/code/00_plot_gene_lengths.R
## This script prepares the three hypothesis trees for the animal tree of life
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# repo_dir                              <- Location of caitlinch/ancient_ILS github repository
# partition_dir                         <- directory to save constrained ML trees estimated from the alignment

repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
partition_dir               <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_partition_files/"




#### 2. Open packages and functions ####
## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))

## Create plot directories
# Check whether plot file exists for repository
repo_output <- paste0(repo_dir, "plots/")
if (dir.exists(repo_output) == FALSE){
  dir.create(repo_output)
}







