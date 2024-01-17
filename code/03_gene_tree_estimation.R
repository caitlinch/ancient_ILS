## caitlinch/ancient_ILS/code/03_gene_tree_estimation.R
# This script estimates maximum likelihood trees under different models of substitution for 14 empirical data sets
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## Specify parameters:
# location                    <- Where the script is being run
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# alignment_dir               <- Directory containing alignments for all data sets
#                                   Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                                   E.g. Cherryh2022.alignment1.aa.fa
# output_dir                  <- Directory for IQ-Tree output (trees and tree mixtures)
# iqtree2                     <- Location of IQ-Tree2 executable
# iqtree2_num_bootstraps      <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree
# iqtree2_num_threads         <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
# astral                      <- Location of ASTRAL executable
# astral_constrained          <- Location of ASTRAL constrained tree version executable

## Specify control parameters (all take logical values TRUE or FALSE):
# prepare.parameters                    <- Create dataframe with parameters for gene tree/ML tree/constrained tree estimation: T/F

location = "dayhoff"
if (location == "local"){
  repo_dir            <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  alignment_dir       <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/"
  output_dir          <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_empirical_tree_estimation/"
  iqtree2             <- "iqtree2"
  iqtree2_num_threads  <- "AUTO"
  astral              <- "/Users/caitlincherryh/Documents/Executables/ASTRAL-5.7.8-master/Astral/astral.5.7.8.jar"
  astral_constrained  <- "/Users/caitlincherryh/Documents/Executables/ASTRAL-Constrained-search-master/Astral2/astral.5.6.9.jar"
  
} else if (location == "dayhoff" | location == "rona" ){
  if (location == "dayhoff"){
    repo_dir <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  } else if (location == "rona"){
    repo_dir <- "/home/caitlin/ancient_ILS/"
  }
  alignment_dir <- paste0(repo_dir, "data_all/")
  output_dir <-  paste0(repo_dir, "output/")
  iqtree2 <- paste0(repo_dir, "iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2")
  iqtree2_num_threads <- 20
  astral <- paste0(repo_dir, "astral/Astral/astral.5.7.8.jar")
  astral_constrained <- paste0(repo_dir, "astral_constrained/Astral/astral.5.6.9.jar")
}

# Set parameters that are identical for all run locations
iqtree2_num_bootstraps <- 1000

# Set control parameters
control_parameters <- list(prepare.parameters = FALSE)



#### 2. Prepare functions, variables and packages ####
# Open packages
library(stringr)

# Source functions and dataset information
source(paste0(repo_dir, "code/func_prepare_trees.R"))
source(paste0(repo_dir, "code/func_constraint_trees.R"))
source(paste0(repo_dir, "code/func_empirical_tree_estimation.R"))
source(paste0(repo_dir, "code/data_dataset_info.R"))

# Prepare directories, if they do not exist
if (dir.exists(paste0(output_dir, "constraint_trees/")) == FALSE){dir.create(paste0(output_dir, "constraint_trees/"))}
if (dir.exists(paste0(output_dir, "tree_estimation/")) == FALSE){dir.create(paste0(output_dir, "tree_estimation/"))}
if (dir.exists(paste0(output_dir, "hypothesis_trees/")) == FALSE){dir.create(paste0(output_dir, "hypothesis_trees/"))}

# Remove the individual dataset lists (only need collated lists) (yes it is a bit cheeky to hard code the removal)
rm(borowiec2015_list, chang2015_list, dunn2008_list, hejnol2009_list, laumer2018_list, laumer2019_list, moroz2014_list, nosenko2013_list, philippe2009_list,
   philippe2011_list, pick2010_list, ryan2013_list, simion2017_list, whelan2015_list, whelan2017_list, models_list, all_taxa, all_models)

# Read in the .tsv file that contains a column for each dataset, with one row per taxa included in that dataset
alignment_taxa_df <- read.table(paste0(repo_dir, "output/dataset_included_taxa.tsv"), header = T)



#### 3. Prepare for tree estimation ####