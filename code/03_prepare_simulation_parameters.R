## caitlinch/ancient_ILS/code/03_prepare_simulation_parameters.R
# This script prepares simulations based on empirical phylogenetic datasets
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
# iqtree2_num_threads         <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
# astral                      <- Location of ASTRAL executable

location = "local"
if (location == "local"){
  repo_dir            <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  output_dir          <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/04_simulation_parameters/"
  iqtree2             <- "iqtree2"
  iqtree2_num_threads  <- "AUTO"
  astral              <- "/Users/caitlincherryh/Documents/Executables/ASTRAL-5.7.8-master/Astral/astral.5.7.8.jar"
  
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
}



#### 2. Prepare functions, variables and packages ####
## Open packages
library(ape)

## Source files
# Open dataset info, move Simion2017 taxa into a separate object and remove unneeded objects
source(paste0(repo_dir, "code/data_dataset_info.R"))
simion2017_clades <- simion2017_list
rm(all_taxa, all_datasets, borowiec2015_list, chang2015_list, dunn2008_list, hejnol2009_list, 
   laumer2018_list, laumer2019_list, matrix_taxa, models_list, moroz2014_list, nosenko2013_list,
   philippe2009_list, philippe2011_list, pick2010_list, ryan2013_list, simion2017_list, 
   whelan2015_list, whelan2017_list)

## Open the Simion 2017 tree files
# List all files
simion_files <- grep("Simion2017", list.files(paste0(repo_dir, "empirical_tree/")), value = T)
# Open gene trees
gene_trees_path <- paste0(repo_dir, "empirical_tree/", grep("gene_trees.treefile", simion_files, value = T))
gene_trees <- read.tree(gene_trees_path)
# Open ASTRAL tree
astral_tree_path <- paste0(repo_dir, "empirical_tree/", grep("ASTRAL_tree.tre", simion_files, value = T))
astral_tree <- read.tree(astral_tree_path)



#### 3. Extract clade monophyly and clade depth ####




#### 4. Extract branch lengths for in and outgroups ####




#### 5. Extract branch lengths leading to outgroups ####




