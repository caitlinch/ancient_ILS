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



location = "local"
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



#### 2. Prepare functions, variables and packages ####
# Source functions and dataset information
source(paste0(repo_dir, "code/func_prepare_trees.R"))
#source(paste0(repo_dir, "code/func_constraint_trees.R"))
#source(paste0(repo_dir, "code/func_empirical_tree_estimation.R"))



#### 3. Extract genes ####
# Create directory for the genes
gene_output_directory <- paste0(alignment_dir, "genes/")
if (dir.exists(gene_output_directory) == FALSE){
  dir.create(gene_output_directory)
}
# Open csv file with alignment information
gene_file_df <- read.csv(paste0(repo_dir, "output/empirical_gene_files.csv"), stringsAsFactors = F)
gene_file_df$gene_directory <- paste0(gene_output_directory, gene_file_df$dataset_id, "/")
# Create a loop to apply to each alignment
for (i in 1:nrow(gene_file_df)){
  # Identify id
  id <- gene_file_df$dataset_id[[i]]
  # Extract alignment and partition for the id
  id_alignment <- paste0(alignment_dir, gene_file_df$alignment_file[i])
  id_partition <- paste0(alignment_dir, gene_file_df$partition_file[i])
  # Create directory to save genes in
  id_directory <- gene_file_df$gene_directory[i]
  if (dir.exists(id_directory) == FALSE){
    dir.create(id_directory)
  }
  # Call the function to separate individual genes
  extract.all.genes(alignment_file = id_alignment, partition_file = id_partition, dataset_id = id, gene_directory = id_directory)
}

alignment_file <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/Philippe2009.Philippe_etal_superalignment_FixedNames.aa.alignment.nex"
partition_file <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/Philippe2009_partitions_formatted.nex"
gene_directory <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/genes/Philippe2009.Philippe_etal_superalignment_FixedNames/"
dataset_id <- "Philippe2009.Philippe_etal_superalignment_FixedNames"

extract.all.genes(alignment_file, partition_file, dataset_id, gene_directory)


