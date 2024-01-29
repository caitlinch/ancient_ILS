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
  alignment_dir       <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/alignments/"
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
# Source required packages
library(ape)

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
# Create directory to save all genes in 
dirs_to_create <- gene_file_df$gene_directory[! dir.exists(gene_file_df$gene_directory)]
if (length(dirs_to_create) > 0){
  lapply(1:length(dirs_to_create), function(x){dir.create(dirs_to_create[x])})
}

# Separate each alignment into individual genes 
# Note: set `create.gene.alignments=FALSE` to extract csv with gene details without extracting genes
id_op <- lapply(1:nrow(gene_file_df), function(i){extract.all.genes(alignment_file = paste0(alignment_dir, gene_file_df$alignment_file[i]),
                                                                    partition_file = paste0(alignment_dir, gene_file_df$partition_file[i]), 
                                                                    dataset_id = gene_file_df$dataset_id[[i]], 
                                                                    gene_directory = gene_file_df$gene_directory[i], 
                                                                    create.gene.alignments = FALSE)} )
gene_df <- as.data.frame(do.call(rbind, id_op))



#### 4. Update constraint trees for each gene ####
# List all constraint trees
constraint_tree_files <- paste0(repo_dir, "constraint_trees/", list.files(paste0(repo_dir, "constraint_trees/")))
# Process constraint trees for the taxa in each gene
ct_op <- lapply(1:nrow(gene_df), trim.constraint.tree.taxa, gene_df = gene_df)

# Works for: Borowiec 2015,

row_id <- 109
