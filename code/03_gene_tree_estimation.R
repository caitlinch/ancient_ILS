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
  repo_dir              <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  alignment_dir         <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/alignments/"
  gene_output_dir       <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/genes/"
  output_dir            <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/05_output_files/"
  iqtree2               <- "iqtree2"
  iqtree2_num_threads   <- "AUTO"
} else if (location == "dayhoff" | location == "rona" ){
  if (location == "dayhoff"){
    repo_dir          <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  } else if (location == "rona"){
    repo_dir          <- "/home/caitlin/ancient_ILS/"
  }
  alignment_dir       <- paste0(repo_dir, "data_all/")
  gene_output_dir     <- paste0(repo_dir, "genes/")
  output_dir          <- paste0(repo_dir, "output/")
  iqtree2             <- paste0(repo_dir, "iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2")
  iqtree2_num_threads <- 20
}

# Set parameters that are identical for all run locations
iqtree2_num_bootstraps <- 1000



#### 2. Prepare functions, variables and packages ####
# Source required packages
library(ape)

# Source functions and dataset information
source(paste0(repo_dir, "code/func_prepare_trees.R"))
source(paste0(repo_dir, "code/func_empirical_tree_estimation.R"))
source(paste0(repo_dir, "code/func_data_analysis.R"))




#### 3. Extract genes ####
gene_df_filepath <- paste0(output_dir, "genes_001_individualGene_files.csv")
# Open gene_df if it exists. If it doesn't, then create it.
if (file.exists(gene_df_filepath) == FALSE){
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
  gene_df <- as.data.frame(do.call(rbind, lapply(1:nrow(gene_file_df), 
                                                 function(i){extract.all.genes(alignment_file = paste0(alignment_dir, gene_file_df$alignment_file[i]),
                                                                               partition_file = paste0(alignment_dir, gene_file_df$partition_file[i]), 
                                                                               dataset_id = gene_file_df$dataset_id[[i]], 
                                                                               gene_directory = gene_file_df$gene_directory[i], 
                                                                               create.gene.alignments = FALSE)} ) ) )
  # Write gene_df to file
  write.csv(gene_df, file = gene_df_filepath, row.names = FALSE)
} else {
  gene_df <- read.csv(gene_df_filepath, stringsAsFactors = FALSE)
}



#### 4. Update constraint trees for each gene ####
constraint_df_filepath <- paste0(output_dir, "genes_002_individualGene_constraintTreeFiles.csv")
# Open constraint_df if it exists. If it doesn't, then create it.
if (file.exists(constraint_df_filepath) == FALSE){
  # List all constraint trees
  constraint_tree_files <- paste0(repo_dir, "constraint_trees/", list.files(paste0(repo_dir, "constraint_trees/")))
  # Process constraint trees for the taxa in each gene
  constraint_df <- as.data.frame(do.call(rbind,  lapply(1:nrow(gene_df), 
                                                        trim.constraint.tree.taxa, 
                                                        gene_df = gene_df, constraint_tree_files = constraint_tree_files) ) )
  # Remove directory path from constraint tree file paths
  constraint_df$constraint_tree_1 <- basename(constraint_df$constraint_tree_1)
  constraint_df$constraint_tree_2 <- basename(constraint_df$constraint_tree_2)
  constraint_df$constraint_tree_3 <- basename(constraint_df$constraint_tree_3)
  ## Add new columns for tree estimation
  constraint_df$initial_model             <- "MFP"
  constraint_df$unconstrained_tree_prefix <- paste0(constraint_df$gene_id, ".MFP")
  # Write the constraint tree to file
  write.csv(constraint_df, file = constraint_df_filepath, row.names = FALSE)
} else {
  constraint_df <- read.csv(constraint_df_filepath, stringsAsFactors = FALSE)
}



#### 5. Estimate unconstrained gene tree ####
## Construct IQ-Tree command line
initial_df_filepath <- paste0(output_dir, "genes_003_individualGene_InitialIQTreeCommand.csv")
# Open initial_df if it exists. If it doesn't, then create it.
if (file.exists(initial_df_filepath) == FALSE){
  # Create IQ-Tree command line
  initial_df <- as.data.frame(do.call(rbind, lapply(1:nrow(constraint_df), estimate.empirical.single.gene.tree.wrapper, 
                                                    dataframe = constraint_df, 
                                                    iqtree2_path = iqtree2, 
                                                    iqtree2_num_threads = iqtree2_num_threads, 
                                                    estimate.trees = FALSE) ) )
  # Write the initial run df to file
  write.csv(initial_df, file = initial_df_filepath, row.names = FALSE)
} else {
  initial_df <- read.csv(initial_df_filepath, stringsAsFactors = FALSE)
}

## Extract models from IQ-Tree initial (unconstrained) run
initial_run_df_filepath <- paste0(output_dir, "genes_004_individualGene_InitialIQTreeResults.csv")
# Open initial_run_df if it exists. If it doesn't, then create it.
if (file.exists(initial_run_df_filepath) == FALSE){
  # Extract output from iqtree files
  initial_run_df <- as.data.frame(do.call(rbind, lapply(1:nrow(initial_df), 
                                                        extract.unconstrained.tree.details, 
                                                        dataframe = initial_df) ) ) 
  # Write the initial run output df to file
  write.csv(initial_run_df, file = initial_run_df_filepath, row.names = FALSE)
} else {
  initial_run_df <- read.csv(initial_df_filepath, stringsAsFactors = FALSE)
}


#### 6. Estimate constrained gene trees and sCF/quartet scores ####
# Estimate gene tree (unconstrained)

# Extract model



#### 7. Calculate AU test for each gene ####







