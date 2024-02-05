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
# gene_output_dir             <- Location for saving genes and running gene-wise analyses
# output_dir                  <- Directory for saving results files (csv and text)
# iqtree2                     <- Location of IQ-Tree2 executable
# iqtree2_num_bootstraps      <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree
# iqtree2_num_threads         <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"

location = "local"
if (location == "local"){
  repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  alignment_dir           <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/alignments/"
  gene_output_dir         <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_02_empirical_genes_initial_tree_estimation/"
  constraint_output_dir   <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_03_empirical_genes_constrained_trees/"
  output_dir              <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/05_output_files/"
  iqtree2                 <- "iqtree2"
  iqtree2_num_threads     <- 3
} else if (location == "dayhoff"){
  repo_dir                <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  alignment_dir           <- paste0(repo_dir, "data_all/")
  gene_output_dir         <- paste0(repo_dir, "genes/")
  constraint_output_dir   <- gene_output_dir
  output_dir              <- paste0(repo_dir, "output/")
  iqtree2                 <- paste0(repo_dir, "iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2")
  iqtree2_num_threads     <- 5
}

# Set parameters that are identical for all run locations
iqtree2_num_bootstraps <- 1000

# Set the Slurm lines for your specific server - so we can automate creaying the slurm files
slurm_start_lines <- c("#!/bin/bash", "#")
slurm_id_line <- "#SBATCH --job-name="
slurm_middle_lines <- c("#SBATCH --output=/mnt/data/dayhoff/home/u5348329/ancient_ILS/slurm_files/%j.%x.out", 
                        "#SBATCH --error=/mnt/data/dayhoff/home/u5348329/ancient_ILS/slurm_files/%j.%x.err",
                        "#SBATCH --partition=Standard",
                        "#",
                        "#SBATCH --time=96:00:00 # 4 days",
                        "#SBATCH --ntasks=1",
                        paste0("#SBATCH --cpus-per-task=", iqtree2_num_threads),
                        "#SBATCH --mem=3000 # 3GB required by IQ-Tree",
                        "#",
                        "#SBATCH --mail-user u5348329@anu.edu.au",
                        "#SBATCH --mail-type TIME_LIMIT,FAIL",
                        "",
                        "# Change to the directory containing each gene",
                        "cd /mnt/data/dayhoff/home/u5348329/ancient_ILS/genes/",
                        "",  
                        "# Run IQ-Tree2 command lines")



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
  if (dir.exists(gene_output_dir) == FALSE){
    dir.create(gene_output_dir)
  }
  # Open csv file with alignment information
  gene_file_df <- read.csv(paste0(repo_dir, "output/empirical_gene_files.csv"), stringsAsFactors = F)
  gene_file_df$gene_directory <- paste0(gene_output_dir, gene_file_df$dataset_id, "/")
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
## Check whether dataframe exists
constraint_df_filepath <- paste0(output_dir, "genes_002_individualGene_constraintTreeFiles.csv")
# Open constraint_df if it exists. If it doesn't, then create it.
if (file.exists(constraint_df_filepath) == FALSE){
  ## Create constraint trees and output in a dataframe
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
## Check whether dataframe exists
initial_df_filepath <- paste0(output_dir, "genes_003_individualGene_InitialIQTreeCommand.csv")
# Open initial_df if it exists. If it doesn't, then create it.
if (file.exists(initial_df_filepath) == FALSE){
  ## Create IQ-Tree command line for estimating unconstrained gene trees
  initial_df <- as.data.frame(do.call(rbind, lapply(1:nrow(constraint_df), 
                                                    estimate.empirical.single.gene.tree.wrapper, 
                                                    dataframe = constraint_df, 
                                                    iqtree2_path = iqtree2, 
                                                    iqtree2_num_threads = iqtree2_num_threads, 
                                                    estimate.trees = FALSE) ) )
  # Write the initial run df to file
  write.csv(initial_df, file = initial_df_filepath, row.names = FALSE)
  # Write just the IQ-Tree command lines to file
  initial_commands_filepath <- paste0(output_dir, "genes_003_individualGene_InitialIQTreeCommand_IQTreeCall.txt")
  write(initial_df$unconstrained_tree_iqtree2_call, file = initial_commands_filepath)
  
  ## Create Slurm files
  # Prepare for extracting specific rows
  filepath_start <- paste0(output_dir, "initial_iqtree_run_")
  max_i     <- 5
  start_seq <- seq(from = 1, to = nrow(initial_df), by = ceiling(nrow(initial_df)/max_i) )
  end_seq   <- c( (start_seq[2:length(start_seq)] - 1), nrow(initial_df))
  for (i in 1:max_i){
    # Extract start and end points for this file
    i_start_row     <- start_seq[i]
    i_end_row       <- end_seq[i]
    # Extract the rows for this file
    i_rows          <- initial_df$unconstrained_tree_iqtree2_call[i_start_row:i_end_row]
    # Make the slurm file
    i_slurm_id_line <- paste0(slurm_id_line, "gene", i)
    i_slurm_txt     <- c(slurm_start_lines,
                         i_slurm_id_line,
                         slurm_middle_lines,
                         i_rows)
    # Save the slurm file
    i_op_file <- paste0(filepath_start, i, ".sh")
    write(i_slurm_txt, i_op_file)
  }
} else {
  initial_df <- read.csv(initial_df_filepath, stringsAsFactors = FALSE)
}

## Extract models from IQ-Tree initial (unconstrained) run
## Check whether dataframe exists
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
  initial_run_df <- read.csv(initial_run_df_filepath, stringsAsFactors = FALSE)
}



#### 6. Estimate constrained gene trees ####
## Construct IQ-Tree command line
constraint_tree_df_filepath <- paste0(output_dir, "genes_005_individualGene_ConstrainedTreeEstimation.csv")
# Open constraint_tree_df if it exists. If it doesn't, then create it.
if (file.exists(constraint_tree_df_filepath) == FALSE){
  # Add columns needed for estimating constrained gene trees
  initial_run_df$CTEN_prefix <- paste0(initial_run_df$gene_id, ".CTEN_tree")
  initial_run_df$PORI_prefix <- paste0(initial_run_df$gene_id, ".PORI_tree")
  initial_run_df$CTEN_PORI_prefix <- paste0(initial_run_df$gene_id, ".CTEN_PORI_tree")
  # Estimate gene trees (constrained)
  # Note: Use model from initial unconstrained run - otherwise cannot compare trees with AU test or MAST
  constraint_tree_df <- as.data.frame(do.call(rbind, lapply(1:nrow(initial_run_df), 
                                                            estimate.constrained.gene.trees.wrapper,
                                                            dataframe = initial_run_df, 
                                                            iqtree2_path = iqtree2, 
                                                            iqtree2_num_threads = iqtree2_num_threads, 
                                                            iqtree_num_bootstraps = iqtree2_num_bootstraps,
                                                            run.UFB = TRUE, 
                                                            estimate.trees = FALSE) ) )
  # Write the constraint tree df to file
  write.csv(constraint_tree_df, file = constraint_tree_df_filepath, row.names = FALSE)
  
  ## Create slurm files (separate all command lines into 5 slurm files)
  # Prepare for extracting specific rows
  filepath_start <- paste0(output_dir, "constraint_tree_run_")
  max_i     <- 10
  start_seq <- seq(from = 1, to = nrow(constraint_tree_df), by = ceiling(nrow(constraint_tree_df)/max_i) )
  end_seq   <- c( (start_seq[2:length(start_seq)] - 1), nrow(constraint_tree_df))
  for (i in 1:max_i){
    # Extract start and end points for this file
    i_start_row     <- start_seq[i]
    i_end_row       <- end_seq[i]
    # Extract the rows for this file
    i_iqtree_calls <- unlist(lapply(i_start_row:i_end_row, function(j){paste(constraint_tree_df$CTEN_iqtree2_call[j], 
                                                                             constraint_tree_df$PORI_iqtree2_call[j], 
                                                                             constraint_tree_df$CTEN_PORI_iqtree2_call[j], 
                                                                             sep = "; ")}))
    # Make the slurm file
    i_slurm_id_line <- paste0(slurm_id_line, "ct_", i)
    i_slurm_txt     <- c(slurm_start_lines,
                         i_slurm_id_line,
                         slurm_middle_lines,
                         i_iqtree_calls)
    # Save the slurm file
    i_op_file <- paste0(filepath_start, i, ".sh")
    write(i_slurm_txt, i_op_file)
  }
} else {
  constraint_tree_df <- read.csv(constraint_tree_df_filepath, stringsAsFactors = FALSE)
}

## Extract results from constrained trees run
## Check whether dataframe exists
constrained_run_df_filepath <- paste0(output_dir, "genes_006_individualGene_ConstrainedTreeResults.csv")
# Open constrained_run_df_ if it exists. If it doesn't, then create it.
if (file.exists(constrained_run_df_filepath) == FALSE){
  # Add constraint tree directory (location of saved constraint tree output files)
  constraint_tree_df$constraint_tree_directory <- paste0(constraint_output_dir, constraint_tree_df$dataset_id, "/")
  # Extract output from iqtree files
  ## NOTE: CURRENTLY ONLY ROWS 1 - 10 COMPLETED FOR TESTING
  constrained_run_df <- as.data.frame(do.call(rbind, lapply(1:nrow(constraint_tree_df), 
                                                            extract.constrained.tree.details, 
                                                            dataframe = constraint_tree_df) ) ) 
  # Write the initial run output df to file
  write.csv(constrained_run_df, file = constrained_run_df_filepath, row.names = FALSE)
} else {
  constrained_run_df <- read.csv(constrained_run_df_filepath, stringsAsFactors = FALSE)
}




#### 7. Estimate sCFs ####
# Calculate sCF: 
# $ iqtree2 -te concat.treefile -s ALN_FILE -m 'best_model' --scfl 100 -pre sCF




#### 8. Calculate AU test for each gene ####
# AU test
# $ iqtree -s gene.fa -te constrained_tree.nex -n 0 -zb 1000 -zw -au -pre AU_test

# MAST
# $ iqtree2 -s gene.fa -m 'best_model+TR' -te constrained_tree.nex -pre MAST







