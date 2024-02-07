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
# iqtree_MAST                 <- Location of IQ-Tree2 release with MAST model
# iqtree2_num_bootstraps      <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree
# iqtree2_num_threads         <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"

location = "local"
if (location == "local"){
  repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  alignment_dir           <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/alignments/"
  gene_output_dir         <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_02_empirical_genes_initial_tree_estimation/"
  constraint_output_dir   <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_03_empirical_genes_constrained_trees/"
  scf_output_dir          <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_04_empirical_genes_scf/"
  au_test_output_dir      <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_05_empirical_genes_AU_tests/"
  mast_output_dir         <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_05_empirical_genes_MAST/"
  output_dir              <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/05_output_files/"
  iqtree2                 <- "iqtree2"
  iqtree_MAST             <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.6.hmmster-MacOSX/bin/iqtree2"
  iqtree2_num_threads     <- 5
} else if (location == "dayhoff"){
  repo_dir                <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  alignment_dir           <- paste0(repo_dir, "data_all/")
  gene_output_dir         <- paste0(repo_dir, "genes/")
  constraint_output_dir   <- gene_output_dir
  scf_output_dir          <- paste0(repo_dir, "gene_scf/")
  au_test_output_dir      <- paste0(repo_dir, "gene_au_test/")
  mast_output_dir         <- paste0(repo_dir, "gene_mast/")
  output_dir              <- paste0(repo_dir, "output/")
  iqtree2                 <- paste0(repo_dir, "iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2")
  iqtree_MAST             <- paste0(repo_dir, "iqtree2/iqtree-2.2.6.hmmster-Linux/bin/iqtree2")
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
                        "#SBATCH --time=192:00:00 # 8 days",
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
source(paste0(repo_dir, "code/data_dataset_info.R"))
rm(borowiec2015_list, chang2015_list, dunn2008_list, hejnol2009_list, laumer2018_list, laumer2019_list,
   moroz2014_list, nosenko2013_list, philippe2009_list, philippe2011_list, pick2010_list, ryan2013_list,
   simion2017_list, whelan2015_list, whelan2017_list, all_models, models_list)



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
  
  ## Create slurm files
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
  ## Extract constraint tree output from iqtree files
  # Add constraint tree directory (location of saved constraint tree output files)
  constraint_tree_df$constraint_tree_directory <- paste0(constraint_output_dir, constraint_tree_df$dataset_id, "/")
  # Check that updated paths are for local run
  constraint_tree_df <- update.directory.paths(any_dataframe = constraint_tree_df, location = "local")
  # Extract output from iqtree files
  constrained_run_df <- as.data.frame(do.call(rbind, lapply(1:nrow(constraint_tree_df), 
                                                            extract.constrained.tree.details, 
                                                            dataframe = constraint_tree_df) ) ) 
  
  ## Separate into missing and present genes
  # Identify missing genes
  missing_gene_rows <- union( which(is.na(constrained_run_df$CTEN_treefile) == TRUE),
                              union( which(is.na(constrained_run_df$PORI_treefile) == TRUE),
                                     which(is.na(constrained_run_df$CTEN_PORI_treefile) == TRUE) ) )
  missing_gene_df <- constrained_run_df[missing_gene_rows, ]
  missing_gene_df_filepath <- paste0(output_dir, "genes_006_individualGene_ConstrainedTreeResults_MissingGenes.csv")
  write.csv(missing_gene_df, file = missing_gene_df_filepath, row.names = FALSE)
  # Remove missing genes
  constrained_run_df <- constrained_run_df[setdiff(1:nrow(constrained_run_df), missing_gene_rows), ]
  
  
  ## Collate constrained trees
  # Create file path for the collated constraint trees
  constrained_run_df$collated_constraint_tree <- paste0(constrained_run_df$gene_id, ".hypothesis.treefile")
  # Combine the constrained trees into one file
  for (i in 1:nrow(constrained_run_df)){
    # Extract constrained tree files
    i_ct_files <- paste0(constrained_run_df$constraint_tree_directory[i], c(constrained_run_df$CTEN_treefile[i],
                                                                            constrained_run_df$PORI_treefile[i],
                                                                            constrained_run_df$CTEN_PORI_treefile[i]))
    i_ct_collated_trees <- c(unlist(lapply(i_ct_files,readLines)), "")
    write(i_ct_collated_trees, file = paste0(constraint_output_dir, constrained_run_df$collated_constraint_tree[i]) )
  }
  
  # Write the initial run output df to file
  write.csv(constrained_run_df, file = constrained_run_df_filepath, row.names = FALSE)
  
  
} else {
  constrained_run_df <- read.csv(constrained_run_df_filepath, stringsAsFactors = FALSE)
}



#### 7. Estimate sCFs ####
## Extimate sCFs
# Check whether dataframe exists
scf_call_df_filepath <- paste0(output_dir, "genes_007_individualGene_sCFCommand.csv")
if (file.exists(scf_call_df_filepath) == FALSE){
  # Update paths for dayhoff
  scf_input_df <- update.directory.paths(any_dataframe = constrained_run_df, location = "dayhoff")
  # Create scf command lines
  scf_call_df <- as.data.frame(do.call(rbind, lapply(1:nrow(scf_input_df), 
                                                     estimate.gene.scf.wrapper, 
                                                     dataframe = scf_input_df,
                                                     iqtree2_path = iqtree2_dayhoff, 
                                                     iqtree2_num_threads = iqtree2_num_threads) ) ) 
  # Write the initial run output df to file
  write.csv(scf_call_df, file = scf_call_df_filepath, row.names = FALSE)
} else {
  scf_call_df <- read.csv(scf_call_df_filepath, stringsAsFactors = FALSE)
}

## Create slurm files
# Prepare for extracting specific rows
filepath_start <- paste0(output_dir, "scf_run_")
max_i     <- 10
start_seq <- seq(from = 1, to = nrow(scf_call_df), by = ceiling(nrow(scf_call_df)/max_i) )
end_seq   <- c( (start_seq[2:length(start_seq)] - 1), nrow(scf_call_df))
for (i in 1:max_i){
  # Extract start and end points for this file
  i_start_row     <- start_seq[i]
  i_end_row       <- end_seq[i]
  # Make the slurm file
  i_slurm_id_line <- paste0(slurm_id_line, "scf_", i)
  # Extract the rows for this file
  i_iqtree_calls <- unlist(lapply(i_start_row:i_end_row, function(j){paste(scf_call_df$CTEN_scf_call[j], 
                                                                           scf_call_df$PORI_scf_call[j], 
                                                                           scf_call_df$CTEN_PORI_scf_call[j], 
                                                                           sep = "; ")}))
  # Save the slurm file
  i_slurm_txt     <- c(slurm_start_lines,
                       i_slurm_id_line,
                       slurm_middle_lines,
                       i_iqtree_calls,
                       "")
  i_op_file <- paste0(filepath_start, i, ".sh")
  write(i_slurm_txt, i_op_file)
}

## Extract sCFs
# Prepare filepath
scf_results_df_filepath <- paste0(output_dir, "genes_008_individualGene_sCFResults.csv")
# Check sCF for key branches
scf_results_df <- as.data.frame(do.call(rbind, lapply(1:length(scf_call_df),
                                                      extract.key.scf,
                                                      dataframe = scf_call_df, 
                                                      all_datasets = all_datasets, 
                                                      matrix_taxa = matrix_taxa)))
# Write SCFs
write.csv(scf_results_df, file = scf_results_df_filepath, row.names = FALSE)



#### 8. Calculate AU test and run MAST for each gene ####
# Check if file exists
topology_test_df_filepath <- paste0(output_dir, "genes_009_individualGene_TreeComparisonCommands.csv")
if (file.exists(topology_test_df_filepath) == FALSE){
  # Update paths for dayhoff
  topology_input_df <- update.directory.paths(any_dataframe = constrained_run_df, location = "dayhoff")
  # Create AU test commands and MAST commands for each gene
  topology_test_df <- as.data.frame(do.call(rbind, lapply(1:nrow(topology_input_df), 
                                                          gene.au.test.mast.command, 
                                                          dataframe = topology_input_df,
                                                          iqtree2_path = iqtree2, 
                                                          iqtree2_MAST_path = iqtree2_MAST,
                                                          iqtree2_num_threads = iqtree2_num_threads) ) ) 
  # Write the initial run output df to file
  write.csv(topology_test_df, file = topology_test_df_filepath, row.names = FALSE)
} else {
  topology_test_df <- read.csv(topology_test_df_filepath, stringsAsFactors = FALSE)
}


## Create slurm files
# Prepare for extracting specific rows
filepath_start_au   <- paste0(output_dir, "au_test_")
filepath_start_mast <- paste0(output_dir, "mast_")
max_i     <- 10
start_seq <- seq(from = 1, to = nrow(topology_test_df), by = ceiling(nrow(topology_test_df)/max_i) )
end_seq   <- c( (start_seq[2:length(start_seq)] - 1), nrow(topology_test_df))
for (i in 1:max_i){
  # Extract start and end points for this file
  i_start_row     <- start_seq[i]
  i_end_row       <- end_seq[i]
  # Make the slurm file for the au tests
  i_slurm_id_line_au  <- paste0(slurm_id_line, "au_test_", i)
  i_rows_au           <- topology_test_df$AU_test_iqtree2_command[i_start_row:i_end_row]
  i_slurm_txt_au      <- c(slurm_start_lines,
                           i_slurm_id_line_au,
                           slurm_middle_lines,
                           i_rows_au,
                           "")
  i_op_file_au        <- paste0(filepath_start_au, i, ".sh")
  write(i_slurm_txt_au, i_op_file_au)
  # Make the slurm file for the MAST run
  i_slurm_id_line_mast  <- paste0(slurm_id_line, "mast_", i)
  i_rows_mast           <- topology_test_df$MAST_iqtree2_command[i_start_row:i_end_row]
  i_slurm_txt_mast      <- c(slurm_start_lines,
                             i_slurm_id_line_mast,
                             slurm_middle_lines,
                             i_rows_mast,
                             "")
  i_op_file_mast        <- paste0(filepath_start_mast, i, ".sh")
  write(i_slurm_txt_mast, i_op_file_mast)
}

## Extract topology test results
# Prepare filepaths
topology_results_df_filepath <- paste0(output_dir, "genes_010_individualGene_TreeComparisonResults.csv")
all_au_test_files <- paste0(update.directory.paths(any_dataframe = constrained_run_df, location = "local")$au_test_directory,
                            topology_test_df$AU_test_prefix, ".iqtree")
# Apply function to extract tree topology details
topology_results_df <- as.data.frame(do.call(rbind,lapply(all_au_test_files, 
                                                          extract.tree.topology.test.results)))
# Save dataframe with AU test results 
write.csv(topology_results_df, file = topology_results_df_filepath, row.names = FALSE)

## Commented out - MAST not good option for single genes - was not able to converge in test runs
# ## Extract MAST results
# # Prepare filepaths
# mast_results_df_filepath <- paste0(output_dir, "genes_010_individualGene_MASTResults.csv")
# all_MAST_files <- paste0(update.directory.paths(any_dataframe = constrained_run_df, location = "local")$mast_directory,
#                          topology_test_df$MAST_prefix, ".iqtree")
# # Apply function to extract tree topology details
# mast_results_df <- as.data.frame(do.call(rbind,lapply(all_MAST_files, 
#                                                       extract.tree.weights,
#                                                       trim.output.columns = FALSE)))
# # Save dataframe with AU test results 
# write(mast_results_df, file = mast_results_df_filepath, row.names = FALSE)







