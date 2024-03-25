## caitlinch/ancient_ILS/code/03_C60_gene_tree_estimation.R
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
  c60_ml_dir              <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_02_empirical_tree_estimation/01_C60_tree_estimation/"
  gene_dir                <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_01_empirical_genes/"
  gene_output_dir         <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_02_C60_empirical_genes_initial_tree_estimation/"
  output_dir              <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/05_output_files/"
  iqtree2                 <- "iqtree2"
  iqtree2_num_threads     <- 5
} else if (location == "dayhoff"){
  repo_dir                <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  alignment_dir           <- paste0(repo_dir, "data_all/")
  gene_output_dir         <- paste0(repo_dir, "genes/")
  output_dir              <- paste0(repo_dir, "output/")
  iqtree2                 <- paste0(repo_dir, "iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2")
  iqtree2_num_threads     <- 30
}

# Set parameters that are identical for all run locations
iqtree2_num_bootstraps  <- 1000
create.slurm.files      <- TRUE

if (create.slurm.files == TRUE){
  # Set the Slurm lines for your specific server - so we can automate creaying the slurm files
  slurm_start_lines <- c("#!/bin/bash", "#")
  slurm_id_line <- "#SBATCH --job-name="
  slurm_middle_lines <- c("#SBATCH --output=/mnt/data/dayhoff/home/u5348329/ancient_ILS/slurm_files/%j.%x.out", 
                          "#SBATCH --error=/mnt/data/dayhoff/home/u5348329/ancient_ILS/slurm_files/%j.%x.err",
                          "#SBATCH --partition=Standard",
                          "#",
                          "#SBATCH --time=168:00:00 # 7 days",
                          "#SBATCH --ntasks=1",
                          paste0("#SBATCH --cpus-per-task=", iqtree2_num_threads),
                          "#SBATCH --mem=5000 # 5GB required by IQ-Tree",
                          "#",
                          "#SBATCH --mail-user u5348329@anu.edu.au",
                          "#SBATCH --mail-type TIME_LIMIT,FAIL",
                          "",
                          "# Change to the directory containing each gene",
                          "cd /mnt/data/dayhoff/home/u5348329/ancient_ILS/C60_gene_trees/",
                          "",  
                          "# Run IQ-Tree2 command lines")
}



#### 2. Prepare functions, variables and packages ####
# Source required packages
library(ape)

# Source functions and dataset information
source(paste0(repo_dir, "code/func_prepare_trees.R"))
source(paste0(repo_dir, "code/func_empirical_tree_estimation.R"))
source(paste0(repo_dir, "code/func_data_analysis.R"))
source(paste0(repo_dir, "code/func_data_processing.R"))
source(paste0(repo_dir, "code/data_dataset_info.R"))
rm(borowiec2015_list, chang2015_list, dunn2008_list, hejnol2009_list, laumer2018_list, laumer2019_list,
   moroz2014_list, nosenko2013_list, philippe2009_list, philippe2011_list, pick2010_list, ryan2013_list,
   simion2017_list, whelan2015_list, whelan2017_list, all_models, models_list, all_taxa)



#### 3. Extract best model from each of the C60 runs ####
# List all iqtree files from C60 ML trees
c60_op_files  <- list.files(c60_ml_dir)
c60_iqtree    <- paste0(c60_ml_dir, grep("\\.iqtree", c60_op_files, value = T))
c60_log       <- paste0(c60_ml_dir, grep("\\.log", c60_op_files, value = T))
# Extract the C60 alias model with rhas for each of the 12 datasets
gene_models <- unlist(lapply( 1:length(c60_log), function(i){ extract.C60.alias(iqtree_log = c60_log[i],
                                                                                iqtree_file = c60_iqtree[i]) }) )



#### 4. Estimate gene trees with C60 models ####
## Create IQ-Tree command lines for C60 models
# Open csv file with alignment information
gene_file_df <- read.csv(paste0(repo_dir, "output/input_gene_files.csv"), stringsAsFactors = F)
gene_file_df <- gene_file_df[which(gene_file_df$dataset != "Simion2017" & gene_file_df$dataset != "Hejnol2009"), ]
# Update the files for Dayhoff paths
gene_file_df <- update.C60.directory.paths(gene_file_df, location = "dayhoff")
# Add model of evolution to the gene_file_df
gene_file_df$model <- c(rep(gene_models[1], length(which(gene_file_df$dataset_id == "Borowiec2015.Best108"))), 
                        rep(gene_models[2], length(which(gene_file_df$dataset_id == "Chang2015.Chang_AA"))),
                        rep(gene_models[3], length(which(gene_file_df$dataset_id == "Dunn2008.Dunn2008_FixedNames"))),
                        rep(gene_models[4], length(which(gene_file_df$dataset_id == "Laumer2018.Tplx_BUSCOeuk"))),
                        rep(gene_models[5], length(which(gene_file_df$dataset_id == "Moroz2014.ED3d"))),
                        rep(gene_models[6], length(which(gene_file_df$dataset_id == "Nosenko2013.nonribosomal_9187_smatrix"))),
                        rep(gene_models[7], length(which(gene_file_df$dataset_id == "Nosenko2013.ribosomal_14615_smatrix"))),
                        rep(gene_models[8], length(which(gene_file_df$dataset_id == "Philippe2009.Philippe_etal_superalignment_FixedNames"))),
                        rep(gene_models[9], length(which(gene_file_df$dataset_id == "Philippe2011.UPDUNN_MB_FixedNames"))),
                        rep(gene_models[10], length(which(gene_file_df$dataset_id == "Ryan2013.REA_EST_includingXenoturbella"))),
                        rep(gene_models[11], length(which(gene_file_df$dataset_id == "Whelan2015.Dataset10"))),
                        rep(gene_models[12], length(which(gene_file_df$dataset_id == "Whelan2017.Metazoa_Choano_RCFV_strict"))) )
# Create output id
gene_file_df$c60_id <- paste0(gene_file_df$gene_id, ".C60")
# Create iqtree2 command lines
gene_file_df$iqtree2_call <- paste0(gene_file_df$iqtree_path, " -s ", gene_file_df$gene_directory, gene_file_df$gene_file, " -m '", gene_file_df$model, "' ",
                                    "-p ", gene_file_df$gene_tree_directory, gene_file_df$c60_id, " -nt ", gene_file_df$iqtree_num_threads)
# Write out df as csv
write.csv(gene_file_df, file = paste0(output_dir, "c60_gene_tree_estimation_input.csv"))

## Create Slurm files
# Prepare for extracting specific rows
filepath_start <- paste0(output_dir, "c60_gene_trees_")
max_i     <- 10
start_seq <- seq(from = 1, to = nrow(gene_file_df), by = ceiling(nrow(gene_file_df)/max_i) )
end_seq   <- c( (start_seq[2:length(start_seq)] - 1), nrow(gene_file_df))
for (i in 1:max_i){
  # Extract start and end points for this file
  i_start_row     <- start_seq[i]
  i_end_row       <- end_seq[i]
  # Extract the rows for this file
  i_rows          <- gene_file_df$unconstrained_tree_iqtree2_call[i_start_row:i_end_row]
  # Make the slurm file
  i_slurm_id_line <- paste0(slurm_id_line, "c60_", i)
  i_slurm_txt     <- c(slurm_start_lines,
                       i_slurm_id_line,
                       slurm_middle_lines,
                       i_rows)
  # Save the slurm file
  i_op_file <- paste0(filepath_start, i, ".sh")
  write(i_slurm_txt, i_op_file)
}



