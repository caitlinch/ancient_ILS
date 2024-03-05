## caitlinch/ancient_ILS/code/02_empirical_concordance_factors.R
# This script estimates site, gene, and quartet concordance factors for empirical gene trees and constrained hypothesis trees.
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

## Specify control parameters (all take logical values TRUE or FALSE):
# create.output.filepaths     <- Open the dataset_df and add output file paths for estimating concordance factors/quartet scores: T/F
# estimate.scf.gcf            <- Run command lines to estimate sCF and gCF in IQ-Tree2: T/F
# estimate.qs                 <- Run command lines to estimate quartet scores in ASTRAL: T/F
# output.command.lines        <- Output dataframe with command lines to estimate quartet scores and concordance factors: T/F

location = "local"
if (location == "local"){
  repo_dir              <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  alignment_dir         <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_01_empirical_data/"
  output_dir            <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_02_empirical_tree_estimation/"
  output_csv_dir        <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/05_output_files/"
  iqtree2               <- "iqtree2"
  iqtree2_num_threads   <- "AUTO"
  astral                <- "/Users/caitlincherryh/Documents/Executables/ASTRAL-5.7.8-master/Astral/astral.5.7.8.jar"
  
} else if (location == "dayhoff"){
  repo_dir              <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  alignment_dir         <- paste0(repo_dir, "data_all/")
  output_dir            <-  paste0(repo_dir, "output/")
  output_csv_dir        <- paste0(repo_dir, "output/")
  iqtree2               <- paste0(repo_dir, "iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2")
  iqtree2_num_threads   <- 20
  astral                <- paste0(repo_dir, "astral/Astral/astral.5.7.8.jar")
}

# Set control parameters
control_parameters      <- list(create.output.filepaths = FALSE,
                                estimate.scf.gcf = FALSE,
                                estimate.qs = FALSE,
                                output.command.lines = FALSE,
                                constrained.concordance.analysis = TRUE)



#### 2. Prepare functions, variables and packages ####
# Source functions and dataset information
source(paste0(repo_dir, "code/func_empirical_tree_estimation.R"))

# Create a new folder for concordance factor and quartet scores
cf_dir <- paste0(output_dir, "concordance_factors/")
if (dir.exists(cf_dir) == FALSE){dir.create(cf_dir)}

# Open the required dataframe with dataset information
precf_dataset_df_file <- paste0(output_dir, "dataset_preparation_concordance_factors.csv")
if (control_parameters$create.output.filepaths == TRUE | file.exists(precf_dataset_df_file) == FALSE){
  # Open the dataframe containing information about each alignment
  all_files <- list.files(output_dir)
  dataset_df_file <- paste0(output_dir, grep("dataset", grep("command_lines", all_files, value = T), value = T))
  dataset_df <- read.csv(dataset_df_file)
  # Add columns to the dataset_df
  dataset_df$prefix_gcf.scf <- paste0(dataset_df$dataset, ".", dataset_df$matrix)
  dataset_df$CTEN_qcf_tree <- paste0(cf_dir, dataset_df$dataset, ".", dataset_df$matrix, ".CTEN.qs.tre")
  dataset_df$CTEN_qcf_log <- paste0(cf_dir, dataset_df$dataset, ".", dataset_df$matrix, ".CTEN.qs.log")
  dataset_df$PORI_qcf_tree <- paste0(cf_dir, dataset_df$dataset, ".", dataset_df$matrix, ".PORI.qs.tre")
  dataset_df$PORI_qcf_log <- paste0(cf_dir, dataset_df$dataset, ".", dataset_df$matrix, ".PORI.qs.log")
  # Output dataset
  write.csv(dataset_df, file = precf_dataset_df_file, row.names = F)
} else {
  dataset_df <- read.csv(precf_dataset_df_file, stringsAsFactors = FALSE)
}



#### 3. Calculate site and gene concordance factors ####
if (control_parameters$estimate.scf.gcf == TRUE){
  # For CTEN-sister
  CTEN_cf_commands <- unlist(lapply(1:nrow(dataset_df), estimate.gcf.scf.wrapper, 
                                    dataframe = dataset_df, constraint_tree_hypothesis = "CTEN", 
                                    iqtree2_path = iqtree2, iqtree2_num_threads = 10, 
                                    estimate.trees = FALSE))
  # For PORI-sister
  PORI_cf_commands <- unlist(lapply(1:nrow(dataset_df), estimate.gcf.scf.wrapper, 
                                    dataframe = dataset_df, constraint_tree_hypothesis = "PORI", 
                                    iqtree2_path = iqtree2, iqtree2_num_threads = 10, 
                                    estimate.trees = FALSE))
  # Add columns to dataset
  dataset_df$CTEN_gCF_command_line <- CTEN_cf_commands[c(T,F)]
  dataset_df$CTEN_sCF_command_line <- CTEN_cf_commands[c(F,T)]
  dataset_df$PORI_gCF_command_line <- PORI_cf_commands[c(T,F)]
  dataset_df$PORI_sCF_command_line <- PORI_cf_commands[c(F,T)]
}



#### 4. Calculate quartet scores ####
if (control_parameters$estimate.qs == TRUE){
  # For CTEN-sister
  dataset_df$CTEN_quartet_score_command_line <- unlist(lapply(1:nrow(dataset_df), estimate.quartet.scores.wrapper, 
                                                              dataframe = dataset_df, constraint_tree_hypothesis = "CTEN",
                                                              astral_path = astral, estimate.trees = FALSE))
  # For PORI-sister
  dataset_df$PORI_quartet_score_command_line <- unlist(lapply(1:nrow(dataset_df), estimate.quartet.scores.wrapper, 
                                                              dataframe = dataset_df, constraint_tree_hypothesis = "PORI",
                                                              astral_path = astral, estimate.trees = FALSE))
}



#### 5. Save dataframe with new command lines  ####
cf_dataset_df_trimmed_file <- paste0(output_dir, "dataset_estimate_concordance_factors_ModelFinder.csv")
cf_dataset_df_file <- paste0(output_dir, "dataset_estimate_concordance_factors.csv")
# Output command lines
if (control_parameters$output.command.lines == TRUE){
  write.csv(dataset_df, file = cf_dataset_df_file, row.names = F)
  write.csv(dataset_df[13:24,], file = cf_dataset_df_trimmed_file, row.names = F)
}



#### 6. Run concordance factors for constrained species trees with constrained and unconstrained gene trees  ####
# prefix = "constrained_cf"
if (control_parameters$constrained.concordance.analysis == TRUE){
  
}

# Open csv file
in_df <- read.csv(paste0(output_csv_dir, "constrained_cf_input.csv"), stringsAsFactors = FALSE)




