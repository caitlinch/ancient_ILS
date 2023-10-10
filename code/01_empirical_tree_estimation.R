## caitlinch/metazoan-mixtures/code/01_estimate_all_trees_parallel.R
# This script estimates maximum likelihood trees under different models of substitution for 14 empirical data sets
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## Specify parameters:
# alignment_dir               <- Directory containing alignments for all data sets
#                                   Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                                   E.g. Cherryh2022.alignment1.aa.fa
# output_dir                  <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# iqtree2                     <- Location of IQ-Tree2 stable release
# iqtree_num_threads          <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
# iqtree_mrate                <- Specify a comma separated list of rate heterogeneity types for model selection in IQ-Tree
#                                   We set iqtree_mrate = "E,I,G,I+G,R,I+R"
#                                   See IQ-Tree documentation for more details (http://www.iqtree.org/doc/Command-Reference)
# ml_tree_bootstraps          <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree
# hypothesis_tree_bootstraps  <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree when estimating constrained maximum likelihood trees
# number_parallel_processes   <- The number of simultaneous processes to run at once using mclapply(). 
#                                   If 1, then all processes will run sequentially

## Specify control parameters (all take logical values TRUE or FALSE):
# estimate.ML.trees             <- TRUE to estimate all maximum likelihood trees (each combination of model and alignment). FALSE to skip.

location = "rona"
if (location == "local"){
  alignment_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/"
  repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  
  iqtree2 <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0-MacOSX/bin/iqtree2"
  
  number_parallel_processes <- 1
  
} else if (location == "soma"){
  alignment_dir <- "/data/caitlin/metazoan-mixtures/data_all/"
  output_dir <- "/data/caitlin/metazoan-mixtures/output/"
  repo_dir <- "/data/caitlin/metazoan-mixtures/"
  
  iqtree2 <- "/data/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  
  number_parallel_processes <- 2
  
} else if (location == "dayhoff"){
  alignment_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/data_all/"
  output_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/output/"
  repo_dir <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/"
  
  iqtree2 <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  
  number_parallel_processes <- 1
  
} else if (location == "laptop"){
  alignment_dir <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/01_alignments/"
  output_dir <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/02_output/"
  repo_dir <- "/Users/caitlin/Repositories/metazoan-mixtures/"
  
  iqtree2 <- "/Users/caitlin/Documents/PhD/Ch03_sponge_mixtures/iqtree-2.2.0-MacOSX/bin/iqtree2"
  
  number_parallel_processes <- 1
} else if (location == "rona"){
  alignment_dir <- "/home/caitlin/metazoan-mixtures/data_all/"
  output_dir <- "/home/caitlin/metazoan-mixtures/output/"
  repo_dir <- "/home/caitlin/metazoan-mixtures/"
  
  iqtree2 <- "/home/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  
  number_parallel_processes <- 1
}

# Set parameters that are identical for all run locations
iqtree_mrate <- "E,I,G,I+G,R,I+R"
iqtree_num_threads <- 15
ml_tree_bootstraps <- 1000
hypothesis_tree_bootstraps <- 1000

# Set control parameters
estimate.ML.trees <- FALSE



#### 2. Prepare functions, variables and packages ####
# Open packages
library(parallel)
library(ape)

# Source functions
source(paste0(repo_dir, "code/func_estimate_trees.R"))
source(paste0(repo_dir, "code/func_data_processing.R"))

# Source information about datasets
source(paste0(repo_dir, "code/data_dataset_info.R"))

# Remove the individual dataset lists (only need collated lists) (yes it is a bit cheeky to hard code the removal)
rm(borowiec2015_list, chang2015_list, dunn2008_list, hejnol2009_list, laumer2018_list, laumer2019_list, moroz2014_list, nosenko2013_list, philippe2009_list,
   philippe2011_list, pick2010_list, ryan2013_list, simion2017_list, whelan2015_list, whelan2017_list, models_list, all_taxa)

# Extract the list of all files from the folder containing alignments/models/partition files
all_files <- list.files(alignment_dir)
if (length(all_files) > 0){
  all_files <- paste0(alignment_dir, all_files)
}
# Extract the list of alignments (i.e. files that contain the word "alignment")
all_alignments <- grep("\\.alignment\\.", all_files, value = T)
# Remove Simion and Hejnol alignments - too large to estimate ML trees with all models
all_alignments <- grep("Hejnol", grep("Simion", all_alignments, value = TRUE, invert = TRUE), value = TRUE, invert = TRUE)

# Sort and categorise models of sequence evolution 
all_models <- sort(unique(unlist(lapply(all_models, sort.model.chunks)))) 
# Identify unique matrices within the models
model_components <- sort(unique(unlist(strsplit(all_models, "\\+"))))
# Remove components that do not correspond to a matrix or model (e.g. remove "I", "F", "G4")
model_components <- model_components[!model_components %in% c("F", "FO", "G", "G4", "I", "R", "R4")]
# Create a vector of models to exclude from analysis
#   PMSF: PMSF is more complex to run than other models - requires separate function/multiple steps (see separate PMSF file)
#   CAT: CAT is not a model in IQ-Tree. Use C10:C60 instead.
#   GTR: GTR is a DNA model
#   F81: F81 is a DNA model
exclude_models <- c("PMSF", "CAT", "GTR", "F81", "C10", "C30", "C40", "C50")
# Remove the models to exclude from the list of models
model_components <- model_components[which(!model_components %in% exclude_models)]
# Add the newer models 
model_components <- c(model_components, "'LG+C20'", "'LG+C60'")
# Sort the mode components
model_components <- sort(model_components)
# Move ModelFinder to the end of the list (will take the longest as has to test all possible models, so run last)
model_components <- model_components[!model_components == "ModelFinder"]
model_components <- c(model_components, "ModelFinder")
# Note: partitioning schemes currently not possible in mixture of trees implementation

# Create output folders
# Create a folder for the ml trees
ml_tree_dir <- paste0(output_dir, "maximum_likelihood_trees/")
if (file.exists(ml_tree_dir) == FALSE){dir.create(ml_tree_dir)}
# Create a folder for the constraint trees
c_tree_dir <- paste0(output_dir, "constraint_trees/")
if (file.exists(c_tree_dir) == FALSE){dir.create(c_tree_dir)}

# Create file paths for output files
txt_op_01_01 <- paste0(output_dir, "01_01_maximum_likelihood_iqtree2_calls.txt")
df_op_01_01 <- paste0(output_dir, "01_01_maximum_likelihood_tree_estimation_parameters.tsv")



#### 3. Estimate maximum likelihood trees for each combination of model and dataset ####
# Estimate ML trees (for each combination of alignment and model)
if (estimate.ML.trees == TRUE){
  # Move to the folder for the maximum likelihood trees
  setwd(ml_tree_dir)
  
  # Create a dataframe of combinations of alignments and models
  ml_tree_df <- expand.grid(dataset = unlist(lapply(strsplit(basename(all_alignments), "\\."), "[[", 1)),
                            model_code = model_components,
                            stringsAsFactors = FALSE)
  
  # Copy the models across to input into IQ-Tree
  ml_tree_df$model_mset <- ml_tree_df$model_code
  
  # Tidy up the model codes for the models with multiple input parameters (e.g. C20+LG)
  ml_tree_df$model_code <- gsub("\\+", "_", gsub("'", "", ml_tree_df$model_code))
  
  # Add other columns
  ml_tree_df$model_m <- NA
  ml_tree_df$model_mrate <- iqtree_mrate
  ml_tree_df$sequence_format = unlist(lapply(strsplit(basename(all_alignments), "\\."), "[[", 3))
  ml_tree_df$matrix_name <- unlist(lapply(strsplit(basename(all_alignments), "\\."), "[[", 2))
  ml_tree_df$prefix <- paste0(ml_tree_df$dataset, ".", ml_tree_df$matrix_name, ".", ml_tree_df$model_code)
  ml_tree_df$iqtree_num_threads <- iqtree_num_threads
  ml_tree_df$iqtree_num_bootstraps <- ml_tree_bootstraps
  ml_tree_df$alignment_file <- all_alignments
  ml_tree_df$iqtree_file <- paste0(ml_tree_df$prefix, ".iqtree")
  ml_tree_df$ml_tree_file <- paste0(ml_tree_df$prefix, ".treefile")
  
  
  # Fix model specification for rows with ModelFinder
  ml_tree_df[ml_tree_df$model_mset == "ModelFinder",]$model_m <- "MFP"
  ml_tree_df$model_mset[which(ml_tree_df$model_code == "ModelFinder")] <- NA
  
  # Sort matrix by dataset and matrix
  ml_tree_df <- ml_tree_df[order(ml_tree_df$dataset, ml_tree_df$matrix_name),]
  
  # Create IQ-Tree2 commands
  ml_tree_df$iqtree2_call <- unlist(lapply(1:nrow(ml_tree_df), ml.iqtree.wrapper, iqtree_path = iqtree2, df = ml_tree_df))
  
  # Save dataframe
  write.table(ml_tree_df, file = df_op_01_01, row.names = FALSE, sep = "\t")
  
  # Save list of iqtree2 commands as text file
  write(ml_tree_df$iqtree2_call, file = txt_op_01_01)
  
  # Run IQ-Tree commands to estimate ML trees for each model/matrix combination
  mclapply(ml_tree_df$iqtree2_call, system, mc.cores = number_parallel_processes)
}


