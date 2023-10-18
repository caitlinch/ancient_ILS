## caitlinch/ancient_ILS/code/01_empirical_tree_estimation.R
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
  astral_constrained  <- "/Users/caitlincherryh/Documents/Executables/ASTRAL-constrained/astral.5.6.9.jar"
  
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
dataset_df_file <- paste0(output_dir, "dataset_parameters.csv")
if (control_parameters$prepare.parameters == TRUE){
  ## Prepare partition files
  # Get all files from the partition folder (and remove any files with "00" as the beginning of the file name)
  all_files <- list.files(alignment_dir)
  # Reformat partition files
  partition_files <- grep("ModelFinder|PMSF", grep("partitions.nex|partitions.txt|smatrix.txt|partitions.part", all_files, value = TRUE), value = TRUE, invert = TRUE)
  reformatted_partition_files <- unlist(lapply(paste0(alignment_dir, partition_files), create.partition.nexus))
  
  ## Assemble nice dataframe of empirical datasets
  # Extract all alignments
  all_als <- grep("ModelFinder|PMSF", grep("\\.nex|\\.phy|\\.phylip|\\.fasta|\\.fa|\\.fas", grep("alignment", all_files, value = T), value = T), value = TRUE, invert = TRUE)
  # Start dataframe
  dataset_df <- data.frame(year = as.numeric(str_extract(all_als, "(\\d)+")),
                           dataset = unlist(lapply(1:length(all_als), function(i){strsplit(all_als, "\\.")[[i]][1]})),
                           matrix = unlist(lapply(1:length(all_als), function(i){strsplit(all_als, "\\.")[[i]][2]})), 
                           alignment_file = all_als,
                           partition_file = basename(reformatted_partition_files) )
  # Copy dataset_df - need one copy for PMSF and one for ModelFinder gene trees
  dataset_df <- rbind(dataset_df, dataset_df)
  # Remove datasets without PMSF best model
  dataset_df <- dataset_df[dataset_df$dataset != "Hejnol2009"  & dataset_df$dataset != "Laumer2019" & dataset_df$dataset != "Simion2017", ]
  rownames(dataset_df) <- 1:nrow(dataset_df)
  # Add model code
  dataset_df$model_code <- c(rep("PMSF", nrow(dataset_df)/2), rep("ModelFinder", nrow(dataset_df)/2))
  # Add PMSF site freq files
  dataset_df$PMSF_sitefreq_file <- c(grep("PMSF", all_files, value = T)[1:3], NA, grep("PMSF", all_files, value = T)[4:11], rep(NA, 12))
  ## TODO: REMOVE HARD CODING BY CREATING PMSF RUNS FOR LAUMER 2018 BUSCO ALIGNMENT
  ##dataset_df$PMSF_sitefreq_file <- grep("PMSF", all_files, value = T)
  ##dataset_df$PMSF_sitefreq_file[which(dataset_df$model_code == "ModelFinder")] <- NA
  
  # Add best model (best models for ModelFinder extracted from ModelFinder runs)
  dataset_df$best_model <- dataset_df$model_code
  dataset_df$best_model[grep("PMSF_LG_C60", dataset_df$PMSF_sitefreq_file)] <- "LG+C60+F+R4"
  dataset_df$best_model[grep("PMSF_C60", dataset_df$PMSF_sitefreq_file)] <- "C60+F+R4"
  dataset_df$best_model[grep("ModelFinder", dataset_df$model_code)] <- c("'Q.insect+R6'", "'Q.insect+R8'", "'LG+R7'", "'LG+F+R7'", 
                                                                         "'LG+F+I+R5'", "'LG+F+I+R6'", "'Q.insect+R7'", "'Q.insect+I+R6'",
                                                                         "'Q.insect+R7'", "'Q.insect+R6'", "'LG+F+R8'", "'Q.insect+R8'")
  # Add columns with output prefixes
  dataset_df$prefix_gene_trees    <- paste0(dataset_df$dataset, ".", dataset_df$matrix, ".", dataset_df$model_code, ".", "gene_trees")
  dataset_df$prefix_ML_tree       <- paste0(dataset_df$dataset, ".", dataset_df$matrix, ".", dataset_df$model_code, ".", "ML_tree")
  dataset_df$prefix_ASTRAL_tree   <- paste0(dataset_df$dataset, ".", dataset_df$matrix, ".", dataset_df$model_code, ".", "ASTRAL_tree")
  dataset_df$prefix_CTEN_ASTRAL_tree     <- paste0(dataset_df$dataset, ".", dataset_df$matrix, ".", dataset_df$model_code, ".", "CTEN_ASTRAL_tree")
  dataset_df$prefix_PORI_ASTRAL_tree     <- paste0(dataset_df$dataset, ".", dataset_df$matrix, ".", dataset_df$model_code, ".", "PORI_ASTRAL_tree")
  dataset_df$prefix_CTEN_ML_tree     <- paste0(dataset_df$dataset, ".", dataset_df$matrix, ".", dataset_df$model_code, ".", "CTEN_ML_tree")
  dataset_df$prefix_PORI_ML_tree     <- paste0(dataset_df$dataset, ".", dataset_df$matrix, ".", dataset_df$model_code, ".", "PORI_ML_tree")
  # Add files for gene and ML tree estimation
  dataset_df$gene_tree_treefile     <- paste0(dataset_df$prefix_gene_trees, ".treefile")
  dataset_df$gene_tree_partition_scheme <- paste0(dataset_df$prefix_gene_trees, ".best_scheme.nex")
  dataset_df$ML_tree_treefile       <- paste0(dataset_df$prefix_ML_tree, ".treefile")
  dataset_df$ASTRAL_tree_treefile   <- paste0(dataset_df$prefix_ASTRAL_tree, ".tre")
  dataset_df$ASTRAL_tree_log        <- paste0(dataset_df$prefix_ASTRAL_tree, ".log")
  # Add files for constrained ML and ASTRAL tree estimation
  dataset_df$CTEN_ML_tree_treefile      <- paste0(dataset_df$prefix_CTEN_ML_tree, ".tre")
  dataset_df$CTEN_ASTRAL_tree_treefile  <- paste0(dataset_df$prefix_CTEN_ASTRAL_tree, ".tre")
  dataset_df$CTEN_ASTRAL_tree_log       <- paste0(dataset_df$prefix_CTEN_ASTRAL_tree, ".log")
  dataset_df$PORI_ML_tree_treefile      <- paste0(dataset_df$prefix_PORI_ML_tree, ".tre")
  dataset_df$PORI_ASTRAL_tree_treefile  <- paste0(dataset_df$prefix_PORI_ASTRAL_tree, ".tre")
  dataset_df$PORI_ASTRAL_tree_log       <- paste0(dataset_df$prefix_PORI_ASTRAL_tree, ".log")
  
  ## Generate constraint trees
  # Create the directory for the constraint trees
  constraint_dir <- paste0(output_dir, "constraint_trees/")
  # Generate new dataframe to create one set of partition trees per dataset (as the dataset_df has duplicates - one for MFP and one for PMSF)
  constraint_tree_df <- dataset_df[1:(nrow(dataset_df)/2), ]
  # Generate constraint tree files
  all_constraint_tree_files <- lapply(1:nrow(constraint_tree_df), constraint.tree.wrapper, 
                                      output_directory = constraint_dir,  dataset_info = all_datasets, 
                                      matrix_taxa_info = matrix_taxa, constraint_tree_df = constraint_tree_df, 
                                      alignment_taxa_df = alignment_taxa_df, force.update.constraint.trees = TRUE)
  dataset_df$constraint_tree_CTEN <- rep(basename(unlist(lapply(all_constraint_tree_files, function(x){x[[1]]}))), 2)
  dataset_df$constraint_tree_PORI <- rep(basename(unlist(lapply(all_constraint_tree_files, function(x){x[[2]]}))), 2)
  dataset_df$constraint_tree_CTEN.PORI <- rep(basename(unlist(lapply(all_constraint_tree_files, function(x){x[[3]]}))), 2)
  
  ## Save dataframe
  write.csv(dataset_df, file = dataset_df_file, row.names = FALSE)
} else {
  dataset_df <- read.csv(dataset_df_file, stringsAsFactors = FALSE)
}



#### 4. Prepare file paths for server runs ####
## Prepare columns for full run by adding full file paths for the run location
# Full paths for data files
dataset_df$alignment_file <- paste0(alignment_dir, basename(dataset_df$alignment_file))
dataset_df$partition_file <- paste0(alignment_dir, basename(dataset_df$partition_file))
dataset_df$PMSF_sitefreq_file[which(is.na(dataset_df$PMSF_sitefreq_file) == FALSE)] <- paste0(alignment_dir, basename(dataset_df$PMSF_sitefreq_file[which(is.na(dataset_df$PMSF_sitefreq_file) == FALSE)]))
# Full paths for ML output files
dataset_df$gene_tree_treefile   <- paste0(output_dir, "tree_estimation/", basename(dataset_df$gene_tree_treefile))
dataset_df$gene_tree_partition_scheme <- gsub("\\.treefile", "\\.best_scheme\\.nex", dataset_df$gene_tree_treefile)
dataset_df$ML_tree_treefile     <- paste0(output_dir, "tree_estimation/", basename(dataset_df$ML_tree_treefile))
dataset_df$ASTRAL_tree_treefile <- paste0(output_dir, "tree_estimation/", basename(dataset_df$ASTRAL_tree_treefile))
dataset_df$ASTRAL_tree_treefile <- paste0(output_dir, "tree_estimation/", basename(dataset_df$ASTRAL_tree_treefile))
# Full paths for constraint tree files
dataset_df$constraint_tree_CTEN <- paste0(output_dir, "constraint_trees/", basename(dataset_df$constraint_tree_CTEN))
dataset_df$constraint_tree_PORI <- paste0(output_dir, "constraint_trees/", basename(dataset_df$constraint_tree_PORI))
dataset_df$constraint_tree_CTEN.PORI <- paste0(output_dir, "constraint_trees/", basename(dataset_df$constraint_tree_CTEN.PORI))
# Full paths for constrained tree estimation output files
dataset_df$CTEN_ML_tree_treefile      <- paste0(output_dir, "hypothesis_trees/", basename(dataset_df$CTEN_ML_tree_treefile))
dataset_df$CTEN_ASTRAL_tree_treefile  <- paste0(output_dir, "hypothesis_trees/", basename(dataset_df$CTEN_ASTRAL_tree_treefile))
dataset_df$CTEN_ASTRAL_tree_log       <- paste0(output_dir, "hypothesis_trees/", basename(dataset_df$CTEN_ASTRAL_tree_log))
dataset_df$PORI_ML_tree_treefile      <- paste0(output_dir, "hypothesis_trees/", basename(dataset_df$PORI_ML_tree_treefile))
dataset_df$PORI_ASTRAL_tree_treefile  <- paste0(output_dir, "hypothesis_trees/", basename(dataset_df$PORI_ASTRAL_tree_treefile))
dataset_df$PORI_ASTRAL_tree_log       <- paste0(output_dir, "hypothesis_trees/", basename(dataset_df$PORI_ASTRAL_tree_log))
# Save with specific parameters for each dataset
dataset_server_path_file <- paste0(output_dir, "dataset_parameters_", location, "_paths.csv")
write.csv(dataset_df, file = dataset_server_path_file, row.names = FALSE)



#### 5. Estimate gene trees (ModelFinder and PMSF) ####
# Create command lines to estimate gene trees
dataset_df$gene_tree_command_line <- unlist(lapply(1:nrow(dataset_df), 
                                                   estimate.empirical.gene.trees.wrapper, 
                                                   dataframe = dataset_df, 
                                                   iqtree2_path = iqtree2, 
                                                   iqtree2_num_threads = iqtree2_num_threads, 
                                                   estimate.trees = FALSE))



#### 6. Estimate partitioned ML and ASTRAL trees (ModelFinder and PMSF) ####
# Create command lines to estimate partitioned ML trees
dataset_df$ML_tree_command_line <- unlist(lapply(1:nrow(dataset_df), 
                                                 estimate.partitioned.ML.tree.wrapper, 
                                                 dataframe = dataset_df, 
                                                 iqtree2_path = iqtree2, 
                                                 iqtree2_num_threads = iqtree2_num_threads,
                                                 iqtree2_num_bootstraps = iqtree2_num_bootstraps,
                                                 use.gene.tree.partitions = TRUE,
                                                 estimate.trees = FALSE))
# Create command lines to estimate ASTRAL trees
dataset_df$ASTRAL_command_line <- unlist(lapply(1:nrow(dataset_df), 
                                                estimate.astral.tree.wrapper, 
                                                dataframe = dataset_df, 
                                                astral_path = astral,
                                                estimate.trees = FALSE))



#### 7. Estimate constrained ML trees (ModelFinder and PMSF) ####
# Create command lines for constrained tree estimation in IQ-Tree
dataset_df$CTEN_ML_command_line <- unlist(lapply(1:nrow(dataset_df), estimate.constrained.partitioned.ML.tree.wrapper, 
                                                 dataframe = dataset_df, iqtree2_path = iqtree2, constraint_tree_hypothesis = "CTEN", 
                                                 iqtree2_num_threads = iqtree2_num_threads, iqtree2_num_bootstraps = iqtree2_num_bootstraps,
                                                 use.gene.tree.partitions = TRUE, estimate.trees = FALSE))
dataset_df$PORI_ML_command_line <- unlist(lapply(1:nrow(dataset_df), estimate.constrained.partitioned.ML.tree.wrapper, 
                                                 dataframe = dataset_df, iqtree2_path = iqtree2, constraint_tree_hypothesis = "PORI", 
                                                 iqtree2_num_threads = iqtree2_num_threads, iqtree2_num_bootstraps = iqtree2_num_bootstraps,
                                                 use.gene.tree.partitions = TRUE, estimate.trees = FALSE))

# Create command lines for constrained tree estimation in ASTRAL
dataset_df$CTEN_ASTRAL_command_line <- unlist(lapply(1:nrow(dataset_df), estimate.constrained.astral.tree.wrapper, 
                                                     dataframe = dataset_df, astral_path = astral_constrained, 
                                                     constraint_tree_hypothesis = "CTEN", estimate.trees = FALSE))
dataset_df$PORI_ASTRAL_command_line <- unlist(lapply(1:nrow(dataset_df), estimate.constrained.astral.tree.wrapper, 
                                                     dataframe = dataset_df, astral_path = astral_constrained, 
                                                     constraint_tree_hypothesis = "PORI", estimate.trees = FALSE))
# Save dataframe
command_line_csv_file <- paste0(output_dir, "dataset_", location, "_command_lines.csv")
write.csv(dataset_df, file = command_line_csv_file)



#### 8. Output text files for each type of command ####
# Gene trees
write(dataset_df$gene_tree_command_line, file = paste0(output_dir, "empirical_commandLines_gene_trees.txt"))
# ML trees
write(dataset_df$ML_tree_command_line, file = paste0(output_dir, "empirical_commandLines_ML_trees.txt"))
# ASTRAL trees
write(dataset_df$ASTRAL_command_line, file = paste0(output_dir, "empirical_commandLines_ASTRAL_trees.txt"))
# ML constrained trees
write(c(dataset_df$CTEN_ML_command_line, dataset_df$PORI_ML_command_line), file = paste0(output_dir, "empirical_commandLines_constrained_ML.txt"))
# ASTRAL constrained trees
write(c(dataset_df$CTEN_ASTRAL_command_line, dataset_df$PORI_ASTRAL_command_line), file = paste0(output_dir, "empirical_commandLines_constrained_ASTRAL.txt"))


