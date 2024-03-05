## caitlinch/ancient_ILS/code/04_data_processing.R
# This script estimates maximum likelihood trees under different models of substitution for 14 empirical data sets
# Caitlin Cherryh 2023



#### 1. Input parameters ####
## Specify parameters:
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# output_dir                  <- Directory for saving results files (csv and text)
# hypothesis_tree_dir         <- Location of constrained ML and ASTRAL trees
# hypothesis_gene_dir         <- Location of constrained ML gene trees

repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_dir              <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/05_output_files/"
hypothesis_tree_dir     <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_02_empirical_tree_estimation/03_hypothesis_trees/"
hypothesis_gene_dir     <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_03_empirical_genes_constrained_trees/"


#### 2. Prepare libraries and functions ####
# Open packages
library(stringr)

# Open functions
source(paste0(repo_dir, "code/func_data_processing.R"))

# Extract list of all output files
all_output_files <- list.files(output_dir)



#### 3. Prepare tree and constrained tree output for figures ####
# Open results file with constrained tree results
raw_tree_df <- read.csv(paste0(output_dir, grep("ConstrainedTreeResults", all_output_files, value = T)), stringsAsFactors = FALSE)
# Add year column
raw_tree_df$year <- as.numeric(str_extract(raw_tree_df$dataset, "(\\d)+"))
# Add column with model code
raw_tree_df$ML_model_code <- unlist(lapply(raw_tree_df$unconstrained_tree_best_model, function(x){strsplit(x, "\\+")[[1]][1]}))
# Trim unnecessary columns
tree_df <- raw_tree_df[, c("dataset", "matrix_name", "dataset_id", "gene_name", "gene_id",
                           "ML_model_code", "unconstrained_tree_best_model", "unconstrained_tree_alisim_model",
                           "unconstrained_tree_logl", "unconstrained_tree_unconstrained_logl", "unconstrained_tree_numFreeParams",
                           "unconstrained_tree_BIC", "unconstrained_tree_length", "unconstrained_tree_sumInternalBL",
                           "CTEN_LogL", "CTEN_Unconstrained_LogL", "CTEN_NumFreeParams",
                           "CTEN_BIC", "CTEN_TreeLength", "CTEN_SumInternalBranchLengths",
                           "PORI_LogL", "PORI_Unconstrained_LogL", "PORI_NumFreeParams",
                           "PORI_BIC", "PORI_TreeLength", "PORI_SumInternalBranchLengths",
                           "CTEN_PORI_LogL", "CTEN_PORI_Unconstrained_LogL", "CTEN_PORI_NumFreeParams",
                           "CTEN_PORI_BIC", "CTEN_PORI_TreeLength", "CTEN_PORI_SumInternalBranchLengths")]
# Rename columns for plotting
names(tree_df) <- c("dataset", "matrix_name", "dataset_id", "gene_name", "gene_id",
                    "ML_model_code", "ML_best_model", "ML_alisim_model",
                    "ML_logl", "ML_unconstrained_logl", "ML_numFreeParams",
                    "ML_BIC", "ML_tree_length", "ML_sumInternalBL",
                    "CTEN_LogL", "CTEN_Unconstrained_LogL", "CTEN_NumFreeParams",
                    "CTEN_BIC", "CTEN_TreeLength", "CTEN_SumInternalBL",
                    "PORI_LogL", "PORI_Unconstrained_LogL", "PORI_NumFreeParams",
                    "PORI_BIC", "PORI_TreeLength", "PORI_SumInternalBL",
                    "CTEN_PORI_LogL", "CTEN_PORI_Unconstrained_LogL", "CTEN_PORI_NumFreeParams",
                    "CTEN_PORI_BIC", "CTEN_PORI_TreeLength", "CTEN_PORI_SumInternalBL")
# Write the output
tree_output_file <- paste0(output_dir, "results_tree_likelihood.csv")
write.csv(tree_df, file = tree_output_file, row.names = FALSE)



#### 4. Prepare scf output for figures ####




#### 5. Prepare AU test output for figures ####
# Open raw AU test output
raw_au_df <- read.csv(paste0(output_dir, grep("TreeComparisonResults", all_output_files, value = T)), stringsAsFactors = FALSE)
# Extract AU results per ID
au_test_df <- as.data.frame(do.call(rbind, lapply(unique(raw_au_df$ID), summarise.AU.test.results, raw_au_df)))
# Remove columns consisting only of NA
au_test_df <- Filter(function(x)!all(is.na(x)), au_test_df)
# Sort output by year
au_test_df <- au_test_df[order(au_test_df$year, au_test_df$dataset, au_test_df$matrix, au_test_df$gene_id),]
# Reorder columns
au_test_df <- au_test_df[, c("gene_id", "dataset", "matrix", "gene", "year",
                             "topology_test", "tree_1", "tree_2", "tree_3")]
# Write the output
au_test_file <- paste0(output_dir, "results_gene_AU_test.csv")
write.csv(au_test_df, file = au_test_file, row.names = FALSE)



#### 5. Prepare ELW output for figures ####
# Open raw AU test output
raw_au_df <- read.csv(paste0(output_dir, grep("TreeComparisonResults", all_output_files, value = T)), stringsAsFactors = FALSE)
# Extract AU results per ID
elw_df <- as.data.frame(do.call(rbind, lapply(unique(raw_au_df$ID), summarise.eLW, raw_au_df)))
# Remove columns consisting only of NA
elw_df <- Filter(function(x)!all(is.na(x)), elw_df)
# Sort output by year
elw_df <- elw_df[order(elw_df$year, elw_df$dataset, elw_df$matrix, elw_df$gene_id),]
# Reorder columns
elw_df <- elw_df[, c("gene_id", "dataset", "matrix", "gene", "year",
                     "topology_test", "tree_1", "tree_2", "tree_3")]
# Write the output
elw_file <- paste0(output_dir, "results_gene_elw.csv")
write.csv(elw_df, file = elw_file, row.names = FALSE)



#### 6. Collate hypothesis trees ####
# List all files in the directory
all_ht_files <- list.files(hypothesis_tree_dir)
all_ht_files <- grep("00_", all_ht_files, value = TRUE, invert = TRUE)
# Create ID for each set of trees
dataset_id <- unique(unlist(lapply(strsplit(all_ht_files, "\\."), function(x){paste0(x[1], ".", x[2])})))
dataset_id <- grep("Simion|Hejnol", dataset_id, value = TRUE, invert = TRUE)
# Call each ID and collate the tree files
for (id in dataset_id){
  # Extract output files for this dataset
  id_files <- grep(id, all_ht_files, value = T)
  # Collate ML files
  ml_tree_files <- grep("\\.treefile", grep("_ML_", id_files, value = T), value = T)
  ml_tree_files_ordered <- paste0(hypothesis_tree_dir, c(grep("\\.CTEN_ML", ml_tree_files, value = T), 
                                                         grep("\\.PORI_ML", ml_tree_files, value = T), 
                                                         grep("\\.CTEN_PORI_ML", ml_tree_files, value = T)))
  ml_trees <- c(unlist(lapply(ml_tree_files_ordered, readLines)), "")
  collated_ml_tree_file <- paste0(hypothesis_tree_dir, id, ".ModelFinder.ML.collated_hypothesis_trees.treefile")
  write(ml_trees, file = collated_ml_tree_file)
  # Collate ASTRAL files
  astral_tree_files <- grep("\\.tre", grep("_ASTRAL_", id_files, value = T), value = T)
  astral_tree_files_ordered <- paste0(hypothesis_tree_dir, c(grep("\\.CTEN_ASTRAL", astral_tree_files, value = T), 
                                                             grep("\\.PORI_ASTRAL", astral_tree_files, value = T), 
                                                             grep("\\.CTEN_PORI_ASTRAL", astral_tree_files, value = T)))  
  astral_trees <- c(unlist(lapply(astral_tree_files_ordered, readLines)), "")
  collated_astral_tree_file <- paste0(hypothesis_tree_dir, id, ".ModelFinder.ASTRAL.collated_hypothesis_trees.treefile")
  write(astral_trees, file = collated_astral_tree_file)
  
}



#### 7. Collate constrained gene trees ####
# List all files within the hypothesis_gene_dir
all_datasets <- list.files(hypothesis_gene_dir)
all_datasets <- grep("_hypothesis_trees|all_collated_constrained_gene_trees", all_datasets, value = TRUE, invert = TRUE)
all_datasets <- grep("ModelFinder.CTEN.gene_trees|ModelFinder.PORI.gene_trees|ModelFinder.CTEN_PORI.gene_trees", all_datasets, value = TRUE, invert = TRUE)
# Collate the unconstrained gene trees from each dataset
for (d in all_datasets){
  # List all files in the d folder
  d_dir       <- paste0(hypothesis_gene_dir, d, "/")
  d_all_files <- list.files(d_dir)
  # Extract list of unconstrained gene trees
  d_MFP_tree_files      <- grep("\\.MFP\\.treefile", d_all_files, value = TRUE)
  # Open list of trees
  d_MFP_trees      <- c(unlist(lapply(paste0(d_dir, d_MFP_tree_files), readLines)), "")
  # Create output files for collated constrained gene trees
  d_MFP_collated_file      <- paste0(hypothesis_gene_dir, d, ".ModelFinder.MFP.gene_trees.treefile")
  # Save list of trees
  write(d_MFP_trees, file = d_MFP_collated_file)
}
# Collate the constrained gene trees from each dataset
for (d in all_datasets){
  # List all files in the d folder
  d_dir       <- paste0(hypothesis_gene_dir, d, "/")
  d_all_files <- list.files(d_dir)
  # Extract list of trees for each topology
  d_CTEN_tree_files      <- grep("\\.CTEN_tree.treefile", d_all_files, value = TRUE)
  d_PORI_tree_files      <- grep("\\.PORI_tree.treefile", d_all_files, value = TRUE)
  d_CTEN_PORI_tree_files <- grep("\\.CTEN_PORI_tree.treefile", d_all_files, value = TRUE)
  # Open list of trees
  d_CTEN_trees      <- c(unlist(lapply(paste0(d_dir, d_CTEN_tree_files), readLines)), "")
  d_PORI_trees      <- c(unlist(lapply(paste0(d_dir, d_PORI_tree_files), readLines)), "")
  d_CTEN_PORI_trees <- c(unlist(lapply(paste0(d_dir, d_CTEN_PORI_tree_files), readLines)), "")
  # Create output files for collated constrained gene trees
  d_CTEN_collated_file      <- paste0(hypothesis_gene_dir, d, ".ModelFinder.CTEN.gene_trees.treefile")
  d_PORI_collated_file      <- paste0(hypothesis_gene_dir, d, ".ModelFinder.PORI.gene_trees.treefile")
  d_CTEN_PORI_collated_file <- paste0(hypothesis_gene_dir, d, ".ModelFinder.CTEN_PORI.gene_trees.treefile")
  # Save list of trees
  write(d_CTEN_trees, file = d_CTEN_collated_file)
  write(d_PORI_trees, file = d_PORI_collated_file)
  write(d_CTEN_PORI_trees, file = d_CTEN_PORI_collated_file)
}




#### 8. Collate best model from each gene tree ####


