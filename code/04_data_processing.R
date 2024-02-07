## caitlinch/ancient_ILS/code/04_data_processing.R
# This script estimates maximum likelihood trees under different models of substitution for 14 empirical data sets
# Caitlin Cherryh 2023



#### 1. Input parameters ####
## Specify parameters:
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# output_dir                  <- Directory for saving results files (csv and text)

repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_dir              <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/05_output_files/"



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

