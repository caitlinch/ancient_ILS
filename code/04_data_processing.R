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
# Open functions
source(paste0(repo_dir, "code/func_data_processing.R"))

# Extract list of all output files
all_output_files <- list.files(output_dir)



#### 3. Prepare tree and constrained tree output for figures ####




#### 4. Prepare scf output for figures ####




#### 5. Prepare AU test output for figures ####
# Open raw AU test output
raw_au_df <- read.csv(paste0(output_dir, grep("TreeComparisonResults", all_output_files, value = T)), stringsAsFactors = FALSE)
# Extract AU results per ID
au_test_df <- as.data.frame(do.call(rbind, lapply(unique(raw_au_df$ID), summarise.AU.test.results, raw_au_df)))
# Remove columns consisting only of NA
au_test_df <- Filter(function(x)!all(is.na(x)), au_test_df)
# Sort output by year
au_test_df <- au_test_df[order(au_test_df$year, au_test_df$dataset, au_test_df$matrix),]
# Reorder columns
au_test_df <- au_test_df[, c("ID", "dataset", "matrix", "gene", "year",
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
elw_df <- elw_df[order(elw_df$year, elw_df$dataset, elw_df$matrix),]
# Reorder columns
elw_df <- elw_df[, c("ID", "dataset", "matrix", "gene", "year",
                     "topology_test", "tree_1", "tree_2", "tree_3")]
# Write the output
elw_file <- paste0(output_dir, "results_gene_elw.csv")
write.csv(elw_df, file = elw_file, row.names = FALSE)

