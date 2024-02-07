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
# Extract list of all output files
all_output_files <- list.files(output_dir)



#### 3. Prepare tree and constrained tree output for figures ####




#### 4. Prepare scf output for figures ####




#### 5. Prepare AU test output for figures ####
raw_au_df <- read.csv(paste0(output_dir, grep("TreeComparisonResults", all_output_files, value = T)), stringsAsFactors = FALSE)



