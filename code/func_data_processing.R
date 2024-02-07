# ancient_ILS/code/func_data_processing.R
## This script includes functions to analyse concordance factors, quartet scores, and phylogenetic trees
# Caitlin Cherryh, 2023

library(stringr)


#### Reformat output dataframes ####
summarise.AU.test.results <- function(dataset_id, au_test_df){
  ## Function to take one dataset and summarise the AU test results in one row
  # Get rows of dataframe that have matching ID
  d_df <- au_test_df[au_test_df$ID == dataset_id,]
  # Extract year of publication
  dataset_year <- as.numeric(str_extract(unique(d_df$dataset), "(\\d)+"))
  # Identify the number of trees
  num_trees <- length(d_df$p_AU)
  # Pad the dataset df with NA (depending on the number of trees)
  au_test_pvalues <- c(d_df$p_AU, rep(NA, (5 - num_trees)) )
  # Create a new row for output
  new_row <- c(d_df$ID[1], d_df$dataset[1], d_df$matrix[1], d_df$gene[1], "AU_test_p_values", 
               au_test_pvalues, dataset_year)
  names(new_row) <- c("gene_id", "dataset", "matrix", "gene", "topology_test", 
                      "tree_1", "tree_2", "tree_3", "tree_4", "tree_5", "year")
  # Return the output
  return(new_row)
}


summarise.eLW <- function(dataset_id, au_test_df){
  ## Function to take one dataset and summarise the AU test results in one row
  # Get rows of dataframe that have matching ID
  d_df <- au_test_df[au_test_df$ID == dataset_id,]
  # Extract year of publication
  dataset_year <- as.numeric(str_extract(unique(d_df$dataset), "(\\d)+"))
  # Identify the number of trees
  num_trees <- length(d_df$c_ELW)
  # Pad the dataset df with NA (depending on the number of trees)
  elw_values <- c(d_df$c_ELW, rep(NA, (5 - num_trees)) )
  # Create a new row for output
  new_row <- c(d_df$ID[1], d_df$dataset[1], d_df$matrix[1], d_df$gene[1], "expected_likelihood_weights", 
               elw_values, dataset_year)
  names(new_row) <- c("gene_id", "dataset", "matrix", "gene", "topology_test", 
                      "tree_1", "tree_2", "tree_3", "tree_4", "tree_5", "year")
  # Return the output
  return(new_row)
}



