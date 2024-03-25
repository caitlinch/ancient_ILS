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



#### Extract C60 model parameters ####
extract.C60.alias <- function(iqtree_log, iqtree_file){
  # Function that will extract the best model of sequence evolution or the model of sequence evolution used,
  #   given a .log file
  
  # Check if the log file exists
  if (file.exists(iqtree_log) == TRUE & file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    log_lines <- readLines(iqtree_log)
    iq_lines  <- readLines(iqtree_file)
    
    # Extract best-fit model
    bfm_line    <- grep("Best-fit model", log_lines, value = T)
    best_model  <- strsplit(bfm_line, " ")[[1]][3]
    best_model_chunks <- strsplit(best_model, "\\+")[[1]]
    
    ## Extract alias model
    alias_lines <- grep("alias", log_lines, value = T)
    best_model_alias_line <- alias_lines[which(duplicated(alias_lines) == TRUE)]
    alias_line_split <- strsplit(best_model_alias_line, " is alias for ")[[1]]
    alias_model <- alias_line_split[2]
    alias_split <- strsplit(alias_model, "\\+")[[1]]
    fmix_term_raw <- grep("FMIX", alias_split, value = T)
    
    ## Extract table with rates and weights for C60 model
    table_start <- intersect(intersect(grep("Component", iq_lines), grep("Rate", iq_lines)), 
                             intersect(grep("Weight", iq_lines), grep("Parameters", iq_lines)))
    table_end   <- which(iq_lines == "")[which(which(iq_lines == "") > table_start)][1] - 1
    table_lines <- iq_lines[table_start:table_end]
    table_split <- strsplit(table_lines, " ")
    table_split_filtered <- lapply(table_split[2:length(table_split)], function(x){x[which(x != "")]})
    table_df    <- as.data.frame(do.call(rbind, table_split_filtered))
    names(table_df) <- table_split[[1]][which(table_split[[1]] != "")]
    # Replace any 0 weights with 0.0001 weights
    replace_weights <- which(as.numeric(table_df$Weight) == 0)
    table_df$Weight[replace_weights] <- "0.0001"
    # Create MIX argument sections
    single_model <- strsplit(table_df[1,]$Component, "\\+")[[1]][1]
    table_df$component_formatted <- c(single_model,
                                      paste0(single_model, "+", gsub("F", "", unlist(lapply(strsplit(table_df[2:nrow(table_df),]$Component, "\\+"), function(x){x[[2]]})))) )
    table_df$component_op <- gsub("1.0000", "1", paste0(table_df$component_formatted, ":", table_df$Rate, ":", table_df$Weight))
    table_df$pi_component <- c("", gsub("F", "", unlist(lapply(strsplit(table_df[2:nrow(table_df),]$Component, "\\+"), function(x){x[[2]]}))))
    table_df$pi_component_weight <- gsub("1.0000", "1", paste0(table_df$pi_component, ":", table_df$Rate, ":", table_df$Weight))
    all_pi_components <- grep("pi", table_df$pi_component_weight, value = T)
    # Paste FMIX arguments together
    fmix_term <- paste0("FMIX{", paste(all_pi_components, collapse = ","), "}")
    
    ## Extract rates and weights
    rhas_lines  <- grep("Site proportion and rates", log_lines, value = TRUE)
    rhas_line   <- rhas_lines[length(rhas_lines)]
    rhas_split  <- strsplit(strsplit(rhas_line, ":")[[1]][2], "\\) \\(")[[1]]
    rhas_clean  <- gsub("\\)", "", gsub("\\(", "", gsub(" ", "", rhas_split)))
    rhas        <- paste(rhas_clean, collapse = ",")
    
    ## Extract aa freqs
    pi_A_check          <- grep("pi\\(A\\)", iq_lines)
    if (identical(pi_A_check, integer(0)) == FALSE){
      pi_V_check        <- grep("pi\\(V\\)", iq_lines)
      aa_table_lines    <- iq_lines[pi_A_check:pi_V_check]
      aa_split          <- strsplit(aa_table_lines, "=")
      aa_str_raw        <- unlist(lapply(aa_split, function(x){x[[2]]}))
      aa_str            <- gsub(" ", "", aa_str_raw)
      aa_op             <- paste(aa_str, collapse = ",")
      aa_freqs          <- paste0("{", aa_op, "}")
    } else {
      aa_freqs <- ""
    }
    
    ## Assemble the model back together
    full_model <- paste0(best_model_chunks[[1]], "+", fmix_term, "+", best_model_chunks[[3]], "+", best_model_chunks[[4]], "{", rhas, "}")
    
  } else if (file.exists(iqtree_log) == FALSE){
    # If the iqtree_log doesn't exist, return NA
    full_model = NA
  } # end if (file.exists(iqtree_log) == TRUE){
  
  # Return C60 parameters in model format
  return(full_model)
}
