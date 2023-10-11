# ancient_ILS/code/func_empirical_tree_estimation.R
## This script includes functions to run tree estimation for empirical datasets 
# Caitlin Cherryh, 2023

## Estimating gene trees
estimate.empirical.gene.trees.wrapper <- function(row_id, dataframe, iqtree2_path, iqtree2_num_threads = "AUTO", estimate.trees = FALSE){
  ## Wrap around estimate.empirical.gene.trees function
  temp_row <- dataframe[row_id, ]
  iqtree_call <- estimate.empirical.gene.trees(alignment_file = temp_row$alignment_file, partition_file = temp_row$partition_file, 
                                               output_prefix = temp_row$prefix_gene_trees, iqtree2_path = iqtree2_path,
                                               iqtree2_num_threads = iqtree2_num_threads, model = temp_row$best_model, 
                                               sitefreqs_file = temp_row$PMSF_sitefreq_file, estimate.trees = estimate.trees)
  return(iqtree_call)
}

estimate.empirical.gene.trees <- function(alignment_file, partition_file, output_prefix, iqtree2_path, iqtree2_num_threads = "AUTO", model = NA, sitefreqs_file = NA, estimate.trees = FALSE){
  ## Estimate a set of gene trees using a partition file
  ##   Command line: iqtree2 -s ALN_FILE -S PARTITION_FILE --prefix loci -T AUTO
  
  # Assemble model command
  if (is.na(model) == FALSE){
    model_call <- paste0(" -m ", model)
  } else if (is.na(model) == TRUE){
    model_call <- ""
  }
  # Collate arguments into command line
  if (is.na(sitefreqs_file) == TRUE){
    # Standard model - included in IQ-Tree
    command_line <- paste0(iqtree2_path, " -s ", alignment_file, " -S ", partition_file, model_call, " --prefix ", output_prefix, " -nt ", iqtree2_num_threads)
  } else {
    # PMSF model - ssfp file is provided
    command_line <- paste0(iqtree2_path, " -s ", alignment_file, " -S ", partition_file, model_call, " -fs ", sitefreqs_file, " --prefix ", output_prefix, " -nt ", iqtree2_num_threads)
  }
  # Call IQ-Tree if desired
  if (estimate.trees == TRUE){
    system(command_line)
  }
  # Return IQ-Tree command line
  return(command_line)
}
