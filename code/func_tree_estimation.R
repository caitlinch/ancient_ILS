# ancient_ILS/code/func_tree_estimation.R
## This script includes functions to run tree estimation for this project
# Caitlin Cherryh, 2023


#### Functions for IQ-Tree2 ####
estimate.one.tree <- function(alignment_path, iqtree2_path, iqtree2_num_threads = "AUTO", iqtree2_num_ufb = 1000,
                              iqtree2_model = NA, use.partitions = FALSE, partition_file = NA, use.model.finder = FALSE){
  # Function to take one simulated alignment and estimate a single tree using the specified model
  
  ## Assemble the iqtree2 command
  # Asesemble the model section of the command
  if (is.na(gene_models) == FALSE){
    # Model provided
    model_call <- paste0("-m ", iqtree2_model)
  } else if (is.na(gene_models) == TRUE & use.model.finder == TRUE){
    # Gene models not provided but use modelfinder
    model_call <- "-m MFP"
  } else if (use.partitions == TRUE & is.na(partition_file) == FALSE & use.model.finder == TRUE){
    # Use partitions flag on and partition file supplied, model finder to be used
    model_call <- paste0("-p ", partition_file, " -m MFP+MERGE")
  } else if (use.partitions == TRUE & is.na(partition_file) == FALSE){
    # Use partitions flag on and partition file supplied. Do not use model finder
    model_call <- paste0("-p ", partition_file)
  }
  
  model_command <- ""
}


