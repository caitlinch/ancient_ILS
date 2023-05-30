# ancient_ILS/code/func_tree_estimation.R
## This script includes functions to run tree estimation for this project
# Caitlin Cherryh, 2023

#### Wrapper function to estimate trees ####
estimate.trees <- function(row_id, df, call.executable.programs = FALSE){
  # Wrapper function to estimate a tree in IQ-Tree and in ASTRAL
  
  # Extract a single row
  df_row <- df[row_id, ]
  
  # Call IQ-Tree to estimate an ML tree
  iqtree2_ML_output <- run.iqtree2(alignment_path = df_row$output_alignment_file, unique_id = df_row$ID, unique.output.path = FALSE,
                                   iqtree2_path = df_row$iqtree2, iqtree2_num_threads = df_row$iqtree2_num_threads, 
                                   iqtree2_num_ufb = df_row$iqtree2_num_ufb, iqtree2_model = NA, 
                                   use.partition.models = TRUE, partition_file = df_row$output_partition_file, use.model.finder = FALSE, 
                                   call.iqtree2 = call.executable.programs)
  # Call IQ-Tree to estimate gene trees
  iqtree2_gt_output <- estimate.gene.trees(alignment_path = df_row$output_alignment_file, partition_path = df_row$output_partition_file, 
                                           unique_id = df_row$ID, unique.output.path = FALSE, iqtree2_path = df_row$iqtree2,
                                           iqtree2_num_threads = df_row$iqtree2_num_threads, iqtree2_num_ufb = NA,
                                           iqtree2_model = NA, use.partition.models = TRUE, call.iqtree2 = call.executable.programs,
                                           include.bootstraps = FALSE)
  # Call ASTRAL
  astral_output <- run.astral(unique_id = df_row$ID, gene_tree_file = iqtree2_gt_output[["iqtree2_gene_tree_treefile"]],
                              output_directory = df_row$output_folder, astral_path = df_row$ASTRAL,
                              call.ASTRAL = call.executable.programs)
  # Assemble all output - ML tree estimation from IQ-Tree2 and summary coalescent tree estimation from ASTRAL
  output_row <- c(as.character(df_row), iqtree2_ML_output, iqtree2_gt_output, astral_output)
  names(output_row) <- c(names(df_row), names(iqtree2_ML_output), names(iqtree2_gt_output), names(astral_output))
  return(output_row)
}


#### Functions for ASTRAL ####
run.astral <- function(unique_id, gene_tree_file, output_directory, astral_path, call.ASTRAL = FALSE){
  # Function to take one simulated alignment and estimate a single tree in ASTRAL
  
  ## Change working directory to be in the output directory (same directory as gene tree file)
  setwd(output_directory)
  
  # Construct the output files
  astral_output_tree <- paste0(output_directory, unique_id, "_ASTRAL_tree.tre")
  astral_output_log <- paste0(output_directory, unique_id, "_ASTRAL_tree.log")
  filepath_status <- c(file.exists(astral_output_tree), file.exists(astral_output_log))
  filepaths_missing <- FALSE %in% filepath_status
  # Construct the astral command
  astral_command <- paste0("java -jar ", astral_path, " -i ", gene_tree_file, " -o ", astral_output_tree, " 2> ", astral_output_log)
  ## If the ASTRAL command has not been run previously AND call.ASTRAL = TRUE, run the iqtree call
  if (filepaths_missing == TRUE){
    if (call.ASTRAL == TRUE){
      system(astral_command)
    }
  }
  
  # Assemble output
  output_vec        <- c(unique_id, astral_command, call.ASTRAL,
                         gene_tree_file, astral_output_tree, astral_output_log)
  names(output_vec) <- c("gene_tree_ID", "ASTRAL_command", "ASTRAL_command_run",
                         "ASTRAL_input_gene_tree_file", "ASTRAL_tree_treefile", "ASTRAL_tree_log_file")
  return(output_vec)
}


#### Functions for IQ-Tree2 ####
estimate.gene.trees <- function(alignment_path, partition_path, unique_id, unique.output.path = TRUE, iqtree2_path, iqtree2_num_threads = "AUTO", iqtree2_num_ufb = 1000,
                                iqtree2_model = NA, use.partition.models = TRUE, call.iqtree2 = FALSE, include.bootstraps = FALSE){
  # Function to take one simulated alignment and estimate a single tree in IQ-Tree2 using the specified model
  
  ## Change working directory to be in the same directory as the alignment 
  setwd(dirname(alignment_path))
  
  # Determine the prefix to give iqtree2 to create unique output
  if (unique.output.path == TRUE){
    gt_prefix = paste0(unique_id, "-", gsub("\\+", "_", gsub("'", "", iqtree2_model)), "_gene_trees")
  } else if (unique.output.path == FALSE){
    gt_prefix = paste0(unique_id, "_gene_trees")
  }
  # Set output prefix call
  output_prefix_call = paste0("-pre ", gt_prefix)
  
  # Assemble iqtree2 file filepaths
  gt_tree_file <- paste0(dirname(alignment_path), "/", gt_prefix, ".treefile")
  gt_iqtree_file <- paste0(dirname(alignment_path), "/", gt_prefix, ".iqtree")
  gt_log_file <- paste0(dirname(alignment_path), "/", gt_prefix, ".log")
  filepath_status <- c(file.exists(gt_tree_file), file.exists(gt_iqtree_file), file.exists(gt_log_file))
  filepaths_missing <- FALSE %in% filepath_status
  
  ## Assemble the iqtree2 command
  # Assemble the alignment call
  alignment_call <- paste0("-s ", alignment_path)
  # Assemble the partition_call
  partition_call <- paste0("-S ", partition_path)
  
  # Assemble the model section of the command
  if (is.na(iqtree2_model) == TRUE){
    if (use.partition.models == TRUE){
      print("Partitions and models provided - use partition file")
      model_call <- ""
      gt_model <- "Partitioned"
    } else if (use.partition.models == FALSE){
      # Model not provided - iqtree2_model is NA
      print("Model not provided for concatenated alignment: use ModelFinder")
      model_call <- "-m MFP"
      gt_model <- "MFP"
    }
  } else if (is.na(iqtree2_model) == FALSE){
    # Model provided - iqtree2_model is NOT NA
    print(paste0("Model provided for concatenated alignment: ", iqtree2_model))
    model_call <- paste0("-m ", iqtree2_model)
    gt_model <- gsub("'", "", iqtree2_model)
  } 
  # Set number of ultrafast bootstraps
  if (include.bootstraps == TRUE){
    # Include ultrafast bootstraps when estimating gene trees
    # Input: include.bootstraps = TRUE
    if (use.partition.models == TRUE){
      # If using a partition model, add the -bsam GENESITE option to resample partitions then resample sites within partitions
      bootstrap_call <- paste0("-bb ", iqtree2_num_ufb, " -bsam GENESITE")
    } else if (use.partition.models == FALSE){
      # If not using a partition model, use default -bsam option (resample alignment sites within partitions)
      bootstrap_call <- paste0("-bb ", iqtree2_num_ufb) 
    }
  } else if (include.bootstraps == FALSE){
    # Do not estimate ultrafast bootstraps
    # Input: include.bootstraps = FALSE
    bootstrap_call <- ""
  }
  # Assemble the number of threads call
  num_threads_call <- paste0("-nt ", iqtree2_num_threads)
  
  # Assemble the whole iqtree2 call
  iqtree2_call <- paste(iqtree2_path, alignment_call, partition_call, model_call, bootstrap_call, num_threads_call, output_prefix_call, sep = " ")
  
  ## If call.iqtree2 = TRUE, run the iqtree call
  if (filepaths_missing == TRUE){
    if (call.iqtree2 == TRUE){
      system(iqtree2_call)
    }
  }
  
  ## Assemble output
  output_vec <- c(alignment_path, partition_path, gt_prefix, gt_model, 
                  iqtree2_call, call.iqtree2, gt_tree_file, gt_iqtree_file, gt_log_file)
  names(output_vec) <- c("gene_tree_alignment_path", "gene_tree_partition_path", "gene_tree_ID", "iqtree2_gene_tree_model",
                         "iqtree2_gene_tree_command", "iqtree2_gene_tree_command_run", "iqtree2_gene_tree_treefile",
                         "iqtree2_gene_tree_iqtree_file", "iqtree2_gene_tree_log_file")
  return(output_vec)
}


run.iqtree2 <- function(alignment_path, unique_id, unique.output.path = TRUE, iqtree2_path, iqtree2_num_threads = "AUTO", iqtree2_num_ufb = 1000,
                        iqtree2_model = NA, use.partition.models = FALSE, partition_file = NA, use.model.finder = FALSE,
                        call.iqtree2 = FALSE){
  # Function to take one simulated alignment and estimate a single tree in IQ-Tree2 using the specified model
  
  ## Change working directory to be in the same directory as the alignment 
  setwd(dirname(alignment_path))
  
  ## Determine the prefix to give iqtree2 to create unique output
  if (unique.output.path == TRUE){
    output_prefix = paste0(unique_id, "-", gsub("\\+", "_", gsub("'", "", iqtree2_model)), "_ML_tree")
  } else if (unique.output.path == FALSE){
    output_prefix = paste0(unique_id, "_ML_tree")
  }
  # Set output prefix call
  output_prefix_call = paste0("-pre ", output_prefix)
  
  ## Assemble iqtree2 file filepaths
  output_tree_file <- paste0(dirname(alignment_path), "/", output_prefix, ".treefile")
  output_iqtree_file <- paste0(dirname(alignment_path), "/", output_prefix, ".iqtree")
  output_log_file <- paste0(dirname(alignment_path), "/", output_prefix, ".log")
  filepath_status <- c(file.exists(output_tree_file), file.exists(output_iqtree_file), file.exists(output_log_file))
  filepaths_missing <- FALSE %in% filepath_status
  
  ## Assemble the iqtree2 command
  # Assemble the alignment call
  alignment_call <- paste0("-s ", alignment_path)
  
  # Assemble the model section of the command
  if (is.na(iqtree2_model) == FALSE){
    # Model provided
    print(paste0("Model provided for concatenated alignment: ", iqtree2_model))
    model_call <- paste0("-m ", iqtree2_model)
    output_model <- gsub("'", "", iqtree2_model)
  } else if (is.na(iqtree2_model) == TRUE & use.model.finder == TRUE & use.partition.models == FALSE){
    # Gene models not provided but use modelfinder
    print("Model not provided for concatenated alignment: use ModelFinder")
    model_call <- "-m MFP"
    output_model <- "MFP"
  } else if (use.partition.models == TRUE & is.na(partition_file) == FALSE & use.model.finder == TRUE){
    # Use partitions flag on and partition file supplied, model finder to be used
    print("Partition file provided for concatenated alignment: use genes from partition file and models from ModelFinder")
    print(paste0("Partition file: ", partition_file))
    model_call <- paste0("-p ", partition_file, " -m MFP+MERGE")
    output_model <- "Partition_MFP+MERGE"
  } else if (use.partition.models == TRUE & is.na(partition_file) == FALSE){
    # Use partitions flag on and partition file supplied. Do not use model finder
    print(paste0("Partition file provided for concatenated alignment: use genes and model from partition file"))
    print(paste0("Partition file: ", partition_file))
    model_call <- paste0("-p ", partition_file)
    output_model <- "Partition"
  } else {
    # Unclear settings. Default to concatenated alignment with modelfinder
    print("Unclear model parameters. Default to concatenated alignment and use ModelFinder")
    model_call <- "-m MFP"
    output_model <- "UNCLEAR_DEFAULT_MFP"
  }
  # Set number of ultrafast bootstraps
  # Set number of ultrafast bootstraps
  if (use.partition.models == TRUE){
    # If using a partition model, add the -bsam GENESITE option to resample partitions then resample sites within partitions
    bootstrap_call <- paste0("-bb ", iqtree2_num_ufb, " -bsam GENESITE")
  } else if (use.partition.models == FALSE){
    # If not using a partition model, use default -bsam option (resample alignment sites within partitions)
    bootstrap_call <- paste0("-bb ", iqtree2_num_ufb) 
  }
  # Assemble the number of threads call
  num_threads_call <- paste0("-nt ", iqtree2_num_threads)
  
  # Assemble the whole iqtree2 call
  iqtree2_call <- paste(iqtree2_path, alignment_call, model_call, bootstrap_call, num_threads_call, output_prefix_call, sep = " ")
  
  ## If the iqtree command has not been run previously AND call.iqtree2 = TRUE, run the iqtree call
  if (filepaths_missing == TRUE){
    if (call.iqtree2 == TRUE){
      system(iqtree2_call)
    }
  }
  
  ## Assemble output
  output_vec <- c(alignment_path, output_prefix, output_model, iqtree2_call, call.iqtree2,
                  output_tree_file, output_iqtree_file, output_log_file)
  names(output_vec) <- c("alignment_path", "alignment_model_ID", "tree_estimation_model", "iqtree2_command", "iqtree2_command_run",
                         "ML_tree_treefile", "ML_tree_iqtree_file", "ML_tree_log_file")
  return(output_vec)
}



