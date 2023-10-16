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

estimate.empirical.gene.trees <- function(alignment_file, partition_file, 
                                          output_prefix, iqtree2_path, 
                                          iqtree2_num_threads = "AUTO", model = "MFP", 
                                          sitefreqs_file = NA, estimate.trees = FALSE){
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
    command_line <- paste0(iqtree2_path, " -s ", alignment_file, " -S ", partition_file, model_call, " -pre ", output_prefix, " -nt ", iqtree2_num_threads)
  } else {
    # PMSF model - ssfp file is provided
    command_line <- paste0(iqtree2_path, " -s ", alignment_file, " -S ", partition_file, model_call, " -fs ", sitefreqs_file, " -pre ", output_prefix, " -nt ", iqtree2_num_threads)
  }
  # Call IQ-Tree if desired
  if (estimate.trees == TRUE){
    system(command_line)
  }
  # Return IQ-Tree command line
  return(command_line)
}




## Estimating partitioned ML tree
estimate.partitioned.ML.tree.wrapper <- function(row_id, dataframe, iqtree2_path, iqtree2_num_threads = "AUTO", 
                                                 iqtree2_num_bootstraps = 1000, use.gene.tree.partitions = FALSE, estimate.trees = FALSE){
  ## Wrap around estimate.partitioned.ML.tree function
  # Extract parameters for a single dataset
  temp_row <- dataframe[row_id, ]
  # Select partition file
  if (use.gene.tree.partitions == FALSE){
    temp_partition_file <- temp_row$partition_file
  } else if (use.gene.tree.partitions == TRUE){
    temp_partition_file = temp_row$gene_tree_partition_scheme
  }
  # Construct IQ-Tree call
  iqtree_call <- estimate.partitioned.ML.tree (alignment_file = temp_row$alignment_file, partition_file = temp_partition_file, 
                                               output_prefix = temp_row$prefix_ML_tree, iqtree2_path = iqtree2_path,
                                               iqtree2_num_threads = iqtree2_num_threads, iqtree2_num_bootstraps = iqtree2_num_bootstraps,
                                               model = temp_row$best_model, sitefreqs_file = temp_row$PMSF_sitefreq_file, 
                                               estimate.trees = estimate.trees)
  return(iqtree_call)
}

estimate.partitioned.ML.tree <- function(alignment_file, partition_file, 
                                         output_prefix, iqtree2_path, 
                                         iqtree2_num_threads = "AUTO", iqtree2_num_bootstraps = 1000, 
                                         model = "MFP", sitefreqs_file = NA, 
                                         estimate.trees = FALSE){
  ## Estimate a partitioned ML tree in IQ-Tree
  
  # Assemble model command
  if (is.na(model) == FALSE){
    model_call <- paste0(" -m ", model)
  } else if (is.na(model) == TRUE){
    model_call <- ""
  }
  # Collate arguments into command line
  if (is.na(sitefreqs_file) == TRUE){
    # Standard model - included in IQ-Tree
    command_line <- paste0(iqtree2_path, " -s ", alignment_file, " -p ", partition_file, model_call, 
                           " -pre ", output_prefix, " -bb ", iqtree2_num_bootstraps, " -nt ", iqtree2_num_threads)
  } else {
    # PMSF model - ssfp file is provided
    command_line <- paste0(iqtree2_path, " -s ", alignment_file, " -p ", partition_file, model_call, " -fs ", sitefreqs_file, 
                           " -pre ", output_prefix, " -bb ", iqtree2_num_bootstraps, " -nt ", iqtree2_num_threads)
  }
  # Call IQ-Tree if desired
  if (estimate.trees == TRUE){
    system(command_line)
  }
  # Return IQ-Tree command line
  return(command_line)
}




## Estimating constrained partitioned ML tree
estimate.constrained.partitioned.ML.tree.wrapper <- function(row_id, dataframe, iqtree2_path, constraint_tree_hypothesis, 
                                                             iqtree2_num_threads = "AUTO", iqtree2_num_bootstraps = 1000, 
                                                             use.gene.tree.partitions = FALSE, estimate.trees = FALSE){
  ## Wrap around estimate.constrained.partitioned.ML.tree function
  # Select dataset to run by row number
  temp_row <- dataframe[row_id, ]
  # Select partition file
  if (use.gene.tree.partitions == FALSE){
    temp_partition_file <- temp_row$partition_file
  } else if (use.gene.tree.partitions == TRUE){
    temp_partition_file = temp_row$gene_tree_partition_scheme
  }
  # Select constraint tree and prefix
  if (constraint_tree_hypothesis == "CTEN"){
    temp_constraint_tree <- temp_row$constraint_tree_CTEN
    temp_prefix <- temp_row$prefix_CTEN_ML_tree
  } else if (constraint_tree_hypothesis == "PORI"){
    temp_constraint_tree <- temp_row$constraint_tree_PORI
    temp_prefix <- temp_row$prefix_PORI_ML_tree
  }
  # Create IQ-Tree command line
  iqtree_call <- estimate.constrained.partitioned.ML.tree (alignment_file = temp_row$alignment_file, partition_file = temp_partition_file, 
                                                           constraint_tree = temp_constraint_tree, output_prefix = temp_prefix, iqtree2_path = iqtree2_path,
                                                           iqtree2_num_threads = iqtree2_num_threads, iqtree2_num_bootstraps = iqtree2_num_bootstraps, 
                                                           model = temp_row$best_model, sitefreqs_file = temp_row$PMSF_sitefreq_file, 
                                                           estimate.trees = estimate.trees)
  return(iqtree_call)
}

estimate.constrained.partitioned.ML.tree <- function(alignment_file, partition_file, constraint_tree, output_prefix, iqtree2_path, 
                                                     iqtree2_num_threads = "AUTO", iqtree2_num_bootstraps = 1000, 
                                                     model = "MFP", sitefreqs_file = NA, estimate.trees = FALSE){
  ## Estimate a constrained partitioned ML tree in IQ-Tree
  
  # Assemble model command
  if (is.na(model) == FALSE){
    model_call <- paste0(" -m ", model)
  } else if (is.na(model) == TRUE){
    model_call <- ""
  }
  # Collate arguments into command line
  if (is.na(sitefreqs_file) == TRUE){
    # Standard model - included in IQ-Tree
    command_line <- paste0(iqtree2_path, " -s ", alignment_file, " -p ", partition_file, model_call, " -g ", constraint_tree, 
                           " -pre ", output_prefix, " -bb ", iqtree2_num_bootstraps,  " -nt ", iqtree2_num_threads)
  } else {
    # PMSF model - ssfp file is provided
    command_line <- paste0(iqtree2_path, " -s ", alignment_file, " -p ", partition_file, model_call, " -g ", constraint_tree, " -fs ", sitefreqs_file, 
                           " -pre ", output_prefix, " -bb ", iqtree2_num_bootstraps, " -nt ", iqtree2_num_threads)
  }
  # Call IQ-Tree if desired
  if (estimate.trees == TRUE){
    system(command_line)
  }
  # Return IQ-Tree command line
  return(command_line)
}




## Estimating ASTRAL tree
estimate.astral.tree.wrapper <- function(row_id, dataframe, astral_path, estimate.trees = FALSE){
  ## Wrap around estimate.astral.tree function
  temp_row <- dataframe[row_id, ]
  astral_call <- estimate.astral.tree(gene_tree_file = temp_row$gene_tree_treefile, astral_tree_file = temp_row$ASTRAL_tree_treefile, 
                                      astral_log_file = temp_row$ASTRAL_tree_log, astral_path = astral_path, estimate.trees = estimate.trees)
  return(astral_call)
}

estimate.astral.tree <- function(gene_tree_file, astral_tree_file, astral_log_file, astral_path, estimate.trees = FALSE){
  ## Estimate a tree in ASTRAL
  
  # Assemble model command
  command_line <- paste0("java -jar ", astral_path, " -i ", gene_tree_file, " -o ", astral_tree_file, " 2> ", astral_log_file)
  # Call ASTRAL if desired
  if (estimate.trees == TRUE){
    system(command_line)
  }
  # Return ASTRAL command line
  return(command_line)
}




## Estimating constrained ASTRAL tree
estimate.constrained.astral.tree.wrapper <- function(row_id, dataframe, astral_path, constraint_tree_hypothesis, estimate.trees = FALSE){
  ## Wrap around estimate.constrained.astral.tree function
  # Select dataset to run by row number
  temp_row <- dataframe[row_id, ]
  # Select constraint tree and prefix
  if (constraint_tree_hypothesis == "CTEN"){
    temp_constraint_tree <- temp_row$constraint_tree_CTEN
    temp_tree_file <- temp_row$CTEN_ASTRAL_tree_treefile
    temp_log_file <- temp_row$CTEN_ASTRAL_tree_log
  } else if (constraint_tree_hypothesis == "PORI"){
    temp_constraint_tree <- temp_row$constraint_tree_PORI
    temp_tree_file <- temp_row$PORI_ASTRAL_tree_treefile
    temp_log_file <- temp_row$PORI_ASTRAL_tree_log
  }
  # Call ASTRAL
  astral_call <- estimate.constrained.astral.tree(gene_tree_file = temp_row$gene_tree_treefile, astral_tree_file = temp_tree_file, 
                                                  astral_log_file = temp_log_file, constraint_tree = temp_constraint_tree,
                                                  astral_path = astral_path, estimate.trees = estimate.trees)
  return(astral_call)
}

estimate.constrained.astral.tree <- function(gene_tree_file, astral_tree_file, astral_log_file, constraint_tree, astral_path, estimate.trees = FALSE){
  ## Estimate a constrained tree in ASTRAL
  
  # Assemble model command
  command_line <- paste0("java -jar ", astral_path, " -i ", gene_tree_file, " -o ", astral_tree_file, " -j ", constraint_tree, " 2> ", astral_log_file)
  # Call ASTRAL if desired
  if (estimate.trees == TRUE){
    system(command_line)
  }
  # Return ASTRAL command line
  return(command_line)
}



