# ancient_ILS/code/func_tree_estimation.R
## This script includes functions to run tree estimation for this project
# Caitlin Cherryh, 2023


#### Functions for IQ-Tree2 ####
estimate.one.tree <- function(alignment_path, unique.output.path = TRUE, iqtree2_path, iqtree2_num_threads = "AUTO", iqtree2_num_ufb = 1000,
                              iqtree2_model = NA, use.partitions = FALSE, partition_file = NA, use.model.finder = FALSE,
                              run.iqtree2 = FALSE){
  # Function to take one simulated alignment and estimate a single tree using the specified model
  
  ## Change working directory to be in the same directory as the alignment 
  setwd(dirname(alignment_path))
  
  ## Assemble the iqtree2 command
  # Assemble the alignment call
  alignment_call <- paste0("-s ", alignment_path)
  
  # Assemble the model section of the command
  if (is.na(iqtree2_model) == FALSE){
    # Model provided
    print(paste0("Model provided for concatenated alignment: ", iqtree2_model))
    model_call <- paste0("-m ", iqtree2_model)
    output_model <- gsub("'", "", iqtree2_model)
  } else if (is.na(iqtree2_model) == TRUE & use.model.finder == TRUE & use.partitions == FALSE){
    # Gene models not provided but use modelfinder
    print("Model not provided for concatenated alignment: use ModelFinder")
    model_call <- "-m MFP"
    output_model <- "MFP"
  } else if (use.partitions == TRUE & is.na(partition_file) == FALSE & use.model.finder == TRUE){
    # Use partitions flag on and partition file supplied, model finder to be used
    print("Partition file provided for concatenated alignment: use genes from partition file and models from ModelFinder")
    print(paste0("Partition file: ", partition_file))
    model_call <- paste0("-p ", partition_file, " -m MFP+MERGE")
    output_model <- "Partition_MFP+MERGE"
  } else if (use.partitions == TRUE & is.na(partition_file) == FALSE){
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
  bootstrap_call <- paste0("-bb ", iqtree2_num_ufb)
  # Assemble the number of threads call
  num_threads_call <- paste0("-nt ", iqtree2_num_threads)
  # Determine the prefix to give iqtree2 to create unique output
  if (unique.output.path == TRUE){
    output_prefix = paste0(gsub("_alignment.fa", "", basename(alignment_path)), "-", gsub("\\+", "_", gsub("'", "", iqtree2_model)))
  } else if (unique.output.path == FALSE){
    output_prefix = gsub("_alignment.fa", "", basename(alignment_path))
  }
  # Set output prefix call
  output_prefix_call = paste0("-pre ", output_prefix)
  # Assemble the whole iqtree2 call
  iqtree2_call <- paste(iqtree2_path, alignment_call, model_call, bootstrap_call, num_threads_call, output_prefix_call, sep = " ")
  
  ## If run.iqtree2 = TRUE, run the iqtree call
  if (run.iqtree2 == TRUE){
    system(iqtree2_call)
  }
  
  ## Return the output files
  # Assemble iqtree2 file filepaths
  output_tree_file <- paste0(dirname(alignment_path), output_prefix, ".treefile")
  output_iqtree_file <- paste0(dirname(alignment_path), output_prefix, ".iqtree")
  output_log_file <- paste0(dirname(alignment_path), output_prefix, ".log")
  
  # Assemble output
  output_vec <- c(alignment_path, output_prefix, output_model, iqtree2_call, run.iqtree2, output_tree_file, output_iqtree_file, output_log_file)
  names(output_vec) <- c("alignment_path", "alignment_model_ID", "tree_estimation_model", "iqtree2_command", "iqtree2_command_run", "ML_tree_treefile", "ML_tree_iqtree_file", "ML_tree_log_file")
  return(output_vec)
}


iqtree2.concordance.factors <- function(alignment_path, iqtree2_path, iqtree2_num_threads = "AUTO", iqtree2_num_ufb = 1000,
                                        iqtree2_model = NA){
  # Function to take a simulated alignment and estimate gCF from it using iqtree2
  
  ## Get the directory and list of files in that directory
  al_dir <- paste0(dirname(alignment_path),"/")
  al_files <- list.files(al_dir)
  # Change to that dorectory
  setwd(al_dir)
  # Extract the unique id for this alignment
  al_id <- tail(unlist(strsplit(al_dir, "/")),1)
  
  ## Create a gene partition file with no models
  # Find and open the alisim partition file
  alisim_partition_file <- paste0(al_dir, grep("log", grep("partition", al_files, value = TRUE), value = TRUE, invert = TRUE))
  # Generate the gcf partition file
  gcf_partition_file <- generate.gcf.partition.file(alisim_partition_file)
  
  ## Inferring species tree
  # Create model call
  if (is.na(iqtree2_model) == TRUE){
    model_call <- " -m MFP+MERGE "
  } else if (is.na(iqtree2_model) == FALSE){
    model_call <- paste0(" -m ", iqtree2_model, " ")
  }
  # Create IQ-Tree call
  species_tree_prefix <- paste0(al_id, "-concat")
  species_tree_call <- paste0(iqtree2_path, " -s ", alignment_path, " -p ", gcf_partition_file, model_call, " --prefix ", species_tree_prefix, " -bb ",  iqtree2_num_ufbm, " -nt ",iqtree2_num_threads)
  system(species_tree_call)
  
  ## Inferring gene/locus trees  
  # Create model call
  if (is.na(iqtree2_model) == TRUE){
    model_call <- " -m MFP+MERGE "
  } else if (is.na(iqtree2_model) == FALSE){
    model_call <- paste0(" -m ", iqtree2_model, " ")
  }
  # Create IQ-Tree call
  gene_tree_prefix <- paste0(al_id, "-gene_trees")
  gene_tree_call <- paste0(iqtree2_path, " -s ", alignment_path, " -S ", gcf_partition_file, model_call, " --prefix ", gene_tree_prefix, " -bb ",  iqtree2_num_ufbm, " -nt ",iqtree2_num_threads)
  system(gene_tree_call)
  
  ## Calculating gene concordance factors
  gcf_tree_prefix <- paste0(al_id, "-concord")
  gcf_call <- paste0(iqtree2_path ," -t ", species_tree_prefix, ".treefile --gcf ", gene_tree_prefix, ".treefile --prefix ", gcf_tree_prefix)
  system(gcf_call)
  
  ## Return output
  
}


#### Functions for partition files ####
generate.gcf.partition.file <- function(partition_file){
  # Quick function to create a gCF partition file from an Alisim partition file 
  
  # Open partition file
  alisim_partitions <- readLines(partition_file)
  # Extract gene names
  gene_names <- grep("gene", unlist(strsplit(grep("charset", alisim_partitions, value = T), " ")), value = T)
  # Make a new charpartition section
  new_charpartition <- paste0("\tcharpartition genes = ", paste(gene_names, collapse = ", "), ";")
  # Replace the charpartition section
  line_ind <- grep("charpartition", alisim_partitions)
  alisim_partitions[line_ind] <- new_charpartition
  # Create file path
  p_dir <- paste0(dirname(partition_file), "/")
  p_id <- tail(unlist(strsplit(p_dir, "/")),1)
  gcf_partitions <- paste0(p_dir, p_id, "_gCF_partitions.nexus")
  # Save the new partition file 
  write(alisim_partitions, file = gcf_partitions)
  # Return the gcf partition file path
  return(gcf_partitions)
}
  
  