# ancient_ILS/code/func_analysis.R
## This script includes functions investigate and analyses the results of the simulations 
# Caitlin Cherryh, 2023


#### Concordance factors ####
extract.input.concordance.factors <- function(alignment_path, iqtree2_path, iqtree2_num_threads = "AUTO", rename.taxa.for.ms = TRUE, renamed_taxa){
  # Function to return concordance factors (calculated in iqtree), for a simulated tree
  #     gCFs estimated from the base tree supplied to ms and from the set of gene trees estimated in AliSim
  
  ## Get the directory and list of files in that directory
  al_dir <- paste0(dirname(alignment_path),"/")
  al_files <- list.files(al_dir)
  # Change to that dorectory
  setwd(al_dir)
  # Extract the unique id for this alignment
  al_id <- tail(unlist(strsplit(al_dir, "/")),1)
  
  ## Identify the files for the starting tree and the gene trees
  tree_file <- paste0(al_dir, grep("branch_lengths_modified", al_files, value = T))
  gene_tree_file <- paste0(al_dir, grep("ms_gene_trees", al_files, value = T))
  
  ## Rename taxa (if necessary)
  # If taxa names are provided and rename.taxa == TRUE, rename the taxa from the tree_file
  if (rename.taxa == TRUE & length(renamed_taxa) > 0){
    # Open the tree
    t <- read.tree(tree_file)
    # Relabel the tips to have the right number
    t$tip.label <- unlist(lapply(t$tip.label, function(x){renamed_taxa[[x]]}))
    # Remove the "t" from the taxa label - ms labels by number only
    t$tip.label <- gsub("t","",t$tip.label)
    # Save the tree
    tree_file_formatted <- gsub(".treefile", "_renamed.treefile", tree_file)
    write.tree(t, file = tree_file_formatted)
  } else {
    tree_file_formatted <- tree_file
  }
  
  ## Calculate gCFS
  actual_gcf_prefix <- paste0(al_id, "-actual")
  actual_gcf_call <- paste0(iqtree2_path ," -t ", tree_file_formatted, " --gcf ", gene_tree_file, " --prefix ", actual_gcf_prefix)
  system(actual_gcf_call)
  
  ## Extract and return gCFs
  # Open the gcf table
  al_files <- paste0(al_dir, list.files(al_dir))
  gcf_files <- grep(actual_gcf_prefix, al_files, value = T)
  gcf_stat_file = grep("cf.stat", gcf_files, value = T)
  gcf_table <- read.table(gcf_stat_file, header = TRUE, sep = "\t")
  # Return the gCF table
  return(gcf_table)
}

extract.output.concordance.factors <- function(alignment_path, iqtree2_path, iqtree2_num_threads = "AUTO", iqtree2_num_ufb = 1000,
                                               iqtree2_model = NA){
  # Function to return concordance factors (calculated in iqtree), for a simulated alignment
  #     gCFs estimated from the ML tree estimated from the partitioned simulated alignment and gene trees estimated from each partition
  
  # Run iqtree to estimate gene concordance factors
  gcf_files <- iqtree2.concordance.factors(alignment_path, iqtree2_path, iqtree2_num_threads, iqtree2_num_ufb, iqtree2_model)
  # Retrieve gcf from output table
  gcf_stat_file <- op_vector[["gCF_table_file"]]
  gcf_table <- read.table(gcf_stat_file, header = TRUE, sep = "\t")
  # Return the gCF table
  return(gcf_table)
  
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
  # Assemble output files
  gcf_tree_file <- paste0(al_dir, gcf_tree_prefix, ".cf.tree")
  gcf_branch_file <- paste0(al_dir, gcf_tree_prefix, ".cf.branch")
  gcf_table_file <- paste0(al_dir, gcf_tree_prefix, ".cf.stat")
  # Assemble output
  op_vector <- c(alignment_path, gcf_tree_file, gcf_branch_file, gcf_table_file)
  names(op_vector) <- c("alignment_path", "gCF_tree_file", "gCF_branch_file", "gCF_table_file")
  return(op_vector)
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
