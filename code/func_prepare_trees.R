# ancient_ILS/code/func_prepare_trees.R
## This script includes functions for formatting and preparing hypothesis trees
# Caitlin Cherryh 2024

library(castor)
library(ape)
library(phytools)
library(phangorn) # read.aa
library(seqinr) # write.fasta

#### Extracting information from phylogenetic trees ####
internal.branch.proportion <- function(tree){
  ## Small function to take any tree and determine the total sum of branch length, the sum of internal branch length, and the ratio between the two
  # Sum length of all branches
  total_branch_length <- sum(tree$edge.length)
  # Identify internal branches and external branches
  num_tips <- Ntip(tree)
  internal_branch_inds <- which(tree$edge[,2] > num_tips)
  external_branch_inds <- which(tree$edge[,2] <= num_tips)
  # Determine ratio of internal to external branches
  sum_internal_branch_lengths <- sum(tree$edge.length[internal_branch_inds])
  percent_internal_branch_length <- (sum_internal_branch_lengths/total_branch_length)*100
  # Collate and return output
  op <- c(total_branch_length, sum_internal_branch_lengths, percent_internal_branch_length)
  return(op)
}


root.to.tip.distances <- function(tree, root.tree = TRUE, outgroup = NA){
  ## Function to take a single taxa, root it, and calculate the root to tip distance
  
  # Root the tree at the outgroup (if requested)
  if (root.tree == TRUE & is.character(outgroup) == TRUE){
    # Check which outgroups are present in the tree
    present_outgroup <- outgroup[which(outgroup %in% gene_trees[[4]]$tip.label)]
    # Root tree
    tree <- root(tree, outgroup = present_outgroup, resolve.root = TRUE)
  }
  # Identify root node
  root_node <- find_root(tree)
  # Calculate the root-to-tip distance for the specified taxa
  r2t_distances <- get_all_distances_to_root(tree)
  # Identify terminal branches
  terminal_branch_edges <- tree$edge[which(tree$edge[ , 2] < (Ntip(tree) + 1) ), ]
  branch_order <- terminal_branch_edges[ , 2]
  # Identify root-to-tip distances for terminal branches
  terminal_r2t_distances <- r2t_distances[branch_order]
  # Reorder tips to match branch_order
  terminal_r2t_tips <- tree$tip.label[branch_order]
  # Return the root-to-tip distance as a named vector
  names(terminal_r2t_distances) <- terminal_r2t_tips
  return(terminal_r2t_distances)
}


add.missing.distances <- function(r2t_d, tree_tips){
  ## Function to take the root to tip distances, and add any missing distances
  
  # Check which tips are missing
  missing_tips <- tree_tips[which(! tree_tips %in% names(r2t_d))]
  # Add all the missing tips as NAs
  new_dists <- c(as.numeric(r2t_d), rep(NA, length(missing_tips)))
  names(new_dists) <- c(names(r2t_d), missing_tips)
  # Sort to match order of the tree tips
  new_dists <- new_dists[tree_tips]
  # Return the distances
  return(new_dists)
}


extract.clade.length <- function(gene_tree, clade_tips, root_tips = "Drosophila_melanogaster"){
  ## Extract the length and depth of a specified clade (given the tips in that clade) from a given gene tree
  ##    This function assumes that each clade is monophyletic - if the clade is not, it returns NA 
  ##        and does not calculate length/depth for that clade
  ##    This function roots each gene tree so between-gene-tree-comparisons are possible
  
  # Root gene tree
  rooted_gene_tree <- root(gene_tree, outgroup = root_tips)
  # Remove tips not in this gene tree
  present_clade_tips <- clade_tips[clade_tips %in% rooted_gene_tree$tip.label]
  # Extract clade if >1 tip is present
  if (length(present_clade_tips) > 1){
    # Identify number of tips present
    num_tips <- length(present_clade_tips)
    # Extract clade
    clade_node <- getMRCA(rooted_gene_tree, present_clade_tips)
    clade <- extract.clade(rooted_gene_tree, clade_node)
    if (Ntip(clade) > length(present_clade_tips)){
      # Tips in clade do not have monophyletic relationship
      # Cannot extract clade details - return NA
      clade_relationship = "paraphyletic"
      clade_length <- NA
      clade_branch_length <- NA
      max_branching_time <- NA
    } else {
      # Tips in clade have monophyletic relationship
      # Can extract clade details - calculate distances
      clade_relationship = "monophyletic"
      # Determine length of clade
      clade_length <- sum(clade$edge.length)
      # Determine depth of clade (won't be proper because tree is not ultrametric)
      ult_clade <- force.ultrametric(clade, method = "extend")
      max_branching_time <- max(branching.times(ult_clade))
      # Determine length of branch leading to clade
      edge_rows <- which(rooted_gene_tree$edge[,2] == clade_node)
      clade_branch_length <- rooted_gene_tree$edge.length[edge_rows]
    }
  } else if (length(present_clade_tips) == 1){
    # Only one tip is present - cannot extract information about a clade
    # Assemble output
    num_tips <- 1
    clade_relationship <- "1_taxa_present"
    clade_length <- NA
    clade_branch_length <- NA
    max_branching_time <- NA
  } else if (length(present_clade_tips) == 0){
    # Only one tip is present - cannot extract information about a clade
    num_tips <- 0
    clade_relationship <- "0_taxa_present"
    clade_length <- NA
    clade_branch_length <- NA
    max_branching_time <- NA
  } else {
    # Insufficient information
    num_tips <- NA
    clade_relationship <- NA
    clade_length <- NA
    clade_branch_length <- NA
    max_branching_time <- NA
  }
  # Assemble output
  op <- c(num_tips, clade_relationship, clade_length, 
          clade_branch_length, max_branching_time, paste(root_tips, collapse = ","))
  names(op) <- c("num_clade_tips", "clade_relationship", "clade_length", 
                 "branch_length_to_clade", "max_branching_time_for_ultrametric_clade", "outgroup")
  # Return results
  return(op)
}


extract.clade.monophyly <- function(gene_tree, clade_tips, drop_tips, remove.specified.tips = FALSE){
  ## Determine whether the set of taxa provided is monophyletic
  
  # Remove tips not in this gene tree
  present_clade_tips <- clade_tips[clade_tips %in% gene_tree$tip.label]
  # If any tips are specified to be removed, remove them
  if ( ((NA %in% drop_tips) == FALSE) & (remove.specified.tips == TRUE) ){
    gene_tree <- drop.tip(gene_tree, tip = drop_tips, trim.internal = TRUE)
  }
  # Check monophyly of clade if more than one tip is present
  if (length(present_clade_tips) > 1){
    # Check relationship between tips
    check_monophyly <- is.monophyletic(gene_tree, tips = present_clade_tips)
    if (check_monophyly == TRUE){
      clade_relationship = "Monophyletic"
    } else if (check_monophyly == FALSE){
      clade_relationship = "Paraphyletic"
    }
  } else if (length(present_clade_tips) == 1){
    # Only one taxa - can't be monophyletic/paraphyletic
    clade_relationship = "One_taxon"
  } else if (length(present_clade_tips) == 0){
    # No taxa - can't be monophyletic/paraphyletic
    clade_relationship = "Zero_taxa"
  }
  # Specify number of tips in this clade for the given gene tree
  num_tips = length(present_clade_tips)
  # Assemble output
  op <- c(num_tips, clade_relationship)
  names(op) <- c("num_clade_tips", "clade_relationship")
  # Return results
  return(op)
}




#### Creating constraint trees ####
format.constraint.tree.clade <- function(clade){
  ## Function to take in a vector of species and return a nicely formatted character object to paste into a constraint tree
  # Check how many taxa are in the clade
  clade_size <- length(clade)
  
  # Format the clade
  if (clade_size == 1){
    clade_formatted = clade
  } else if (clade_size > 1){
    clade_formatted = paste0("(", paste(clade, collapse = ", "), ")")
  } else if (clade_size < 1){
    clade_formatted <- ""
  }
  
  # Return the nicely formatted clade
  return(clade_formatted)
}




#### Reformatting partition lines ####
partition.one.model.line <- function(model_line){
  ## Small function to take one line from the models file and return in partition format
  # Break line into two chunks
  line_break <- strsplit(model_line, "\\=")
  model_chunk <- line_break[[1]][1]
  genes_chunk <- line_break[[1]][2]
  # Extract the model
  model <- strsplit(model_chunk, ",")[[1]][1]
  name_subset <- strsplit(model_chunk, ",")[[1]][2]
  # Return the extracted sections
  return(c(model, name_subset, genes_chunk))
}




#### Reading gene lengths from partition files ####
create.partition.nexus <- function(partition_file, return.charpartition = FALSE){
  ## Take partition format and reformat into a nice nexus partition file
  
  # Check file type and extract gene ranges
  if (grepl("partitions.nex|partition.nex|partitions.nexus|partition.nexus", basename(partition_file)) == TRUE){
    if (grepl("Philippe2009|Philippe2011", basename(partition_file)) == TRUE){
      gene_df <- gene.lengths.nexus.philippe(partition_file, return.lengths = FALSE)
    } else {
      gene_df <- gene.lengths.nexus(partition_file, return.lengths = FALSE)
    }
  } else if (grepl("partitions.txt|partition.txt|partitions.raxml|partition.raxml|partitions.part|partition.part", basename(partition_file)) == TRUE){
    gene_df <- gene.lengths.raxml(partition_file, return.lengths = FALSE)
  } else if (grepl("smatrix.txt|partition_smatrix", basename(partition_file)) == TRUE){
    gene_df <- gene.lengths.smatrix(partition_file, return.lengths = FALSE)
  } 
  # Assemble sections for partition file
  new_p_start <- c("#nexus", "begin sets;")
  new_p_end <- c("end;", "")
  new_p_charsets <- paste0("\tcharset ", gene_df$gene_name, " = ", gene_df$gene_start, " - ", gene_df$gene_end, ";")
  if (return.charpartition == FALSE){
    # Return charsets only
    # Assemble partition file
    p_text <- c(new_p_start, new_p_charsets, new_p_end)
  } else if (return.charpartition == TRUE){
    # Return all charsets with one charpartition "all", which includes all charsets
    # Create charpartition
    new_p_charpartition <- paste0("\tcharpartition all = ", paste(gene_df$gene_name, collapse = ","), ";")
    # Assemble partition file
    p_text <- c(new_p_start, new_p_charsets, new_p_charpartition, new_p_end)
  }
  # Create new output path file name
  # Split path at "."
  p_path_split <- strsplit(basename(partition_file), "\\.")[[1]]
  new_p_path <- paste0(dirname(partition_file), "/", paste(p_path_split[1:(length(p_path_split) - 1)], collapse = "."), "_formatted.nex")
  new_p_path <- gsub("_smatrix", "", new_p_path)
  # Write out partitions in nexus format to the directory containing the partition file
  write(p_text, file = new_p_path)
  # Return the new partition file
  return(new_p_path)
}


gene.lengths.nexus <- function(partition_file, return.lengths = TRUE){
  ## Extract all gene lengths from a nexus partition file
  # Open partition file
  lines <- readLines(partition_file)
  # Extract all lines with a charset
  charset_lines <- grep("charset", lines, ignore.case = TRUE, value = TRUE)
  # Split the charset lines at the "="
  gene_lines <- unlist(lapply(strsplit(charset_lines, "="), function(x){x[2]}))
  # Split the genes into chunks by breaking at the commas ","
  gene_chunks <- unlist(strsplit(gene_lines, ","))
  # Format the gene chunks nicely
  gene_chunks_nospace <- gsub(" ", "", gene_chunks)
  gene_chunks_noend <- gsub(";", "", gene_chunks_nospace)
  # Get start and end of each gene
  gene_start <- as.numeric(unlist(lapply(strsplit(gene_chunks_noend, "-"), function(x){x[1]})))
  gene_end <- as.numeric(unlist(lapply(strsplit(gene_chunks_noend, "-"), function(x){x[2]})))
  # Create the gene names
  gene_names <- sprintf("gene_%04d", 1:length(gene_start))
  # Create a nice little dataframe for the genes
  gene_df <- data.frame(gene_name = gene_names, gene_range = gene_chunks_noend, gene_start = gene_start, gene_end = gene_end)
  # Calculate the gene length
  gene_df$gene_length <- gene_end - (gene_start - 1) # subtract one from gene_start to count the starting site in the gene length
  # Set output
  if (return.lengths == TRUE){
    output <- gene_df$gene_length
  } else {
    output <- gene_df
  }
  # Return the gene length or gene dataframe
  return(output)
}


gene.lengths.nexus.philippe <- function(partition_file, return.lengths = TRUE){
  ## Extract all gene lengths from a nexus partition file
  # Open partition file
  lines <- readLines(partition_file)
  # Extraction depends on dataset
  if (grepl("Philippe2009", partition_file) == TRUE){
    # Extract gene start, end, and range
    gene_lines <- grep("charset", lines, value = T)
    gene_line_split <- strsplit(gene_lines, "=")[[1]][2]
    gene_line_split <- gsub(";", "", gene_line_split)
    gene_range <- gsub(" ", "", strsplit(gene_line_split, ",")[[1]])
    gene_start <- as.numeric(unlist(lapply(strsplit(gene_range, "-"), function(x){x[1]})))
    gene_end <- as.numeric(unlist(lapply(strsplit(gene_range, "-"), function(x){x[2]})))
  } else if (grepl("Philippe2011", partition_file) == TRUE){
    # Extract gene start, end, and range
    gene_lines <- grep("partition part", lines, value = T)
    gene_line_split <- strsplit(gene_lines, ":")[[1]][2]
    gene_line_split <- gsub(";", "", gene_line_split)
    gene_range <- gsub(" ", "", strsplit(gene_line_split, ",")[[1]])
    gene_start <- as.numeric(unlist(lapply(strsplit(gene_range, "-"), function(x){x[1]})))
    gene_end <- as.numeric(unlist(lapply(strsplit(gene_range, "-"), function(x){x[2]})))
  }
  # Create the gene names
  gene_names <- sprintf("gene_%04d", 1:length(gene_start))
  # Create a nice little dataframe for the genes
  gene_df <- data.frame(gene_name = gene_names, gene_range = gene_range, gene_start = gene_start, gene_end = gene_end)
  # Calculate the gene length
  gene_df$gene_length <- gene_end - (gene_start - 1) # subtract one from gene_start to count the starting site in the gene length
  # Set output
  if (return.lengths == TRUE){
    output <- gene_df$gene_length
  } else {
    output <- gene_df
  }
  # Return the gene length or gene dataframe
  return(output)
}


gene.lengths.raxml <- function(partition_file, return.lengths = TRUE){
  ## Extract all gene lengths from a RAxML partition file 
  # Open partition file
  lines <- readLines(partition_file)
  # Different processing for different datasets (subtle differences in partition file depending on original dataset)
  if (grepl("Borowiec2015|Simion2017", basename(partition_file)) == TRUE){
    # Extract all lines with an equals sign (these lines will define a gene)
    eq_lines <- grep("\\=", lines, ignore.case = TRUE, value = TRUE)
    # Split the charset lines at the "="
    gene_lines <- unlist(lapply(strsplit(eq_lines, "="), function(x){x[2]}))
    # Split the genes into chunks by breaking at the commas ","
    gene_chunks <- unlist(strsplit(gene_lines, ","))
    # Format the gene chunks nicely
    gene_range <- gsub(" ", "", gene_chunks)
    # Remove semicolons for Simion 2017
    if (grepl("Simion2017", basename(partition_file)) == TRUE){
      gene_range <- gsub(";", "", gene_range)
    }
    # Get start and end of each gene
    gene_start <- as.numeric(unlist(lapply(strsplit(gene_range, "-"), function(x){x[1]})))
    gene_end <- as.numeric(unlist(lapply(strsplit(gene_range, "-"), function(x){x[2]})))
  } else if (grepl("Chang2015", basename(partition_file)) == TRUE){
    # Extract all lines with an equals sign (these lines will define a gene)
    eq_lines <- grep("\\=", lines, ignore.case = TRUE, value = TRUE)
    # Split the charset lines at the "="
    gene_lines <- unlist(lapply(strsplit(eq_lines, "="), function(x){x[2]}))
    # Split the genes into chunks by breaking at the commas ","
    gene_chunks <-  unlist(strsplit(gene_lines, ","))
    # Format the gene chunks nicely
    gene_chunks_nospace <- gsub(" ", "", gene_chunks)
    # Get start and end of each gene
    gene_start_unordered <- as.numeric(unlist(lapply(strsplit(gene_chunks_nospace, "-"), function(x){x[1]})))
    gene_end_unordered <- as.numeric(unlist(lapply(strsplit(gene_chunks_nospace, "-"), function(x){x[2]})))
    # Reorder gene chunks
    order_df <- data.frame(gene_range = gene_chunks_nospace, gene_start = gene_start_unordered, gene_end = gene_end_unordered)
    order_df <- order_df[order(order_df$gene_start), ]
    # Extract gene start, end and range
    gene_start <- order_df$gene_start
    gene_end <- order_df$gene_end
    gene_range <- order_df$gene_range
  } else if (grepl("Laumer2018", basename(partition_file)) == TRUE){
    # Extract all lines with an equals sign (these lines will define a gene)
    dash_lines <- grep("\\-", lines, ignore.case = TRUE, value = TRUE)
    # Split the genes into chunks by breaking at the commas ","
    gene_chunks <-  unlist(lapply(strsplit(dash_lines, ","), function(x){x[2]}))
    # Format the gene chunks nicely
    gene_range <- gsub(";", "", gsub(" ", "", gene_chunks))
    # Get start and end of each gene
    gene_start <- as.numeric(unlist(lapply(strsplit(gene_range, "-"), function(x){x[1]})))
    gene_end <- as.numeric(unlist(lapply(strsplit(gene_range, "-"), function(x){x[2]})))
  }
  # Create the gene names
  gene_names <- sprintf("gene_%04d", 1:length(gene_start))
  # Create a nice little dataframe for the genes
  gene_df <- data.frame(gene_name = gene_names, gene_range = gene_range, gene_start = gene_start, gene_end = gene_end)
  # Calculate the gene length
  gene_df$gene_length <- gene_end - (gene_start - 1) # subtract one from gene_start to count the starting site in the gene length
  # Set output
  if (return.lengths == TRUE){
    output <- gene_df$gene_length
  } else {
    output <- gene_df
  }
  # Return the gene length or gene dataframe
  return(output)
}


gene.lengths.smatrix <- function(partition_file, return.lengths = TRUE){
  ## Extract all gene lengths from an smatrix.txt file 
  # Open file
  lines <- readLines(partition_file)
  # Extract all lines with arrow "=>" (these lines will define a gene)
  arrow_lines <- grep("\\=>", lines, value = TRUE)
  # Split the charset lines at the tab "\t"
  gene_chunks <- unlist(lapply(strsplit(arrow_lines, "\t"), function(x){x[1]}))
  # Format the gene chunks nicely
  gene_chunks_nospace <- gsub(" ", "", gene_chunks)
  # Get start and end of each gene
  gene_start <- as.numeric(unlist(lapply(strsplit(gene_chunks_nospace, "\\=>"), function(x){x[1]})))
  gene_end <- as.numeric(unlist(lapply(strsplit(gene_chunks_nospace, "\\=>"), function(x){x[2]})))
  # Create the gene names
  gene_names <- sprintf("gene_%04d", 1:length(gene_start))
  # Create a nice little dataframe for the genes
  gene_df <- data.frame(gene_name = gsub(" ", "", gene_names), gene_range = gene_chunks_nospace, gene_start = gene_start, gene_end = gene_end)
  # Calculate the gene length
  gene_df$gene_length <- gene_end - (gene_start - 1) # subtract one from gene_start to count the starting site in the gene length
  # Set output
  if (return.lengths == TRUE){
    output <- gene_df$gene_length
  } else {
    output <- gene_df
  }
  # Return the gene length or gene dataframe
  return(output)
}


redo.gene.names <- function(partition_file){
  ## Function to take a partition file, update all the names, and return it
  
  ## Open partition file
  lines <- readLines(partition_file)
  
  ## Create new charsets charsets
  # Extract all lines with a charset
  charset_lines <- grep("charset", lines, ignore.case = TRUE, value = TRUE)
  # Split the charset lines at the "="
  gene_lines <- unlist(lapply(strsplit(charset_lines, "="), function(x){x[2]}))
  # Create the gene names
  new_gene_names <- sprintf("gene_%04d", 1:length(gene_lines))
  # Identify old gene names
  old_gene_name_chunks <- unlist(lapply(strsplit(charset_lines, "="), function(x){x[1]}))
  old_gene_name_raw <- unlist(lapply(strsplit(old_gene_name_chunks, "charset"), function(x){x[2]}))
  old_gene_names <- gsub(" ", "", old_gene_name_raw)
  # Create new charsets
  new_charsets <- paste0('\tcharset ', new_gene_names, " =", gene_lines)
  
  ## Create new charpartition
  # Extract models
  charpartition_lines <- grep(":", lines, value = T)
  model_chunk_raw <- unlist(lapply(strsplit(charpartition_lines,":"), function(x){x[1]}))
  models <- gsub(" ", "", model_chunk_raw)
  # New charpartition sections
  cpart_lines <- paste0("\t", models, ": ", new_gene_names)
  # New charpartition
  new_charpartition <- paste0("charpartition mymodels =\n", paste(cpart_lines, collapse = ",\n"), ";")
  
  ## Complete partition file 
  # Create whole partition file
  new_p_file <- c("#nexus", "begin sets;", new_charsets, new_charpartition, "end;", "")
  # Write out updated partition file to the input partition file path
  write(new_p_file, file = partition_file)
  # Return the partition file
  return(partition_file)
}




#### Extract genes ####
extract.all.genes <- function(alignment_file, partition_file, dataset_id, gene_directory, create.gene.alignments = TRUE){
  ## Given a partition file and an alignment, extract and output all genes as fasta files in the specified location
  
  ## Identify genes
  # Open partition file
  partition_lines <- readLines(partition_file)
  # Extract gene start and end
  charsets <- grep("charset", partition_lines, ignore.case = TRUE, value = TRUE)
  charset_split <- strsplit(charsets, "=")
  gene_locations_chunk <- gsub(";", "", unlist(lapply(charset_split, function(x){x[[2]]})))
  gene_start <- unlist(lapply(strsplit(gene_locations_chunk, "-"), function(x){x[[1]]}))
  gene_end <- unlist(lapply(strsplit(gene_locations_chunk, "-"), function(x){x[[2]]}))
  # Extract gene names
  gene_name_chunk <- unlist(lapply(charset_split, function(x){x[[1]]}))
  gene_names_raw <- unlist(lapply(strsplit(gene_name_chunk, " "), function(x){x[[2]]}))
  gene_names <- gsub(" ", "", gene_names_raw)
  # Turn gene locations into a table
  gene_table <- data.frame(name = gene_names, 
                           start = as.numeric(gsub(" ", "", gene_start)), 
                           end = as.numeric(gsub(" ", "", gene_end)) )
  
  ## Prepare alignment
  # Identify alignment type
  suffix <- tail(strsplit(basename(alignment_file), "\\.")[[1]], 1)
  # Open alignment
  if (suffix == "phy" | suffix == "phylip"){
    # Read alignment in as matrix
    al <- read.aa(alignment_file)
    # Convert to AAbin object
    al <- as.AAbin(al)
  } else if (suffix == "fasta" | suffix == "fa" | suffix == "fas"){
    # Read alignment in as AAbin list object
    al <- read.FASTA(file = alignment_file, type = "AA")
    # Convert AAbin list to AAbin matrix 
    al <- as.matrix(al)
  } else if (suffix == "nex" | suffix == "nexus"){
    # Read alignment in as list 
    al <- read.nexus.data(alignment_file)
    # Convert alignment to AAbin matrix
    al <- as.matrix(as.AAbin(al))
  }
  
  ## Split alignment (default - does not need specification in function call)
  if (create.gene.alignments == TRUE){
    # Use lapply to call extract.one.gene function for each row in the gene_table
    lapply_output <- lapply(1:nrow(gene_table), function(i){extract.one.gene(gene_details = gene_table[i, ], 
                                                                             al = al, 
                                                                             dataset_id = dataset_id, 
                                                                             gene_directory = gene_directory)})
  }
  
  ## Return gene details
  # Rename table rows
  names(gene_table) <- c("gene_name", "gene_start", "gene_end")
  # Add extra columns with details about the alignments
  gene_table$alignment_file <- basename(alignment_file)
  gene_table$partition_file <- basename(partition_file)
  gene_table$gene_directory <- gene_directory
  gene_table$dataset_id <- dataset_id
  gene_table$dataset <- strsplit(dataset_id, "\\.")[[1]][1]
  gene_table$matrix_name <- strsplit(dataset_id, "\\.")[[1]][2]
  gene_table$gene_id <- paste0(dataset_id, ".", gene_table$gene_name)
  gene_table$gene_file <- paste0(dataset_id, ".", gene_table$gene_name, ".fa")
  # Rearrange gene columns
  gene_table <- gene_table[ , c("dataset", "matrix_name", "dataset_id", "gene_name",
                                "gene_id", "gene_start", "gene_end", "gene_file",
                                "alignment_file", "partition_file", "gene_directory") ]
  # Return the table
  return(gene_table)
}



extract.one.gene <- function(gene_details, al, dataset_id, gene_directory){
  # Extract a single gene from a larger alignment
  # gene_details is a single row from a dataframe with the columns "name, "start" and "stop"
  # al is an alignment in matrix format
  
  ## Extract gene
  # Filter alignment matrix using gene start and end
  gene_matrix <- al[ , gene_details$start:gene_details$end ]
  # Convert gene to AAbin
  gene_al <- as.list(gene_matrix)
  
  ## Check for missing taxa
  # Identify the number of unique characters in the sequence for each taxa
  taxa_unique_chars <- unlist(lapply(1:length(gene_al), function(i){length(unique(as.character(gene_al)[[i]]))}))
  # Identify which taxa have only one unique character
  check_taxa_1 <- which(taxa_unique_chars == 1)
  check_taxa_2 <- which(taxa_unique_chars == 2)
  # Remove any taxa with no sequence for this gene, if required
  # Check all taxa with one unique character
  if (length(check_taxa_1) > 0){
    # Trimming required
    # Identify the unique character for the taxa with only one character in the gene
    check_chars_1 <- unlist(lapply(check_taxa_1, function(i){unique(as.character(gene_al)[[i]])}))
    names(check_chars_1) <- check_taxa_1
    # Identify tips to remove and tips to keep
    tips_to_remove <- as.numeric(names(check_chars_1 == "-" | check_chars_1 == "?"))
    tips_to_keep <- setdiff(1:length(gene_al), tips_to_remove)
    # Remove empty taxa
    gene_al_trimmed <- gene_al[tips_to_keep]
  } else {
    # No trimming required
    gene_al_trimmed <- gene_al
  }
  # Check all taxa with two unique characters
  if (length(check_taxa_2) > 0){
    # Trimming required
    # Output gene details for manual check
    print("2 missing chars")
    print(dataset_id)
    print(gene_details)
    # Identify the unique characters for each taxa with 2 characters
    check_chars_2 <- lapply(check_taxa_2, function(i){unique(as.character(gene_al)[[i]])})
    # Identify any options that contain both "?" and "-"
    tips_to_remove <- c()
    for (i in 1:length(check_chars_2)){
      if (setequal(check_chars_2[[i]], c("?", "-"))){
        tips_to_remove <- c(tips_to_remove, i)
      }
    }
    # Identify tips to keep
    tips_to_keep <- setdiff(1:length(gene_al_trimmed), tips_to_remove)
    # Remove empty taxa
    gene_al_complete <- gene_al_trimmed[tips_to_keep]
  } else {
    # No trimming required
    gene_al_complete <- gene_al_trimmed
  }
  
  ## Output gene
  # Assemble output file
  gene_output_file <- paste0(gene_directory, dataset_id, ".", gene_details$name, ".fa")
  # Save gene as a fasta file
  write.fasta(sequences = gene_al_complete, names = names(gene_al_complete), file.out = gene_output_file)
}




#### Trim constraint trees ####
## Remove taxa from the constraint trees that aren't present in each gene

trim.constraint.tree.taxa <- function(row_id, gene_df, constraint_tree_files){
  ## Remove unneeded taxa from all constraint trees and save each updated tree
  
  # Extract row for this rownumber
  gene_row <- gene_df[row_id, ]
  # Identify relevant constraint tree files
  gene_ct_files <- grep(gene_row$dataset_id, constraint_tree_files, value = T)
  # Open alignment and identify which taxa are present
  gene_al_file <- paste0(gene_row$gene_directory, gene_row$gene_file)
  gene_taxa <- names(read.fasta(gene_al_file))
  # Open trees and remove unneeded tips
  gene_id <- gene_row$gene_id
  trimmed_tree_output_directory <- paste0(dirname(gene_al_file), "/")
  # Trim each of the three constraint trees and output each row with the new information
  gene_row$constraint_tree_1 <- trim.one.constraint.tree(gene_id = gene_id, 
                                                         gene_taxa = gene_taxa, 
                                                         trimmed_tree_output_directory = trimmed_tree_output_directory, 
                                                         constraint_tree_file = gene_ct_files[1])
  gene_row$constraint_tree_2 <- trim.one.constraint.tree(gene_id = gene_id, 
                                                         gene_taxa = gene_taxa, 
                                                         trimmed_tree_output_directory = trimmed_tree_output_directory, 
                                                         constraint_tree_file = gene_ct_files[2])
  gene_row$constraint_tree_3 <- trim.one.constraint.tree(gene_id = gene_id, 
                                                         gene_taxa = gene_taxa, 
                                                         trimmed_tree_output_directory = trimmed_tree_output_directory, 
                                                         constraint_tree_file = gene_ct_files[3])
  # Return the updated row
  return(gene_row)
}

trim.one.constraint.tree <- function(gene_id, gene_taxa, trimmed_tree_output_directory, constraint_tree_file){
  ## Function to remove unneeded taxa from one constraint tree
  
  # Open constraint tree
  t <- read.tree(constraint_tree_file)
  # Remove any taxa from the gene taxa list that don't appear in the tree 
  #   (i.e., Trichoplax which is unconstrained)
  keep_taxa <- intersect(gene_taxa, t$tip.label)
  # Remove missing tips
  trimmed_t <- keep.tip(t, keep_taxa)
  # Create new filename for outputting trimmed constraint tree
  constraint_tree_split <- strsplit(basename(constraint_tree_file), "\\.")[[1]]
  trimmed_tree_output_file <- paste0(trimmed_tree_output_directory, gene_id, ".", constraint_tree_split[3], ".nex")
  # Write the trimmed tree
  write.tree(trimmed_t, file = trimmed_tree_output_file, append = FALSE)
  # Return the trimmed tree filepath
  return(trimmed_tree_output_file)
}




#### Checking simulation completion ####
check.rep <- function(dir_name){
  # Check that the ASTRAL tree, gene trees and IQ-Tree ML tree all exist
  
  # Assemble output file names
  id <- basename(dir_name)
  astral_tree_path <- paste0(id, "_ASTRAL_tree.tre")
  gene_iqtree_path <- paste0(id, "_gene_trees.iqtree")
  # List files in directory
  dir_files <- list.files(dir_name)
  # Check for file presence
  a_tree_exists <- TRUE %in% grepl(astral_tree_path, dir_files)
  g_tree_exists <- TRUE %in% grepl(gene_iqtree_path, dir_files)
  # Assemble output row
  op <- c(id, a_tree_exists, g_tree_exists, dir_name)
  names(op) <- c("id", "astral_tree_complete", "gene_trees_complete", "directory_path")
  # Return output vector
  return(op)
}



