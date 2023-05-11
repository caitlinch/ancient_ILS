# ancient_ILS/code/func_prepare_trees.R
## This script includes functions for formatting and preparing hypothesis trees
# Caitlin Cherryh, 2023

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
gene.lengths.nexus <- function(partition_file){
  ## Extract all gene lengths from a nexus partition file
  # Open partition file
  lines <- readLines(partition_file)
  # Extract all lines with a charset
  charset_lines <- grep("charset", lines, ignore.case = TRUE, value = TRUE)
  # Split the charset lines at the "="
  gene_chunks <- unlist(lapply(strsplit(charset_lines, "="), function(x){x[2]}))
  # Format the gene chunks nicely
  gene_chunks_nospace <- gsub(" ", "", gene_chunks)
  gene_chunks_noend <- gsub(";", "", gene_chunks_nospace)
  # Get start and end of each gene
  gene_start <- as.numeric(unlist(lapply(strsplit(gene_chunks_noend, "-"), function(x){x[1]})))
  gene_end <- as.numeric(unlist(lapply(strsplit(gene_chunks_noend, "-"), function(x){x[2]})))
  # Create a nice little dataframe for the genes
  gene_df <- data.frame(gene_range = gene_chunks_noend, gene_start = gene_start, gene_end = gene_end)
  # Calculate the gene length
  gene_df$gene_length <- gene_end - (gene_start - 1) # subtract one from gene_start to count the starting site in the gene length
  # Return the gene length
  return(gene_df$gene_length)
}


gene.lengths.raxml <- function(partition_file){
  ## Extract all gene lengths from a RAxML partition file 
  # Open partition file
  lines <- readLines(partition_file)
  # Extract all lines with an equals sign (these lines will define a gene)
  eq_lines <- grep("\\=", lines, ignore.case = TRUE, value = TRUE)
  # Split the charset lines at the "="
  gene_chunks <- unlist(lapply(strsplit(eq_lines, "="), function(x){x[2]}))
  # Format the gene chunks nicely
  gene_chunks_nospace <- gsub(" ", "", gene_chunks)
  # Get start and end of each gene
  gene_start <- as.numeric(unlist(lapply(strsplit(gene_chunks_nospace, "-"), function(x){x[1]})))
  gene_end <- as.numeric(unlist(lapply(strsplit(gene_chunks_nospace, "-"), function(x){x[2]})))
  # Create a nice little dataframe for the genes
  gene_df <- data.frame(gene_range = gene_chunks_nospace, gene_start = gene_start, gene_end = gene_end)
  # Calculate the gene length
  gene_df$gene_length <- gene_end - (gene_start - 1) # subtract one from gene_start to count the starting site in the gene length
  # Return the gene length
  return(gene_df$gene_length)
}


gene.lengths.smatrix <- function(partition_file){
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
  # Create a nice little dataframe for the genes
  gene_df <- data.frame(gene_range = gene_chunks_nospace, gene_start = gene_start, gene_end = gene_end)
  # Calculate the gene length
  gene_df$gene_length <- gene_end - (gene_start - 1) # subtract one from gene_start to count the starting site in the gene length
  # Return the gene length
  return(gene_df$gene_length)
}

