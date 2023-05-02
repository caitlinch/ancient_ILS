# ancient_ILS/code/func_prepare_trees.R
## This script includes functions for formatting and preparing hypothesis trees
# Caitlin Cherryh, 2023

format.constraint.tree.clade <- function(clade){
  # Function to take in a vector of species and return a nicely formatted character object to paste into a constraint tree
  
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
  # Small function to take one line from the models file and return in partition format
  
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
