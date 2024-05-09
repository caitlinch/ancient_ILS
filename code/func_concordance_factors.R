# ancient_ILS/code/func_concordance_factors.R
## This script includes functions to extract concordance factors for different datasets
# Caitlin Cherryh, 2023

library(ape)

#### Utility functions ####
# Small function to create a table of nodes and node labels for any "phylo" object
select.tip.or.node <- function(element, tree) {
  # Function selecting the tip or node name corresponding to the edge row
  ifelse(element < Ntip(tree)+1,
         tree$tip.label[element],
         tree$node.label[element-Ntip(tree)])
}



#### Extract qCF ####
extract.qcf.wrapper <- function(i, qcf_df, 
                                matrix_taxa = matrix_taxa, all_datasets = all_datasets, alignment_taxa_df = alignment_taxa_df){
  ## Extract the single row from the qCF dataframe
  row <- qcf_df[i, ]
  
  ## Extract the relevant list of taxa for this dataframe
  # First, check whether this matrix is included in the keys of the matrix_taxa list
  row_dataset <- row$dataset
  row_key <- paste0(row$dataset, ".", row$matrix_name, ".", "aa")
  list_keys <- names(matrix_taxa)
  # Check if row_key in list_key
  if ((row_key %in% list_keys) == FALSE){
    # If row key is not in list key, then all taxa for this dataset have the same names
    # Extract the object containing those taxa names
    constraint_clades <- all_datasets[[row_dataset]]
  } else if ((row_key %in% list_keys) == TRUE){
    # First, identify the list of taxa in this matrix
    keep_taxa <- matrix_taxa[[row_key]]
    # Secondly, extract the taxa clades for this dataset
    dataset_taxa_clades <- all_datasets[[row_dataset]]
    # Make a copy of the clades object
    constraint_clades <- dataset_taxa_clades
    # Lastly, remove any taxa that is NOT in the keep_taxa from the constraint clades
    #   i.e., remove any taxa from this dataset that are NOT present in this matrix
    #   (as some datasets have multiple matrices, with different taxon sampling or different taxon naming conventions)
    constraint_clades$Bilateria <- dataset_taxa_clades$Bilateria[which(dataset_taxa_clades$Bilateria %in% keep_taxa)]
    constraint_clades$Cnidaria <- dataset_taxa_clades$Cnidaria[which(dataset_taxa_clades$Cnidaria %in% keep_taxa)]
    constraint_clades$Placozoa <- dataset_taxa_clades$Placozoa[which(dataset_taxa_clades$Placozoa %in% keep_taxa)]
    constraint_clades$Porifera <- dataset_taxa_clades$Porifera[which(dataset_taxa_clades$Porifera %in% keep_taxa)]
    constraint_clades$Ctenophora <- dataset_taxa_clades$Ctenophora[which(dataset_taxa_clades$Ctenophora %in% keep_taxa)]
    constraint_clades$Outgroup <- dataset_taxa_clades$Outgroup[which(dataset_taxa_clades$Outgroup %in% keep_taxa)]
  }
  
  ## Remove any taxa from the constraint_clades that are not included in the ML tree for the alignment
  # Extract the column of tip labels from the relevant unique_id column of the alignment_taxa_df
  al_col_key <- paste0(row$dataset, ".", row$matrix_name)
  tree_tips_raw <- alignment_taxa_df[[c(al_col_key)]]
  tree_tips_cleaned <- na.omit(tree_tips_raw)
  tree_tips <- as.character(tree_tips_cleaned)
  # Check each of the clades and remove any tips not in the list of tree tips
  constraint_clades$Bilateria <- constraint_clades$Bilateria[(constraint_clades$Bilateria %in% tree_tips)]
  constraint_clades$Cnidaria <- constraint_clades$Cnidaria[(constraint_clades$Cnidaria %in% tree_tips)]
  constraint_clades$Placozoa <- constraint_clades$Placozoa[(constraint_clades$Placozoa %in% tree_tips)]
  constraint_clades$Porifera <- constraint_clades$Porifera[(constraint_clades$Porifera %in% tree_tips)]
  constraint_clades$Ctenophora <- constraint_clades$Ctenophora[(constraint_clades$Ctenophora %in% tree_tips)]
  constraint_clades$Outgroup <- constraint_clades$Outgroup[(constraint_clades$Outgroup %in% tree_tips)]
  
  ## Extract qCF depending on topology
  qcf_extracted <- extract.qcf(dataset = row_dataset, matrix_name = row$matrix_name, topology = row$tree_topology,
                               tree_file = row$qcf_tree_file, constraint_clades = constraint_clades)
  
  ## Return the qCF values for this tree
  return(qcf_extracted)
}



extract.qcf <- function(dataset, matrix_name, topology, 
                        tree_file, constraint_clades){
  # Function to extract qCF for CTEN tree or PORI tree
  
  ## Open tree with qCF annotation
  q_tree <- read.tree(tree_file)
  # Root tree at outgroup
  q_rooted <- root(q_tree, constraint_clades$Outgroup)
  
  ## Make an edge table with nodes and node labels
  q_edge_table <- data.frame(
    "parent" = q_rooted$edge[,1],
    "par.name" = sapply(q_rooted$edge[,1],
                        select.tip.or.node,
                        tree = q_rooted),
    "child" = q_rooted$edge[,2],
    "chi.name" = sapply(q_rooted$edge[,2],
                        select.tip.or.node,
                        tree = q_rooted)
  )
  
  ## Branches to extract length and node values:
  ## All animals (CTEN+PORI+CNID+BILAT+PLAC)
  # Identify tips in this group
  met_taxa <- c(constraint_clades$Bilateria, constraint_clades$Cnidaria, constraint_clades$Placozoa, 
                constraint_clades$Porifera, constraint_clades$Ctenophora)
  # Extract MRCA
  met_cn <- getMRCA(q_rooted, met_taxa) # child node
  met_pn <- q_rooted$edge[which(q_rooted$edge[,2] == met_cn), 1] # parent_node
  # Extract branch id
  met_branch_id <- which(q_rooted$edge[,1] == met_pn & q_rooted$edge[,2] == met_cn)
  # Extract branch length
  met_branch_length <- q_rooted$edge.length[met_branch_id]
  # Extract node value 
  met_node_value <- q_edge_table[which(q_edge_table$parent == met_pn & q_edge_table$child == met_cn), ]$chi.name
  
  ## Key branch (leading to ALL OTHER ANIMALS aka PLAC+CNID+BILAT)
  # Identify tips in this group - do not include PLAC when identifying MRCA, 
  #     as sometimes PLAC placement is sister to PORI which will result in 
  #     extracting the wrong branch
  if (topology == "CTEN"){
    key_taxa <- c(constraint_clades$Porifera, constraint_clades$Cnidaria, constraint_clades$Bilateria)
  } else if (topology == "PORI"){
    key_taxa <- c(constraint_clades$Ctenophora, constraint_clades$Cnidaria, constraint_clades$Bilateria)
  } else if (topology == "CTEN_PORI" | topology == "CTEN.PORI"){
    key_taxa <- c(constraint_clades$Ctenophora, constraint_clades$Porifera)
  }
  # If there are taxa within the key branch, extract the branch length and node values
  if (length(key_taxa > 1)){
    # Extract MRCA
    key_cn <- getMRCA(q_rooted, key_taxa) # child node
    key_pn <- q_rooted$edge[which(q_rooted$edge[,2] == key_cn), 1] # parent node
    # Extract branch id
    key_branch_id <- which(q_rooted$edge[,1] == key_pn & q_rooted$edge[,2] == key_cn)
    # Extract branch length
    key_branch_length <- q_rooted$edge.length[key_branch_id]
    # Extract node value
    key_node_value <- q_edge_table[which(q_edge_table$parent == key_pn & q_edge_table$child == key_cn), ]$chi.name
  } else {
    # Assign NA if only 1 tip
    key_branch_length  <- NA
    key_node_value     <- NA
  }
  
  ## CTEN (leading to CTEN branch)
  # Identify tips in this group
  cten_taxa <- c(constraint_clades$Ctenophora)
  if (length(cten_taxa) > 1){
    # Extract MRCA
    cten_cn <- getMRCA(q_rooted, cten_taxa) # child node
    cten_pn <- q_rooted$edge[which(q_rooted$edge[,2] == cten_cn), 1] # parent node
    # Extract branch id
    cten_branch_id <- which(q_rooted$edge[,1] == cten_pn & q_rooted$edge[,2] == cten_cn)
    # Extract branch length
    cten_branch_length <- q_rooted$edge.length[cten_branch_id]
    # Extract node value
    cten_node_value <- q_edge_table[which(q_edge_table$parent == cten_pn & q_edge_table$child == cten_cn), ]$chi.name
  } else {
    # Assign NA if only 1 tip
    cten_branch_length  <- NA
    cten_node_value     <- NA
  }
  
  ## PORI (leading to PORI branch)
  # Identify tips in this group
  pori_taxa <- c(constraint_clades$Porifera)
  if (length(pori_taxa) > 1){
    # Extract MRCA
    pori_cn <- getMRCA(q_rooted, pori_taxa) # child node
    pori_pn <- q_rooted$edge[which(q_rooted$edge[,2] == pori_cn), 1] # parent node
    # Extract branch id
    pori_branch_id <- which(q_rooted$edge[,1] == pori_pn & q_rooted$edge[,2] == pori_cn)
    # Extract branch length
    pori_branch_length <- q_rooted$edge.length[pori_branch_id]
    # Extract node value
    pori_node_value <- q_edge_table[which(q_edge_table$parent == pori_pn & q_edge_table$child == pori_cn), ]$chi.name
  } else {
    # Assign NA if only 1 tip
    pori_branch_length  <- NA
    pori_node_value     <- NA
  }
  
  # Assemble output vector
  qcf_output <- c(met_branch_length, met_node_value, key_branch_length, key_node_value, 
                  cten_branch_length, cten_node_value, pori_branch_length, pori_node_value)
  names(qcf_output) <- c("MET_branch_length", "MET_node_value", "KEY_branch_length", "KEY_branch_value",
                         "CTEN_branch_length", "CTEN_node_value", "PORI_branch_length", "PORI_branch_value")
  return(qcf_output)
}





