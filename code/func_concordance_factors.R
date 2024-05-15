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



#### Extract gCF ####
extract.gcf.wrapper <- function(i, gcf_df, 
                                matrix_taxa = matrix_taxa, all_datasets = all_datasets, alignment_taxa_df = alignment_taxa_df){
  ## Extract the single row from the qCF dataframe
  row <- gcf_df[i, ]
  
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
  gcf_extracted <- extract.gcf(dataset = row_dataset, matrix_name = row$matrix_name, topology = row$tree_topology,
                               tree_file = row$gcf_branch_files, table_file = row$gcf_stat_files,
                               constraint_clades = constraint_clades)
  
  ## Return the qCF values for this tree along with the tree parameters
  return(gcf_extracted)
}



extract.gcf <- function(dataset, matrix_name, topology, 
                        tree_file, table_file, 
                        constraint_clades){
  # Function to extract gCF for any constrained tree topology (CTEN, PORI, CTEN_PORI)
  
  ## Open tree with gCF annotation
  g_tree <- read.tree(tree_file)
  # Root tree at outgroup
  g_rooted <- root(g_tree, constraint_clades$Outgroup)
  
  ## Make an edge table with nodes and node labels
  g_edge_table <- data.frame(
    "parent" = g_rooted$edge[,1],
    "par.name" = sapply(g_rooted$edge[,1],
                        select.tip.or.node,
                        tree = g_rooted),
    "child" = g_rooted$edge[,2],
    "chi.name" = sapply(g_rooted$edge[,2],
                        select.tip.or.node,
                        tree = g_rooted)
  )
  
  ## Open the stat table
  g_table <- read.table(table_file, header = T, sep = "", fill = T)
  # Check columns and correct if necessary
  if (length(which(is.na(g_table$Length))) == nrow(g_table)){
    # Check columns to see if the "Label" column is missing (in which case, the "Length" column will be all NA)
    g_table <- g_table[ , 1:11]
    names(g_table) <- c("ID", "gCF", "gCF_N", "gDF1", "gDF1_N", "gDF2", "gDF2_N", "gDFP", "gDFP_N", "gN", "Length")
  } else {
    # Remove "Label" column
    g_table <- g_table[ , c(1:10,12)]
    names(g_table) <- c("ID", "gCF", "gCF_N", "gDF1", "gDF1_N", "gDF2", "gDF2_N", "gDFP", "gDFP_N", "gN", "Length")
  }
  
  ## Branches to extract length and node values:
  ## All animals (CTEN+PORI+CNID+BILAT+PLAC)
  # Identify tips in this group
  met_taxa <- c(constraint_clades$Bilateria, constraint_clades$Cnidaria, constraint_clades$Placozoa, 
                constraint_clades$Porifera, constraint_clades$Ctenophora)
  # Extract MRCA
  met_cn <- getMRCA(g_rooted, met_taxa) # child node
  met_pn <- g_rooted$edge[which(g_rooted$edge[,2] == met_cn), 1] # parent_node
  # Potential branch ids
  potential_met_branch_id <- as.numeric(c(g_edge_table[which(g_edge_table$parent == met_pn & g_edge_table$child == met_cn), ]$par.name,
                                          g_edge_table[which(g_edge_table$parent == met_pn & g_edge_table$child == met_cn), ]$chi.name))
  # Potential branch branch lengths
  potential_met_branch_lengths <- g_table$Length[c(which(g_table$ID == potential_met_branch_id[1]),
                                                   which(g_table$ID == potential_met_branch_id[2]))]
  # Length of actual branch
  actual_met_branch          <- which(g_rooted$edge[,1] == met_pn & g_rooted$edge[,2] == met_cn)
  actual_met_branch_length  <- g_rooted$edge.length[actual_met_branch]
  # Identify the met_table_id (g_table "ID" column value) by length
  met_table_id <- potential_met_branch_id[which(round(potential_met_branch_lengths, digits = 5) == round(actual_met_branch_length, digits = 5))]
  # Check again with 1 less digit if no identical branch lengths
  if (length(met_table_id) == 0){
    met_table_id <- potential_met_branch_id[which(round(potential_met_branch_lengths, digits = 4) == round(actual_met_branch_length, digits = 4))]
  }
  # Extract the row from the stat table for this met_table_id
  met_values <- g_table[which(g_table$ID == met_table_id), ]
  names(met_values) <- paste0("MET_", names(met_values))
  
  ## Key branch (leading to ALL OTHER ANIMALS aka PLAC+CNID+BILAT)
  # Identify tips in this group - do not include PLAC when identifying MRCA, 
  #     as sometimes PLAC placement is sister to PORI which will result in 
  #     extracting the wrong branch
  if (topology == "CTEN"){
    # Extract taxa
    key_taxa <- c(constraint_clades$Porifera, constraint_clades$Cnidaria, constraint_clades$Bilateria)
  } else if (topology == "PORI"){
    # Extract taxa
    key_taxa <- c(constraint_clades$Ctenophora, constraint_clades$Cnidaria, constraint_clades$Bilateria)
  } else if (topology == "CTEN_PORI" | topology == "CTEN.PORI"){
    # Extract taxa
    key_taxa <- c(constraint_clades$Ctenophora, constraint_clades$Porifera)
  }
  # Extract MRCA
  key_cn <- getMRCA(g_rooted, key_taxa) # child node
  key_pn <- g_rooted$edge[which(g_rooted$edge[,2] == key_cn), 1] # parent node
  # Potential branch ids
  potential_key_branch_id <- as.numeric(c(g_edge_table[which(g_edge_table$parent == key_pn & g_edge_table$child == key_cn), ]$par.name,
                                          g_edge_table[which(g_edge_table$parent == key_pn & g_edge_table$child == key_cn), ]$chi.name))
  # Potential branch branch lengths
  potential_key_branch_lengths <- g_table$Length[c(which(g_table$ID == potential_key_branch_id[1]),
                                                   which(g_table$ID == potential_key_branch_id[2]))]
  # Length of actual branch
  actual_key_branch          <- which(g_rooted$edge[,1] == key_pn & g_rooted$edge[,2] == key_cn)
  actual_key_branch_length  <- g_rooted$edge.length[actual_key_branch]
  # Identify the key_table_id (g_table "ID" column value) by length
  key_table_id <- potential_key_branch_id[which(round(potential_key_branch_lengths, digits = 5) == round(actual_key_branch_length, digits = 5))]
  # Check again with 1 less digit if no identical branch lengths
  if (length(key_table_id) == 0){
    key_table_id <- potential_key_branch_id[which(round(potential_key_branch_lengths, digits = 4) == round(actual_key_branch_length, digits = 4))]
  }
  # Extract the row from the stat table for this key_table_id
  key_values <- g_table[which(g_table$ID == key_table_id), ]
  names(key_values) <- paste0("KEY_", names(key_values))
  
  ## CTEN (leading to CTEN branch)
  # Identify tips in this group
  cten_taxa <- c(constraint_clades$Ctenophora)
  # Check more than one tip present
  if (length(cten_taxa) > 1){
    # Check monophyly
    if (is.monophyletic(g_rooted, cten_taxa) == TRUE){
      # Extract MRCA
      cten_cn <- getMRCA(g_rooted, cten_taxa) # child node
      cten_pn <- g_rooted$edge[which(g_rooted$edge[,2] == cten_cn), 1] # parent node
      # Potential branch ids
      potential_cten_branch_id <- as.numeric(c(g_edge_table[which(g_edge_table$parent == cten_pn & g_edge_table$child == cten_cn), ]$par.name,
                                               g_edge_table[which(g_edge_table$parent == cten_pn & g_edge_table$child == cten_cn), ]$chi.name))
      # Potential branch branch lengths
      potential_cten_branch_lengths <- g_table$Length[c(which(g_table$ID == potential_cten_branch_id[1]),
                                                        which(g_table$ID == potential_cten_branch_id[2]))]
      # Length of actual branch
      actual_cten_branch          <- which(g_rooted$edge[,1] == cten_pn & g_rooted$edge[,2] == cten_cn)
      actual_cten_branch_length  <- g_rooted$edge.length[actual_cten_branch]
      # Identify the cten_table_id (g_table "ID" column value) by length
      cten_table_id <- potential_cten_branch_id[which(round(potential_cten_branch_lengths, digits = 5) == round(actual_cten_branch_length, digits = 5))]
      # Check again with 1 less digit if no identical branch lengths
      if (length(cten_table_id) == 0){
        cten_table_id <- potential_cten_branch_id[which(round(potential_cten_branch_lengths, digits = 4) == round(actual_cten_branch_length, digits = 4))]
      }
      # Extract the row from the stat table for this cten_table_id
      cten_values <- g_table[which(g_table$ID == cten_table_id), ]
      names(cten_values) <- paste0("CTEN_", names(cten_values))
      cten_clade_monophyly <- "Monophyletic"
    } else {
      # Assign NA if paraphyletic
      cten_values <- rep(NA, 11)
      names(cten_values) <- paste0("CTEN_", names(g_table))
      cten_clade_monophyly <- "Paraphyletic"
    }
  } else {
    # Assign NA if only 1 tip
    cten_values <- rep(NA, 11)
    names(cten_values) <- paste0("CTEN_", names(g_table))
    cten_clade_monophyly <- "Single_tip"
  }
  
  ## PORI (leading to PORI branch)
  # Identify tips in this group
  pori_taxa <- c(constraint_clades$Porifera)
  # Check more than one tip present
  if (length(pori_taxa) > 1){
    # Check monophyly
    if (is.monophyletic(g_rooted, pori_taxa) == TRUE){
      # Extract MRCA
      pori_cn <- getMRCA(g_rooted, pori_taxa) # child node
      pori_pn <- g_rooted$edge[which(g_rooted$edge[,2] == pori_cn), 1] # parent node
      # Potential branch ids
      potential_pori_branch_id <- as.numeric(c(g_edge_table[which(g_edge_table$parent == pori_pn & g_edge_table$child == pori_cn), ]$par.name,
                                               g_edge_table[which(g_edge_table$parent == pori_pn & g_edge_table$child == pori_cn), ]$chi.name))
      # Potential branch branch lengths
      potential_pori_branch_lengths <- g_table$Length[c(which(g_table$ID == potential_pori_branch_id[1]),
                                                        which(g_table$ID == potential_pori_branch_id[2]))]
      # Length of actual branch
      actual_pori_branch          <- which(g_rooted$edge[,1] == pori_pn & g_rooted$edge[,2] == pori_cn)
      actual_pori_branch_length  <- g_rooted$edge.length[actual_pori_branch]
      # Identify the pori_table_id (g_table "ID" column value) by length
      pori_table_id <- potential_pori_branch_id[which(round(potential_pori_branch_lengths, digits = 5) == round(actual_pori_branch_length, digits = 5))]
      # Check again with 1 less digit if no identical branch lengths
      if (length(pori_table_id) == 0){
        pori_table_id <- potential_pori_branch_id[which(round(potential_pori_branch_lengths, digits = 4) == round(actual_pori_branch_length, digits = 4))]
      }
      # Extract the row from the stat table for this pori_table_id
      pori_values <- g_table[which(g_table$ID == pori_table_id), ]
      names(pori_values) <- paste0("PORI_", names(pori_values))
      pori_clade_monophyly <- "Monophyletic"
    } else {
      # Assign NA if paraphyletic
      pori_values <- rep(NA, 11)
      names(pori_values) <- paste0("PORI_", names(g_table))
      pori_clade_monophyly <- "Paraphyletic"
    }
  } else {
    # Assign NA if only 1 tip
    pori_values <- rep(NA, 11)
    names(pori_values) <- paste0("PORI_", names(g_table))
    pori_clade_monophyly <- "Single_tip"
  }
  
  # Assemble output vector
  gcf_output <- as.character(c(met_values, key_values, cten_values, cten_clade_monophyly, pori_values, pori_clade_monophyly))
  names(gcf_output) <- c(names(met_values), names(key_values), names(cten_values), "CTEN_monophyly", names(pori_values), "PORI_monophyly")
  return(gcf_output)
}



#### Reformat gCF dataframe ####
reformat.gCF.df <- function(input_df){
  ## Reformat the gCF dataframe into wide format
  
  # Identify rows for each topology
  cten_rows <- which(input_df$tree_topology == "CTEN")
  pori_rows <- which(input_df$tree_topology == "PORI")
  ctenpori_rows <- which(input_df$tree_topology == "CTEN_PORI")
  # Create new dataframe with columns for each topology
  new_df <- data.frame(id = unlist(lapply(strsplit(gcf_collated_df$id[c(T,F,F)], "\\."), function(x){paste0(x[[1]], ".", x[[2]], ".", x[[3]])})),
                       dataset = input_df$dataset[c(T,F,F)],
                       matrix_name = input_df$matrix_name[c(T,F,F)],
                       dataset_id = input_df$dataset_id[c(T,F,F)],
                       model = input_df$model[c(T,F,F)],
                       CTEN.MET_gCF = input_df$MET_gCF[cten_rows],
                       PORI.MET_gCF = input_df$MET_gCF[pori_rows],
                       CTENPORI.MET_gCF = input_df$MET_gCF[ctenpori_rows],
                       CTEN.MET_gCF_N = input_df$MET_gCF_N[cten_rows],
                       PORI.MET_gCF_N = input_df$MET_gCF_N[pori_rows],
                       CTENPORI.MET_gCF_N = input_df$MET_gCF_N[ctenpori_rows],
                       CTEN.MET_gDFP = input_df$MET_gDFP[cten_rows],
                       PORI.MET_gDFP = input_df$MET_gDFP[pori_rows],
                       CTENPORI.MET_gDFP = input_df$MET_gDFP[ctenpori_rows],
                       CTEN.MET_gDFP_N = input_df$MET_gDFP_N[cten_rows],
                       PORI.MET_gDFP_N = input_df$MET_gDFP_N[pori_rows],
                       CTENPORI.MET_gDFP_N = input_df$MET_gDFP_N[ctenpori_rows],
                       CTEN.MET_Length = input_df$MET_Length[cten_rows],
                       PORI.MET_Length = input_df$MET_Length[pori_rows],
                       CTENPORI.MET_Length = input_df$MET_Length[ctenpori_rows],
                       CTEN.KEY_gCF = input_df$KEY_gCF[cten_rows],
                       PORI.KEY_gCF = input_df$KEY_gCF[pori_rows],
                       CTENPORI.KEY_gCF = input_df$KEY_gCF[ctenpori_rows],
                       CTEN.KEY_gCF_N = input_df$KEY_gCF_N[cten_rows],
                       PORI.KEY_gCF_N = input_df$KEY_gCF_N[pori_rows],
                       CTENPORI.KEY_gCF_N = input_df$KEY_gCF_N[ctenpori_rows],
                       CTEN.KEY_gDFP = input_df$KEY_gDFP[cten_rows],
                       PORI.KEY_gDFP = input_df$KEY_gDFP[pori_rows],
                       CTENPORI.KEY_gDFP = input_df$KEY_gDFP[ctenpori_rows],
                       CTEN.KEY_gDFP_N = input_df$KEY_gDFP_N[cten_rows],
                       PORI.KEY_gDFP_N = input_df$KEY_gDFP_N[pori_rows],
                       CTENPORI.KEY_gDFP_N = input_df$KEY_gDFP_N[ctenpori_rows],
                       CTEN.KEY_Length = input_df$KEY_Length[cten_rows],
                       PORI.KEY_Length = input_df$KEY_Length[pori_rows],
                       CTENPORI.KEY_Length = input_df$KEY_Length[ctenpori_rows],
                       CTEN.CTEN_gCF = input_df$CTEN_gCF[cten_rows],
                       PORI.CTEN_gCF = input_df$CTEN_gCF[pori_rows],
                       CTENPORI.CTEN_gCF = input_df$CTEN_gCF[ctenpori_rows],
                       CTEN.CTEN_gCF_N = input_df$CTEN_gCF_N[cten_rows],
                       PORI.CTEN_gCF_N = input_df$CTEN_gCF_N[pori_rows],
                       CTENPORI.CTEN_gCF_N = input_df$CTEN_gCF_N[ctenpori_rows],
                       CTEN.CTEN_gDFP = input_df$CTEN_gDFP[cten_rows],
                       PORI.CTEN_gDFP = input_df$CTEN_gDFP[pori_rows],
                       CTENPORI.CTEN_gDFP = input_df$CTEN_gDFP[ctenpori_rows],
                       CTEN.CTEN_gDFP_N = input_df$CTEN_gDFP_N[cten_rows],
                       PORI.CTEN_gDFP_N = input_df$CTEN_gDFP_N[pori_rows],
                       CTENPORI.CTEN_gDFP_N = input_df$CTEN_gDFP_N[ctenpori_rows],
                       CTEN.CTEN_Length = input_df$CTEN_Length[cten_rows],
                       PORI.CTEN_Length = input_df$CTEN_Length[pori_rows],
                       CTENPORI.CTEN_Length = input_df$CTEN_Length[ctenpori_rows],
                       CTEN.PORI_gCF = input_df$PORI_gCF[cten_rows],
                       PORI.PORI_gCF = input_df$PORI_gCF[pori_rows],
                       CTENPORI.PORI_gCF = input_df$PORI_gCF[ctenpori_rows],
                       CTEN.PORI_gCF_N = input_df$PORI_gCF_N[cten_rows],
                       PORI.PORI_gCF_N = input_df$PORI_gCF_N[pori_rows],
                       CTENPORI.PORI_gCF_N = input_df$PORI_gCF_N[ctenpori_rows],
                       CTEN.PORI_gDFP = input_df$PORI_gDFP[cten_rows],
                       PORI.PORI_gDFP = input_df$PORI_gDFP[pori_rows],
                       CTENPORI.PORI_gDFP = input_df$PORI_gDFP[ctenpori_rows],
                       CTEN.PORI_gDFP_N = input_df$PORI_gDFP_N[cten_rows],
                       PORI.PORI_gDFP_N = input_df$PORI_gDFP_N[pori_rows],
                       CTENPORI.PORI_gDFP_N = input_df$PORI_gDFP_N[ctenpori_rows],
                       CTEN.PORI_Length = input_df$PORI_Length[cten_rows],
                       PORI.PORI_Length = input_df$PORI_Length[pori_rows],
                       CTENPORI.PORI_Length = input_df$PORI_Length[ctenpori_rows] )
  # Return the reformatted dataframe
  return(new_df)
}



#### Extract qCF ####
extract.qcf.wrapper <- function(i, qcf_df, 
                                matrix_taxa = matrix_taxa, all_datasets = all_datasets, alignment_taxa_df = alignment_taxa_df){
  ## Extract the single row from the qCF dataframe
  # Extract the relevant rows from the qcf_df
  i_split     <- strsplit(i, "\\.")[[1]]
  i_dataset   <- i_split[[1]]
  i_matrix    <- i_split[[2]]
  i_model     <- i_split[[3]]
  i_rows      <- qcf_df[which(qcf_df$dataset == i_dataset & qcf_df$matrix_name == i_matrix & qcf_df$model == i_model), ]
  
  ## Extract the relevant list of taxa for this dataframe
  # First, check whether this matrix is included in the keys of the matrix_taxa list
  row_key   <- paste0(i_dataset, ".", i_matrix, ".", "aa")
  list_keys <- names(matrix_taxa)
  # Check if row_key in list_key
  if ((row_key %in% list_keys) == FALSE){
    # If row key is not in list key, then all taxa for this dataset have the same names
    # Extract the object containing those taxa names
    constraint_clades   <- all_datasets[[i_dataset]]
  } else if ((row_key %in% list_keys) == TRUE){
    # First, identify the list of taxa in this matrix
    keep_taxa           <- matrix_taxa[[row_key]]
    # Secondly, extract the taxa clades for this dataset
    dataset_taxa_clades <- all_datasets[[i_dataset]]
    # Make a copy of the clades object
    constraint_clades   <- dataset_taxa_clades
    # Lastly, remove any taxa that is NOT in the keep_taxa from the constraint clades
    #   i.e., remove any taxa from this dataset that are NOT present in this matrix
    #   (as some datasets have multiple matrices, with different taxon sampling or different taxon naming conventions)
    constraint_clades$Bilateria   <- dataset_taxa_clades$Bilateria[which(dataset_taxa_clades$Bilateria %in% keep_taxa)]
    constraint_clades$Cnidaria    <- dataset_taxa_clades$Cnidaria[which(dataset_taxa_clades$Cnidaria %in% keep_taxa)]
    constraint_clades$Placozoa    <- dataset_taxa_clades$Placozoa[which(dataset_taxa_clades$Placozoa %in% keep_taxa)]
    constraint_clades$Porifera    <- dataset_taxa_clades$Porifera[which(dataset_taxa_clades$Porifera %in% keep_taxa)]
    constraint_clades$Ctenophora  <- dataset_taxa_clades$Ctenophora[which(dataset_taxa_clades$Ctenophora %in% keep_taxa)]
    constraint_clades$Outgroup    <- dataset_taxa_clades$Outgroup[which(dataset_taxa_clades$Outgroup %in% keep_taxa)]
  }
  
  ## Remove any taxa from the constraint_clades that are not included in the ML tree for the alignment
  # Extract the column of tip labels from the relevant unique_id column of the alignment_taxa_df
  al_col_key        <- paste0(i_dataset, ".", i_matrix)
  tree_tips_raw     <- alignment_taxa_df[[c(al_col_key)]]
  tree_tips_cleaned <- na.omit(tree_tips_raw)
  tree_tips         <- as.character(tree_tips_cleaned)
  # Check each of the clades and remove any tips not in the list of tree tips
  constraint_clades$Bilateria   <- constraint_clades$Bilateria[(constraint_clades$Bilateria %in% tree_tips)]
  constraint_clades$Cnidaria    <- constraint_clades$Cnidaria[(constraint_clades$Cnidaria %in% tree_tips)]
  constraint_clades$Placozoa    <- constraint_clades$Placozoa[(constraint_clades$Placozoa %in% tree_tips)]
  constraint_clades$Porifera    <- constraint_clades$Porifera[(constraint_clades$Porifera %in% tree_tips)]
  constraint_clades$Ctenophora  <- constraint_clades$Ctenophora[(constraint_clades$Ctenophora %in% tree_tips)]
  constraint_clades$Outgroup    <- constraint_clades$Outgroup[(constraint_clades$Outgroup %in% tree_tips)]
  
  ## Extract CTEN_PORI clade values
  ctenpori_qcf <- extract.CTENPORI.qcf(dataset = i_dataset, matrix_name = i_matrix, topology = "CTEN_PORI", 
                                       tree_file = i_rows[which(i_rows$tree_topology == "CTEN_PORI"), "qcf_tree_file"], 
                                       constraint_clades)
  # Extract QC and EN values
  ctenpori_node_split <- strsplit(ctenpori_qcf[["KEY_node_value"]], "\\;")[[1]]
  ctenpori_qc <- gsub("QC=", "", grep("QC", ctenpori_node_split, value = T))
  ctenpori_en <- gsub("EN=", "", grep("EN", ctenpori_node_split, value = T))
  
  ## Extract qCF depending on topology
  cten_qcf <- extract.qcf(dataset = i_dataset, matrix_name = i_matrix, topology = "CTEN",
                          tree_file = i_rows[which(i_rows$tree_topology == "CTEN"), "qcf_tree_file"], 
                          ctenpori_qc = ctenpori_qc, ctenpori_en = ctenpori_en,
                          constraint_clades = constraint_clades)
  pori_qcf <- extract.qcf(dataset = i_dataset, matrix_name = i_matrix, topology = "PORI",
                          tree_file = i_rows[which(i_rows$tree_topology == "PORI"), "qcf_tree_file"], 
                          ctenpori_qc = ctenpori_qc, ctenpori_en = ctenpori_en,
                          constraint_clades = constraint_clades)
  
  ## Return the qCF values for this tree along with the tree parameters
  return(qcf_extracted)
}



extract.qcf <- function(dataset, matrix_name, topology, 
                        ctenpori_qc, ctenpori_en, 
                        tree_file, constraint_clades){
  # Function to extract qCF for CTEN or PORI constrained tree topologies
  
  ## Open tree with qCF annotation
  q_tree <- read.tree(tree_file)
  # Root tree at outgroup
  q_rooted <- root(q_tree, constraint_clades$Outgroup)
  # Drop Placozoa taxa
  if (length(constraint_clades$Placozoa) > 0){
    q_rooted <- drop.tip(q_rooted, constraint_clades$Placozoa)
  }
  
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
  
  
  ## Key branch (leading to ALL OTHER ANIMALS aka PLAC+CNID+BILAT)
  # Identify tips in this group - do not include PLAC when identifying MRCA, 
  #     as sometimes PLAC placement is sister to PORI which will result in 
  #     extracting the wrong branch
  if (topology == "CTEN" | topology == "PORI"){
    # Specify key taxa
    if (topology == "CTEN"){
      key_taxa <- c(constraint_clades$Porifera, constraint_clades$Cnidaria, constraint_clades$Bilateria)
    } else if (topology == "PORI"){
      key_taxa <- c(constraint_clades$Ctenophora, constraint_clades$Cnidaria, constraint_clades$Bilateria)
    } 
    # Extract node value, if there's more than 1 taxon in the key_taxa clade
    if (length(key_taxa) > 1){
      # Extract MRCA
      key_cn <- getMRCA(q_rooted, key_taxa) # child node
      key_pn <- q_rooted$edge[which(q_rooted$edge[,2] == key_cn), 1] # parent node
      # Extract the node that's identical to the key_check_node (i.e., the node that connects to the MET node)
      if (key_cn == key_check_node){
        # If the key_check_node and the child node are identical, then we want to extract values from the parent node
        # Extract child annotation from tree
        key_pn_lab <- q_rooted$node.label[(key_pn-Ntip(q_rooted))]
        # Clean string
        key_node_value <- gsub("\\[|\\]|'", "",  key_pn_lab)
      } else if (key_pn == key_check_node){
        # If the key_check_node and the parent node are identical, then we want to extract values from the child node
        # Extract parent annotation from tree
        key_cn_lab <- q_rooted$node.label[(key_cn-Ntip(q_rooted))]
        # Clean string
        key_node_value <- gsub("\\[|\\]|'", "",  key_cn_lab)
      }
      # Extract branch length
      key_branch_length <- q_rooted$edge.length[which(q_rooted$edge[,1] == key_pn & q_rooted$edge[,2] == key_cn)]
    } else {
      # Assign NA if only 1 tip
      key_branch_length  <- NA
      key_node_value     <- NA
    }
  } else if (topology == "CTEN_PORI" | topology == "CTEN.PORI"){
    # Specify key taxa
    key_taxa <- c(constraint_clades$Ctenophora, constraint_clades$Porifera)
    # Extract node value, if there's more than 1 taxon in the key_taxa clade
    if (length(key_taxa) > 1){
      # Extract MRCA
      key_cn <- getMRCA(q_rooted, key_taxa) # child node
      key_pn <- q_rooted$edge[which(q_rooted$edge[,2] == key_cn), 1] # parent node
      # Extract branch length
      key_branch_length <- q_rooted$edge.length[which(q_rooted$edge[,1] == key_pn & q_rooted$edge[,2] == key_cn)]
      # Identify correct node for the CTEN-PORI clade
      #     If the key_pn == key_check_node, then the key_cn is the correct node for the CTEN-PORI clade
      if (key_pn == key_check_node){
        # Extract child annotation from tree
        key_cn_lab <- q_rooted$node.label[(key_cn-Ntip(q_rooted))]
        # Clean string
        key_node_value <- gsub("\\[|\\]|'", "",  key_cn_lab)
      } else {
        # Extract the edges with the same key_pn (KEY parent node)
        key_edge_nodes  <- q_rooted$edge[c(which(q_rooted$edge[,1] == key_pn), which(q_rooted$edge[,2] == key_pn)), ]
        # Identify the possible child nodes
        test_edges <- c(which(q_rooted$edge[,1] %in% unique(as.numeric(key_edge_nodes))), which(q_rooted$edge[,2] %in% unique(as.numeric(key_edge_nodes))))
        test_nodes <- q_rooted$edge[test_edges,] # Extract all nodes connected to these nodes
        key_check_nodes <- sort(unique(as.numeric(test_nodes)), decreasing = TRUE) # Remove duplicates and order in decreasing order
        key_check_nodes <- key_check_nodes[which(key_check_nodes > Ntip(q_rooted))] # Remove external branches. Keep internal nodes.
        # Iterate and test the nodes
        key_node <- c()
        for (i in 1:length(key_check_nodes)){
          # Extract the tips in the clade defined by this node
          i_node          <- key_check_nodes[i]
          i_clade         <- extract.clade(q_rooted, i_node)
          i_tips          <- i_clade$tip.label
          i_tips_noPlac   <- i_tips[ which( ! (i_tips %in% constraint_clades$Placozoa) ) ]
          # Check whether the tips in this clade are the same as the tips in the CTEN_PORI clade
          if (setequal(i_tips_noPlac, key_taxa) == TRUE){
            key_node <- c(i_node)
          } # end: if (setequal(i_tips_noPlac, key_taxa) == TRUE){
        } # end: for (i in 1:length(key_check_nodes)){
        # Extract node annotation from tree
        key_lab <- q_rooted$node.label[(key_node-Ntip(q_rooted))]
        # Clean string
        key_node_value <- gsub("\\[|\\]|'", "",  key_lab)
      } # end: if (key_pn == key_check_node){ 
    } else {
      # Assign NA if only 1 tip
      key_branch_length  <- NA
      key_node_value     <- NA
    } # end: if (length(key_taxa) > 1){
  } # end: if (topology == "CTEN" | topology == "PORI")
  
  ## CTEN (leading to CTEN branch)
  # Identify tips in this group
  cten_taxa <- c(constraint_clades$Ctenophora)
  # Check that more than one taxa is present
  if (length(cten_taxa) > 1){
    # Check monophyly
    if (is.monophyletic(q_rooted, cten_taxa) == TRUE){
      # Extract MRCA
      cten_cn <- getMRCA(q_rooted, cten_taxa) # child node
      cten_pn <- q_rooted$edge[which(q_rooted$edge[,2] == cten_cn), 1] # parent node
      # Extract the q_edge_table row
      cten_q_edge_row <- q_edge_table[which(q_edge_table$parent == cten_pn & q_edge_table$child == cten_cn),]
      # Extract parent and child nodes
      clade_cten_cn <- extract.clade(q_rooted, cten_cn)
      clade_cten_pn <- extract.clade(q_rooted, cten_pn)
      # Check which clade has the right tips and use that node
      if (setequal(clade_cten_cn$tip.label, cten_taxa) == TRUE){
        # Extract child annotation from tree
        cten_cn_lab <- q_rooted$node.label[(cten_cn-Ntip(q_rooted))]
        # Clean string
        cten_node_value <- gsub("\\[|\\]|'", "",  cten_cn_lab)
      } else if (setequal(clade_cten_pn$tip.label, cten_taxa) == TRUE){
        # Extract parent annotation from tree
        cten_pn_lab <- q_rooted$node.label[(cten_pn-Ntip(q_rooted))]
        # Clean string
        cten_node_value <- gsub("\\[|\\]|'", "",  cten_pn_lab)
      } # END: if (setequal(clade_cten_cn$tip.label, cten_taxa) == TRUE){
      # Extract branch length
      cten_branch_length <- q_rooted$edge.length[which(q_rooted$edge[,1] == cten_pn & q_rooted$edge[,2] == cten_cn)]
      # Add description of clade
      cten_clade_monophyly    <- "Monophyletic"
    } else {
      # Assign NA if not monophyletic
      cten_branch_length      <- NA
      cten_node_value         <- NA
      cten_clade_monophyly    <- "Paraphyletic"
    } # END: if (is.monophyletic(q_rooted, cten_taxa) == TRUE){
  } else {
    # Assign NA if only 1 tip
    cten_branch_length  <- NA
    cten_node_value     <- NA
    cten_clade_monophyly    <- "Single_tip"
  }
  
  ## PORI (leading to PORI branch)
  # Identify tips in this group
  pori_taxa <- c(constraint_clades$Porifera)
  # Check number of taxa > 1
  if (length(pori_taxa) > 1){
    # Check monophyly
    if (is.monophyletic(q_rooted, pori_taxa) == TRUE){
      # Extract MRCA
      pori_cn <- getMRCA(q_rooted, pori_taxa) # child node
      pori_pn <- q_rooted$edge[which(q_rooted$edge[,2] == pori_cn), 1] # parent node
      # Extract the q_edge_table row
      pori_q_edge_row <- q_edge_table[which(q_edge_table$parent == pori_pn & q_edge_table$child == pori_cn),]
      # Extract parent and child nodes
      clade_pori_cn <- extract.clade(q_rooted, pori_cn)
      clade_pori_pn <- extract.clade(q_rooted, pori_pn)
      # Check which clade has the right tips and use that node
      if (setequal(clade_pori_cn$tip.label, pori_taxa) == TRUE){
        # Extract child annotation from tree
        pori_cn_lab <- q_rooted$node.label[(pori_cn-Ntip(q_rooted))]
        # Clean string
        pori_node_value <- gsub("\\[|\\]|'", "",  pori_cn_lab)
      } else if (setequal(clade_pori_pn$tip.label, pori_taxa) == TRUE){
        # Extract parent annotation from tree
        pori_pn_lab <- q_rooted$node.label[(pori_pn-Ntip(q_rooted))]
        # Clean string
        pori_node_value <- gsub("\\[|\\]|'", "",  pori_pn_lab)
      } # END: if (setequal(clade_pori_cn$tip.label, pori_taxa) == TRUE){
      # Extract branch length
      pori_branch_length <- q_rooted$edge.length[which(q_rooted$edge[,1] == pori_pn & q_rooted$edge[,2] == pori_cn)]
      # Add description of clade
      pori_clade_monophyly    <- "Monophyletic"
    } else {
      # Assign NA if not monophyletic
      pori_branch_length      <- NA
      pori_node_value         <- NA
      pori_clade_monophyly    <- "Paraphyletic"
    } # END: if (is.monophyletic(q_rooted, pori_taxa) == TRUE){
  } else {
    # Assign NA if only 1 tip
    pori_branch_length      <- NA
    pori_node_value         <- NA
    pori_clade_monophyly    <- "Single_tip"
  } # END:   if (length(pori_taxa) > 1){
  
  # Assemble output vector
  qcf_output <- c(met_branch_length, met_node_value, 
                  key_branch_length, key_node_value, 
                  cten_branch_length, cten_clade_monophyly, cten_node_value, 
                  pori_branch_length, pori_clade_monophyly,  pori_node_value)
  names(qcf_output) <- c("MET_branch_length", "MET_node_value", 
                         "KEY_branch_length", "KEY_node_value",
                         "CTEN_branch_length", "CTEN_monophyly", "CTEN_node_value", 
                         "PORI_branch_length", "PORI_monophyly", "PORI_node_value")
  return(qcf_output)
}



extract.CTENPORI.qcf <- function(dataset, matrix_name, topology, 
                                 tree_file, constraint_clades){
  # Function to extract qCF for any constrained tree topology (CTEN, PORI, CTEN_PORI)
  
  ## Open tree with qCF annotation
  q_tree <- read.tree(tree_file)
  # Root tree at outgroup
  q_rooted <- root(q_tree, constraint_clades$Outgroup)
  # Drop Placozoa taxa
  if (length(constraint_clades$Placozoa) > 0){
    q_rooted <- drop.tip(q_rooted, constraint_clades$Placozoa)
  }
  
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
  
  ## For the CTEN_PORI clade, identify key_check_node within the METAZOA clade
  key_check_node <- getMRCA(q_rooted, c(constraint_clades$Porifera, constraint_clades$Cnidaria))

  ## Key branch (leading to ALL OTHER ANIMALS aka PLAC+CNID+BILAT)
  # Specify key taxa
  key_taxa <- c(constraint_clades$Ctenophora, constraint_clades$Porifera)
  # Extract node value, if there's more than 1 taxon in the key_taxa clade
  if (length(key_taxa) > 1){
    # Extract MRCA
    key_cn <- getMRCA(q_rooted, key_taxa) # child node
    key_pn <- q_rooted$edge[which(q_rooted$edge[,2] == key_cn), 1] # parent node
    # Extract branch length
    key_branch_length <- q_rooted$edge.length[which(q_rooted$edge[,1] == key_pn & q_rooted$edge[,2] == key_cn)]
    # Identify correct node for the CTEN-PORI clade
    #     If the key_pn == key_check_node, then the key_cn is the correct node for the CTEN-PORI clade
    if (key_pn == key_check_node){
      # Extract child annotation from tree
      key_cn_lab <- q_rooted$node.label[(key_cn-Ntip(q_rooted))]
      # Clean string
      key_node_value <- gsub("\\[|\\]|'", "",  key_cn_lab)
    } else {
      # Extract the edges with the same key_pn (KEY parent node)
      key_edge_nodes  <- q_rooted$edge[c(which(q_rooted$edge[,1] == key_pn), which(q_rooted$edge[,2] == key_pn)), ]
      # Identify the possible child nodes
      test_edges <- c(which(q_rooted$edge[,1] %in% unique(as.numeric(key_edge_nodes))), which(q_rooted$edge[,2] %in% unique(as.numeric(key_edge_nodes))))
      test_nodes <- q_rooted$edge[test_edges,] # Extract all nodes connected to these nodes
      key_check_nodes <- sort(unique(as.numeric(test_nodes)), decreasing = TRUE) # Remove duplicates and order in decreasing order
      key_check_nodes <- key_check_nodes[which(key_check_nodes > Ntip(q_rooted))] # Remove external branches. Keep internal nodes.
      # Iterate and test the nodes
      key_node <- c()
      for (i in 1:length(key_check_nodes)){
        # Extract the tips in the clade defined by this node
        i_node          <- key_check_nodes[i]
        i_clade         <- extract.clade(q_rooted, i_node)
        i_tips          <- i_clade$tip.label
        i_tips_noPlac   <- i_tips[ which( ! (i_tips %in% constraint_clades$Placozoa) ) ]
        # Check whether the tips in this clade are the same as the tips in the CTEN_PORI clade
        if (setequal(i_tips_noPlac, key_taxa) == TRUE){
          key_node <- c(i_node)
        } # end: if (setequal(i_tips_noPlac, key_taxa) == TRUE){
      } # end: for (i in 1:length(key_check_nodes)){
      # Extract node annotation from tree
      key_lab <- q_rooted$node.label[(key_node-Ntip(q_rooted))]
      # Clean string
      key_node_value <- gsub("\\[|\\]|'", "",  key_lab)
    } # end: if (key_pn == key_check_node){ 
  } else {
    # Assign NA if only 1 tip
    key_branch_length  <- NA
    key_node_value     <- NA
  } # end: if (length(key_taxa) > 1){
  
  ## CTEN (leading to CTEN branch)
  # Identify tips in this group
  cten_taxa <- c(constraint_clades$Ctenophora)
  # Check that more than one taxa is present
  if (length(cten_taxa) > 1){
    # Check monophyly
    if (is.monophyletic(q_rooted, cten_taxa) == TRUE){
      # Extract MRCA
      cten_cn <- getMRCA(q_rooted, cten_taxa) # child node
      cten_pn <- q_rooted$edge[which(q_rooted$edge[,2] == cten_cn), 1] # parent node
      # Extract the q_edge_table row
      cten_q_edge_row <- q_edge_table[which(q_edge_table$parent == cten_pn & q_edge_table$child == cten_cn),]
      # Extract parent and child nodes
      clade_cten_cn <- extract.clade(q_rooted, cten_cn)
      clade_cten_pn <- extract.clade(q_rooted, cten_pn)
      # Check which clade has the right tips and use that node
      if (setequal(clade_cten_cn$tip.label, cten_taxa) == TRUE){
        # Extract child annotation from tree
        cten_cn_lab <- q_rooted$node.label[(cten_cn-Ntip(q_rooted))]
        # Clean string
        cten_node_value <- gsub("\\[|\\]|'", "",  cten_cn_lab)
      } else if (setequal(clade_cten_pn$tip.label, cten_taxa) == TRUE){
        # Extract parent annotation from tree
        cten_pn_lab <- q_rooted$node.label[(cten_pn-Ntip(q_rooted))]
        # Clean string
        cten_node_value <- gsub("\\[|\\]|'", "",  cten_pn_lab)
      } # END: if (setequal(clade_cten_cn$tip.label, cten_taxa) == TRUE){
      # Extract branch length
      cten_branch_length <- q_rooted$edge.length[which(q_rooted$edge[,1] == cten_pn & q_rooted$edge[,2] == cten_cn)]
      # Add description of clade
      cten_clade_monophyly    <- "Monophyletic"
    } else {
      # Assign NA if not monophyletic
      cten_branch_length      <- NA
      cten_node_value         <- NA
      cten_clade_monophyly    <- "Paraphyletic"
    } # END: if (is.monophyletic(q_rooted, cten_taxa) == TRUE){
  } else {
    # Assign NA if only 1 tip
    cten_branch_length  <- NA
    cten_node_value     <- NA
    cten_clade_monophyly    <- "Single_tip"
  }
  
  ## PORI (leading to PORI branch)
  # Identify tips in this group
  pori_taxa <- c(constraint_clades$Porifera)
  # Check number of taxa > 1
  if (length(pori_taxa) > 1){
    # Check monophyly
    if (is.monophyletic(q_rooted, pori_taxa) == TRUE){
      # Extract MRCA
      pori_cn <- getMRCA(q_rooted, pori_taxa) # child node
      pori_pn <- q_rooted$edge[which(q_rooted$edge[,2] == pori_cn), 1] # parent node
      # Extract the q_edge_table row
      pori_q_edge_row <- q_edge_table[which(q_edge_table$parent == pori_pn & q_edge_table$child == pori_cn),]
      # Extract parent and child nodes
      clade_pori_cn <- extract.clade(q_rooted, pori_cn)
      clade_pori_pn <- extract.clade(q_rooted, pori_pn)
      # Check which clade has the right tips and use that node
      if (setequal(clade_pori_cn$tip.label, pori_taxa) == TRUE){
        # Extract child annotation from tree
        pori_cn_lab <- q_rooted$node.label[(pori_cn-Ntip(q_rooted))]
        # Clean string
        pori_node_value <- gsub("\\[|\\]|'", "",  pori_cn_lab)
      } else if (setequal(clade_pori_pn$tip.label, pori_taxa) == TRUE){
        # Extract parent annotation from tree
        pori_pn_lab <- q_rooted$node.label[(pori_pn-Ntip(q_rooted))]
        # Clean string
        pori_node_value <- gsub("\\[|\\]|'", "",  pori_pn_lab)
      } # END: if (setequal(clade_pori_cn$tip.label, pori_taxa) == TRUE){
      # Extract branch length
      pori_branch_length <- q_rooted$edge.length[which(q_rooted$edge[,1] == pori_pn & q_rooted$edge[,2] == pori_cn)]
      # Add description of clade
      pori_clade_monophyly    <- "Monophyletic"
    } else {
      # Assign NA if not monophyletic
      pori_branch_length      <- NA
      pori_node_value         <- NA
      pori_clade_monophyly    <- "Paraphyletic"
    } # END: if (is.monophyletic(q_rooted, pori_taxa) == TRUE){
  } else {
    # Assign NA if only 1 tip
    pori_branch_length      <- NA
    pori_node_value         <- NA
    pori_clade_monophyly    <- "Single_tip"
  } # END:   if (length(pori_taxa) > 1){
  
  # Assemble output vector
  qcf_output <- c(key_branch_length, key_node_value, 
                  cten_branch_length, cten_clade_monophyly, cten_node_value, 
                  pori_branch_length, pori_clade_monophyly,  pori_node_value)
  names(qcf_output) <- c("KEY_branch_length", "KEY_node_value",
                         "CTEN_branch_length", "CTEN_monophyly", "CTEN_node_value", 
                         "PORI_branch_length", "PORI_monophyly", "PORI_node_value")
  return(qcf_output)
}



#### Reformat qCF dataframe ####
reformat.qCF.df <- function(input_df){
  ## Reformat the qCF dataframe into wide format
  
  # Identify rows for each topology
  cten_rows <- which(input_df$tree_topology == "CTEN")
  pori_rows <- which(input_df$tree_topology == "PORI")
  ctenpori_rows <- which(input_df$tree_topology == "CTEN_PORI")
  # Create new dataframe with columns for each topology
  new_df <- data.frame(id = unlist(lapply(strsplit(input_df$id[c(T,F,F)], "\\."), function(x){paste0(x[[1]], ".", x[[2]], ".", x[[3]])})),
                       dataset = input_df$dataset[c(T,F,F)],
                       matrix_name = input_df$matrix_name[c(T,F,F)],
                       dataset_id = input_df$dataset_id[c(T,F,F)],
                       model = input_df$model[c(T,F,F)],
                       CTEN.KEY_quartet_support = input_df$KEY_quartet_support[cten_rows],
                       PORI.KEY_quartet_support = input_df$KEY_quartet_support[pori_rows],
                       CTENPORI.KEY_quartet_support = input_df$KEY_quartet_support[ctenpori_rows],
                       CTEN.KEY_quartet_tree_freq = input_df$KEY_quartet_tree_freq[cten_rows],
                       PORI.KEY_quartet_tree_freq = input_df$KEY_quartet_tree_freq[pori_rows],
                       CTENPORI.KEY_quartet_tree_freq = input_df$KEY_quartet_tree_freq[ctenpori_rows],
                       CTEN.KEY_lpp = input_df$KEY_lpp[cten_rows],
                       PORI.KEY_lpp = input_df$KEY_lpp[pori_rows],
                       CTENPORI.KEY_lpp = input_df$KEY_lpp[ctenpori_rows],
                       CTEN.KEY_QC = input_df$KEY_QC[cten_rows],
                       PORI.KEY_QC = input_df$KEY_QC[pori_rows],
                       CTENPORI.KEY_QC = input_df$KEY_QC[ctenpori_rows],
                       CTEN.KEY_EN = input_df$KEY_EN[cten_rows],
                       PORI.KEY_EN = input_df$KEY_EN[pori_rows],
                       CTENPORI.KEY_EN = input_df$KEY_EN[ctenpori_rows],
                       CTEN.KEY_branch_length = input_df$KEY_branch_length[cten_rows],
                       PORI.KEY_branch_length = input_df$KEY_branch_length[pori_rows],
                       CTENPORI.KEY_branch_length = input_df$KEY_branch_length[ctenpori_rows],
                       CTEN.CTEN_quartet_support = input_df$CTEN_quartet_support[cten_rows],
                       PORI.CTEN_quartet_support = input_df$CTEN_quartet_support[pori_rows],
                       CTENPORI.CTEN_quartet_support = input_df$CTEN_quartet_support[ctenpori_rows],
                       CTEN.CTEN_quartet_tree_freq = input_df$CTEN_quartet_tree_freq[cten_rows],
                       PORI.CTEN_quartet_tree_freq = input_df$CTEN_quartet_tree_freq[pori_rows],
                       CTENPORI.CTEN_quartet_tree_freq = input_df$CTEN_quartet_tree_freq[ctenpori_rows],
                       CTEN.CTEN_lpp = input_df$CTEN_lpp[cten_rows],
                       PORI.CTEN_lpp = input_df$CTEN_lpp[pori_rows],
                       CTENPORI.CTEN_lpp = input_df$CTEN_lpp[ctenpori_rows],
                       CTEN.CTEN_QC = input_df$CTEN_QC[cten_rows],
                       PORI.CTEN_QC = input_df$CTEN_QC[pori_rows],
                       CTENPORI.CTEN_QC = input_df$CTEN_QC[ctenpori_rows],
                       CTEN.CTEN_EN = input_df$CTEN_EN[cten_rows],
                       PORI.CTEN_EN = input_df$CTEN_EN[pori_rows],
                       CTENPORI.CTEN_EN = input_df$CTEN_EN[ctenpori_rows],
                       CTEN.CTEN_branch_length = input_df$CTEN_branch_length[cten_rows],
                       PORI.CTEN_branch_length = input_df$CTEN_branch_length[pori_rows],
                       CTENPORI.CTEN_branch_length = input_df$CTEN_branch_length[ctenpori_rows],
                       CTEN.PORI_quartet_support = input_df$PORI_quartet_support[cten_rows],
                       PORI.PORI_quartet_support = input_df$PORI_quartet_support[pori_rows],
                       CTENPORI.PORI_quartet_support = input_df$PORI_quartet_support[ctenpori_rows],
                       CTEN.PORI_quartet_tree_freq = input_df$PORI_quartet_tree_freq[cten_rows],
                       PORI.PORI_quartet_tree_freq = input_df$PORI_quartet_tree_freq[pori_rows],
                       CTENPORI.PORI_quartet_tree_freq = input_df$PORI_quartet_tree_freq[ctenpori_rows],
                       CTEN.PORI_lpp = input_df$PORI_lpp[cten_rows],
                       PORI.PORI_lpp = input_df$PORI_lpp[pori_rows],
                       CTENPORI.PORI_lpp = input_df$PORI_lpp[ctenpori_rows],
                       CTEN.PORI_QC = input_df$PORI_QC[cten_rows],
                       PORI.PORI_QC = input_df$PORI_QC[pori_rows],
                       CTENPORI.PORI_QC = input_df$PORI_QC[ctenpori_rows],
                       CTEN.PORI_EN = input_df$PORI_EN[cten_rows],
                       PORI.PORI_EN = input_df$PORI_EN[pori_rows],
                       CTENPORI.PORI_EN = input_df$PORI_EN[ctenpori_rows],
                       CTEN.PORI_branch_length = input_df$PORI_branch_length[cten_rows],
                       PORI.PORI_branch_length = input_df$PORI_branch_length[pori_rows],
                       CTENPORI.PORI_branch_length = input_df$PORI_branch_length[ctenpori_rows] )
  # Return the reformatted dataframe
  return(new_df)
}


