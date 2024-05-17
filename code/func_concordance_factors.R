# ancient_ILS/code/func_concordance_factors.R
## This script includes functions to extract concordance factors for different datasets
# Caitlin Cherryh, 2023

library(ape)
library(dplyr)



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
  
  ## Extract the list of taxa in this dataset
  constraint_clades <- set.taxa(dataset_name = row_dataset, matrix_name = row$matrix_name,
                                all_datasets = all_datasets, matrix_taxa = matrix_taxa, alignment_taxa_df = alignment_taxa_df)
  
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
extract.qcf.wrapper <- function(analysis_id, qcf_df, 
                                matrix_taxa = matrix_taxa, all_datasets = all_datasets, alignment_taxa_df = alignment_taxa_df){
  ## Extract the single row from the qCF dataframe
  # Extract the relevant rows from the qcf_df
  id_split     <- strsplit(analysis_id, "\\.")[[1]]
  id_dataset   <- id_split[[1]]
  id_matrix    <- id_split[[2]]
  id_model     <- id_split[[3]]
  id_rows      <- qcf_df[which(qcf_df$dataset == id_dataset & qcf_df$matrix_name == id_matrix & qcf_df$model == id_model), ]
  
  ## Extract the list of taxa in this dataset
  constraint_clades <- set.taxa(dataset_name = id_dataset, matrix_name = id_matrix,
                                all_datasets = all_datasets, matrix_taxa = matrix_taxa, alignment_taxa_df = alignment_taxa_df)
  
  ## Extract CTEN_PORI clade values
  key_qcf <- extract.key.qcf(dataset = id_dataset, matrix_name = id_matrix, model = id_model, 
                             cten_tree_file = id_rows[which(id_rows$tree_topology == "CTEN"), "qcf_tree_file"], 
                             pori_tree_file = id_rows[which(id_rows$tree_topology == "PORI"), "qcf_tree_file"], 
                             ctenpori_tree_file = id_rows[which(id_rows$tree_topology == "CTEN_PORI"), "qcf_tree_file"], 
                             constraint_clades = constraint_clades)
  
  ## Extract qCF depending on topology
  cten_qcf          <- extract.CTEN.PORI.clade.qcf(dataset = id_dataset, matrix_name = id_matrix, topology = "CTEN",
                                                   tree_file = id_rows[which(id_rows$tree_topology == "CTEN"), "qcf_tree_file"], 
                                                   constraint_clades = constraint_clades)
  cten_qcf_values     <- format.CTEN.PORI.clade.qcf(cten_qcf)
  pori_qcf            <- extract.CTEN.PORI.clade.qcf(dataset = id_dataset, matrix_name = id_matrix, topology = "PORI",
                                                     tree_file = id_rows[which(id_rows$tree_topology == "PORI"), "qcf_tree_file"], 
                                                     constraint_clades = constraint_clades)
  pori_qcf_values     <- format.CTEN.PORI.clade.qcf(pori_qcf)
  ctenpori_qcf        <- extract.CTEN.PORI.clade.qcf(dataset = id_dataset, matrix_name = id_matrix, topology = "CTEN_PORI",
                                                     tree_file = id_rows[which(id_rows$tree_topology == "CTEN_PORI"), "qcf_tree_file"], 
                                                     constraint_clades = constraint_clades)
  ctenpori_qcf_values <- format.CTEN.PORI.clade.qcf(ctenpori_qcf)
  # Collate into a separate dataframe
  all_clade_qcf_df <- as.data.frame(matrix(data = c(cten_qcf_values, pori_qcf_values, ctenpori_qcf_values), nrow = 3, ncol = 14, byrow = TRUE))
  colnames(all_clade_qcf_df) <- names(ctenpori_qcf_values)
  
  ## Collate the KEY node dataframe with the CTEN/PORI node dataframe
  collated_output_df <- cbind(key_qcf, all_clade_qcf_df)
  
  ## Return the qCF values for this tree along with the tree parameters
  return(collated_output_df)
}



extract.key.qcf <- function(dataset, matrix_name, model,
                            cten_tree_file, pori_tree_file, ctenpori_tree_file, 
                            constraint_clades){
  # Function to extract the qCF for the key branches of the three topologies
  
  ## Open each tree and turn the nodes into a dataframe
  # Extract CTEN node labels
  cten_tree <- root(read.tree(cten_tree_file), constraint_clades$Outgroup)
  cten_output <- extract.node.details(in_tree = cten_tree, topology = "CTEN", constraint_clades = constraint_clades)
  cten_node_df <- cten_output$node_df[which(cten_output$node_df$node %in% cten_output$possible_nodes), ]
  cten_node_df$topology <- "CTEN"
  # Extract PORI node labels
  pori_tree <- root(read.tree(pori_tree_file), constraint_clades$Outgroup)
  pori_output <- extract.node.details(in_tree = pori_tree, topology = "PORI", constraint_clades = constraint_clades)
  pori_node_df <- pori_output$node_df[which(pori_output$node_df$node %in% pori_output$possible_nodes), ]
  pori_node_df$topology <- "PORI"
  # Extract CTEN_PORI node labels
  ctenpori_tree <- root(read.tree(ctenpori_tree_file), constraint_clades$Outgroup)
  ctenpori_output <- extract.node.details(in_tree = ctenpori_tree, topology = "CTEN_PORI", constraint_clades = constraint_clades)
  ctenpori_node_df <- ctenpori_output$node_df[which(ctenpori_output$node_df$node %in% ctenpori_output$possible_nodes), ]
  ctenpori_node_df$topology <- "CTEN_PORI"
  # Bind dataframes
  nodes_all           <- bind_rows(cten_node_df, pori_node_df, ctenpori_node_df)
  nodes_all$QC_EN_id  <- paste0(nodes_all$QC, ":", nodes_all$EN)
  # Extract any rows that have 3 identical QC and EN values
  candidate_QC_EN_id <- names(which(sort(table(nodes_all$QC_EN_id)) == 3))
  if (length(candidate_QC_EN_id) > 1){
    # Check all nodes with 3 matching QC/EN values to find the correct one by summing q1 values
    for (j in candidate_QC_EN_id){
      j_df <- nodes_all[which(nodes_all$QC_EN_id == j), ]
      j_qcf_sum <- sum(j_df$q1)
      if (round(j_qcf_sum, 2) == 1.00){
        j_node <- j
      }
    }
    # Extract the correct node
    key_node_df <- nodes_all[which(nodes_all$QC_EN_id == j_node), ]
  } else if (length(candidate_QC_EN_id) == 1){
    # Extract the correct node
    key_node_df <- nodes_all[which(nodes_all$QC_EN_id == candidate_QC_EN_id), ]
  } else if (length(candidate_QC_EN_id) == 0){
    # Print error value as no correct node found
    print(paste0("ERROR: ", dataset, " - ", matrix_name, " - ", model, " has no matching QC/EN nodes across all 3 topologies"))
  }
  # Format the dataframe
  key_node_output_df               <- key_node_df[c(which(key_node_df$topology == "CTEN"), 
                                                    which(key_node_df$topology == "PORI"), 
                                                    which(key_node_df$topology == "CTEN_PORI")), ]
  key_node_output_df$dataset        <- dataset
  key_node_output_df$matrix_name    <- matrix_name
  key_node_output_df$dataset_id     <- paste0(dataset, ".", matrix_name)
  key_node_output_df$model          <- model
  key_node_output_df$id             <- paste0(dataset, ".", matrix_name, ".", model, ".",  key_node_df$topolog)
  key_node_output_df$branch_length  <- c(cten_output$key_branch_length, pori_output$key_branch_length, ctenpori_output$key_branch_length)
  key_node_output_df                <- key_node_output_df[ , c("dataset", "matrix_name", "dataset_id", "model", "id", "topology", 
                                                               "node", "q1", "f1", "pp1", "QC", "EN", "branch_length")]
  names(key_node_output_df)         <- c("dataset", "matrix_name", "dataset_id", "model", "id", "topology",
                                         "KEY_node", "KEY_q1", "KEY_f1", "KEY_pp1", "KEY_QC", "KEY_EN", "KEY_branch_length")
  # Return the dataframe
  return(key_node_output_df)
}



extract.node.details <- function(in_tree, topology, constraint_clades){
  # Extract all nodes and node labels in the tree
  node_df <- data.frame("node" = c(in_tree$edge[,1], in_tree$edge[,2]), 
                        "node.lab" = c(sapply(in_tree$edge[,1],
                                              select.tip.or.node,
                                              tree = in_tree),
                                       sapply(in_tree$edge[,2],
                                              select.tip.or.node,
                                              tree = in_tree) ) )
  node_df <- unique(node_df) # Remove duplicate rows
  node_df <- node_df[grep("q1=", node_df$node.lab),] # Get internal nodes
  node_df$node.lab <- gsub("\\[|\\]|'", "",  node_df$node.lab) # Remove extra punctuation
  # Break out values into separate columns
  split_labs <- strsplit(node_df$node.lab, ";")
  node_df$q1   <- as.numeric(unlist(lapply(strsplit(unlist(lapply(split_labs, function(x){x[[1]]})), "\\="), function(x){x[[2]]})))
  node_df$f1   <- as.numeric(unlist(lapply(strsplit(unlist(lapply(split_labs, function(x){x[[4]]})), "\\="), function(x){x[[2]]})))
  node_df$pp1  <- as.numeric(unlist(lapply(strsplit(unlist(lapply(split_labs, function(x){x[[7]]})), "\\="), function(x){x[[2]]})))
  node_df$QC   <- as.numeric(unlist(lapply(strsplit(unlist(lapply(split_labs, function(x){x[[10]]})), "\\="), function(x){x[[2]]})))
  node_df$EN   <- as.numeric(unlist(lapply(strsplit(unlist(lapply(split_labs, function(x){x[[11]]})), "\\="), function(x){x[[2]]})))
  # Extract all edges 
  edge_df <- data.frame("parent" = in_tree$edge[,1],
                        "par.name" = sapply(in_tree$edge[,1],
                                            select.tip.or.node,
                                            tree = in_tree),
                        "child" = in_tree$edge[,2],
                        "chi.name" = sapply(in_tree$edge[,2],
                                            select.tip.or.node,
                                            tree = in_tree) )
  edge_df <- edge_df[sort(c(grep("q1=", edge_df$chi.name), which(edge_df$chi.name == ""))),] # Get internal parent nodes 
  edge_df$par.name <- gsub("\\[|\\]|'", "",  edge_df$par.name) # Remove extra punctuation in parent node label
  edge_df$chi.name <- gsub("\\[|\\]|'", "",  edge_df$chi.name) # Remove extra punctuation in child node label
  
  # Test for which nodes to exclude
  if (topology == "CTEN"){
    key_taxa <- c(constraint_clades$Porifera, constraint_clades$Cnidaria, constraint_clades$Bilateria)
  } else if (topology == "PORI"){
    key_taxa <- c(constraint_clades$Ctenophora, constraint_clades$Cnidaria, constraint_clades$Bilateria)
  } else if (topology == "CTEN_PORI"){
    key_taxa <- c(constraint_clades$Ctenophora, constraint_clades$Porifera)
  }
  # Get the MRCA for the key taxa
  key_node <- getMRCA(in_tree, key_taxa)
  # Get nodes attached to the key node
  key_branch <- which(in_tree$edge[,2] == key_node)
  key_branch_length <- in_tree$edge.length[key_branch]
  # Extract possible nodes
  possible_nodes <- unique(as.numeric(in_tree$edge[which(in_tree$edge[,1] == key_node | in_tree$edge[,2] == key_node),]))
  # Construct output
  output <- list("node_df" = node_df,
                 "edge_df" = edge_df,
                 "possible_nodes" = possible_nodes,
                 "key_branch_length" = key_branch_length)
  # Return the node df
  return(output)
}


format.CTEN.PORI.clade.qcf <- function(output_string){
  # Reformat extract.CTEN.PORI.clade.qcf output
  
  # Split up the vector
  cten_clade_lab  <- output_string[["CTEN_node_value"]]
  cten_split      <- strsplit(cten_clade_lab, ";")[[1]]
  pori_clade_lab  <- output_string[["PORI_node_value"]]
  pori_split      <- strsplit(pori_clade_lab, ";")[[1]]
  # Create the output vectors for CTEN clade, if there is more than one taxon in the clade
  if (is.na(output_string[["CTEN_node_value"]]) == TRUE){
    cten_output_vector <- c(rep(NA, 5), output_string[["CTEN_branch_length"]], output_string[["CTEN_monophyly"]])
  } else {
    cten_output_vector <- c(gsub("q1=", "", cten_split[[1]]), gsub("f1=", "", cten_split[[4]]), gsub("pp1=", "", cten_split[[7]]),
                            gsub("QC=", "", cten_split[[10]]), gsub("EN=", "", cten_split[[11]]), 
                            output_string[["CTEN_branch_length"]], output_string[["CTEN_monophyly"]])
  }
  # Create the output vectors for PORI cladess if there is more than one taxon in the clade
  if (is.na(output_string[["PORI_node_value"]]) == TRUE){
    pori_output_vector <- c(rep(NA, 5), output_string[["PORI_branch_length"]], output_string[["PORI_monophyly"]])
  } else {
    pori_output_vector <- c(gsub("q1=", "", pori_split[[1]]), gsub("f1=", "", pori_split[[4]]), gsub("pp1=", "", pori_split[[7]]),
                            gsub("QC=", "", pori_split[[10]]), gsub("EN=", "", pori_split[[11]]), 
                            output_string[["PORI_branch_length"]], output_string[["PORI_monophyly"]])
  }
  # Assemble the output vector
  node_value_names <- c("q1", "f1", "pp1", "QC", "EN", "branch_length", "monophyly")
  output_vector <- c(cten_output_vector, pori_output_vector)
  names(output_vector) <- paste0(rep(c("CTEN_", "PORI_"), each = 7), rep(node_value_names, 2))
  # Return the output vector
  return(output_vector)
}



extract.CTEN.PORI.clade.qcf <- function(dataset, matrix_name, topology, model, 
                                        tree_file, constraint_clades){
  # Function to extract qCF for the CTEN and PORI clades for any constrained tree topology (CTEN, PORI, CTEN_PORI)
  
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
  qcf_output <- c(cten_branch_length, cten_clade_monophyly, cten_node_value, 
                  pori_branch_length, pori_clade_monophyly,  pori_node_value)
  names(qcf_output) <- c("CTEN_branch_length", "CTEN_monophyly", "CTEN_node_value", 
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



#### Extract list of taxa within the dataset ####
set.taxa <- function(dataset_name, matrix_name,
                     all_datasets, matrix_taxa, alignment_taxa_df){
  ## Identify the taxa in each clade using details about the dataset
  
  ## Extract the relevant list of taxa for this dataframe
  # First, check whether this matrix is included in the keys of the matrix_taxa list
  row_key   <- paste0(dataset_name, ".", matrix_name, ".", "aa")
  # Get the list of possible keys
  list_keys <- names(matrix_taxa)
  # Check if row_key in list_key
  if ((row_key %in% list_keys) == FALSE){
    # If row key is not in list key, then all taxa for this dataset have the same names
    # Extract the object containing those taxa names
    constraint_clades   <- all_datasets[[dataset_name]]
  } else if ((row_key %in% list_keys) == TRUE){
    # First, identify the list of taxa in this matrix
    keep_taxa           <- matrix_taxa[[row_key]]
    # Secondly, extract the taxa clades for this dataset
    dataset_taxa_clades <- all_datasets[[dataset_name]]
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
  al_col_key        <- paste0(dataset_name, ".", matrix_name)
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
  
  ## Return the constraint clades
  return(constraint_clades)
}



#### Remove Placozoa taxa ####
tree.remove.Plac <- function(tree_file, 
                             all_datasets, matrix_taxa, alignment_taxa_df){
  # Remove Plac taxa from a single ML tree
  
  ## Extract details about the alignment
  file_split      <- strsplit(basename(tree_file), "\\.")[[1]]
  file_dataset    <- file_split[[1]]
  file_matrix     <- file_split[[2]]
  
  ## Extract tips within the tree
  constraint_clades <- set.taxa(dataset_name = file_dataset, matrix_name = file_matrix,
                                all_datasets = all_datasets, matrix_taxa = matrix_taxa, alignment_taxa_df = alignment_taxa_df)
  
  ## Remove Plac tips from tree
  full_tree             <- read.tree(tree_file)
  if (length(constraint_clades$Placozoa) > 0){
    # Trim internal branches in drop.tip function - otherwise will get NA nodes when multiple PLAC taxa are removed
    noPlac_tree         <- drop.tip(full_tree, constraint_clades$Placozoa, trim.internal = TRUE)
  } else {
    noPlac_tree         <- full_tree
  }
  
  ## Remove node labels in case they get distorted by moving branches around
  noPlac_tree$node.label <- NULL
  
  ## Root tree and save
  noPlac_tree_rooted      <- root(noPlac_tree, constraint_clades$Outgroup)
  noPlac_tree_rooted_file_path   <- gsub("\\.treefile", ".noPlac_rooted.treefile", tree_file)
  write.tree(noPlac_tree_rooted, file = noPlac_tree_rooted_file_path)
  
  ## Save and return unrooted noPlac tree
  new_file_path   <- gsub("\\.treefile", ".noPlac.treefile", tree_file)
  write.tree(noPlac_tree, file = new_file_path)
  return(new_file_path)
}



gene.trees.remove.Plac <- function(gene_tree_file, 
                                   all_datasets, matrix_taxa, alignment_taxa_df){
  # Remove Plac taxa from a single ML tree
  
  ## Extract details about the alignment
  file_split      <- strsplit(basename(gene_tree_file), "\\.")[[1]]
  file_dataset    <- file_split[[1]]
  file_matrix     <- file_split[[2]]
  
  ## Extract tips within the tree
  constraint_clades <- set.taxa(dataset_name = file_dataset, matrix_name = file_matrix,
                                all_datasets = all_datasets, matrix_taxa = matrix_taxa, alignment_taxa_df = alignment_taxa_df)
  
  ## Remove Plac tips from tree
  full_multiphylo         <- read.tree(gene_tree_file)
  # Trim internal branches in drop.tip function - otherwise will get NA nodes when multiple PLAC taxa are removed
  full_multiphylo_noPlac  <- lapply(1:length(full_multiphylo), 
                                    function(x){noPlac.gene.trees(full_multiphylo[[x]], constraint_clades)})
  # Change class from list (result of lapply) to multiPhylo
  class(full_multiphylo_noPlac) <- "multiPhylo"
  
  ## Save updated tree
  new_file_path   <- gsub("\\.treefile", ".noPlac.treefile", gene_tree_file)
  write.tree(full_multiphylo_noPlac, file = new_file_path)
  return(new_file_path)
}


noPlac.gene.trees <- function(gene_tree, constraint_clades){
  # Remove PLAC from gene trees
  
  # Check if any PLAC present in gene tree
  if (length(constraint_clades$Placozoa) > 0){
    # Identify which PLAC taxa present in gene tree
    gene_tree_taxa <- gene_tree$tip.label
    plac_to_remove <- constraint_clades$Placozoa[which(constraint_clades$Placozoa %in% gene_tree_taxa)]
    if (length(plac_to_remove) > 0){
      # Remove Placazoa taxa
      # Trim internal branches in drop.tip function - otherwise will get NA nodes when multiple PLAC taxa are removed
      noPlac_gene_tree <- drop.tip(gene_tree, plac_to_remove, trim.internal = TRUE)
    } else {
      # No Placazoa taxa in this gene tree- return unchanged
      noPlac_gene_tree = gene_tree
    }
  } else {
    # No Placozoa in dataset - return unchanged
    noPlac_gene_tree = gene_tree
  }
  
  # Return the gene tree
  return(noPlac_gene_tree)
}




