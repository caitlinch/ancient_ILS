## caitlinch/ancient_ILS/code/func_constraint_trees.R
# This script estimates maximum likelihood trees under different models of substitution for 14 empirical data sets
# Caitlin Cherryh 2024


constraint.tree.wrapper <- function(i, output_directory, dataset_info, matrix_taxa_info, 
                                    constraint_tree_df, alignment_taxa_df, 
                                    force.update.constraint.trees = TRUE){
  
  # Function to take a row from the ML output dataframe and create the constraint trees plus parameters to estimate the hypothesis trees
  
  ## Extract the relevant row
  row <- constraint_tree_df[i, ]
  
  ## Extract the relevant list of taxa for this dataframe
  # First, check whether this matrix is included in the keys of the matrix_taxa list
  row_dataset <- row$dataset
  row_key <- paste0(row$dataset, ".", row$matrix, ".", "aa")
  list_keys <- names(matrix_taxa_info)
  # Check if row_key in list_key
  if ((row_key %in% list_keys) == FALSE){
    # If row key is not in list key, then all taxa for this dataset have the same names
    # Extract the object containing those taxa names
    constraint_clades <- dataset_info[[row_dataset]]
  } else if ((row_key %in% list_keys) == TRUE){
    # First, identify the list of taxa in this matrix
    keep_taxa <- matrix_taxa_info[[row_key]]
    # Secondly, extract the taxa clades for this dataset
    dataset_taxa_clades <- dataset_info[[row_dataset]]
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
    constraint_clades$Sponges_Calcarea <- dataset_taxa_clades$Sponges_Calcarea[which(dataset_taxa_clades$Sponges_Calcarea %in% keep_taxa)]
    constraint_clades$Sponges_Homoscleromorpha <- dataset_taxa_clades$Sponges_Homoscleromorpha[which(dataset_taxa_clades$Sponges_Homoscleromorpha %in% keep_taxa)]
    constraint_clades$Sponges_Hexactinellida <- dataset_taxa_clades$Sponges_Hexactinellida[which(dataset_taxa_clades$Sponges_Hexactinellida %in% keep_taxa)]
    constraint_clades$Sponges_Demospongiae <- dataset_taxa_clades$Sponges_Demospongiae[which(dataset_taxa_clades$Sponges_Demospongiae %in% keep_taxa)]
  }
  
  ## Remove any taxa from the constraint_clades that are not included in the ML tree for the alignment
  # Extract the column of tip labels from the relevant unique_id column of the alignment_taxa_df
  al_col_key <- paste0(row$dataset, ".", row$matrix)
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
  constraint_clades$Sponges_Calcarea <- constraint_clades$Sponges_Calcarea[(constraint_clades$Sponges_Calcarea %in% tree_tips)]
  constraint_clades$Sponges_Homoscleromorpha <- constraint_clades$Sponges_Homoscleromorpha[(constraint_clades$Sponges_Homoscleromorpha %in% tree_tips)]
  constraint_clades$Sponges_Hexactinellida <- constraint_clades$Sponges_Hexactinellida[(constraint_clades$Sponges_Hexactinellida %in% tree_tips)]
  constraint_clades$Sponges_Demospongiae <- constraint_clades$Sponges_Demospongiae[(constraint_clades$Sponges_Demospongiae %in% tree_tips)]
  
  ## Create the constraint tree dataframe by calling the create.constraint.trees function
  constraint_files <- create.constraint.trees(dataset = row$dataset, 
                                              matrix_name = row$matrix,
                                              dataset_constraint_tree_dir = output_directory,
                                              constraint_clades = constraint_clades,
                                              force.update.constraint.trees = TRUE)
  ## Return the output files
  return(constraint_files)
} # end function



create.constraint.trees <- function(dataset, matrix_name, dataset_constraint_tree_dir, 
                                    constraint_clades, force.update.constraint.trees = TRUE){
  # Function to create the constraint trees and constraint tree information data frame, for a given dataset and model
  # Does not include Placozoa - Placozoa is not relevant to the question I am investigating (and is only ~1 taxa in most data sets)
  
  ### Construct constraint trees
  # Generate file names for all 5 constraint trees
  constraint_tree_1_file_name <- paste0(dataset_constraint_tree_dir, dataset, ".", matrix_name, ".constraint_tree_", "1", ".nex")
  constraint_tree_2_file_name <- paste0(dataset_constraint_tree_dir, dataset, ".", matrix_name, ".constraint_tree_", "2", ".nex")
  constraint_tree_3_file_name <- paste0(dataset_constraint_tree_dir, dataset, ".", matrix_name, ".constraint_tree_", "3", ".nex")
  
  # Split the taxa into clades
  outgroup_taxa = constraint_clades$Outgroup
  ctenophora_taxa = constraint_clades$Ctenophora
  porifera_taxa = constraint_clades$Porifera
  cnidaria_taxa = constraint_clades$Cnidaria
  bilateria_taxa = constraint_clades$Bilateria
  
  # Format the outgroup clades nicely
  outgroup_taxa_formatted <- format.constraint.tree.clade(outgroup_taxa)
  ctenophora_taxa_formatted <- format.constraint.tree.clade(ctenophora_taxa)
  porifera_taxa_formatted <- format.constraint.tree.clade(porifera_taxa)
  ctenophora_porifera_taxa_formatted <- format.constraint.tree.clade(c(ctenophora_taxa, porifera_taxa))
  cnidaria_taxa_formatted <- format.constraint.tree.clade(cnidaria_taxa)
  bilateria_taxa_formatted <- format.constraint.tree.clade(bilateria_taxa)
  cnidaria_bilateria_taxa_formatted <- format.constraint.tree.clade(c(cnidaria_taxa, bilateria_taxa))
  
  # Only need to create one set of hypothesis trees per dataset/matrix combination - create constraint tree if it doesn't exist
  ### Hypothesis 1: Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (porifera_taxa, (cnidaria_taxa, bilateria_taxa))));
  if ((file.exists(constraint_tree_1_file_name) == FALSE) | (force.update.constraint.trees == TRUE)){
    ## Construct constraint tree
    constraint_tree_1 <- paste0("(", outgroup_taxa_formatted, ", (", ctenophora_taxa_formatted, ", (", porifera_taxa_formatted, ", ", cnidaria_bilateria_taxa_formatted, ")));")
    write(constraint_tree_1, file = constraint_tree_1_file_name)
  } else if (file.exists(constraint_tree_1_file_name) == TRUE){ 
    constraint_tree_1 <- readLines(constraint_tree_1_file_name)
  }
  
  ### Hypothesis 2: Porifera-sister
  # Tree: (outgroup_taxa, (porifera_taxa, (ctenophora_taxa, (cnidaria_taxa, bilateria_taxa))));
  if ((file.exists(constraint_tree_2_file_name) == FALSE) | (force.update.constraint.trees == TRUE)){
    ## Construct constraint tree
    constraint_tree_2 <- paste0("(", outgroup_taxa_formatted, ", (", porifera_taxa_formatted, ", (", ctenophora_taxa_formatted, ", ", cnidaria_bilateria_taxa_formatted, ")));")
    write(constraint_tree_2, file = constraint_tree_2_file_name)
  } else if (file.exists(constraint_tree_2_file_name) == TRUE){ 
    constraint_tree_2 <- readLines(constraint_tree_2_file_name)
  }
  
  ### Hypothesis 3: (Ctenophore+Porifera)-sister
  # Tree: (outgroup_taxa, ((porifera_taxa, ctenophora_taxa), (cnidaria_taxa, bilateria_taxa)));
  if ((file.exists(constraint_tree_3_file_name) == FALSE) | (force.update.constraint.trees == TRUE)){
    ## Construct constraint tree
    constraint_tree_3 <- paste0("(", outgroup_taxa_formatted, ", ((", ctenophora_taxa_formatted, ", ", porifera_taxa_formatted, "), ", cnidaria_bilateria_taxa_formatted, "));")
    write(constraint_tree_3, file = constraint_tree_3_file_name)
  } else if (file.exists(constraint_tree_3_file_name) == TRUE){ 
    constraint_tree_3 <- readLines(constraint_tree_3_file_name)
  }
  
  # Output constraint tree files
  cf_files <- c(constraint_tree_1_file_name, constraint_tree_2_file_name, constraint_tree_3_file_name)
  return(cf_files)
} # end function



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


