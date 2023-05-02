# ancient_ILS/code/func_simulations.R
## This script includes functions to run the simulations for this project
# Caitlin Cherryh, 2023

library(ape)
library(phytools)

run.one.simulation <- function(sim_row){
  # Function to take one row from the simulation parameters dataframe and run start to finish
  
  # Find the nodes for the branch depending on which tree is being used
  if (basename(sim_row$hypothesis_tree_file) == "Whelan2017_hypothesis_tree_1_Cten.treefile"){
    a_start = 91
    a_end = 90
    b_start = 91
    b_end = 92
  } else if (basename(sim_row$hypothesis_tree_file) == "Whelan2017_hypothesis_tree_2_Pori.treefile"){
    a_start = NA
    a_end = NA
    b_start = NA
    b_end = NA
  } else if (basename(sim_row$hypothesis_tree_file) == "Whelan2017_hypothesis_tree_3_CtenPori.treefile"){
    a_start = NA
    a_end = NA
    b_start = NA
    b_end = NA
  }
  
  # Create the folder for this replicate
  if (dir.exists(sim_row$output_folder) == FALSE){
    dir.create(sim_row$output_folder)
  }
  
  # Copy the hypothesis tree to the folder
  hyp_tree_file = paste0(sim_row$output_folder, sim_row$ID, ".treefile")
  file.copy(from = sim_row$hypothesis_tree_file, to = hyp_tree_file)
  
  # Open the hypothesis tree
  tree <- read.tree(hyp_tree_file)
  # Root the hypothesis tree
  rooted_tree <- root(tree, outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"))
  # Make the tree ultrametric
  rooted_tree <- force.ultrametric(rooted_tree, method = "extend")
  # Scale rooted tree to be tree age
  rooted_tree$edge.length <- rooted_tree$edge.length * (sim_row$tree_length / max(branching.times(rooted_tree)))
  # Find the nodes for the branch depending on which tree is being used
  if (basename(sim_row$hypothesis_tree_file) == "Whelan2017_hypothesis_tree_1_Cten.treefile"){
    # Extract the nodes by checking for monophyletic clades
    a_end <- getMRCA(rooted_tree, c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster","Daphnia_pulex",
                                          "Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera", 
                                          "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona","Craseo_lathetica",
                                          "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla", "Trichoplax_adhaerens", "Cliona_varians", "Sycon_coactum", "Sycon_ciliatum",
                                          "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum", "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux",
                                          "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica", "Petrosia_ficiformis", "Spongilla_lacustris", 
                                          "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", "Crella_elegans", "Kirkpatrickia_variolosa"))
    a_start <- rooted_tree$edge[which(rooted_tree$edge[,2] == a_end),1]
    b_end <- getMRCA(rooted_tree, c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus",
                                          "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                                          "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                                          "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                                          "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                                          "Ctenophora_sp_Florida_USA"))
    b_start <- rooted_tree$edge[which(rooted_tree$edge[,2] == b_end),1]
  } else if (basename(sim_row$hypothesis_tree_file) == "Whelan2017_hypothesis_tree_2_Pori.treefile"){
    a_start = NA
    a_end = NA
    b_start = NA
    b_end = NA
  } else if (basename(sim_row$hypothesis_tree_file) == "Whelan2017_hypothesis_tree_3_CtenPori.treefile"){
    a_start = NA
    a_end = NA
    b_start = NA
    b_end = NA
  }
  # Identify branch a
  branch_a <- which(rooted_tree$edge[,1] == a_start & rooted_tree$edge[,2] == a_end)
  # Identify branch b
  branch_b <- which(rooted_tree$edge[,1] == b_start & rooted_tree$edge[,2] == b_end)
  # Modify branch b length
  if (sim_row$simulation_number == "sim1" | sim_row$simulation_number == "sim3"){
    # LBA simulation - vary branch b (branch that leads to Ctenophore clade)
    branch_b_value <- rooted_tree$edge.length[branch_b]
    # Multiply existing branch to be that percent of the current tree height
    # e.g. if branch_b_percent_height is 30% and the tree is 1.0 sub/site long, then branch b should be 0.33 sub/site long
    new_branch_b_value <- max(branching.times(rooted_tree)) * (sim_row$branch_b_percent_height/100)
    rooted_tree$edge.length[branch_b] <- new_branch_b_value
  }
  # Modify branch a length
  if (sim_row$simulation_number == "sim2" | sim_row$simulation_number == "sim3"){
    # ILS simulation - vary branch a (branch that allows more time for the two species to differentiate)
    branch_a_value <- rooted_tree$edge.length[branch_a]
    # Multiply existing branch to be that percent of the current tree height
    # e.g. if branch_a_percent_height is 50% and the tree is 1.0 sub/site long, then branch a should be 0.5 sub/site long
    new_branch_a_value <- max(branching.times(rooted_tree)) * (sim_row$branch_a_percent_height/100)
    rooted_tree$edge.length[branch_a] <- new_branch_a_value
  }
  # Scale rooted tree to be tree age (reset impacts from modifying branch lengths)
  rooted_tree$edge.length <- rooted_tree$edge.length * (sim_row$tree_length / max(branching.times(rooted_tree)))
  # Save the tree
  write.tree(rooted_tree, file = hyp_tree_file)
  
  # Write the tree in ms command line format
  ms_files <- ms.generate.trees(unique_id = sim_row$ID, base_tree = rooted_tree, ntaxa = sim_row$num_taxa, 
                                ntrees = sim_row$num_genes, tree_depth = sim_row$tree_length, 
                                output_directory = sim_row$output_folder, ms_path = sim_row$ms)
  
}



#### Functions for ms ####
ms.generate.trees <- function(unique_id, base_tree, ntaxa, ntrees, tree_depth, output_directory, ms_path = "ms"){
  ## Randomly generate a tree with n taxa; format into an ms command and run ms; generate and save the resulting gene trees
  
  ## Generate file paths using the unique id
  t_path <- paste0(output_directory, unique_id, "_starting_tree.txt")
  ms_op_path <- paste0(output_directory, unique_id, "_ms_output.txt")
  ms_gene_trees_path <- paste0(output_directory, unique_id, "_ms_gene_trees.txt")
  
  ## Convert the tree into the format for ms
  # Calculate times for ms -ej commands by finding coalescence times (coalescent intervals found using ape::coalescent.intervals)
  ms_coal_ints <- calculate.ms.coalescent.times(base_tree$Nnode, coalescent.intervals(base_tree))
  # Determine the nodes that lead to non-terminal branches {e.g. which(node.depth(t) != 1) }
  nodes <- (ntaxa+1):(ntaxa+base_tree$Nnode)
  # Extract information about all clades from tree
  node_df <- do.call(rbind.data.frame, lapply(nodes, extract.clade.from.node, tree = base_tree, coalescent_times = ms_coal_ints))
  names(node_df) <- c("node", "tip_names", "tip_numbers", "ms_tip_order", "ntips", "ndepth", "coalescence_time", "removed_taxa", "ms_input")
  # Format coalescences for ms input
  node_df <- determine.coalescence.taxa(node_df)
  # Create a new column containing -ej event for each row
  node_df$ej <- paste0("-ej ", node_df$coalescence_time, " ", node_df$ms_input)
  # No recombination event is present. Do not add any extra splitting (-es) or joining (-ej) events
  
  ## Generate gene trees in ms
  # Paste together all the -ej coalescence events for this tree
  all_ej <- paste(node_df$ej, collapse = " ")
  # Construct the ms command line using the -ej events
  coal_call <- paste0(ms_path, " ", ntaxa, " ", ntrees, " -T -I ", ntaxa," ", paste(rep(1, ntaxa), collapse = " "), " ", all_ej)
  # Write all output to file
  write(ms_op, file = ms_op_path)
  # Call ms
  ms_op <- system(coal_call, intern = TRUE)
  
  ## Format and save gene trees
  # Remove non-gene tree lines from the ms output
  ms_txt <- ms_op[3:length(ms_op)] # Remove first two lines (ms command and random seeds lines)
  ms_txt <- ms_txt[which(ms_txt != "")] # Remove empty lines
  ms_txt <- ms_txt[grep("//", ms_txt, invert = TRUE)] # Remove separation lines between gene trees ("//")
  write(ms_txt, file = ms_gene_trees_path)
  
  ## Return file paths for output 
  output_files <- c(t_path, ms_op_path, ms_gene_trees_path)
  names(output_files) <- c("starting_tree_file", "ms_output_file", "ms_gene_tree_file")
  return(output_files)
}



extract.clade.from.node <- function(node, tree, coalescent_times){
  ## Small function to take a node, extract the clade from that node, and return the number and names of taxa in that node
  
  # Extract clade
  clade <- extract.clade(tree, node)
  # Extract information about clade
  tip_names <- clade$tip.label
  tip_numbers <- gsub("t", "", tip_names)
  tip_order <- tip_numbers[order(as.numeric(tip_numbers), decreasing = TRUE)]
  ntips <- length(clade$tip.label)
  # Determine depth of this node (how many species does this node contain)
  ndepth <- node.depth(tree)[node]
  # Determine which coalescent time is associated with this node
  n_tree_tips <- length(tree$tip.label)
  coal_index <- node - n_tree_tips # node numbering for non-trivial tips starts at n_tree_tips+1
  coal_time <- coalescent_times[coal_index]
  # Determine which taxa to remove
  if (ndepth == 2){
    removed_taxa = tip_order[1]
    ms_input = paste0(tip_order[1], " ", tip_order[2])
  } else {
    removed_taxa = NA
    ms_input = NA
  }
  # Assemble results into a vector
  o <- c(node = node, tip_names = paste(tip_names, collapse = ","), tip_numbers = paste(tip_numbers, collapse = ","), 
         ms_tip_order = paste(tip_order, collapse = ","), ntips = ntips, ndepth = ndepth, coalescence_time = coal_time,
         removed_taxa = removed_taxa, ms_input = ms_input)
  # Return vector
  return(o)
}


calculate.ms.coalescent.times <- function(number_of_nodes, coalescent_intervals){
  ## Small function to take a number of nodes and determine all the coalescent times needed to run ms
  
  ints <- coalescent_intervals$interval.length
  # The interval length is the length between two coalescent events: to find the time for e.g. the second event, add the first and second interval together
  # The last interval should be the same as the total depth
  times_vec <- c()
  for (i in 1:number_of_nodes){
    temp_time <- sum(ints[1:i])
    times_vec <- c(times_vec, temp_time)
  }
  
  # Round to 6dp (to allow for values as small as 0.000001)
  times_vec <- round(times_vec, digits = 6)
  # Reverse vector so that the longest time aligns with the deepest node
  times_vec <- rev(times_vec)
  # Return coalescent times
  return(times_vec)
}


determine.coalescence.taxa <- function(node_dataframe){
  ## Take the node dataframe and work out which taxa will be coalescing into which (essential for the ms command line)
  
  # Convert columns used for filtering into numeric so they can be ordered properly
  node_dataframe$node <- as.numeric(node_dataframe$node)
  node_dataframe$ntips <- as.numeric(node_dataframe$ntips)
  node_dataframe$ndepth <- as.numeric(node_dataframe$ndepth)
  node_dataframe$coalescence_time <- as.numeric(node_dataframe$coalescence_time)
  # Order dataframe by node depth value
  node_dataframe <- node_dataframe[order(node_dataframe$ndepth),]
  # Make a list of all the taxa to remove
  removed_taxa <- as.numeric(node_dataframe$removed_taxa[!is.na(node_dataframe$removed_taxa)])
  # Iterate through the dataframe row by row to check 
  rows_to_process = which(is.na(node_dataframe$ms_input))
  for (i in rows_to_process){
    # Extract row
    row <- node_dataframe[i,]
    # Identify taxa in row
    row_tips <- as.numeric(unlist(strsplit(row$ms_tip_order, ",")))
    # Remove any taxa in removed_taxa
    keep_tips <- setdiff(row_tips, removed_taxa)
    ordered_keep_tips <- keep_tips[order(keep_tips, decreasing = TRUE)]
    row_ms_input <- paste0(ordered_keep_tips[1], " ", ordered_keep_tips[2])
    row_removed_tips <- c(removed_taxa, ordered_keep_tips[1])
    # Attach results back to row
    row$ms_input <- row_ms_input
    row$removed_taxa <- paste(row_removed_tips, collapse = ",")
    # Attach row back to dataframe
    node_dataframe[i,] <- row
    # Add removed taxa to list of removed taxa and remove duplicates
    removed_taxa <- unique(c(removed_taxa, row_removed_tips))
    # Repeat process with next row: the list of removed taxa will continue to grow as the node depth increases
  }
  
  # Each column should now have a nicely formatted ms_input column consisting of two lineages, with the smallest number lineage second
  #     meaning the lineages in the first population will be moved into the second subpopulation 
  #     (in the forward direction, this is equivalent to a population splitting)
  # Now, reorder the rows from longest to shortest coalescence time (in ms, the events should be put in from longest to shortest coalescent time,
  #     with the root on the left of the command and the tips on the right of the command)
  node_dataframe <- node_dataframe[order(node_dataframe$coalescence_time, decreasing = TRUE),]
  # Return the dataframe
  return(node_dataframe)
}


