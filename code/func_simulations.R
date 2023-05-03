# ancient_ILS/code/func_simulations.R
## This script includes functions to run the simulations for this project
# Caitlin Cherryh, 2023

library(ape)
library(phytools)

#### Functions to run the simulation pipeline ####
generate.one.alignment <- function(sim_row, renamed_taxa, partition_path, gene_models, rerun = FALSE){
  # Function to take one row from the simulation parameters dataframe and run start to finish
  
  ## Create the folder for this replicate
  if (dir.exists(sim_row$output_folder) == FALSE){
    dir.create(sim_row$output_folder)
  }
  
  ## Create output file paths
  hyp_tree_file = paste0(sim_row$output_folder, sim_row$ID, "_hypothesis_tree.treefile")
  sim_row_rooted_tree_file = paste0(sim_row$output_folder, sim_row$ID, "_branch_lengths_modified.treefile")
  sim_row_partition_file <- paste0(sim_row$output_folder, sim_row$ID, "_partitions.nexus")
  sim_row_output_alignment_file <- paste0(sim_row$output_folder, sim_row$ID, "_alignment")
  check_ms_gene_tree_file <- paste0(sim_row$output_folder, sim_row$ID, "_ms_gene_trees.txt")
  
  # Check if all the above files exist
  all_files_exist <- file.exists(hyp_tree_file) & file.exists(sim_row_rooted_tree_file) & 
    file.exists(sim_row_partition_file) & file.exists(paste0(sim_row_output_alignment_file, ".fa")) & 
    file.exists(check_ms_gene_tree_file)
  
  if (all_files_exist == FALSE | rerun == TRUE){
    ## Modify the branch lengths of the starting tree based on the simulation parameters
    # Copy the hypothesis tree to the folder
    file.copy(from = sim_row$hypothesis_tree_file, to = hyp_tree_file, overwrite = TRUE)
    # Open the hypothesis tree
    tree <- read.tree(hyp_tree_file)
    # Modify branch lengths for this simulation replicate
    rooted_tree <- manipulate.branch.lengths(starting_tree = tree, parameters_row = sim_row)
    # Save the tree
    write.tree(rooted_tree, file = sim_row_rooted_tree_file)
    
    ## Write the tree in ms command line format and generate gene trees in ms
    ms_files <- ms.generate.trees(unique_id = sim_row$ID, base_tree = rooted_tree, ntaxa = sim_row$num_taxa, 
                                  ntrees = sim_row$num_genes, output_directory = sim_row$output_folder, 
                                  ms_path = sim_row$ms, renamed_taxa = renamed_taxa)
    # Extract gene tree file
    sim_row_gene_tree_file <- ms_files[["ms_gene_tree_file"]]
    
    ## Generate DNA data using Alisim in IQ-Tree
    # Extract gene start and end points from partition file
    gene_ranges <- extract.genes.from.partition.file(partition_path, return.dataframe = FALSE)
    # Create partition file
    partition.gene.trees(num_trees = sim_row$num_genes, gene_ranges = gene_ranges, sequence_type = "AA", 
                         models = gene_models, rescaled_tree_lengths = sim_row$tree_length, 
                         output_filepath = sim_row_partition_file)
    # Generate alignments along gene trees
    alisim.topology.unlinked.partition.model(iqtree_path = sim_row$iqtree2, output_alignment_path = sim_row_output_alignment_file,
                                             partition_file_path = sim_row_partition_file, trees_path = sim_row_gene_tree_file, 
                                             output_format = "fasta", sequence_type = "AA")
    
    # Append information to the simulation row
    sim_row$output_base_tree_file <- sim_row_rooted_tree_file
    sim_row$output_gene_tree_file <- sim_row_gene_tree_file
    sim_row$output_partition_file <- sim_row_partition_file
    sim_row$output_alignment_file <- paste0(sim_row_output_alignment_file, ".fa")
  } else {
    # If all output files already exist, append information to the simulation row
    sim_row$output_base_tree_file <- sim_row_rooted_tree_file
    sim_row$output_gene_tree_file <- check_ms_gene_tree_file
    sim_row$output_partition_file <- sim_row_partition_file
    sim_row$output_alignment_file <- paste0(sim_row_output_alignment_file, ".fa")
  }
  
  # Return the updated simulation row
  return(sim_row)
}



#### Functions to modify trees and manipulate branch lengths
manipulate.branch.lengths <- function(starting_tree, parameters_row){
  # Root the hypothesis tree
  rooted_tree <- root(starting_tree, outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"),
                      resolve.root = TRUE)
  # Make the tree ultrametric
  rooted_tree <- force.ultrametric(rooted_tree, method = "extend")
  # Scale rooted tree to be tree age
  rooted_tree$edge.length <- rooted_tree$edge.length * (parameters_row$tree_length / max(branching.times(rooted_tree)))
  # Find the nodes for the branch depending on which tree is being used
  if (basename(parameters_row$hypothesis_tree_file) == "Whelan2017_hypothesis_tree_1_Cten.treefile"){
    # Extract the nodes by checking for monophyletic clades
    a_end <- getMRCA(rooted_tree, c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster","Daphnia_pulex",
                                    "Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera", 
                                    "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona","Craseo_lathetica",
                                    "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla", "Cliona_varians", "Sycon_coactum", "Sycon_ciliatum",
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
  } else if (basename(parameters_row$hypothesis_tree_file) == "Whelan2017_hypothesis_tree_2_Pori.treefile"){
    a_end <- getMRCA(rooted_tree, c("Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum",
                                    "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica",
                                    "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", 
                                    "Crella_elegans", "Kirkpatrickia_variolosa"))
    a_start <- rooted_tree$edge[which(rooted_tree$edge[,2] == a_end),1]
    b_end <- getMRCA(rooted_tree, c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus",
                                    "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                                    "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                                    "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                                    "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                                    "Ctenophora_sp_Florida_USA"))
    b_start <- rooted_tree$edge[which(rooted_tree$edge[,2] == b_end),1]
  } else if (basename(parameters_row$hypothesis_tree_file) == "Whelan2017_hypothesis_tree_3_CtenPori.treefile"){
    a_end <- getMRCA(rooted_tree, c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus",
                                    "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                                    "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                                    "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                                    "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                                    "Ctenophora_sp_Florida_USA", "Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum",
                                    "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica",
                                    "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", 
                                    "Crella_elegans", "Kirkpatrickia_variolosa"))
    a_start <- rooted_tree$edge[which(rooted_tree$edge[,2] == a_end),1]
    b_end <- getMRCA(rooted_tree, c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus",
                                    "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                                    "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                                    "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                                    "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                                    "Ctenophora_sp_Florida_USA"))
    b_start <- rooted_tree$edge[which(rooted_tree$edge[,2] == b_end),1]
  }
  # Identify branch a
  branch_a <- which(rooted_tree$edge[,1] == a_start & rooted_tree$edge[,2] == a_end)
  # Identify branch b
  branch_b <- which(rooted_tree$edge[,1] == b_start & rooted_tree$edge[,2] == b_end)
  # Modify branch b length
  if (parameters_row$simulation_number == "sim1" | parameters_row$simulation_number == "sim3"){
    # LBA simulation - vary branch b (branch that leads to Ctenophore clade)
    branch_b_value <- rooted_tree$edge.length[branch_b]
    # Multiply existing branch to be that percent of the current tree height
    # e.g. if branch_b_percent_height is 30% and the tree is 1.0 sub/site long, then branch b should be 0.33 sub/site long
    new_branch_b_value <- max(branching.times(rooted_tree)) * (parameters_row$branch_b_percent_height/100)
    rooted_tree$edge.length[branch_b] <- new_branch_b_value
  }
  # Modify branch a length
  if (parameters_row$simulation_number == "sim2" | parameters_row$simulation_number == "sim3"){
    # ILS simulation - vary branch a (branch that allows more time for the two species to differentiate)
    branch_a_value <- rooted_tree$edge.length[branch_a]
    # Multiply existing branch to be that percent of the current tree height
    # e.g. if branch_a_percent_height is 50% and the tree is 1.0 sub/site long, then branch a should be 0.5 sub/site long
    new_branch_a_value <- max(branching.times(rooted_tree)) * (parameters_row$branch_a_percent_height/100)
    rooted_tree$edge.length[branch_a] <- new_branch_a_value
  }
  # Scale rooted tree to be tree age (reset impacts from modifying branch lengths)
  rooted_tree$edge.length <- rooted_tree$edge.length * (parameters_row$tree_length / max(branching.times(rooted_tree)))
  # Make the tree ultrametric
  rooted_tree <- force.ultrametric(rooted_tree, method = "extend")
  
  # Return the manipulated and modified tree
  return(rooted_tree)
}



#### Functions for ms ####
ms.generate.trees <- function(unique_id, base_tree, ntaxa, ntrees, output_directory, ms_path = "ms", renamed_taxa){
  ## Randomly generate a tree with n taxa; format into an ms command and run ms; generate and save the resulting gene trees
  
  ## Generate file paths using the unique id
  t_path <- paste0(output_directory, unique_id, "_starting_tree.txt")
  ms_op_path <- paste0(output_directory, unique_id, "_ms_output.txt")
  ms_gene_trees_path <- paste0(output_directory, unique_id, "_ms_gene_trees.txt")
  
  ## Rename taxa to short versions
  base_tree$tip.label <- unlist(lapply(base_tree$tip.label, function(x){renamed_taxa[[x]]}))
  
  ## Save the input tree
  write.tree(base_tree, file = t_path)
  
  ## Convert the tree into the format for ms
  # Calculate times for ms -ej commands by finding coalescence times (coalescent intervals found using ape::coalescent.intervals)
  ms_coal_ints <- calculate.ms.coalescent.times(base_tree$Nnode, coalescent.intervals(base_tree))
  # Determine the nodes that lead to non-terminal branches {e.g. which(node.depth(t) != 1) }
  nodes <- (ntaxa+1):(ntaxa+base_tree$Nnode)
  # Extract information about all clades from tree
  node_df <- do.call(rbind.data.frame, lapply(nodes, find.branching.times, tree = base_tree))
  names(node_df) <- c("node", "tip_names", "tip_numbers", "ms_tip_order", "ntips", "ndepth", "max_branching_time", "removed_taxa", "ms_input")
  # Order by descending branching time
  node_df <- node_df[order(node_df$max_branching_time, decreasing = TRUE),]
  # Add coalescent times in descending order and reorder columns
  node_df$coalescence_time <- ms_coal_ints
  node_df <- node_df[,c("node", "tip_names", "tip_numbers", "ms_tip_order", "ntips", "ndepth", "max_branching_time", "coalescence_time", "removed_taxa", "ms_input")]
  # Format coalescences for ms input
  node_df <- determine.coalescence.taxa(node_df)
  # Determine which taxa have not yet coalesced
  root_taxa <- select.noncoalesced.taxa(node_df)
  if (length(root_taxa) == 1){
    # Update the coalescence time of the root taxa
    node_df$ms_input_1 <- as.numeric(unlist(lapply(strsplit(node_df$ms_input, " "), function(x){x[1]})))
    node_df$ms_input_2 <- as.numeric(unlist(lapply(strsplit(node_df$ms_input, " "), function(x){x[2]})))
    root_df <- node_df[which(node_df$ms_input_2 == root_taxa),]
    root_row <- root_df[which(root_df$coalescence_time == max(root_df$coalescence_time)),]
    # Check whether the coalescence time is smaller or equal to any other coalescence time
    time_check <- root_row$coalescence_time > max(node_df$coalescence_time)
    if (time_check == FALSE){
      # Slightly increase the coalescent time/max branching time for the root
      new_coal_time <- round(as.numeric(root_row$coalescence_time) + 0.02, digits = 6)
      new_branch_time <- round(as.numeric(root_row$max_branching_time) + 0.02, digits = 6)
      # Determine which row to update in the dataframe
      row_id <- which(node_df$ms_input_1 == root_row$ms_input_1 & node_df$ms_input_2 == root_row$ms_input_2 & 
                        node_df$coalescence_time == root_row$coalescence_time)
      # Update the row
      node_df$coalescence_time[row_id] <- new_coal_time
      node_df$max_branching_time[row_id] <- new_branch_time
    }
  } else if (length(root_taxa) == 2){
    # Make a new row
    new_row <- rep(NA, 10)
    names(new_row) <- c("node", "tip_names", "tip_numbers", "ms_tip_order", "ntips", "ndepth", "max_branching_time", "coalescence_time", "removed_taxa", "ms_input")
    # Create a new coalescence event
    # Specify the two species
    new_row["ms_input"] <- paste0(max(as.numeric(root_taxa)), " ", min(as.numeric(root_taxa)))
    new_row["max_branching_time"] <- round(as.numeric(max(node_df$max_branching_time)) + 0.02, digits = 6)
    new_row["coalescence_time"] <- round(as.numeric(max(node_df$coalescence_time)) + 0.02, digits = 6)
    node_df <- rbind(node_df, new_row)
  } else {
    break
  }
  
  ## Generate gene trees in ms
  # Sort all rows by branching times
  node_df <- node_df[order(node_df$max_branching_time, decreasing = TRUE),]
  # Create a new column containing -ej event for each row
  node_df$ej <- paste0("-ej ", node_df$coalescence_time, " ", node_df$ms_input)
  # No recombination event is present. Do not add any extra splitting (-es) or joining (-ej) events
  # Paste together all the -ej coalescence events for this tree
  all_ej <- paste(node_df$ej, collapse = " ")
  # Construct the ms command line using the -ej events
  coal_call <- paste0(ms_path, " ", ntaxa, " ", ntrees, " -T -I ", ntaxa," ", paste(rep(1, ntaxa), collapse = " "), " ", all_ej)
  # Call ms
  ms_op <- system(coal_call, intern = TRUE)
  # Write all output to file
  write(ms_op, file = ms_op_path)
  
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



extract.clade.from.node <- function(node, tree){
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
  # Determine the maximum branching time associated with this node
  max_branching_time <- max(branching.times(clade))
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
         ms_tip_order = paste(tip_order, collapse = ","), ntips = ntips, ndepth = ndepth, max_branching_time = max_branching_time,
         removed_taxa = removed_taxa, ms_input = ms_input)
  # Return vector
  return(o)
}


find.branching.times <- function(node, tree){
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
  # Determine the maximum branching time associated with this node
  max_branching_time <- max(branching.times(clade))
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
         ms_tip_order = paste(tip_order, collapse = ","), ntips = ntips, ndepth = ndepth, max_branching_time = max_branching_time,
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


select.noncoalesced.taxa <- function(df){
  # Function to add root by coalescing the two remaining lineages together
  
  # Split the ms input column into two separate columns (one for the first number and one for the second number)
  df$ms_input_1 <- as.numeric(unlist(lapply(strsplit(df$ms_input, " "), function(x){x[1]})))
  df$ms_input_2 <- as.numeric(unlist(lapply(strsplit(df$ms_input, " "), function(x){x[2]})))
  # Identify all taxa
  all_taxa <- unique(sort(c(df$ms_input_1, df$ms_input_2)))
  # Identify which taxa are in both (means they have coalesced)
  coalesced_taxa <- sort(intersect(df$ms_input_1, df$ms_input_2))
  remaining_taxa <- setdiff(all_taxa, coalesced_taxa)
  # Check whether the remaining taxa were coalesced at the first node, meaning they are not present in both columns
  coal_check_list <- lapply(remaining_taxa, check.coalesced, coalesced_taxa, df)
  coal_check_df <- as.data.frame(do.call(rbind, coal_check_list))
  # Reduce to only FALSE rows
  noncoal_df <- coal_check_df[coal_check_df$Coalesced == "FALSE",]
  # Return the non-coalesced taxa
  noncoal_taxa <- noncoal_df$Taxon
  return(noncoal_taxa)
}


check.coalesced <- function(test_taxon, coalesced_taxa, df){
  # Quick function to check whether a single taxa coalesced at the first node
  
  # Get the row of the dataframe where this taxon coalesced
  check_rows <- which(df$ms_input_1 == test_taxon)
  # Check how many rows exist
  num_rows <- length(check_rows)
  if (num_rows == 0){
    # There are no rows with this taxon - it does not coalesce into any other taxon
    coalesce = FALSE
    coalesced_into = NA
  } else {
    # There are rows with this taxon - it does coalesce. Determine the coalescence.
    # Extract the relevant rows
    df_row <- df[check_rows,]
    # Check whether the taxon coalesced
    # Conditions:
    #       - The test_taxon is listed in the "removed_taxa" column 
    #         AND the node depth is 2 (meaning this is the first coalescent event for both taxa)
    #       OR
    #       - The taxon to be removed (first in the ms statement) is identical to the test_taxon, 
    #         AND the test_taxon coalesces into a taxon that is in the coalesced_taxa (i.e. a taxon that is coalesced and coalesces)
    #       OR
    #       - The taxon to be removed (first in the ms statement) is identical to the test_taxon,
    #         AND the taxa that the test_taxon coalesces into is present in the tip_numbers for this row 
    if ( (grepl(test_taxon, df_row$removed_taxa) == TRUE & as.numeric(df_row$ndepth) == 2) | 
         (df_row$ms_input_1 == test_taxon & df_row$ms_input_2 %in% coalesced_taxa) | 
         (df_row$ms_input_1 == test_taxon & grepl(df_row$ms_input_2, df_row$tip_numbers) == TRUE) ){
      coalesce = TRUE
      coalesced_into = df_row$ms_input_2
    } else {
      coalesce = FALSE
      coalesced_into = NA
    }
  }
  # Construct output
  op <- c(as.character(test_taxon), coalesce, coalesced_into)
  names(op) <- c("Taxon", "Coalesced", "Coalesced_into")
  return(op)
}



#### Functions for partition files ####
extract.genes.from.partition.file <- function(partition_path, return.dataframe = FALSE){
  # Small function to take a partition file and return a list of gene lengths
  
  # Open the partition file
  p_lines <- readLines(partition_path)
  charset_lines <- grep("charset", p_lines, value = TRUE)
  # Get only the gene parts of the charset lines
  gene_subsets <- unlist(lapply(strsplit(charset_lines,"="), function(x){x[2]}))
  # Remove empty spaces and semi-colons
  gene_subsets <- gsub(" ", "", gene_subsets)
  gene_subsets <- gsub(";", "", gene_subsets)
  # Split into individual gene chunks
  gene_ranges <- unlist(strsplit(gene_subsets,","))
  # Get start and end of each gene
  gene_start <- as.numeric(unlist(lapply(strsplit(gene_ranges, "-"), function(x){x[1]})))
  gene_end <- as.numeric(unlist(lapply(strsplit(gene_ranges, "-"), function(x){x[2]})))
  # Assemble gene start and end points into a df
  gene_df <- data.frame(gene_range = gene_ranges, gene_start = gene_start, gene_end = gene_end)
  # Calculate the gene length for each gene
  gene_df$gene_length <- gene_end - (gene_start - 1) # subtract one from gene_start to count the starting site in the gene length
  # Sort the loci in order
  gene_df <- gene_df[order(gene_df$gene_start, decreasing = FALSE),]
  rownames(gene_df) <- 1:nrow(gene_df)
  # Prepare output
  if (return.dataframe == TRUE){
    # Return the dataframe
    output <- gene_df
  } else if (return.dataframe == FALSE){
    # Return just the gene lengths
    output <- gene_df$gene_range
  }
  # Return output
  return(output)
}

partition.gene.trees <- function(num_trees, gene_ranges, sequence_type, models = NA, rescaled_tree_lengths = NA, output_filepath){
  # This function generates a charpartition file for a set of genes extracted from a partition file
  
  # Create gene names
  gene_names <- paste0("gene_", 1:length(gene_ranges))
  # Generate charsets
  csets <- paste0("\t charset ", gene_names, " = AA, ", gene_ranges, ";")
  
  # Generate model arguments
  if (!is.na(models) == TRUE){
    # If models have been provided, generate a charpartition line for the partition file
    # "models" should be a list of models with one model corresponding to one gene, or 1 model to be applied to all genes
    if (is.na(rescaled_tree_lengths) == TRUE){
      # If no rescaled_tree_lengths are provided, create the chpartitions by directly assembling the models and gene names
      gene_models <- paste0(models, ":", gene_names)
    } else if (is.na(rescaled_tree_lengths) == FALSE){
      # If rescaled_gene_lengths are provided, add them to the gene_models vector. 
      # This allows different genes to have different evolutionary rates (by specifying gene-specific tree lengths within the charpartition command)
      gene_models <- paste0(models, ":", gene_names, "{", rescaled_tree_lengths, "}")
    }
    # Assemble all sections of the charpartition command
    cpart <- paste0("\tcharpartition genes = ", paste(gene_models, collapse = ", "), ";")
  } else{
    # If no models provided, do not create a charpartition line for the partition file
    cpart <- NULL
  }
  
  # Assemble components into partition file
  c_file <- c("#nexus", "begin sets;", csets, cpart, "end;")
  
  # Write partition file
  write(c_file, file = output_filepath)
  
  # Print completion statement
  print("Partition file complete")
}



#### Functions for IQ-Tree and Alisim ####
alisim.topology.unlinked.partition.model <- function(iqtree_path, output_alignment_path, partition_file_path, trees_path, 
                                                     output_format = "fasta", sequence_type = "DNA"){
  # This function uses the topology-unlinked partition model in Alisim to generate a sequence alignment
  #     containing multiple concatenated genes, each with its own tree topology and branch lengths
  
  # Assemble function call 
  function_call <- paste0(iqtree_path, " --alisim ", output_alignment_path, " -Q ", partition_file_path, 
                          " -t ", trees_path, " --seqtype ", sequence_type, " -af ", output_format)
  # Invoke the OS command and call IQ-Tree
  system(function_call)
  
  # Print completion statement
  print("Alisim (IQ-Tree2) run complete")
}


