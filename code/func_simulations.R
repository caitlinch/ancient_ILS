# ancient_ILS/code/func_simulations.R
## This script includes functions to run the simulations for this project
# Caitlin Cherryh, 2023

library(ape)
library(phytools)



#### Functions to run the simulation pipeline ####
generate.one.alignment.wrapper <- function(row_id, sim_df, converted_taxa_names, rerun = FALSE){
  ## Wrapper to easily run generate.one.alignment via dataframe rows
  
  # Apply generate.one.alignment to a single row
  id_row <- sim_df[row_id,]
  output <- generate.one.alignment(sim_row = id_row, converted_taxa_names = converted_taxa_names, rerun = rerun)
  # Return output
  return(output)
}



generate.one.alignment <- function(sim_row, converted_taxa_names, rerun = FALSE){
  ## Function to take one row from the simulation parameters dataframe and run start to finish
  
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
                                  ms_path = sim_row$ms, rename.taxa = TRUE, converted_taxa_names = converted_taxa_names,
                                  preferred_time_difference = sim_row$minimum_coalescent_time_difference)
    # Extract gene tree file
    sim_row_gene_tree_file <- ms_files[["ms_gene_tree_file"]]
    # Open gene trees
    gene_trees <- read.tree(sim_row_gene_tree_file)
    for (i in 1:length(gene_trees)){
      # Select one gene tree
      gt <- gene_trees[[i]]
      # Scale gene tree to match ML tree depth (convert gene trees into subs/site)
      gt$edge.length <- gt$edge.length * (sim_row$ML_tree_depth / max(branching.times(gt)))
      # Replace gene tree
      gene_trees[[i]] <- gt
    }
    
    ## Generate DNA data using Alisim in IQ-Tree
    # Create ranges for gene partitions of the form X-Y (where x and y are numbers)
    num_sites <- sim_row$gene_length * sim_row$num_genes
    gene_end <- seq(sim_row$gene_length, sim_row$gene_length * sim_row$num_genes, sim_row$gene_length)
    gene_start <- c(1, gene_end[1:(length(gene_end) - 1)]+1)
    gene_ranges <- paste0(gene_start, "-", gene_end)
    # Create partition file
    partition.gene.trees(num_trees = sim_row$num_genes, gene_ranges = gene_ranges, sequence_type = "AA", 
                         models = sim_row$alisim_gene_models, rescaled_tree_lengths = NA, 
                         output_filepath = sim_row_partition_file)
    # Generate alignments along gene trees
    alisim.topology.unlinked.partition.model(iqtree_path = sim_row$iqtree2, output_alignment_path = sim_row_output_alignment_file,
                                             partition_file_path = sim_row_partition_file, trees_path = sim_row_gene_tree_file, 
                                             output_format = "fasta", sequence_type = "AA")
  }
  
  # Calculate the ratio of internal branches
  rooted_tree <- read.tree(sim_row_rooted_tree_file)
  tree_length_params <- tree.length.ratio(rooted_tree)
  
  # Construct output
  sim_row_output_vector <- c(as.character(sim_row), as.character(round(tree_length_params, digits = 5)), 
                             sim_row_rooted_tree_file, check_ms_gene_tree_file, 
                             sim_row_partition_file, paste0(sim_row_output_alignment_file, ".fa"))
  names(sim_row_output_vector) <- c(names(sim_row), names(tree_length_params), 
                                    "output_base_tree_file",  "output_gene_tree_file", 
                                    "output_partition_file", "output_alignment_file")

  # Return the updated simulation row
  return(sim_row_output_vector)
}






#### Functions to modify trees and manipulate branch lengths
manipulate.branch.lengths <- function(starting_tree, parameters_row, change.internal.branch.length.percentage = FALSE){
  ## Function to manipulate branch lengths and scale tree
  
  ## Root the tree at the outgroup and scale the tree length 
  # Root the hypothesis tree
  rooted_tree <- root(starting_tree, outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"),
                      resolve.root = TRUE)
  # Make all terminal branch lengths equal to 0.01
  rooted_tree$edge.length[which(is.nan(rooted_tree$edge.length))] <- 0.01
  ## Scale rooted tree to be tree age
  rooted_tree$edge.length <- rooted_tree$edge.length * (as.numeric(parameters_row$ASTRAL_tree_depth) / max(branching.times(rooted_tree)))
  # Make a copy of the tree for modifying the branch lengths
  bl_tree <- rooted_tree
  
  ## Identify the two branches to vary
  # Find the nodes for the branch depending on which tree is being used
  if (grepl("_Cten.tre", basename(parameters_row$hypothesis_tree_file)) == TRUE){
    # Extract the nodes by checking for monophyletic clades
    a_end <- getMRCA(bl_tree, c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster","Daphnia_pulex",
                                "Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera", 
                                "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona","Craseo_lathetica",
                                "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla", "Cliona_varians", "Sycon_coactum", "Sycon_ciliatum",
                                "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum", "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux",
                                "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica", "Petrosia_ficiformis", "Spongilla_lacustris", 
                                "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", "Crella_elegans", "Kirkpatrickia_variolosa"))
    a_start <- bl_tree$edge[which(bl_tree$edge[,2] == a_end),1]
    b_end <- getMRCA(bl_tree, c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus",
                                "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                                "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                                "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                                "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                                "Ctenophora_sp_Florida_USA"))
    b_start <- bl_tree$edge[which(bl_tree$edge[,2] == b_end),1]
  } else if (grepl("_Pori.tre", basename(parameters_row$hypothesis_tree_file)) == TRUE){
    a_end <- getMRCA(bl_tree, c("Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum",
                                "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica",
                                "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", 
                                "Crella_elegans", "Kirkpatrickia_variolosa"))
    a_start <- bl_tree$edge[which(bl_tree$edge[,2] == a_end),1]
    b_end <- getMRCA(bl_tree, c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus",
                                "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                                "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                                "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                                "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                                "Ctenophora_sp_Florida_USA"))
    b_start <- bl_tree$edge[which(bl_tree$edge[,2] == b_end),1]
  } else if (grepl("_CtenPori.tre", basename(parameters_row$hypothesis_tree_file)) == TRUE){
    a_end <- getMRCA(bl_tree, c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus",
                                "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                                "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                                "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                                "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                                "Ctenophora_sp_Florida_USA", "Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum",
                                "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica",
                                "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", 
                                "Crella_elegans", "Kirkpatrickia_variolosa"))
    a_start <- bl_tree$edge[which(bl_tree$edge[,2] == a_end),1]
    b_end <- getMRCA(bl_tree, c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus",
                                "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                                "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                                "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                                "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                                "Ctenophora_sp_Florida_USA"))
    b_start <- bl_tree$edge[which(bl_tree$edge[,2] == b_end),1]
  }
  
  ## Modify branch a length
  # Identify branch a
  branch_a <- which(bl_tree$edge[,1] == a_start & bl_tree$edge[,2] == a_end)
  # ILS simulation - vary branch a (branch that allows more time for the two species to differentiate)
  bl_tree$edge.length[branch_a] <- as.numeric(parameters_row$branch_a_length)
  
  ## Modify branch b length
  # Identify branch b
  branch_b <- which(bl_tree$edge[,1] == b_start & bl_tree$edge[,2] == b_end)
  # LBA simulation - vary branch b (branch that leads to Ctenophore clade)
  bl_tree$edge.length[branch_b] <- as.numeric(parameters_row$branch_b_length)
  
  ## Modify the other branch lengths
  # Branch c - branch leading to Bilateria/Cnidaria
  c_end <- getMRCA(bl_tree, c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster","Daphnia_pulex",
                              "Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera", 
                              "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona","Craseo_lathetica",
                              "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla"))
  c_start <- bl_tree$edge[which(bl_tree$edge[,2] == c_end),1]
  branch_c <- which(bl_tree$edge[,1] == c_start & bl_tree$edge[,2] == c_end)
  bl_tree$edge.length[branch_c] <- as.numeric(parameters_row$branch_c_length)
  # Branch leading to Cnidaria
  cnidaria_end <- getMRCA(bl_tree, c("Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera", 
                                     "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona","Craseo_lathetica",
                                     "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla"))
  cnidaria_start <- bl_tree$edge[which(bl_tree$edge[,2] == cnidaria_end),1]
  branch_cnidaria <- which(bl_tree$edge[,1] == cnidaria_start & bl_tree$edge[,2] == cnidaria_end)
  bl_tree$edge.length[branch_cnidaria] <- as.numeric(parameters_row$branch_cnidaria_length)
  # Branch leading to Bilateria
  bilateria_end <- getMRCA(bl_tree, c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster","Daphnia_pulex"))
  bilateria_start <- bl_tree$edge[which(bl_tree$edge[,2] == bilateria_end),1]
  branch_bilateria <- which(bl_tree$edge[,1] == bilateria_start & bl_tree$edge[,2] == bilateria_end)
  bl_tree$edge.length[branch_bilateria] <- as.numeric(parameters_row$branch_bilateria_length)
  # Branch leading to Porifera
  porifera_end <- getMRCA(bl_tree, c("Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum",
                                     "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica",
                                     "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", 
                                     "Crella_elegans", "Kirkpatrickia_variolosa"))
  porifera_start <- bl_tree$edge[which(bl_tree$edge[,2] == porifera_end),1]
  branch_porifera <- which(bl_tree$edge[,1] == porifera_start & bl_tree$edge[,2] == porifera_end)
  bl_tree$edge.length[branch_porifera] <- as.numeric(parameters_row$branch_porifera_length)
  # Branch leading to all animals
  all_end <- getMRCA(bl_tree, c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster","Daphnia_pulex",
                                "Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera", 
                                "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona","Craseo_lathetica",
                                "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla",
                                "Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum",
                                "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica",
                                "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", 
                                "Crella_elegans", "Kirkpatrickia_variolosa",
                                "Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus",
                                "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                                "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                                "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                                "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                                "Ctenophora_sp_Florida_USA"))
  all_start <- bl_tree$edge[which(bl_tree$edge[,2] == all_end),1]
  branch_all <- which(bl_tree$edge[,1] == all_start & bl_tree$edge[,2] == all_end)
  bl_tree$edge.length[branch_all] <- as.numeric(parameters_row$branch_all_animals_length)
  # Branch leading to outgroup
  outgroup_end <- getMRCA(bl_tree, c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"))
  outgroup_start <- bl_tree$edge[which(bl_tree$edge[,2] == outgroup_end),1]
  branch_outgroup <- which(bl_tree$edge[,1] == outgroup_start & bl_tree$edge[,2] == outgroup_end)
  bl_tree$edge.length[branch_outgroup] <- as.numeric(parameters_row$branch_outgroup_length)
  
  ## Force the tree to be ultrametric
  # Copy the tree across
  um_tree <- bl_tree
  # Make the tree ultrametric
  um_tree <- force.ultrametric(um_tree, method = "extend")
  
  ## Manipulate the tree depth and length of terminal branches 
  if (change.internal.branch.length.percentage == TRUE){
    # Identify all terminal branches
    terminal_branch_ids <- which(um_tree$edge[,2] <= Ntip(um_tree))
    internal_branch_ids <- which(um_tree$edge[,2] > Ntip(um_tree))
    # Determine how long the sum of all terminal branches would have to be to match the parameters_row$proportion_internal_branches
    current_tree_length <- sum(um_tree$edge.length)
    current_external_length <- sum(um_tree$edge.length[terminal_branch_ids])
    current_internal_length <- sum(um_tree$edge.length[internal_branch_ids])
    # Get minimum and maximum terminal and internal branch lengths
    min_terminal_bl <- min(um_tree$edge.length[terminal_branch_ids])
    max_terminal_bl <- max(um_tree$edge.length[terminal_branch_ids])
    min_internal_bl <- min(um_tree$edge.length[internal_branch_ids])
    max_internal_bl <- max(um_tree$edge.length[internal_branch_ids])
    # Check whether to extend the internal or external branches
    final_percent_internal_branches <- parameters_row$proportion_internal_branches * 100
    um_tree_percent_internal_branches <- current_internal_length/current_tree_length * 100
    if (um_tree_percent_internal_branches < final_percent_internal_branches){
      # The internal branch length is shorter than it should be.
      # Reduce the terminal branches
      final_tree_length <- current_external_length/(100-final_percent_internal_branches)*100
      final_internal_branch_length <- final_tree_length * (final_percent_internal_branches/100)
      # Determine how much each terminal branch must be extended to reach that final tree length
      extend_internal_branches_by <- final_internal_branch_length/Ntip(um_tree)
      # Reduce terminal branches
      um_tree$edge.length[terminal_branch_ids] <- um_tree$edge.length[terminal_branch_ids] + extend_terminal_branches_by
    } else if (um_tree_percent_internal_branches > final_percent_internal_branches){
      # The internal branch length is longer than it should be.
      # Extend the terminal branches
      final_tree_length <- current_internal_length/final_percent_internal_branches*100
      final_terminal_branch_length <- final_tree_length * ((100-final_percent_internal_branches)/100)
      # Determine how much each terminal branch must be extended to reach that final tree length
      extend_terminal_branches_by <- final_terminal_branch_length/Ntip(um_tree)
      # Extend terminal branches
      um_tree$edge.length[terminal_branch_ids] <- um_tree$edge.length[terminal_branch_ids] + extend_terminal_branches_by
    }
  }
  
  ## Return the manipulated and modified tree
  return(um_tree)
}






#### Functions for ms ####
ms.generate.trees <- function(unique_id, base_tree, ntaxa, ntrees, output_directory, ms_path = "ms", 
                              rename.taxa = FALSE, converted_taxa_names, preferred_time_difference = 0.001){
  ## Take a given tree; format into an ms command and run ms; generate and save the resulting gene trees
  # Works for random coalescent trees and trees generated from Metazoan datasets (although needs occasional 
  #     manual tweaking within the function when simultaneous coalescent events involving the same taxa occur)
  
  ## Generate file paths using the unique id
  t_path <- paste0(output_directory, unique_id, "_starting_tree.txt")
  ms_op_path <- paste0(output_directory, unique_id, "_ms_output.txt")
  ms_gene_trees_path <- paste0(output_directory, unique_id, "_ms_gene_trees.txt")
  
  ## Rename taxa to short versions
  if (rename.taxa == TRUE){
    base_tree$tip.label <- unlist(lapply(base_tree$tip.label, function(x){converted_taxa_names[[x]]}))
  }
  
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
  node_df$max_branching_time <- as.numeric(node_df$max_branching_time)
  node_df <- node_df[order(node_df$max_branching_time, decreasing = TRUE),]
  # Add coalescent times in descending order and reorder columns
  node_df$coalescence_time <- ms_coal_ints
  node_df <- node_df[,c("node", "tip_names", "tip_numbers", "ms_tip_order", "ntips", "ndepth", "max_branching_time", "coalescence_time", "removed_taxa", "ms_input")]
  # Order by descending coalescence time
  node_df <- node_df[order(node_df$coalescence_time, decreasing = TRUE),]
  # Format coalescences for ms input
  node_df <- determine.coalescence.taxa(node_df)
  # Check for duplicated times
  node_df <- check.duplicated.coalescent.times(node_df, preferred_time_difference)
  # Determine which taxa have not yet coalesced
  root_taxa <- select.noncoalesced.taxa(node_df)
  node_df <- check.root.coalesced(root_taxa, node_df)
  
  ## Generate gene trees in ms
  # Sort all rows by coalescence time (put in correct order for ms call)
  node_df <- node_df[order(node_df$coalescence_time, decreasing = TRUE),]
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



check.duplicated.coalescent.times <- function(node_df, preferred_time_difference = 0.001){
  ## Small function to check for duplicate coalescent events and space them out so no events occur at the same time
  
  # Identify all duplicate times
  if (length(unique(node_df$coalescence_time)) != length(node_df$coalescence_time)){
    ## Work out which rows have multiple values
    # Determine the matching rows
    duplicate_rows <- which(duplicated(node_df$coalescence_time))
    # For each duplicated row, get all rows with the same time
    rows_equal_to_duplicate_rows <- which(node_df$coalescence_time %in% node_df$coalescence_time[duplicate_rows])
    unduplicate_rows <- setdiff(rows_equal_to_duplicate_rows, duplicate_rows)
    
    ## Find smallest difference between coalescent times
    all_diff_times <- unlist(lapply(1:length(node_df$coalescence_time), function(i){node_df$coalescence_time[i]-node_df$coalescence_time[i+1]}))
    diff_times <- all_diff_times[which(all_diff_times != 0)]
    min_diff_time <- min(diff_times, na.rm = TRUE)
    
    ## For each of the unduplicated rows, need to update coalescent times so that there is no duplicate
    for (u in unduplicate_rows){
      # Get which rows are identical to this row
      u_dupes <- which(node_df$coalescence_time == node_df$coalescence_time[u])
      # Check whether those rows have the same taxa as the u row
      dupes_ms_input <- node_df$ms_input[u_dupes]
      dupe_time_df <- data.frame(row = u_dupes,
                                 coalescence_time = node_df$coalescence_time[u],
                                 taxa_1 = unlist(lapply(1:length(dupes_ms_input), function(x){strsplit(dupes_ms_input, " ")[[x]][1]})),
                                 taxa_2 = unlist(lapply(1:length(dupes_ms_input), function(x){strsplit(dupes_ms_input, " ")[[x]][2]})))
      # Add min taxa and max taxa
      dupe_time_df$min_taxa <- lapply(1:length(dupes_ms_input), function(x){min(as.numeric(strsplit(dupes_ms_input, " ")[[x]]))})
      dupe_time_df$max_taxa <- lapply(1:length(dupes_ms_input), function(x){max(as.numeric(strsplit(dupes_ms_input, " ")[[x]]))})
      # Find the next time (the constraint for updating the identical coalescent time)
      next_coalescence_time <- node_df$coalescence_time[max(u_dupes)+1]
      # Find which rows have the same taxa involved AND the same coalescence time
      check_dupe_rows_1 <- u_dupes[duplicated(dupe_time_df$taxa_2) | duplicated(dupe_time_df$taxa_2)]
      check_dupe_rows_2 <- dupe_time_df$row[which(dupe_time_df$min_taxa %in% dupe_time_df$max_taxa)]
      # If there are any rows with the same taxa and the same coalescent units, process them further: 
      if ( (identical(check_dupe_rows_1, integer(0)) == FALSE) | (identical(check_dupe_rows_2, integer(0)) == FALSE) ){
        # Construct the rows to update
        if (identical(check_dupe_rows_1, integer(0)) == FALSE & identical(check_dupe_rows_2, integer(0)) == TRUE){
          row_to_update_id <- check_dupe_rows_1
          dupe_time_update_row_id <- which((duplicated(dupe_time_df$taxa_2) | duplicated(dupe_time_df$taxa_2)))
        } else if (identical(check_dupe_rows_1, integer(0)) == TRUE & identical(check_dupe_rows_2, integer(0)) == FALSE){
          row_to_update_id <- check_dupe_rows_2
          dupe_time_update_row_id <- which(dupe_time_df$row == row_to_update_id)
        } else if (identical(check_dupe_rows_1, integer(0)) == FALSE & identical(check_dupe_rows_2, integer(0)) == FALSE){
          row_to_update_id <- unique(c(check_dupe_rows_1, check_dupe_rows_2))
          dupe_time_update_row_id <- unique(c(which((duplicated(dupe_time_df$taxa_2) | duplicated(dupe_time_df$taxa_2))),
                                              which(dupe_time_df$row == row_to_update_id) ) )
        }
        # Get times to update
        dupe_coalescent_time <- dupe_time_df$coalescence_time[dupe_time_update_row_id]
        # Determine how many coalescent times are needed
        num_coal_times_needed <- length(row_to_update_id)
        # Check whether the preferred time difference is too big
        test_times <- unlist(lapply(num_coal_times_needed, function(x){dupe_coalescent_time - (x * preferred_time_difference)}))
        check_test_times <- identical(which(test_times < next_coalescence_time), integer(0))
        # If the test times are larger than the next coalescent time, use the test times as the new coalescent times
        # Otherwise if the coalescent times are very close together, make the new coalescent times equally spoaced between coalescent times
        if (check_test_times == TRUE){
          # Make the new coalescent times spaced by the preferred time difference
          new_coalescent_times <- test_times
        } else if (check_test_times == FALSE){
          # Make the new coalescent times equally spaced between the existing coalescent times
          # Create the required number of new times
          new_coalescent_times <- seq(from = next_coalescence_time, to = dupe_coalescent_time, by = ((dupe_coalescent_time-next_coalescence_time)/(num_coal_times_needed+1)))
          # Remove existing coalescent times (remove first and last times, which are the to and from times fed to seq.along)
          new_coalescent_times <- new_coalescent_times[2:(length(new_coalescent_times)-1)]
        }
        # Update any rows with the same taxa involved AND the same coalescence time
        for (j in 1:length(row_to_update_id)){
          # Get correct row and new coalescent time
          n_row_id <- row_to_update_id[j]
          n_new_coal_time <- new_coalescent_times[j]
          # Update coalescent time in original node_df
          node_df$coalescence_time[n_row_id] <- n_new_coal_time
        } # end for (j in 1:length(row_to_update_id)){ -- update coalescent times for events with duplicate times + taxa
      } # end if ( (identical(check_dupe_rows_1, integer(0)) == FALSE) | (identical(check_dupe_rows_2, integer(0)) == FALSE) ){ -- check whether any events have duplicate times + taxa 
    } # end for (u in unduplicate_rows){ -- process rows with duplicate times
  } # end if (length(unique(node_df$coalescence_time)) != length(node_df$coalescence_time)){ -- check whether any events have duplicate times
  
  # Return the node_df with updated coalescent times
  return(node_df)
}


get.conflicting.coalescent.event <- function(row_number, node_df){
  ## Identify coalescent events occuring at same time with same taxa
  
  
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
  node_dataframe <- node_dataframe[order(node_dataframe$ndepth, node_dataframe$coalescence_time),]
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
  ## Function to add root by coalescing the two remaining lineages together
  
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
  ## Quick function to check whether a single taxa coalesced at the first node
  
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



check.root.coalesced <- function(root_taxa, node_df){
  ## Function to add a final coalescence event, if there is not one 
  
  if (length(root_taxa) == 1){
    # Update the coalescence time of the root taxa
    node_df$ms_input_1 <- as.numeric(unlist(lapply(strsplit(node_df$ms_input, " "), function(x){x[1]})))
    node_df$ms_input_2 <- as.numeric(unlist(lapply(strsplit(node_df$ms_input, " "), function(x){x[2]})))
    root_df <- node_df[which(node_df$ms_input_2 == root_taxa),]
    root_row <- root_df[which(root_df$coalescence_time == max(root_df$coalescence_time)),]
    # Check whether the coalescence time is smaller or equal to any other coalescence time
    time_check <- root_row$coalescence_time > max(node_df$coalescence_time)
    if ( (TRUE %in% time_check) == FALSE) {
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
  }
  # return the node_df
  return(node_df)
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
  ## This function generates a charpartition file for a set of genes extracted from a partition file
  
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
  ## This function uses the topology-unlinked partition model in Alisim to generate a sequence alignment
  #     containing multiple concatenated genes, each with its own tree topology and branch lengths
  
  # Assemble function call 
  function_call <- paste0(iqtree_path, " --alisim ", output_alignment_path, " -Q ", partition_file_path, 
                          " -t ", trees_path, " --seqtype ", sequence_type, " -af ", output_format)
  # Invoke the OS command and call IQ-Tree
  system(function_call)
  
  # Print completion statement
  print("Alisim (IQ-Tree2) run complete")
}





##### Calculate tree length ####
tree.length.ratio <- function(tree){
  ## Quick function to calculate the tree length and percentage of internal branches for any tree
  terminal_branch_ids <- which(tree$edge[,2] <= Ntip(tree))
  internal_branch_ids <- which(tree$edge[,2] > Ntip(tree))
  # Determine how long the sum of all terminal branches would have to be to match the parameters_row$proportion_internal_branches
  current_tree_length <- sum(tree$edge.length)
  current_external_length <- sum(tree$edge.length[terminal_branch_ids])
  current_internal_length <- sum(tree$edge.length[internal_branch_ids])
  # Calculate percentage
  percentage_internal_length <- current_internal_length / current_tree_length * 100
  # Assemble output
  op_vec <- c(current_tree_length, current_internal_length, percentage_internal_length)
  names(op_vec) <- c("total_tree_length", "sum_internal_branch_lengths", "percentage_internal_branch_lengths")
  # Return output
  return(op_vec)
}



