# ancient_ILS/code/func_analysis.R
## This script includes functions investigate and analyses the results of the simulations 
# Caitlin Cherryh, 2023

library(ape)
library(phangorn)

#### Analysis wrapper function ####
analysis.wrapper <- function(row_id, df, hypothesis_tree_dir, renamed_taxa){
  # Wrapper function to calculate:
  #   - [x] actual and estimated gcfs (IQ-Tree)
  #   - [ ] actual and estimated qcfs (ASTRAL)
  #   - [x] hypothesis tree distances for both ML and ASTRAL trees
  #   - [ ] hypothesis tests in IQ-Tree to see if any hypothesis can be ignored for this alignment
  
  # Open the row of interest
  df_row <- df[row_id, ]
  
  # Calculate the differences between the three trees
  iqtree2_tree_diffs <- calculate.distance.between.three.trees(tree_path = df_row$ML_tree_treefile, hypothesis_tree_dir, tree_type = "ML", 
                                                               rename.hypothesis.tree.tips = TRUE, renamed_taxa = renamed_taxa)
  astral_tree_diffs <- calculate.distance.between.three.trees(tree_path = df_row$ASTRAL_tree_treefile, hypothesis_tree_dir, tree_type = "ASTRAL", 
                                                              rename.hypothesis.tree.tips = TRUE, renamed_taxa = renamed_taxa)
  
  # Calculate the gCFs using IQ-Tree
  iqtree2_gcfs <- gcf.wrapper(alignment_path = df_row$output_alignment_file, iqtree2_path = df_row$iqtree2, iqtree2_model = NA,
                              iqtree2_num_threads = df_row$iqtree2_num_threads, rename.taxa.for.ms = TRUE, renamed_taxa = renamed_taxa)
  
  # Calculate the qCFs using ASTRAL
  astral_qcfs <- qcf.wrapper(ID = df_row$ID, starting_tree = df_row$output_base_tree_file, ms_gene_trees = df_row$output_gene_tree_file,
                             ASTRAL_tree = df_row$ASTRAL_tree_treefile, ML_gene_trees = df_row$iqtree2_gene_tree_treefile, 
                             ASTRAL_path = df_row$ASTRAL, call.astral = TRUE)
  
  # Perform hypothesis tests in IQ-Tree
  hypothesis_tests <- ""
  
  # Trim unwanted columns
  trimmed_df_row <- df_row[, c("dataset", "dataset_type", "ID", "output_folder", "simulation_number", "simulation_type",
                               "hypothesis_tree", "replicates", "num_taxa", "num_genes", "gene_length", "num_sites",
                               "ML_tree_estimation_models", "branch_a_length", "branch_b_length", "branch_c_length",
                               "branch_all_animals_length", "branch_bilateria_length", "branch_cnidaria_length",
                               "branch_outgroup_length", "branch_porifera_length", "total_tree_length",
                               "sum_internal_branch_lengths", "percentage_internal_branch_lengths")]
  # Assemble output vector
  analysis_output         <- c(as.character(trimmed_df_row), iqtree2_tree_diffs, astral_tree_diffs,
                               iqtree2_gcfs, astral_qcfs, hypothesis_tests)
  names(analysis_output)  <- c(names(trimmed_df_row), names(iqtree2_tree_diffs), names(astral_tree_diffs),
                               names(iqtree2_gcfs), names(astral_qcfs), names(hypothesis_tests))
}







#### Perform hypothesis tests in IQ-Tree2 ####
perform.hypothesis.tests <- function(ID, alignment_path, hypothesis_tree_file, iqtree2_path){
  # Function to perform hypothesis tests in IQ-Tree2 on the simulated alignment with either 2 or 3 hypothesis trees
  
  ## Change directories to the same directory as the alignment
  
  ## Assemble the IQ-Tree command and call IQ-Tree
  
  ## Extract output from the IQ-Tree file
  
  ## Format output
  
  ## Return output
  
}






#### Calculate distance between trees ####
calculate.distance.between.three.trees <- function(tree_path, hypothesis_tree_dir, tree_type, rename.hypothesis.tree.tips = FALSE, renamed_taxa){
  # Function to calculate the distance between a tree and three hypothesis trees
  # tree_type = "ASTRAL" | tree_type = "IQ-Tree2"
  
  ## Open the hypothesis trees
  # Extract the tree files
  all_hyp_dir_files <- list.files(hypothesis_tree_dir)
  # Extract only either the ASTRAL or IQ-Tree tree files
  if (tree_type == "ASTRAL" | tree_type == "astral"){
    tree_files <- grep("ASTRAL", all_hyp_dir_files,value = T)
  } else if (tree_type == "IQ-Tree2" | tree_type == "IQ-Tree" | tree_type == "ML" | tree_type == "iqtree" | tree_type == "iqtree2"){
    tree_files <- grep("treefile", all_hyp_dir_files,value = T)    
  }
  # Extract the three hypothesis tree files
  h1_path <- paste0(hypothesis_tree_dir, grep("hypothesis_tree_1", tree_files, value = T))
  h2_path <- paste0(hypothesis_tree_dir, grep("hypothesis_tree_2", tree_files, value = T))
  h3_path <- paste0(hypothesis_tree_dir, grep("hypothesis_tree_3", tree_files, value = T))
  # Open the tree files
  h1_tree <- read.tree(file = h1_path)
  h2_tree <- read.tree(file = h2_path)
  h3_tree <- read.tree(file = h3_path)
  
  ## Prepare the hypothesis trees for comparison with the test tree
  # Reroot the trees
  h1_tree <- root(h1_tree, outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"),
                  resolve.root = TRUE)
  h2_tree <- root(h2_tree, outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"),
                  resolve.root = TRUE)
  h3_tree <- root(h3_tree, outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"),
                  resolve.root = TRUE)
  if (rename.hypothesis.tree.tips == TRUE){
    # Rename taxa to short versions (numbers only)
    h1_tree$tip.label <- unlist(lapply(h1_tree$tip.label, function(x){renamed_taxa[[x]]}))
    h1_tree$tip.label <- gsub("t", "", h1_tree$tip.label)
    h2_tree$tip.label <- unlist(lapply(h2_tree$tip.label, function(x){renamed_taxa[[x]]}))
    h2_tree$tip.label <- gsub("t", "", h2_tree$tip.label)
    h3_tree$tip.label <- unlist(lapply(h3_tree$tip.label, function(x){renamed_taxa[[x]]}))
    h3_tree$tip.label <- gsub("t", "", h3_tree$tip.label)
  }
  
  ## Open the tree of interest
  t <- read.tree(file = tree_path)
  
  ## Calculate RF and wRF for the trees
  output_vec <- c(tree_path, RF.dist(t, h1_tree), wRF.dist(t, h1_tree), RF.dist(t, h2_tree), wRF.dist(t, h2_tree), RF.dist(t, h3_tree), wRF.dist(t, h3_tree))
  names(output_vec) <- c(paste0(tree_type, "_tree_path"), paste0(tree_type, "_h1_RF_dist"), paste0(tree_type, "_h1_wRF_dist"), paste0(tree_type, "_h2_RF_dist"),
                         paste0(tree_type, "_h2_wRF_dist"), paste0(tree_type, "_h3_RF_dist"), paste0(tree_type, "_h3_wRF_dist") )
  
  ## Return the RF/wRF distances
  return(output_vec)
}






#### Quartet concordance factors ####
qcf.wrapper <- function(ID, starting_tree, ms_gene_trees, ASTRAL_tree, ML_gene_trees, ASTRAL_path, call.astral = TRUE){
  ## Function to determine expected and estimated qCF in ASTRAL and return summary statistics
  
  ## Identify directory
  qcf_dir <- paste0(dirame(starting_tree), "/")
  
  ## Calculate actual quartet concordance factors
  expected_qcf_paths <- qcf.call(output_id = paste0(ID, "_expected_qcfs"), output_directory = qcf_dir, 
                                 tree_path = starting_tree, gene_trees_path = ms_gene_trees,
                                 ASTRAL_path = ASTRAL_path, call.astral = TRUE)
  
  ## Calculate estimated quartet concordance factors
  estimated_qcf_paths <- qcf.call(output_id = paste0(ID, "_estimated_qcfs"), output_directory = qcf_dir, 
                                  tree_path = ASTRAL_tree, gene_trees_path = ML_gene_trees,
                                  ASTRAL_path = ASTRAL_path, call.astral = TRUE)
  
  ## Extract relevant qCF values
  expected_qcf_values <-  extract.qcf.values(qcf_tree = expected_qcf_paths[["qcf_output_tree"]])
  estimated_qcf_values <-  extract.qcf.values(qcf_tree = estimated_qcf_paths[["qcf_output_tree"]])
  
  ## Assemble output
  qcf_output <- c(expected_qcf_paths, expected_qcf_values,
                  estimated_qcf_paths, estimated_qcf_values)
  names(qcf_output) <- c(paste0("expected_", names(expected_qcf_paths)), paste0("expected_", names(expected_qcf_values)),
                         paste0("estimated_", names(estimated_qcf_paths)), paste0("estimated_", names(estimated_qcf_values)))
  
  ## Return output
  return(qcf_output)
}



qcf.call <- function(output_id, output_directory, tree_path, gene_trees_path, ASTRAL_path, call.astral = TRUE){
  ## Function to call ASTRAL and calculate quartet concordance factors
  
  ## Check whether the tree tips match
  test_tree <- read.tree(tree_path)
  tree_tips <- test_tree$tip.label
  if (TRUE %in% grepl("t", tree_tips)){
    # Replace tree tips without "t"
    new_tree_tips <- gsub("t", "", tree_tips)
    # Replace tree tips in original tree
    test_tree$tip.label <- new_tree_tips
    # Create new file path
    split_path <- unlist(strsplit(basename(tree_path), "\\."))
    relabelled_tree_path <- paste0(dirname(tree_path), paste(head(split_path, (length(split_path) - 1)), collapse = "."), "_relabelled.", tail(split_path, 1))
    # Save the relabelled tree
    write.tree(test_tree, file = relabelled_tree_path)
  } else {
    relabelled_tree_path <- tree_path
  }
  
  # Assemble ASTRAL command to calculate quartet concordance factors
  output_tree <- paste0(output_directory, output_id, ".tre")
  output_log <- paste0(output_directory, output_id, ".log")
  qcf_command <- paste0("java -jar ", ASTRAL_path, " -q ", tree_path, " -i ", gene_trees_path, " -o ", output_tree, " 2> ", output_log)
  # if call.astral == TRUE, run ASTRAL to calculate quartet concordance factors
  if (call.astral == TRUE){
    system(qcf_command)
  }
  # Assemble output vector
  qcf_output <- c(output_id, qcf_command, call.astral, output_tree, output_log)
  names(qcf_output) <- c("qcf_ID", "astral_qCF_command", "astral_qCF_run", "qcf_output_tree", "qcf_output_log")
  # Return the output
  return(qcf_output)
}



extract.qcf.values <- function(qcf_tree){
  ## Function to open a qCF tree, extract the values of interest and return some summary statistics
}





#### Gene concordance factors wrapper ####
gcf.wrapper <- function(alignment_path, iqtree2_path, iqtree2_model = NA, iqtree2_num_threads = "AUTO", rename.taxa.for.ms = TRUE, renamed_taxa){
  # Function to calculate the estimated and empirical gCF and return relevant gCFs
  
  ## Calculate the exact/input/actual gCFs using IQ-Tree2
  expected_gcfs <- extract.input.concordance.factors(alignment_path, iqtree2_path, iqtree2_num_threads, rename.taxa.for.ms = TRUE, renamed_taxa)
  
  ## Calculate the expected gCFs using IQ-Tree2
  estimated_gcfs <- extract.output.concordance.factors(alignment_path, iqtree2_path, iqtree2_num_threads, iqtree2_model)
  
  ## Extract information from the actual gCFS
  # Extract the gcf tables
  a_gcfs <- expected_gcfs$gcf_table 
  # Read in the branch labelled tree
  a_tree <- read.tree(file = expected_gcfs$gcf_branch_file)
  # Find the "a" branch and the "b" branch
  if (grepl("h1", alignment_path) == TRUE){
    # Extract the nodes by checking for monophyletic clades
    a_end <- getMRCA(a_tree, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                               "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
                               "39", "40"))
    a_start <- a_tree$edge[which(a_tree$edge[,2] == a_end),1]
    b_end <- getMRCA(a_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                               "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
    b_start <- a_tree$edge[which(a_tree$edge[,2] == b_end),1]
  } else if (grepl("h2", alignment_path) == TRUE){
    a_end <- getMRCA(a_tree, c("22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40"))
    a_start <- a_tree$edge[which(a_tree$edge[,2] == a_end),1]
    b_end <- getMRCA(a_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                               "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
    b_start <- a_tree$edge[which(a_tree$edge[,2] == b_end),1]
  } else if (grepl("h3", alignment_path) == TRUE){
    a_end <- getMRCA(a_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                               "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "22", "23", "24", "25", "26", "27",
                               "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40"))
    a_start <- a_tree$edge[which(a_tree$edge[,2] == a_end),1]
    b_end <- getMRCA(a_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                               "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
    b_start <- a_tree$edge[which(a_tree$edge[,2] == b_end),1]
  }
  # Identify branch "a" and branch "b"
  branch_a <- which(a_tree$edge[,1] == a_start & a_tree$edge[,2] == a_end)
  branch_b <- which(a_tree$edge[,1] == b_start & a_tree$edge[,2] == b_end)
  # Identify which labels belong to branch "a" and branch "b" for the gCFs
  branch_a_label = a_end
  branch_b_label = b_end
  # Extract information from those two rows
  a_gcfs_output_names <- c(paste0("actual_branch_a_", c("ID", "gCF", "gCF_N", "gDF1", "gDF1_N", "gDF2", "gDF2_N", "gDFP", "gDFP_N", "gN", "Label", "Length")),
                           paste0("actual_branch_b_", c("ID", "gCF", "gCF_N", "gDF1", "gDF1_N", "gDF2", "gDF2_N", "gDFP", "gDFP_N", "gN", "Label", "Length")),
                           "actual_mean_gCF")
  a_gcf_output <- c(as.character(a_gcfs[which(a_gcfs$ID == branch_a_label),]), as.character(a_gcfs[which(a_gcfs$ID == branch_b_label),]),
                    round(mean(a_gcfs$gCF), digits = 2) )
  names(a_gcf_output) <- a_gcfs_output_names
  
  ## Extract information from the estimated gCFS
  # Extract the gcf tables
  e_gcfs <- estimated_gcfs$gcf_table 
  # Read in the branch labelled tree
  e_tree <- read.tree(file = estimated_gcfs$gcf_branch_file)
  # Find the "a" branch and the "b" branch
  if (grepl("h1", alignment_path) == TRUE){
    # Extract the nodes by checking for monophyletic clades
    a_end <- getMRCA(e_tree, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                               "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
                               "39", "40"))
    a_start <- e_tree$edge[which(e_tree$edge[,2] == a_end),1]
    b_end <- getMRCA(e_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                               "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
    b_start <- e_tree$edge[which(e_tree$edge[,2] == b_end),1]
  } else if (grepl("h2", alignment_path) == TRUE){
    a_end <- getMRCA(e_tree, c("22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40"))
    a_start <- e_tree$edge[which(e_tree$edge[,2] == a_end),1]
    b_end <- getMRCA(e_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                               "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
    b_start <- e_tree$edge[which(e_tree$edge[,2] == b_end),1]
  } else if (grepl("h3", alignment_path) == TRUE){
    a_end <- getMRCA(e_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                               "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "22", "23", "24", "25", "26", "27",
                               "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40"))
    a_start <- e_tree$edge[which(e_tree$edge[,2] == a_end),1]
    b_end <- getMRCA(e_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                               "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
    b_start <- e_tree$edge[which(e_tree$edge[,2] == b_end),1]
  }
  # Identify branch "a" and branch "b"
  branch_a <- which(e_tree$edge[,1] == a_start & e_tree$edge[,2] == a_end)
  branch_b <- which(e_tree$edge[,1] == b_start & e_tree$edge[,2] == b_end)
  # Identify which labels belong to branch "a" and branch "b" for the gCFs
  branch_a_label = a_end
  branch_b_label = b_end
  # Extract information from those two rows
  e_gcfs_output_names <- c(paste0("estimated_branch_a_", c("ID", "gCF", "gCF_N", "gDF1", "gDF1_N", "gDF2", "gDF2_N", "gDFP", "gDFP_N", "gN", "Label", "Length")),
                           paste0("estimated_branch_b_", c("ID", "gCF", "gCF_N", "gDF1", "gDF1_N", "gDF2", "gDF2_N", "gDFP", "gDFP_N", "gN", "Label", "Length")),
                           "estimated_mean_gCF")
  e_gcf_output <- c(as.character(e_gcfs[which(e_gcfs$ID == branch_a_label),]), as.character(e_gcfs[which(e_gcfs$ID == branch_b_label),]),
                    round(mean(e_gcfs$gCF), digits = 2) )
  names(e_gcf_output) <- e_gcfs_output_names
  
  ## Collate gCF output
  collated_output <- c(alignment_path, a_gcf_output, e_gcf_output)
  collated_output_names <- c("alignment_path", a_gcfs_output_names, e_gcfs_output_names)
  names(collated_output) <- collated_output_names
  
  # Return collated output (both actual/expected and estimated gCFS)
  return(collated_output)
}






#### Calculate gene concordance factors ####
extract.input.concordance.factors <- function(alignment_path, iqtree2_path, iqtree2_num_threads = "AUTO", rename.taxa.for.ms = TRUE, renamed_taxa){
  # Function to return concordance factors (calculated in iqtree), for a simulated tree
  #     gCFs estimated from the base tree supplied to ms and from the set of gene trees estimated in AliSim
  
  ## Get the directory and list of files in that directory
  al_dir <- paste0(dirname(alignment_path),"/")
  al_files <- list.files(al_dir)
  # Change to that dorectory
  setwd(al_dir)
  # Extract the unique id for this alignment
  al_id <- tail(unlist(strsplit(al_dir, "/")),1)
  
  ## Prepare prefix for gCF run in IQ-Tree
  actual_gcf_prefix <- paste0(al_id, "-actual")
  
  ## Check if the gCF has already been calculated
  check_gcf_files <- grep(actual_gcf_prefix, al_files, value = T)
  
  ## If the file length is 0, then continue
  if (length(check_gcf_files) == 0){
    ## Identify the files for the starting tree and the gene trees
    tree_file <- paste0(al_dir, grep("branch_lengths_modified", al_files, value = T))
    gene_tree_file <- paste0(al_dir, grep("ms_gene_trees", al_files, value = T))
    
    ## Rename taxa (if necessary)
    # If taxa names are provided and rename.taxa.for.ms == TRUE, rename the taxa from the tree_file
    if (rename.taxa.for.ms == TRUE & length(renamed_taxa) > 0){
      # Open the tree
      t <- read.tree(tree_file)
      # Relabel the tips to have the right number
      t$tip.label <- unlist(lapply(t$tip.label, function(x){renamed_taxa[[x]]}))
      # Remove the "t" from the taxa label - ms labels by number only
      t$tip.label <- gsub("t","",t$tip.label)
      # Save the tree
      tree_file_formatted <- gsub(".treefile", "_renamed.treefile", tree_file)
      write.tree(t, file = tree_file_formatted)
    } else {
      tree_file_formatted <- tree_file
    }
    
    ## Calculate gCFS
    actual_gcf_call <- paste0(iqtree2_path ," -t ", tree_file_formatted, " --gcf ", gene_tree_file, " --prefix ", actual_gcf_prefix)
    system(actual_gcf_call)
  }
  
  ## Extract and return gCFs
  # Open the gcf table
  al_files <- paste0(al_dir, list.files(al_dir))
  gcf_files <- grep(actual_gcf_prefix, al_files, value = T)
  gcf_stat_file = grep("cf.stat", gcf_files, value = T)
  gcf_branch_file = grep("cf.branch", gcf_files, value = T)
  gcf_tree_file = grep("cf.tree.nex", grep("cf.tree", gcf_files, value = T), value = T, invert = T)
  gcf_table <- read.table(gcf_stat_file, header = TRUE, sep = "\t")
  # Return the gCF table and file names
  output_gcf_list <- list("gcf_stat_file" = gcf_stat_file, "gcf_branch_file" = gcf_branch_file, "gcf_tree_file" = gcf_tree_file, 
                          "gcf_table" = gcf_table)
  return(output_gcf_list)
}



extract.output.concordance.factors <- function(alignment_path, iqtree2_path, iqtree2_num_threads = "AUTO", iqtree2_model = NA){
  # Function to return concordance factors (calculated in iqtree), for a simulated alignment
  #     gCFs estimated from the ML tree estimated from the partitioned simulated alignment and gene trees estimated from each partition
  
  # Run iqtree to estimate gene concordance factors
  gcf_files <- iqtree2.concordance.factors(alignment_path, iqtree2_path, iqtree2_num_threads, iqtree2_num_ufb, iqtree2_model)
  # Retrieve gcf from output table
  gcf_stat_file <- op_vector[["gCF_table_file"]]
  gcf_table <- read.table(gcf_stat_file, header = TRUE, sep = "\t")
  # Return the gCF table
  # Return the gCF table and file names
  output_gcf_list <- list("gcf_stat_file" = gcf_files[4], "gcf_branch_file" = gcf_files[3], "gcf_tree_file" = gcf_files[2], 
                          "gcf_table" = gcf_table)
  return(output_gcf_list)
}



iqtree2.concordance.factors <- function(alignment_path, iqtree2_path, iqtree2_num_threads = "AUTO", iqtree2_model = NA){
  # Function to take a simulated alignment and estimate gCF from it using iqtree2
  
  ## Get the directory and list of files in that directory
  al_dir <- paste0(dirname(alignment_path),"/")
  al_files <- list.files(al_dir)
  # Change to that directory
  setwd(al_dir)
  # Extract the unique id for this alignment
  al_id <- tail(unlist(strsplit(al_dir, "/")),1)
  
  ## Construct a prefix and check whether the gCF has already been calculate
  check_gcf_tree_prefix <- paste0(al_id, "-concord")
  check_files <- grep(check_gcf_tree_prefix, al_files, value = TRUE)
  
  ## If the estimated gCFs have not been calculated, calculate them now
  if (length(check_files) == 0){
    ## Inferring species tree
    # Check for presence of species tree
    species_tree_file_check <- paste0(dirname(alignment_path), "/", df_row$ID,"_ML_tree.treefile")
    species_tree_iqfile_check <- paste0(dirname(alignment_path), "/", df_row$ID,"_ML_tree.iqtree")
    if (file.exists(species_tree_file_check) == FALSE | file.exists(species_tree_iqfile_check) == FALSE){
      ## If the output files don't exist, call IQ-Tree to estimate species tree
      # Create a gene partition file with no models
      # Find and open the alisim partition file
      alisim_partition_file <- paste0(al_dir, grep("log", grep("partition", al_files, value = TRUE), value = TRUE, invert = TRUE))
      # Generate the gcf partition file
      gcf_partition_file <- generate.gcf.partition.file(alisim_partition_file)
      
      # Create model call
      if (is.na(iqtree2_model) == TRUE){
        model_call <- " -m MFP+MERGE "
      } else if (is.na(iqtree2_model) == FALSE){
        model_call <- paste0(" -m ", iqtree2_model, " ")
      }
      # Create IQ-Tree call
      species_tree_prefix <- paste0(al_id, "-concat")
      species_tree_call <- paste0(iqtree2_path, " -s ", alignment_path, " -p ", gcf_partition_file, model_call, " --prefix ", species_tree_prefix, " -nt ",iqtree2_num_threads)
      system(species_tree_call)
    }
    
    ## Inferring gene/locus trees  
    # Check for presence of gene tree file
    gene_tree_file_check <- paste0(dirname(alignment_path), "/", df_row$ID,"_gene_trees.treefile")
    gene_tree_iqfile_check <- paste0(dirname(alignment_path), "/", df_row$ID,"_gene_trees.iqtree")
    if (file.exists(gene_tree_file_check) == FALSE | file.exists(gene_tree_iqfile_check) == FALSE){
      ## If the output files don't exist, call IQ-Tree to estimate gene trees
      # Create model call
      if (is.na(iqtree2_model) == TRUE){
        model_call <- " -m MFP+MERGE "
      } else if (is.na(iqtree2_model) == FALSE){
        model_call <- paste0(" -m ", iqtree2_model, " ")
      }
      # Create IQ-Tree call
      gene_tree_prefix <- paste0(al_id, "-gene_trees")
      gene_tree_call <- paste0(iqtree2_path, " -s ", alignment_path, " -S ", gcf_partition_file, model_call, " --prefix ", gene_tree_prefix, " -nt ",iqtree2_num_threads)
      system(gene_tree_call)
    }
    
    ## Calculating gene concordance factors
    gcf_tree_prefix <- paste0(al_id, "-concord")
    species_tree_file <- paste0(dirname(alignment_path), "/", df_row$ID,"_ML_tree.treefile")
    gene_tree_file <- paste0(dirname(alignment_path), "/", df_row$ID,"_gene_trees.treefile")
    gcf_call <- paste0(iqtree2_path ," -t ", species_tree_file, " --gcf ", gene_tree_file, " --prefix ", gcf_tree_prefix)
    system(gcf_call)
  }
  
  ## Return output
  # Assemble output files
  gcf_tree_prefix <- paste0(al_id, "-concord")
  gcf_tree_file <- paste0(al_dir, gcf_tree_prefix, ".cf.tree")
  gcf_branch_file <- paste0(al_dir, gcf_tree_prefix, ".cf.branch")
  gcf_table_file <- paste0(al_dir, gcf_tree_prefix, ".cf.stat")
  # Assemble output
  op_vector <- c(alignment_path, gcf_tree_file, gcf_branch_file, gcf_table_file)
  names(op_vector) <- c("alignment_path", "gCF_tree_file", "gCF_branch_file", "gCF_table_file")
  return(op_vector)
}






#### Functions for partition files ####
generate.gcf.partition.file <- function(partition_file){
  # Quick function to create a gCF partition file from an Alisim partition file 
  
  # Open partition file
  alisim_partitions <- readLines(partition_file)
  # Extract gene names
  gene_names <- grep("gene", unlist(strsplit(grep("charset", alisim_partitions, value = T), " ")), value = T)
  # Make a new charpartition section
  new_charpartition <- paste0("\tcharpartition genes = ", paste(gene_names, collapse = ", "), ";")
  # Replace the charpartition section
  line_ind <- grep("charpartition", alisim_partitions)
  alisim_partitions[line_ind] <- new_charpartition
  # Create file path
  p_dir <- paste0(dirname(partition_file), "/")
  p_id <- tail(unlist(strsplit(p_dir, "/")),1)
  gcf_partitions <- paste0(p_dir, p_id, "_gCF_partitions.nexus")
  # Save the new partition file 
  write(alisim_partitions, file = gcf_partitions)
  # Return the gcf partition file path
  return(gcf_partitions)
}






#### Functions for manipulating trees ####
check.tip.labels <- function(tree_path){
  # Function to ensure tip labels in gene trees and trees match before calculating qCF values
  
  # Open tree and get tip labels  
  test_tree <- read.tree(tree_path)
  tree_tips <- test_tree$tip.label
  # Check if tip labels have a "t" in them - e.g. "t72"
  if (TRUE %in% grepl("t", tree_tips)){
    # Replace tree tips without "t"
    new_tree_tips <- gsub("t", "", tree_tips)
    # Replace tree tips in original tree
    test_tree$tip.label <- new_tree_tips
    # Create new file path
    split_path <- unlist(strsplit(basename(tree_path), "\\."))
    relabelled_tree_path <- paste0(dirname(tree_path), 
                                   paste(head(split_path, (length(split_path) - 1)), collapse = "."), 
                                   "_relabelled.", 
                                   tail(split_path, 1))
    # Save the relabelled tree
    write.tree(test_tree, file = relabelled_tree_path)
  } else {
    relabelled_tree_path <- tree_path
  }
  # Return the relabelled tree path and use that to calculate qcfs
  return(relabelled_tree_path)
}

