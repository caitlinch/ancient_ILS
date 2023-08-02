# ancient_ILS/code/func_analysis.R
## This script includes functions investigate and analyses the results of the simulations 
# Caitlin Cherryh, 2023

library(ape)
library(phangorn)

#### Analysis wrapper function ####
analysis.wrapper <- function(row_id, df, ASTRAL_path, hypothesis_tree_dir, renamed_taxa){
  # Wrapper function to calculate:
  #   - [x] actual and estimated qcfs (ASTRAL)
  #   - [x] hypothesis tree distances for ASTRAL tree
  #   - [x] hypothesis test in IQ-Tree to see if any hypothesis can be ignored for this alignment
  
  ## Open the row of interest
  df_row <- df[row_id, ]
  
  ## Calculate the differences between the three trees
  astral_tree_diffs <- calculate.distance.between.three.trees(tree_path = df_row$ASTRAL_tree_treefile, hypothesis_tree_dir, tree_type = "ASTRAL", 
                                                              rename.hypothesis.tree.tips = TRUE, renamed_taxa = renamed_taxa)
  
  ## Calculate the qCFs using ASTRAL
  astral_qcfs <- qcf.wrapper(ID = df_row$ID, starting_tree = df_row$output_base_tree_file, ms_gene_trees = df_row$output_gene_tree_file,
                             ASTRAL_tree = df_row$ASTRAL_tree_treefile, ML_gene_trees = df_row$iqtree2_gene_tree_treefile, 
                             ASTRAL_path = ASTRAL_path, call.astral = TRUE, renamed_taxa)
  
  ## Calculate qcf for branch leading to Ctenophora+Porifera clade for both ms and iqtree gene trees
  # Identify correct file path for hypothesis tree
  all_hyp_paths <- paste0(hypothesis_tree_dir, grep("numeric|relabelled", list.files(hypothesis_tree_dir), value = T, invert = T))
  Cten_hyp_tree_path <- grep("ASTRAL", grep("_Cten.tre", all_hyp_paths, value = T), value = T)
  Pori_hyp_tree_path <- grep("ASTRAL", grep("_Pori.tre", all_hyp_paths, value = T), value = T)
  CtenPori_hyp_tree_path <- grep("ASTRAL", grep("_CtenPori.tre", all_hyp_paths, value = T), value = T)
  # Calculate for actual (i.e. ms) gene trees
  actual_hyp1_qcf <- check.qcf.Cten(output_id = "actual_testHyp1_Cten", 
                                    gene_trees_path = df_row$output_gene_tree_file,
                                    Cten_hypothesis_tree_path = Cten_hyp_tree_path, 
                                    ASTRAL_path = ASTRAL_path, 
                                    call.astral = TRUE, 
                                    renamed_taxa = renamed_taxa)
  actual_hyp2_qcf <- check.qcf.Pori(output_id = "actual_testHyp2_Pori", 
                                    gene_trees_path = df_row$output_gene_tree_file,
                                    Pori_hypothesis_tree_path = Pori_hyp_tree_path, 
                                    ASTRAL_path = ASTRAL_path, 
                                    call.astral = TRUE, 
                                    renamed_taxa = renamed_taxa)
  actual_hyp3_qcf <- check.qcf.CtenPori(output_id = "actual_testHyp3_CtenPori", 
                                        gene_trees_path = df_row$output_gene_tree_file, 
                                        CtenPori_hypothesis_tree_path = CtenPori_hyp_tree_path, 
                                        ASTRAL_path = ASTRAL_path, 
                                        call.astral = TRUE, 
                                        renamed_taxa = renamed_taxa)
  actual_clade_qcf <- check.qcf.clades(output_id = "actual_clade",
                                       gene_trees_path = df_row$output_gene_tree_file,
                                       Cten_hypothesis_tree_path = Cten_hyp_tree_path,
                                       ASTRAL_path = ASTRAL_path, 
                                       call.astral = TRUE, 
                                       renamed_taxa = renamed_taxa)
  # Calculate for estimated (i.e. iqtree2) gene trees
  estimated_hyp1_qcf <- check.qcf.Cten(output_id = "estimated_testHyp1_Cten", 
                                       gene_trees_path = df_row$iqtree2_gene_tree_treefile, 
                                       Cten_hypothesis_tree_path = Cten_hyp_tree_path, 
                                       ASTRAL_path = ASTRAL_path, 
                                       call.astral = TRUE, 
                                       renamed_taxa = renamed_taxa)
  estimated_hyp2_qcf <- check.qcf.Pori(output_id = "estimated_testHyp2_Pori", 
                                       gene_trees_path = df_row$iqtree2_gene_tree_treefile, 
                                       Pori_hypothesis_tree_path = Pori_hyp_tree_path, 
                                       ASTRAL_path = ASTRAL_path, 
                                       call.astral = TRUE, 
                                       renamed_taxa = renamed_taxa)
  estimated_hyp3_qcf <- check.qcf.CtenPori(output_id = "estimated_testHyp3_CtenPori", 
                                           gene_trees_path = df_row$iqtree2_gene_tree_treefile, 
                                           CtenPori_hypothesis_tree_path = CtenPori_hyp_tree_path, 
                                           ASTRAL_path = ASTRAL_path, 
                                           call.astral = TRUE, 
                                           renamed_taxa = renamed_taxa)
  estimated_clade_qcf <- check.qcf.clades(output_id = "estimated_clade",
                                          gene_trees_path = df_row$iqtree2_gene_tree_treefile,
                                          Cten_hypothesis_tree_path = Cten_hyp_tree_path,
                                          ASTRAL_path = ASTRAL_path, 
                                          call.astral = TRUE, 
                                          renamed_taxa = renamed_taxa)
  
  ## Get lengths and maximum branching time for starting tree
  # Root trees at Choanoflagellates
  start_tree <- read.tree(df_row$output_base_tree_file)
  astral_tree <- read.tree(df_row$ASTRAL_tree_treefile)
  tree_length <- c(extract.tree.length(start_tree), 
                   extract.tree.depth(start_tree,  c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), root.tree = TRUE),
                   extract.tree.length(astral_tree)[2])
  names(tree_length) <- c(paste0("actual_", names(tree_length)[1:5]), paste0("estimated_", names(tree_length)[6]))
  
  ## Output results
  # Trim unwanted columns
  trimmed_df_row <- df_row[, c("dataset", "dataset_type", "ID", "simulation_number", "simulation_type",
                               "hypothesis_tree", "replicates", "num_taxa", "num_genes", "gene_length", "num_sites",
                               "ML_tree_estimation_models", "branch_a_length", "branch_b_length", "branch_c_length",
                               "branch_all_animals_length", "branch_bilateria_length", "branch_cnidaria_length",
                               "branch_outgroup_length", "branch_porifera_length", "total_tree_length",
                               "sum_internal_branch_lengths", "percentage_internal_branch_lengths")]
  # Assemble output vector
  analysis_output         <- c(as.character(trimmed_df_row), astral_tree_diffs,  astral_qcfs,
                               actual_hyp1_qcf, actual_hyp2_qcf, actual_hyp3_qcf, actual_clade_qcf,
                               estimated_hyp1_qcf, estimated_hyp2_qcf, estimated_hyp3_qcf, estimated_clade_qcf,
                               tree_length)
  names(analysis_output)  <- c(names(trimmed_df_row), names(astral_tree_diffs),  names(astral_qcfs),
                               names(actual_hyp1_qcf), names(actual_hyp2_qcf), names(actual_hyp3_qcf), names(actual_clade_qcf),
                               names(estimated_hyp1_qcf), names(estimated_hyp2_qcf), names(estimated_hyp3_qcf), names(estimated_clade_qcf),
                               names(tree_length))
  ## Return the output
  return(analysis_output)
}



topology.tests.wrapper <- function(row_id, df, hypothesis_tree_dir, test.three.hypothesis.trees = TRUE){
  # Wrapper function to calculate:
  #   - [x] topology tests in IQ-Tree
  
  ## Open the row of interest
  df_row <- df[row_id, ]
  
  # Extract all files from folder
  all_files <- list.files(hypothesis_tree_dir)
  # Extract hypothesis tree files
  if (test.three.hypothesis.trees == TRUE){
    # Test all three hypothesis trees
    astral_hyp_trees <- paste0(hypothesis_tree_dir, grep("ASTRAL_hypothesis_trees_relabelled_tbl", all_files, value = TRUE))
  } else {
    # Test only two hypothesis trees
    astral_hyp_trees <- paste0(hypothesis_tree_dir, grep("ASTRAL_hypothesis_trees_2tree_relabelled_tbl", all_files, value = TRUE))
  }
  # Perform hypothesis tests
  astral_hyp_tests <- perform.hypothesis.tests(ID = paste0(df_row$ID, "_ASTRAL_hyps"), alignment_path = df_row$output_alignment_file, hypothesis_tree_file = astral_hyp_trees,
                                               iqtree2_path = df_row$iqtree2, iqtree2_num_threads = df_row$iqtree2_num_threads, iqtree2_model = df_row$ML_tree_estimation_models)
  # Extract hypothesis test results
  hypothesis_tests <- c(astral_hyp_tests[["t1_logL"]], astral_hyp_tests[["t2_logL"]], astral_hyp_tests[["t3_logL"]],
                        astral_hyp_tests[["t1_deltaL"]], astral_hyp_tests[["t2_deltaL"]], astral_hyp_tests[["t3_deltaL"]],
                        astral_hyp_tests[["t1_bpRELL"]], astral_hyp_tests[["t2_bpRELL"]], astral_hyp_tests[["t3_bpRELL"]],
                        astral_hyp_tests[["t1_cELW"]], astral_hyp_tests[["t2_cELW"]], astral_hyp_tests[["t3_cELW"]],
                        astral_hyp_tests[["t1_pAU"]], astral_hyp_tests[["t2_pAU"]], astral_hyp_tests[["t3_pAU"]])
  trimmed_hypothesis_test_names <- c("t1_logL", "t2_logL", "t3_logL", "t1_deltaL", "t2_deltaL", "t3_deltaL",
                                     "t1_bpRELL", "t2_bpRELL", "t3_bpRELL", "t1_cELW", "t2_cELW", "t3_cELW",
                                     "t1_pAU", "t2_pAU", "t3_pAU")
  names(hypothesis_tests) <- c(paste0("ASTRAL_hyp_", trimmed_hypothesis_test_names))
  
  ## Trim unwanted columns
  trimmed_df_row <- df_row[, c("dataset", "dataset_type", "ID", "output_folder", "simulation_number", "simulation_type",
                               "hypothesis_tree", "replicates", "num_taxa", "num_genes", "gene_length", "num_sites",
                               "ML_tree_estimation_models", "branch_a_length", "branch_b_length", "branch_c_length",
                               "branch_all_animals_length", "branch_bilateria_length", "branch_cnidaria_length",
                               "branch_outgroup_length", "branch_porifera_length", "total_tree_length",
                               "sum_internal_branch_lengths", "percentage_internal_branch_lengths")]
  ## Assemble output vector
  analysis_output         <- c(as.character(trimmed_df_row), hypothesis_tests)
  names(analysis_output)  <- c(names(trimmed_df_row), names(hypothesis_tests))
  
  ## Return the output
  return(analysis_output)
}



#### Calculate actual qCF ####
calculate.actual.qCF <- function(row_id, df, call.ASTRAL = TRUE, renamed_taxa){
  # Wrapper function to calculate actual qCFS 
  
  ## Open the row of interest
  df_row <- df[row_id, ]
  
  ## Calculate the qCFs using ASTRAL
  ## Identify directory
  qcf_dir <- paste0(dirname(df_row$output_base_tree_file), "/")
  ## Calculate actual quartet concordance factors
  actual_qcf_paths <- qcf.call(output_id = paste0(df_row$ID, "_actual_qcfs"), output_directory = qcf_dir, 
                                 tree_path = df_row$output_base_tree_file, gene_trees_path = df_row$output_gene_tree_file,
                                 ASTRAL_path = df_row$ASTRAL, call.astral = call.ASTRAL, 
                                 rename.tree.tips = TRUE, renamed_taxa)
  ## Extract relevant qCF values
  actual_qcf_values <-  extract.qcf.values(qcf_tree_path = actual_qcf_paths[["qcf_output_tree"]], 
                                             qcf_log_path = actual_qcf_paths[["qcf_output_log"]])
  ## Assemble output
  qcf_output <- c(actual_qcf_paths, actual_qcf_values)
  names(qcf_output) <- c(paste0("actual_", names(actual_qcf_paths)), paste0("actual_", names(actual_qcf_values)))
  
  ## Trim unwanted columns
  trimmed_df_row <- df_row[, c("dataset", "dataset_type", "ID", "output_folder", "simulation_number", "simulation_type",
                               "hypothesis_tree", "replicates", "num_taxa", "num_genes", "gene_length", "num_sites",
                               "ML_tree_estimation_models", "branch_a_length", "branch_b_length", "branch_c_length",
                               "branch_all_animals_length", "branch_bilateria_length", "branch_cnidaria_length",
                               "branch_outgroup_length", "branch_porifera_length", "total_tree_length",
                               "sum_internal_branch_lengths", "percentage_internal_branch_lengths")]
  ## Assemble output vector
  analysis_output         <- c(as.character(trimmed_df_row), qcf_output)
  names(analysis_output)  <- c(names(trimmed_df_row), names(qcf_output))
  
  ## Return the output
  return(analysis_output)
}





#### Perform hypothesis tests in IQ-Tree2 ####
perform.hypothesis.tests <- function(ID, alignment_path, hypothesis_tree_file, iqtree2_path, iqtree2_num_threads, iqtree2_model = NA){
  # Function to perform hypothesis tests in IQ-Tree2 on the simulated alignment with either 2 or 3 hypothesis trees
  
  ## Change directories to the same directory as the alignment
  alignment_directory <- paste0(dirname(alignment_path), "/")
  setwd(alignment_directory)
  
  ## Assemble the IQ-Tree command and call IQ-Tree
  # Specify model, if one is provided
  if (is.na(iqtree2_model) == TRUE){
    model_code = "MFP"
  } else {
    model_code <- iqtree2_model
  }
  test_command <- paste0(iqtree2_path, " -s ", alignment_path, " -z ", hypothesis_tree_file, " -m ", model_code,
                         " -n 0 -zb 10000 -zw -au -nt ", iqtree2_num_threads, " -pre ", ID, " -safe")
  system(test_command)
  
  ## Extract output from the IQ-Tree file
  # List all files in the directory
  all_files <- list.files(alignment_directory)
  # Find .iqtree file
  htest_files <- grep(ID, all_files, value = TRUE)
  htest_iqfile <- paste0(alignment_directory, grep("\\.iqtree", htest_files, value = TRUE))
  htest_trees <- paste0(alignment_directory, grep("\\.trees", htest_files, value = TRUE))
  # Open .iqtree file
  iq_lines <- readLines(htest_iqfile)
  # Find headings for table
  header_ind <- intersect(intersect(grep("Tree", iq_lines), grep("logL", iq_lines)), 
                          intersect(grep("p-AU", iq_lines), grep("bp-RELL", iq_lines)) )
  # Split lines for trees into components 
  t1_line <- iq_lines[header_ind+2]
  t2_line <- iq_lines[header_ind+3]
  t3_line <- iq_lines[header_ind+4]
  t1_split_raw <- unlist(strsplit(t1_line, " "))
  t2_split_raw <- unlist(strsplit(t2_line, " "))
  t3_split_raw <- unlist(strsplit(t3_line, " "))
  # Remove strings - "", "-", "_", "+"
  t1_split <- gsub("+","_", t1_split_raw, fixed = TRUE)
  t1_split <- t1_split[which(t1_split != "")]
  t1_split <- t1_split[which(t1_split != "_")]
  t1_split <- t1_split[which(t1_split != "-")]
  t2_split <- gsub("+","_", t2_split_raw, fixed = TRUE)
  t2_split <- t2_split[which(t2_split != "")]
  t2_split <- t2_split[which(t2_split != "_")]
  t2_split <- t2_split[which(t2_split != "-")]
  t3_split <- gsub("+","_", t3_split_raw, fixed = TRUE)
  t3_split <- t3_split[which(t3_split != "")]
  t3_split <- t3_split[which(t3_split != "_")]
  t3_split <- t3_split[which(t3_split != "-")]
  
  ## Format output
  # Assemble output
  test_output <- c(t1_split[2], t2_split[2], t3_split[2],
                   t1_split[3], t2_split[3], t3_split[3],
                   t1_split[4], t2_split[4], t3_split[4],
                   t1_split[5], t2_split[5], t3_split[5],
                   t1_split[6], t2_split[6], t3_split[6],
                   t1_split[7], t2_split[7], t3_split[7],
                   t1_split[8], t2_split[8], t3_split[8],
                   t1_split[9], t2_split[9], t3_split[9],
                   t1_split[10], t2_split[10], t3_split[10],
                   htest_trees)
  output_names <- c("t1_logL", "t2_logL", "t3_logL",
                    "t1_deltaL", "t2_deltaL", "t3_deltaL",
                    "t1_bpRELL", "t2_bpRELL", "t3_bpRELL",
                    "t1_pKH", "t2_pKH", "t3_pKH",
                    "t1_pSH", "t2_pSH", "t3_pSH",
                    "t1_pWKH", "t2_pWKH", "t3_pWKH",
                    "t1_pWSH", "t2_pWSH", "t3_pWSH",
                    "t1_cELW", "t2_cELW", "t3_cELW",
                    "t1_pAU", "t2_pAU", "t3_pAU",
                    "trees")   
  names(test_output) <- output_names
  ## Return output
  return(test_output)
}






#### Calculate distance between trees ####
calculate.distance.between.three.trees <- function(tree_path, hypothesis_tree_dir, tree_type, rename.hypothesis.tree.tips = FALSE, renamed_taxa){
  # Function to calculate the distance between a tree and three hypothesis trees
  # tree_type = "ASTRAL" | tree_type = "IQ-Tree2"
  
  ## Open the hypothesis trees
  # Extract the tree files
  all_hyp_dir_files <- list.files(hypothesis_tree_dir)
  # Remove relabelled trees
  all_hyp_dir_files <- grep("numericTipLabels|relabelled", all_hyp_dir_files, value = T, invert = T)
  # Extract only either the ASTRAL or IQ-Tree tree files
  if (tree_type == "ASTRAL" | tree_type == "astral"){
    tree_files <- grep("\\.tre", grep("ASTRAL", all_hyp_dir_files,value = T), value = T)
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
  # Reroot the trees (only at one outgroup taxon, in case outgroup is paraphyletic)
  h1_tree <- root(h1_tree, outgroup = c("Monosiga_ovata"),
                  resolve.root = TRUE)
  h2_tree <- root(h2_tree, outgroup = c("Monosiga_ovata"),
                  resolve.root = TRUE)
  h3_tree <- root(h3_tree, outgroup = c("Monosiga_ovata"),
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
  # Root at Monosiga ovata
  t <- root(t, outgroup = gsub("t", "", renamed_taxa["Monosiga_ovata"]),
            resolve.root = TRUE)
  
  ## Calculate RF and wRF for the trees
  output_vec <- c(tree_path, 
                  RF.dist(t, h1_tree, normalize = TRUE), round(wRF.dist(t, h1_tree, normalize = TRUE), digits = 2), 
                  RF.dist(t, h2_tree, normalize = TRUE), round(wRF.dist(t, h2_tree, normalize = TRUE), digits = 2), 
                  RF.dist(t, h3_tree, normalize = TRUE), round(wRF.dist(t, h3_tree, normalize = TRUE), digits = 2))
  names(output_vec) <- c(paste0(tree_type, "_tree_path"), 
                         paste0(tree_type, "_h1_norm_RF_dist"), paste0(tree_type, "_h1_norm_wRF_dist"), 
                         paste0(tree_type, "_h2_norm_RF_dist"), paste0(tree_type, "_h2_norm_wRF_dist"), 
                         paste0(tree_type, "_h3_norm_RF_dist"), paste0(tree_type, "_h3_norm_wRF_dist") )
  
  ## Return the RF/wRF distances
  return(output_vec)
}






#### Quartet concordance factors ####
qcf.wrapper <- function(ID, starting_tree, ms_gene_trees, ASTRAL_tree, ML_gene_trees, ASTRAL_path, call.astral = TRUE, renamed_taxa){
  ## Function to determine actual and estimated qCF in ASTRAL and return summary statistics
  
  ## Identify directory
  qcf_dir <- paste0(dirname(starting_tree), "/")
  
  ## Calculate actual quartet concordance factors
  actual_qcf_paths <- qcf.call(output_id = paste0(ID, "_actual_qcfs"), output_directory = qcf_dir, 
                               tree_path = starting_tree, gene_trees_path = ms_gene_trees,
                               ASTRAL_path = ASTRAL_path, call.astral, 
                               rename.tree.tips = TRUE, renamed_taxa)
  
  ## Calculate estimated quartet concordance factors
  estimated_qcf_paths <- qcf.call(output_id = paste0(ID, "_estimated_qcfs"), output_directory = qcf_dir, 
                                  tree_path = ASTRAL_tree, gene_trees_path = ML_gene_trees,
                                  ASTRAL_path = ASTRAL_path, call.astral,
                                  rename.tree.tips = FALSE, renamed_taxa)
  
  ## Extract relevant qCF values
  actual_qcf_values <-  extract.qcf.values(qcf_tree_path = actual_qcf_paths[["qcf_output_tree"]], 
                                           qcf_log_path = actual_qcf_paths[["qcf_output_log"]])
  estimated_qcf_values <-  extract.qcf.values(qcf_tree_path = estimated_qcf_paths[["qcf_output_tree"]], 
                                              qcf_log_path = estimated_qcf_paths[["qcf_output_log"]])
  
  ## Assemble output
  qcf_output <- c(actual_qcf_values, estimated_qcf_values)
  names(qcf_output) <- c(paste0("actual_", names(actual_qcf_values)) , paste0("estimated_", names(estimated_qcf_values)) )
  
  ## Return output
  return(qcf_output)
}



check.qcf.Cten <- function(output_id, gene_trees_path, Cten_hypothesis_tree_path, ASTRAL_path, call.astral = TRUE, renamed_taxa){
  # Function to check the qCF for the Ctenophora+Porifera branch by using the Ctenophora+Porifera tree
  
  ## Prepare tips for hypothesis
  # Identify tip labels for hypothesis branches
  BilatCnidPori_tips <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                          "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", 
                          "39", "40")
  Outgroup_tips = c("71", "72", "73", "74", "75")
  
  ## Calculate qCF
  # Set output directory - save in same file as gene_trees_path
  qcf_dir = paste0(dirname(gene_trees_path), "/")
  # Call ASTRAL with the gene tree file as the gene trees and the CtenPori hypothesis tree as the species tree
  qcf_paths <- qcf.call(output_id = output_id, output_directory = qcf_dir, 
                        tree_path = Cten_hypothesis_tree_path, gene_trees_path = gene_trees_path,
                        ASTRAL_path = ASTRAL_path, call.astral = call.astral, 
                        rename.tree.tips = TRUE, renamed_taxa)
  # Identify and open ASTRAL tree with calculated qcf values
  qcf_tree <- read.tree(qcf_paths[["qcf_output_tree"]])
  # Root tree at outgroup
  qcf_tree <- root(qcf_tree, Outgroup_tips, resolve.root = TRUE)
  
  ## BRANCH A (ILS - all other metazoans)
  # Identify start and end of branch
  a_branch_tipwise <- getMRCA(qcf_tree, BilatCnidPori_tips)
  a_branch_rootwise <- qcf_tree$edge[which(qcf_tree$edge[,2] == a_branch_tipwise), 1]
  # Extract qcf from the branch leading to BilatCnidCten clade
  a_qcf_value <- qcf_tree$node.label[ ( a_branch_rootwise - Ntip(qcf_tree) ) ]
  
  ## Assemble the output
  qcf_op <- c(a_qcf_value)
  names(qcf_op) <- paste0(output_id, "_branch_a_qcf_value")
  # Return output
  return(qcf_op)
}



check.qcf.Pori <- function(output_id, gene_trees_path, Pori_hypothesis_tree_path, ASTRAL_path, call.astral = TRUE, renamed_taxa){
  # Function to check the qCF for the Ctenophora+Porifera branch by using the Ctenophora+Porifera tree
  
  ## Prepare tips for hypothesis
  # Identify tip labels for hypothesis branches
  BilatCnidCten_tips <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                          "21","41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57",
                          "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70")
  Outgroup_tips = c("71", "72", "73", "74", "75")
  
  ## Calculate qCF
  # Set output directory - save in same file as gene_trees_path
  qcf_dir = paste0(dirname(gene_trees_path), "/")
  # Call ASTRAL with the gene tree file as the gene trees and the CtenPori hypothesis tree as the species tree
  qcf_paths <- qcf.call(output_id = output_id, output_directory = qcf_dir, 
                        tree_path = Pori_hypothesis_tree_path, gene_trees_path = gene_trees_path,
                        ASTRAL_path = ASTRAL_path, call.astral = call.astral, 
                        rename.tree.tips = TRUE, renamed_taxa)
  # Identify and open ASTRAL tree with calculated qcf values
  qcf_tree <- read.tree(qcf_paths[["qcf_output_tree"]])
  # Root tree at outgroup
  qcf_tree <- root(qcf_tree, Outgroup_tips, resolve.root = TRUE)
  
  ## BRANCH A (ILS - all other metazoans)
  # Identify start and end of branch
  a_branch_tipwise <- getMRCA(qcf_tree, BilatCnidCten_tips)
  a_branch_rootwise <- qcf_tree$edge[which(qcf_tree$edge[,2] == a_branch_tipwise), 1]
  # Extract qcf from the branch leading to BilatCnidCten clade
  a_qcf_value <- qcf_tree$node.label[ ( a_branch_rootwise - Ntip(qcf_tree) ) ]
  
  ## Assemble the output
  qcf_op <- c(a_qcf_value)
  names(qcf_op) <- paste0(output_id, "_branch_a_qcf_value")
  # Return output
  return(qcf_op)
}



check.qcf.CtenPori <- function(output_id, gene_trees_path, CtenPori_hypothesis_tree_path, ASTRAL_path, call.astral = TRUE, renamed_taxa){
  # Function to check the qCF for the Ctenophora+Porifera branch by using the Ctenophora+Porifera tree
  
  ## Prepare tips for hypothesis
  # Identify tip labels for the possible clade made up of Ctenophora and Porifera
  CtenPori_tips = c("22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",
                    "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59",
                    "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70")
  Outgroup_tips = c("71", "72", "73", "74", "75")
  
  ## Calculate qCF
  # Set output directory - save in same file as gene_trees_path
  qcf_dir = paste0(dirname(gene_trees_path), "/")
  # Call ASTRAL with the gene tree file as the gene trees and the CtenPori hypothesis tree as the species tree
  qcf_paths <- qcf.call(output_id = output_id, output_directory = qcf_dir, 
                        tree_path = CtenPori_hypothesis_tree_path, gene_trees_path = gene_trees_path,
                        ASTRAL_path = ASTRAL_path, call.astral = call.astral, 
                        rename.tree.tips = TRUE, renamed_taxa)
  # Identify and open ASTRAL tree with calculated qcf values
  qcf_tree <- read.tree(qcf_paths[["qcf_output_tree"]])
  # Root tree at outgroup
  qcf_tree <- root(qcf_tree, Outgroup_tips, resolve.root = TRUE)
  
  ## Calculate number of gene trees with Ctenophora and Porifera as a monophyletic clade
  # Identify start and end of branch
  CtenPori_branch_end <- getMRCA(qcf_tree, CtenPori_tips)
  CtenPori_branch_start <- qcf_tree$edge[which(qcf_tree$edge[,2] == CtenPori_branch_end), 1]
  # Extract qcf from the branch leading to CtenPori clade
  CtenPori_qcf_value <- qcf_tree$node.label[ ( CtenPori_branch_end - Ntip(qcf_tree) ) ]
  
  ## Assemble the output
  qcf_op <- c(CtenPori_qcf_value)
  names(qcf_op) <- paste0(output_id, "_qcf_value")
  # Return output
  return(qcf_op)
}



check.qcf.clades <- function(output_id, gene_trees_path, Cten_hypothesis_tree_path, ASTRAL_path, call.astral = TRUE, renamed_taxa){
  # Function to check the qCF for the Ctenophora+Porifera branch by using the Ctenophora+Porifera tree
  
  ## Prepare tips for hypothesis
  # Identify tip labels for the possible clade made up of Ctenophora and Porifera
  Bilateria_tips = c("1", "2", "3", "4", "5", "6")
  Cnidaria_tips = c("7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21")
  Porifera_tips = c("22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40")
  Ctenophora_tips = c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59",
                      "60", "61", "62", "63", "64","65", "66", "67", "68", "69", "70")
  Outgroup_tips = c("71", "72", "73", "74", "75")
  BilatCnid_tips = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21")
  Animals_tips = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21",
                   "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",
                   "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59",
                   "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70")
  
  ## Calculate qCF
  # Set output directory - save in same file as gene_trees_path
  qcf_dir = paste0(dirname(gene_trees_path), "/")
  # Call ASTRAL with the gene tree file as the gene trees and the CtenPori hypothesis tree as the species tree
  qcf_paths <- qcf.call(output_id = output_id, output_directory = qcf_dir, 
                        tree_path = Cten_hypothesis_tree_path, gene_trees_path = gene_trees_path,
                        ASTRAL_path = ASTRAL_path, call.astral = call.astral, 
                        rename.tree.tips = TRUE, renamed_taxa)
  # Identify and open ASTRAL tree with calculated qcf values
  qcf_tree <- read.tree(qcf_paths[["qcf_output_tree"]])
  # Root tree at outgroup
  qcf_tree <- root(qcf_tree, Outgroup_tips, resolve.root = TRUE)
  
  ## Calculate qCF value for each clade
  clade_qcfs <- c(extract.clade.qcf.value(qcf_tree, Bilateria_tips),
                  extract.clade.qcf.value(qcf_tree, Cnidaria_tips),
                  extract.clade.qcf.value(qcf_tree, Porifera_tips),
                  extract.clade.qcf.value(qcf_tree, Ctenophora_tips),
                  extract.clade.qcf.value(qcf_tree, Outgroup_tips),
                  extract.clade.qcf.value(qcf_tree, BilatCnid_tips),
                  extract.clade.qcf.value(qcf_tree, Animals_tips))
  
  ## Assemble the output
  names(clade_qcfs) <- paste0(output_id, "_", c("BILAT", "CNID", "PORI", "CTEN", "CHOANO", "BILAT_CTEN", "ANIMALS"), "_qcf_value")
  # Return output
  return(clade_qcfs)
}



extract.clade.qcf.value <- function(qcf_tree, tips){
  # Quick function to take a monophyletic clade of specified tips and return the qcf value for that branch
  
  # Find start and end of branch
  branch_end <- getMRCA(qcf_tree, tips)
  branch_start <- qcf_tree$edge[which(qcf_tree$edge[,2] == branch_end), 1]
  # Extract qcf from the branch leading to the tips clade
  qcf_value <- qcf_tree$node.label[ ( branch_end - Ntip(qcf_tree) ) ]
  # Return the extracted qCF value
  return(qcf_value)
}



update.tree.tips <- function(tree_path, renamed_taxa){
  # Open tree
  og_tree <- read.tree(tree_path)
  # Update tip labels
  new_tip_labels <- unlist(lapply(og_tree$tip.label, function(i){renamed_taxa[[i]]}))
  # Replace og tip labels with new tip labels
  numericTipLabels_tree <- og_tree
  numericTipLabels_tree$tip.label <- new_tip_labels
  # Create file path to save new tree
  split_path <- unlist(strsplit(basename(tree_path), "\\."))
  numericTipLabels_tree_path <- paste0(dirname(tree_path), "/",
                                       paste(head(split_path, (length(split_path) - 1)), collapse = "."), 
                                       "_numericTipLabels", ".",
                                       tail(split_path, 1))
  # Save new tree
  write.tree(numericTipLabels_tree, numericTipLabels_tree_path)
  # Update tree path file
  tree_path <- numericTipLabels_tree_path
  # Remove any "t" in tip labels
  input_tree_path <- check.tip.labels(tree_path)
  # Return new path
  return(input_tree_path)
}



qcf.call <- function(output_id, output_directory, tree_path, gene_trees_path, ASTRAL_path, call.astral = TRUE, rename.tree.tips = TRUE, renamed_taxa, redo = FALSE){
  ## Function to call ASTRAL and calculate quartet concordance factors
  
  ## Check whether the tree tips match
  # Replace tip labels with numeric labels
  if (rename.tree.tips == TRUE){
    tree_path <- update.tree.tips(tree_path, renamed_taxa)
  }
  # Remove any "t" in tip labels
  input_tree_path <- check.tip.labels(tree_path)
  
  # Assemble ASTRAL command to calculate quartet concordance factors
  output_tree <- paste0(output_directory, output_id, ".tre")
  output_log <- paste0(output_directory, output_id, ".log")
  qcf_command <- paste0("java -jar ", ASTRAL_path, " -q ", input_tree_path, " -i ", gene_trees_path, " -o ", output_tree, " 2> ", output_log)
  # if call.astral == TRUE, run ASTRAL to calculate quartet concordance factors
  if ((call.astral == TRUE) & (file.exists(output_tree) == FALSE) & (file.exists(output_log) == FALSE)){
    system(qcf_command)
  } else if (redo == TRUE){
    system(qcf_command)
  }
  # Assemble output vector
  qcf_output <- c(output_id, qcf_command, call.astral, output_tree, output_log)
  names(qcf_output) <- c("qcf_ID", "astral_qCF_command", "astral_qCF_run", "qcf_output_tree", "qcf_output_log")
  # Return the output
  return(qcf_output)
}



extract.qcf.values <- function(qcf_tree_path, qcf_log_path){
  ## Function to open an ASTRAL output qCF tree and log file, extract the values of interest and return some summary statistics
  
  ## Open the log file and extract some values of interest
  log_lines <- readLines(qcf_log_path)
  num_quartet_trees_ind <- grep("Number of quartet trees in the gene trees:", log_lines)
  num_quartet_trees <- gsub(" ", "", unlist(strsplit(log_lines[num_quartet_trees_ind], ":"))[2])
  final_quartet_score_ind <- grep("Final quartet score is:", log_lines)
  final_quartet_score <- gsub(" ", "", unlist(strsplit(log_lines[final_quartet_score_ind], ":"))[2])
  normalised_quartet_score_ind <- grep("Final normalized quartet score is:", log_lines)
  normalised_quartet_score <- gsub(" ", "", unlist(strsplit(log_lines[normalised_quartet_score_ind], ":"))[2])
  
  ## Open the tree and extract quartet concordance factors for specific branches
  # Read in the qcf output tree
  qcf_tree <- read.tree(file = qcf_tree_path)
  # Extract the node labels (i.e. qCF values) from the qCF tree
  node_labs <- qcf_tree$node.label
  # Extract summary statistics (mean, median, minimum, maximum)
  qcf_summary_stats <- c(mean(as.numeric(node_labs), na.rm = TRUE), median(as.numeric(node_labs), na.rm = TRUE),
                         min(as.numeric(node_labs), na.rm = TRUE), max(as.numeric(node_labs), na.rm = TRUE))
  names(qcf_summary_stats) <- c("qcf_mean", "qcf_median", "qcf_min", "qcf_max")
  
  ## Assemble the output
  qcf_extracted <- c(num_quartet_trees, final_quartet_score, normalised_quartet_score, qcf_summary_stats)
  names(qcf_extracted) <- c("num_quartet_trees", "final_quartet_score", "final_normalised_quartet_score", names(qcf_summary_stats))
  
  ## Return the qCF values 
  return(qcf_extracted)
}



qcf.branch.values <- function(qcf_tree_path){
  ## Function to test for monophyly of a single clade, extract the values of interest and return some summary statistics
  
  # Read in the qcf output tree
  qcf_tree <- read.tree(file = qcf_tree_path)
  # Extract the node labels (i.e. qCF values) from the qCF tree
  node_labs <- qcf_tree$node.label
  # Find the "a" branch and the "b" branch
  if (grepl("h1", basename(qcf_tree_path)) == TRUE){
    # Extract the nodes by checking for monophyletic clades
    a_end <- getMRCA(qcf_tree, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                 "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
                                 "39", "40"))
    a_start <- qcf_tree$edge[which(qcf_tree$edge[,2] == a_end),1]
    branch_a_mono <- is.monophyletic(qcf_tree, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                                 "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
                                                 "39", "40"))
    b_end <- getMRCA(qcf_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                                 "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
    b_start <- qcf_tree$edge[which(qcf_tree$edge[,2] == b_end),1]
    branch_b_mono <- is.monophyletic(qcf_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                                                 "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
  } else if (grepl("h2", basename(qcf_tree_path)) == TRUE){
    a_end <- getMRCA(qcf_tree, c("22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40"))
    a_start <- qcf_tree$edge[which(qcf_tree$edge[,2] == a_end),1]
    branch_a_mono <- is.monophyletic(qcf_tree, c("22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40"))
    b_end <- getMRCA(qcf_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                                 "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
    b_start <- qcf_tree$edge[which(qcf_tree$edge[,2] == b_end),1]
    branch_b_mono <- is.monophyletic(qcf_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                                                 "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
  } else if (grepl("h3", basename(qcf_tree_path)) == TRUE){
    a_end <- getMRCA(qcf_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                                 "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "22", "23", "24", "25", "26", "27",
                                 "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40"))
    a_start <- qcf_tree$edge[which(qcf_tree$edge[,2] == a_end),1]
    branch_a_mono <- is.monophyletic(qcf_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                                                 "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "22", "23", "24", "25", "26", "27",
                                                 "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40"))
    b_end <- getMRCA(qcf_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                                 "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
    b_start <- qcf_tree$edge[which(qcf_tree$edge[,2] == b_end),1]
    branch_b_mono <- is.monophyletic(qcf_tree, c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58",
                                                 "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
  }
  # Extract information from those two branches
  qcf_branch_a <- node_labs[ ( a_end-Ntip(qcf_tree) ) ]
  qcf_branch_a_mono <- branch_a_mono
  qcf_branch_b <- node_labs[ ( b_end-Ntip(qcf_tree) ) ]
  qcf_branch_b_mono <- branch_b_mono
  
  # Assemble the output
  qcf_branches <- c(qcf_branch_a_mono, qcf_branch_a, qcf_branch_b_mono, qcf_branch_b)
  names(qcf_branches) <- c("qcf_branch_a_monophyletic", "qcf_branch_a_value", "qcf_branch_b_monophyletic", "qcf_branch_b_value")
  # Return the output
  return(qcf_branches)
}



qcf.clade.values <- function(clade_of_interest, qcf_tree){
  ## Function to test for monophyly of a single clade, extract the values of interest and return some summary statistics
  
  # Update clade of interest to include the "_numbers" suffix
  clade_numbers <- paste0(clade_of_interest, "_numbers")
  # Specify taxa names for each clade (both species names and converted to numbers)
  whelan2017_list <- list("Bilateria" = c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster","Daphnia_pulex"),
                          "Cnidaria" = c("Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera", 
                                         "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona","Craseo_lathetica",
                                         "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla"),
                          "Porifera" = c("Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum",
                                         "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica",
                                         "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", 
                                         "Crella_elegans", "Kirkpatrickia_variolosa"),
                          "Ctenophora" = c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus",
                                           "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                                           "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                                           "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                                           "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                                           "Ctenophora_sp_Florida_USA"),
                          "Outgroup" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"),
                          "Bilateria_numbers" = c("1", "2", "3", "4", "5", "6"),
                          "Cnidaria_numbers" = c("7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"),
                          "Porifera_numbers" = c("22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40"),
                          "Ctenophora_numbers" = c("41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59",
                                                   "60", "61", "62", "63", "64","65", "66", "67", "68", "69", "70"),
                          "Outgroup_numbers" = c("71", "72", "73", "74", "75"),
                          "BilatCnid_numbers" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"),
                          "BilatCnidPori_numbers" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                                      "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", 
                                                      "39", "40"),
                          "BilatCnidCten_numbers" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                                      "21","41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57",
                                                      "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"),
                          "CtenPori_numbers" = c("22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",
                                                 "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59",
                                                 "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"),
                          "Animals_numbers" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21",
                                                "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40",
                                                "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59",
                                                "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70"))
  # Identify start and end of branch
  clade_branch_end <- getMRCA( qcf_tree, whelan2017_list[[clade_numbers]] )
  clade_branch_start <- qcf_tree$edge[which(qcf_tree$edge[,2] == clade_branch_end), 1]
  clade_branch_mono <- is.monophyletic( qcf_tree, whelan2017_list[[clade_numbers]] )
  # Extract the node labels (i.e. qCF values) from the qCF tree
  node_labs <- qcf_tree$node.label
  # Extract information from those two branches
  if (clade_branch_mono == TRUE){
    qcf_branch_clade <- node_labs[ ( clade_branch_end - Ntip(qcf_tree) ) ]
  } else {
    qcf_branch_clade = NA
  }
  # Assemble the output
  qcf_clade <- c(clade_branch_mono, qcf_branch_clade)
  names(qcf_clade) <- paste0(clade_of_interest, c("_monophyletic", "_qcf_value"))
  # Return the output
  return(qcf_clade)
}







#### Gene concordance factors wrapper ####
gcf.wrapper <- function(alignment_path, iqtree2_path, iqtree2_model = NA, iqtree2_num_threads = "AUTO", rename.taxa.for.ms = TRUE, renamed_taxa){
  # Function to calculate the estimated and empirical gCF and return relevant gCFs
  
  ## Calculate the actual gCFs using IQ-Tree2
  actual_gcfs <- extract.input.concordance.factors(alignment_path, iqtree2_path, iqtree2_num_threads, rename.taxa.for.ms = TRUE, renamed_taxa)
  
  ## Calculate the estimated gCFs using IQ-Tree2
  estimated_gcfs <- extract.output.concordance.factors(alignment_path, iqtree2_path, iqtree2_num_threads, iqtree2_model)
  
  ## Extract information from the actual gCFS
  # Extract the gcf tables
  a_gcfs <- actual_gcfs$gcf_table 
  # Read in the branch labelled tree
  a_tree <- read.tree(file = actual_gcfs$gcf_branch_file)
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
  
  # Return collated output (both actual and estimated gCFS)
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
    tree_file <- paste0(al_dir, grep("numericTipLabels", grep("renamed", grep("branch_lengths_modified", al_files, value = T), value = T, invert = T), value = T, invert = T))
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
  gcf_files <- iqtree2.concordance.factors(alignment_path, iqtree2_path, iqtree2_num_threads, iqtree2_model)
  # Retrieve gcf from output table
  gcf_stat_file <- gcf_files[["gCF_table_file"]]
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
    species_tree_file_check <- paste0(dirname(alignment_path), "/", al_id,"_ML_tree.treefile")
    species_tree_iqfile_check <- paste0(dirname(alignment_path), "/", al_id,"_ML_tree.iqtree")
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
    gene_tree_file_check <- paste0(dirname(alignment_path), "/", al_id,"_gene_trees.treefile")
    gene_tree_iqfile_check <- paste0(dirname(alignment_path), "/", al_id,"_gene_trees.iqtree")
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
    species_tree_file <- paste0(dirname(alignment_path), "/", al_id,"_ML_tree.treefile")
    gene_tree_file <- paste0(dirname(alignment_path), "/", al_id,"_gene_trees.treefile")
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
  
  # Create new file path
  split_path <- unlist(strsplit(basename(tree_path), "\\."))
  relabelled_tree_path <- paste0(dirname(tree_path), "/",
                                 paste(head(split_path, (length(split_path) - 1)), collapse = "."), 
                                 "_relabelled", ".",
                                 tail(split_path, 1))
  # If the file path already exists, skip the rest of this section. 
  # Otherwise, open the tree and edit the tip labels
  if (file.exists(relabelled_tree_path) == FALSE){
    # Open tree and get tip labels  
    test_tree <- read.tree(tree_path)
    tree_tips <- test_tree$tip.label
    # Check if tip labels have a "t" in them - e.g. "t72"
    if (TRUE %in% grepl("t", tree_tips)){
      # Replace tree tips without "t"
      new_tree_tips <- gsub("t", "", tree_tips)
      # Replace tree tips in original tree
      test_tree$tip.label <- new_tree_tips
      # Save the relabelled tree
      write.tree(test_tree, file = relabelled_tree_path)
    } else {
      relabelled_tree_path <- tree_path
    }
  }
  # Return the relabelled tree path and use that to calculate qcfs
  return(relabelled_tree_path)
}


extract.tree.length <- function(tree){
  ## Quick function to calculate the tree length and percentage of internal branches for any tree
  
  # Identify terminal branches and internal branches
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


extract.tree.depth <- function(tree, outgroup, root.tree = TRUE){
  ## Quick function to calculate the depth of any tree
  
  # Check if outgroup is monophyletic
  check_monophyly <- is.monophyletic(phy = tree, tips = outgroup)
  # Calculate tree depth
  if (root.tree == TRUE){
    # Root tree if requested AND if outgroup is monophyletic
    # Otherwise return NaN
    if (check_monophyly == TRUE){
      # Root tree
      rooted_tree <- root(tree, outgroup)
      # Find maximum branching time for the tree
      tree_depth <- max(branching.times(rooted_tree))
    } else {
      # Outgroup is not monophyletic - return NaN
      tree_depth = NaN
    }
  } else {
    # Find maximum branching time for the tree
    tree_depth <- max(branching.times(tree))
  }
  
  # Assemble output
  op_vec <- c(as.character(check_monophyly), tree_depth)
  names(op_vec) <- c("outgroup_monophyletic", "tree_depth")
  # Return output
  return(op_vec)
}

