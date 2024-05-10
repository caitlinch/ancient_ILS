## caitlinch/ancient_ILS/code/03_empirical_concordance_factors.R
# This script estimates gene, and quartet concordance factors for empirical gene trees and constrained hypothesis trees.
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## Specify parameters:
# location                    <- Where the script is being run
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# output_csv_dir              <- Location of csv input/output/results files
# cf_dir                      <- Location of concordance factor output files
# iqtree2_server              <- Location of IQ-Tree2 executable
# iqtree2_num_threads         <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
# astral_server               <- Location of ASTRAL executable

repo_dir              <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_csv_dir        <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/07_output_files/"
cf_dir                <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/06_cf_analyses/"

# Server file paths to run qCF/gCF analyses
iqtree2_server        <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2"
iqtree2_num_threads   <- 30
astral_server         <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/astral/Astral/astral.5.7.8.jar"

# Control commands
control <- list(run.cf.analyses = FALSE,
                extract.gcf = FALSE,
                extract.qcf = FALSE)



#### 2. Prepare functions, variables and packages ####
# Read in alignment_taxa_df, which contains the list of taxa in the ML tree for each dataset
alignment_taxa_df <- read.table(paste0(repo_dir, "output/dataset_included_taxa.tsv"), header = T)

# Source functions and dataset information
source(paste0(repo_dir, "code/func_concordance_factors.R"))
source(paste0(repo_dir, "code/data_dataset_info.R"))

# Remove unneeded dataset information
rm(all_taxa, all_models, models_list, borowiec2015_list, chang2015_list, dunn2008_list, hejnol2009_list, 
   laumer2018_list, laumer2019_list, moroz2014_list, nosenko2013_list, philippe2009_list, philippe2011_list,
   pick2010_list, ryan2013_list, simion2017_list, whelan2015_list, whelan2017_list)



#### 3. Create command lines for calculating gCF and qCF ####
if (control$run.cf.analyses == TRUE){
  # Open the required dataframe with dataset information
  input_df <- read.csv(paste0(output_csv_dir, "cf_analysis_input_paths.csv"), stringsAsFactors = FALSE)
  
  # Add C60 gCF commands
  # $ iqtree2 -t concat.treefile --gcf loci.treefile --prefix concord
  input_df$c60_cten_gcf_command <- paste0(iqtree2_server, " -te ", input_df$C60_CTEN_tree,
                                          " --gcf ", input_df$c60_gene_trees, 
                                          " -pre ", input_df$c60_gcf_output_dir, input_df$dataset_id, ".C60.CTEN.gcf",
                                          " -nt ", iqtree2_num_threads)
  input_df$c60_cten_gcf_prefix <- paste0(input_df$c60_gcf_output_dir, input_df$dataset_id, ".C60.CTEN.gcf")
  input_df$c60_pori_gcf_command <- paste0(iqtree2_server, " -te ", input_df$C60_PORI_tree,
                                          " --gcf ", input_df$c60_gene_trees, 
                                          " -pre ", input_df$c60_gcf_output_dir, input_df$dataset_id, ".C60.PORI.gcf",
                                          " -nt ", iqtree2_num_threads)
  input_df$c60_pori_gcf_prefix <- paste0(input_df$c60_gcf_output_dir, input_df$dataset_id, ".C60.PORI.gcf")
  input_df$c60_cten_pori_gcf_command <- paste0(iqtree2_server, " -te ", input_df$C60_CTEN_PORI_tree,
                                               " --gcf ", input_df$c60_gene_trees, 
                                               " -pre ", input_df$c60_gcf_output_dir, input_df$dataset_id, ".C60.CTEN_PORI.gcf",
                                               " -nt ", iqtree2_num_threads)
  input_df$c60_cten_pori_gcf_prefix <- paste0(input_df$c60_gcf_output_dir, input_df$dataset_id, ".C60.CTEN_PORI.gcf")
  
  # Add C60 qCF commands
  # $ java -jar astral.5.7.8.jar -q test_data/1kp.tre -i test_data/1KP-genetrees.tre -t 2 -o test_data/1kp-scored-t2.tre
  input_df$c60_cten_qcf_command <- paste0("java -jar ", astral_server, 
                                          " -q ", input_df$C60_CTEN_tree, 
                                          " -i ", input_df$c60_gene_trees, " -t 2 ",
                                          " -o ", input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN.qcf.tre", 
                                          " 2> ", input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN.qcf.log")
  input_df$c60_CTEN_qcf_tree <- paste0(input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN.qcf.tre")
  input_df$c60_CTEN_qcf_log <- paste0(input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN.qcf.log")
  input_df$c60_pori_qcf_command <- paste0("java -jar ", astral_server, 
                                          " -q ", input_df$C60_PORI_tree, 
                                          " -i ", input_df$c60_gene_trees, " -t 2 ",
                                          " -o ", input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.PORI.qcf.tre", 
                                          " 2> ", input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.PORI.qcf.log")
  input_df$c60_PORI_qcf_tree <- paste0(input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN.qcf.tre")
  input_df$c60_PORI_qcf_log <- paste0(input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN.qcf.log")
  input_df$c60_cten_pori_qcf_command <- paste0("java -jar ", astral_server, 
                                               " -q ", input_df$C60_CTEN_PORI_tree, 
                                               " -i ", input_df$c60_gene_trees, " -t 2 ",
                                               " -o ", input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN_PORI.qcf.tre", 
                                               " 2> ", input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN_PORI.qcf.log")
  input_df$c60_CTEN_PORI_qcf_tree <- paste0(input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN_PORI.qcf.tre")
  input_df$c60_CTEN_PORI_qcf_log <- paste0(input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN_PORI.qcf.log")
  
  # Add Partitioned gCF commands
  # $ iqtree2 -t concat.treefile --gcf loci.treefile --prefix concord
  input_df$partition_cten_gcf_command <- paste0(iqtree2_server, " -te ", input_df$partition_CTEN_tree,
                                                " --gcf ", input_df$mfp_gene_trees, 
                                                " -pre ", input_df$partition_gcf_output_dir, input_df$dataset_id, ".Partition.CTEN.gcf",
                                                " -nt ", iqtree2_num_threads)
  input_df$partition_cten_gcf_prefix <- paste0(input_df$partition_gcf_output_dir, input_df$dataset_id, ".Partition.CTEN.gcf")
  input_df$partition_pori_gcf_command <- paste0(iqtree2_server, " -te ", input_df$partition_PORI_tree,
                                                " --gcf ", input_df$mfp_gene_trees, 
                                                " -pre ", input_df$partition_gcf_output_dir, input_df$dataset_id, ".Partition.PORI.gcf",
                                                " -nt ", iqtree2_num_threads)
  input_df$partition_pori_gcf_prefix <- paste0(input_df$partition_gcf_output_dir, input_df$dataset_id, ".Partition.PORI.gcf")
  input_df$partition_cten_pori_gcf_command <- paste0(iqtree2_server, " -te ", input_df$partition_CTEN_PORI_tree,
                                                     " --gcf ", input_df$mfp_gene_trees, 
                                                     " -pre ", input_df$partition_gcf_output_dir, input_df$dataset_id, ".Partition.CTEN_PORI.gcf",
                                                     " -nt ", iqtree2_num_threads)
  input_df$partition_pori_gcf_prefix <- paste0(input_df$partition_gcf_output_dir, input_df$dataset_id, ".Partition.CTEN_PORI.gcf")
  
  # Add Partitioned qCF commands
  # $ java -jar astral.5.7.8.jar -q test_data/1kp.tre -i test_data/1KP-genetrees.tre -t 2 -o test_data/1kp-scored-t2.tre
  input_df$partition_cten_qcf_command <- paste0("java -jar ", astral_server, 
                                                " -q ", input_df$partition_CTEN_tree, 
                                                " -i ", input_df$mfp_gene_trees, " -t 2 ",
                                                " -o ", input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN.qcf.tre", 
                                                " 2> ", input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN.qcf.log")
  input_df$partition_CTEN_qcf_tree <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN.qcf.tre")
  input_df$partition_CTEN_qcf_log <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN.qcf.log")
  input_df$partition_pori_qcf_command <- paste0("java -jar ", astral_server, 
                                                " -q ", input_df$partition_PORI_tree, 
                                                " -i ", input_df$mfp_gene_trees, " -t 2 ",
                                                " -o ", input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.PORI.qcf.tre", 
                                                " 2> ", input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.PORI.qcf.log")
  input_df$partition_PORI_qcf_tree <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.PORI.qcf.tre")
  input_df$partition_PORI_qcf_log <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.PORI.qcf.log")
  input_df$partition_cten_pori_qcf_command <- paste0("java -jar ", astral_server, 
                                                     " -q ", input_df$partition_CTEN_PORI_tree, 
                                                     " -i ", input_df$mfp_gene_trees, " -t 2 ",
                                                     " -o ", input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN_PORI.qcf.tre", 
                                                     " 2> ", input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN_PORI.qcf.log")
  input_df$partition_CTEN_PORI_qcf_tree <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN_PORI.qcf.tre")
  input_df$partition_CTEN_PORI_qcf_log <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN_PORI.qcf.log")
  
  # Write out command lines to text files
  write(c(input_df$c60_cten_gcf_command, input_df$c60_pori_gcf_command, input_df$c60_cten_pori_gcf_command), 
        file = paste0(output_csv_dir, "c60_gcf_commands.txt"))
  write(c(input_df$c60_cten_qcf_command, input_df$c60_pori_qcf_command, input_df$c60_cten_pori_qcf_command), 
        file = paste0(output_csv_dir, "c60_qcf_commands.txt"))
  write(c(input_df$partition_cten_gcf_command, input_df$partition_pori_gcf_command, input_df$partition_cten_pori_gcf_command), 
        file = paste0(output_csv_dir, "partition_gcf_commands.txt"))
  write(c(input_df$partition_cten_qcf_command, input_df$partition_pori_qcf_command, input_df$partition_cten_pori_qcf_command), 
        file = paste0(output_csv_dir, "partition_qcf_commands.txt"))
}



#### 4. Extract gCF values from key clades  ####
if (control$extract.qcf == TRUE){
  ## Get gCF output
  # Specify gCF parameters df
  gcf_df_file <- paste0(output_csv_dir, "gCF_tree_files.csv")
  # Open gCF parameters dataframe
  if (file.exists(gcf_df_file)){
    # Extract files from cf_analyses folder
    all_files <- list.files(cf_dir, recursive = TRUE)
    # Extract only gCF files
    gcf_files <- grep("gcf", all_files, value = T)
    # Extract cf.branch and cf.stat files
    cf_stat_files <- grep("cf.stat", gcf_files, value = T)
    cf_branch_files <- grep("cf.branch", gcf_files, value = T)
    # Extract gCF files into a dataframe
    gcf_df <- data.frame(gcf_stat_files = paste0(cf_dir, cf_stat_files),
                         gcf_branch_files = paste0(cf_dir, cf_branch_files))
    # Add other required columns to dataframe
    gcf_df$id             <- gsub(".gcf.cf.stat", "", basename(gcf_df$gcf_stat_files))
    split_id              <- strsplit(gcf_df$id, "\\.")
    gcf_df$dataset        <- unlist(lapply(split_id, function(x){x[[1]]}))
    gcf_df$matrix_name    <- unlist(lapply(split_id, function(x){x[[2]]}))
    gcf_df$dataset_id     <- paste0(gcf_df$dataset, ".", gcf_df$matrix_name)
    gcf_df$model          <- unlist(lapply(split_id, function(x){x[[3]]}))
    gcf_df$tree_topology  <- unlist(lapply(split_id, function(x){x[[4]]}))
    # Rearrange order of columns
    gcf_df <- gcf_df[, c("id", "dataset", "matrix_name", "dataset_id", "model", "tree_topology", "gcf_branch_files", "gcf_stat_files")]
    # Write qCF_df
    write.csv(gcf_df, file = gcf_df_file, row.names = FALSE)
  } else {
    gcf_df <- read.csv(gcf_df_file, stringsAsFactors = TRUE)
  }
  
  ## Extract gCF values
  gcf_output_list <- lapply(1:nrow(gcf_df), extract.gcf.wrapper, gcf_df = gcf_df, 
                            matrix_taxa = matrix_taxa, all_datasets = all_datasets, 
                            alignment_taxa_df = alignment_taxa_df)
  
}

# test 11, 12, 10
i = 10

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
  # Extract branch_id from the node.label 
  met_branch_id <- as.numeric(g_edge_table[which(g_edge_table$parent == met_pn & g_edge_table$child == met_cn), ]$par.name)
  # Extract the row from the stat table for this branch_id
  met_values <- g_table[which(g_table$ID == met_branch_id), ]
  names(met_values) <- paste0("MET_", names(met_values))
  
  ## Key branch (leading to ALL OTHER ANIMALS aka PLAC+CNID+BILAT)
  # Identify tips in this group - do not include PLAC when identifying MRCA, 
  #     as sometimes PLAC placement is sister to PORI which will result in 
  #     extracting the wrong branch
  if (topology == "CTEN" | topology == "PORI"){
    if (topology == "CTEN"){
      # Extract taxa
      key_taxa <- c(constraint_clades$Porifera, constraint_clades$Cnidaria, constraint_clades$Bilateria)
    } else if (topology == "PORI"){
      # Extract taxa
      key_taxa <- c(constraint_clades$Ctenophora, constraint_clades$Cnidaria, constraint_clades$Bilateria)
    }
    # Extract MRCA
    key_cn <- getMRCA(g_rooted, key_taxa) # child node
    key_pn <- g_rooted$edge[which(g_rooted$edge[,2] == key_cn), 1] # parent node
    # Extract branch_id from the node.label 
    key_branch_id <- as.numeric(g_edge_table[which(g_edge_table$parent == key_pn & g_edge_table$child == key_cn), ]$par.name)
    # Extract the row from the stat table for this branch_id
    key_values <- g_table[which(g_table$ID == key_branch_id), ]
    names(key_values) <- paste0("KEY_", names(key_values))
  } else if (topology == "CTEN_PORI" | topology == "CTEN.PORI"){
    # Extract taxa
    key_taxa <- c(constraint_clades$Ctenophora, constraint_clades$Porifera)
    # Extract MRCA
    key_cn <- getMRCA(g_rooted, key_taxa) # child node
    key_pn <- g_rooted$edge[which(g_rooted$edge[,2] == key_cn), 1] # parent node
    # Extract branch_id from the node.label 
    key_branch_id <- as.numeric(g_edge_table[which(g_edge_table$parent == key_pn & g_edge_table$child == key_cn), ]$chi.name)
    # Extract the row from the stat table for this branch_id
    key_values <- g_table[which(g_table$ID == key_branch_id), ]
    names(key_values) <- paste0("KEY_", names(key_values))
  }
  
  ## CTEN (leading to CTEN branch)
  # Identify tips in this group
  cten_taxa <- c(constraint_clades$Ctenophora)
  if (length(cten_taxa) > 1){
    # Extract MRCA
    cten_cn <- getMRCA(g_rooted, cten_taxa) # child node
    cten_pn <- g_rooted$edge[which(g_rooted$edge[,2] == cten_cn), 1] # parent node
    # Extract branch_id from the node.label 
    cten_branch_id <- as.numeric(g_edge_table[which(g_edge_table$parent == cten_pn & g_edge_table$child == cten_cn), ]$chi.name)
    # Extract the row from the stat table for this branch_id
    cten_values <- g_table[which(g_table$ID == cten_branch_id), ]
    names(cten_values) <- paste0("CTEN_", names(cten_values))
  } else {
    # Assign NA if only 1 tip
    cten_values <- rep(NA, 11)
    names(cten_values) <- paste0("CTEN_", names(g_table))
  }
  
  ## PORI (leading to PORI branch)
  # Identify tips in this group
  pori_taxa <- c(constraint_clades$Porifera)
  if (length(pori_taxa) > 1){
    # Extract MRCA
    pori_cn <- getMRCA(g_rooted, pori_taxa) # child node
    pori_pn <- g_rooted$edge[which(g_rooted$edge[,2] == pori_cn), 1] # parent node
    # Extract branch_id from the node.label 
    pori_branch_id <- as.numeric(g_edge_table[which(g_edge_table$parent == pori_pn & g_edge_table$child == pori_cn), ]$chi.name)
    # Extract the row from the stat table for this branch_id
    pori_values <- g_table[which(g_table$ID == pori_branch_id), ]
    names(pori_values) <- paste0("PORI_", names(pori_values))
  } else {
    # Assign NA if only 1 tip
    pori_values <- rep(NA, 11)
    names(pori_values) <- paste0("PORI_", names(g_table))
  }
  
  # Assemble output vector
  gcf_output <- as.character(c(met_values, key_values, cten_values, pori_values))
  names(gcf_output) <- c(names(met_values), names(key_values), names(cten_values), names(pori_values))
  return(gcf_output)
}





#### 6. Extract qCF values from key clades  ####
if (control$extract.qcf == TRUE){
  ## Get qCF output
  # Specify qCF parameters df
  qcf_df_file <- paste0(output_csv_dir, "qCF_tree_files.csv")
  # Open qCF parameters dataframe
  if (file.exists(qcf_df_file)){
    # Extract files from cf_analyses folder
    all_files <- list.files(cf_dir, recursive = TRUE)
    # Extract only gCF files
    qcf_tree_files <- grep("\\.tre", grep("qcf", all_files, value = T), value = T)
    # Extract qCF files into a dataframe
    qcf_df <- data.frame(qcf_tree_file = paste0(cf_dir, qcf_tree_files))
    # Add other required columns to dataframe
    qcf_df$id             <- gsub(".qcf.tre", "", basename(qcf_df$qcf_tree_file))
    split_id              <- strsplit(qcf_df$id, "\\.")
    qcf_df$dataset        <- unlist(lapply(split_id, function(x){x[[1]]}))
    qcf_df$matrix_name    <- unlist(lapply(split_id, function(x){x[[2]]}))
    qcf_df$dataset_id     <- paste0(qcf_df$dataset, ".", qcf_df$matrix_name)
    qcf_df$model          <- unlist(lapply(split_id, function(x){x[[3]]}))
    qcf_df$tree_topology  <- unlist(lapply(split_id, function(x){x[[4]]}))
    # Rearrange order of columns
    qcf_df <- qcf_df[, c("id", "dataset", "matrix_name", "dataset_id", "model", "tree_topology", "qcf_tree_file")]
    # Write qCF_df
    write.csv(qcf_df, file = qcf_df_file, row.names = FALSE)
  } else {
    qcf_df <- read.csv(qcf_df_file, stringsAsFactors = TRUE)
  }
  
  ## Extract qCF values
  qcf_output_list <- lapply(1:nrow(qcf_df), extract.qcf.wrapper, qcf_df = qcf_df, 
                            matrix_taxa = matrix_taxa, all_datasets = all_datasets, 
                            alignment_taxa_df = alignment_taxa_df)
  
  ## Format and output qCF dataframe
  qcf_output_df <- as.data.frame(do.call(rbind, qcf_output_list), stringsAsFactors = FALSE)
  # Break down MET node values column
  qcf_output_df$MET_quartet_support <- gsub("q1=", "", unlist(lapply(strsplit(qcf_output_df$MET_node_value, ";"), function(x){x[[1]]})))
  qcf_output_df$MET_quartet_tree_freq <- gsub("f1=", "", unlist(lapply(strsplit(qcf_output_df$MET_node_value, ";"), function(x){x[[4]]})))
  qcf_output_df$MET_lpp <- gsub("pp1=", "", unlist(lapply(strsplit(qcf_output_df$MET_node_value, ";"), function(x){x[[7]]})))
  qcf_output_df$MET_QC <- gsub("QC=", "", unlist(lapply(strsplit(qcf_output_df$MET_node_value, ";"), function(x){x[[10]]})))
  qcf_output_df$MET_EN <- gsub("EN=", "", unlist(lapply(strsplit(qcf_output_df$MET_node_value, ";"), function(x){x[[11]]})))
  # Break down KEY node values column
  qcf_output_df$KEY_quartet_support <- gsub("q1=", "", unlist(lapply(strsplit(qcf_output_df$KEY_node_value, ";"), function(x){x[[1]]})))
  qcf_output_df$KEY_quartet_tree_freq <- gsub("f1=", "", unlist(lapply(strsplit(qcf_output_df$KEY_node_value, ";"), function(x){x[[4]]})))
  qcf_output_df$KEY_lpp <- gsub("pp1=", "", unlist(lapply(strsplit(qcf_output_df$KEY_node_value, ";"), function(x){x[[7]]})))
  qcf_output_df$KEY_QC <- gsub("QC=", "", unlist(lapply(strsplit(qcf_output_df$KEY_node_value, ";"), function(x){x[[10]]})))
  qcf_output_df$KEY_EN <- gsub("EN=", "", unlist(lapply(strsplit(qcf_output_df$KEY_node_value, ";"), function(x){x[[11]]})))
  # Break down CTEN node values column
  qcf_output_df$CTEN_quartet_support <- gsub("q1=", "", unlist(lapply(strsplit(qcf_output_df$CTEN_node_value, ";"), function(x){x[[1]]})))
  qcf_output_df$CTEN_quartet_tree_freq <- gsub("f1=", "", unlist(lapply(strsplit(qcf_output_df$CTEN_node_value, ";"), function(x){ifelse((length(x)>1), x[[4]], NA)})))
  qcf_output_df$CTEN_lpp <- gsub("pp1=", "", unlist(lapply(strsplit(qcf_output_df$CTEN_node_value, ";"), function(x){ifelse((length(x)>1), x[[7]], NA)})))
  qcf_output_df$CTEN_QC <- gsub("QC=", "", unlist(lapply(strsplit(qcf_output_df$CTEN_node_value, ";"), function(x){ifelse((length(x)>1), x[[10]], NA)})))
  qcf_output_df$CTEN_EN <- gsub("EN=", "", unlist(lapply(strsplit(qcf_output_df$CTEN_node_value, ";"), function(x){ifelse((length(x)>1), x[[11]], NA)})))
  # Break down PORI node values column
  qcf_output_df$PORI_quartet_support <- gsub("q1=", "", unlist(lapply(strsplit(qcf_output_df$PORI_node_value, ";"), function(x){x[[1]]})))
  qcf_output_df$PORI_quartet_tree_freq <- gsub("f1=", "", unlist(lapply(strsplit(qcf_output_df$PORI_node_value, ";"), function(x){ifelse((length(x)>1), x[[4]], NA)})))
  qcf_output_df$PORI_lpp <- gsub("pp1=", "", unlist(lapply(strsplit(qcf_output_df$PORI_node_value, ";"), function(x){ifelse((length(x)>1), x[[7]], NA)})))
  qcf_output_df$PORI_QC <- gsub("QC=", "", unlist(lapply(strsplit(qcf_output_df$PORI_node_value, ";"), function(x){ifelse((length(x)>1), x[[10]], NA)})))
  qcf_output_df$PORI_EN <- gsub("EN=", "", unlist(lapply(strsplit(qcf_output_df$PORI_node_value, ";"), function(x){ifelse((length(x)>1), x[[11]], NA)})))
  # Cbind to the parameters dataframe
  qcf_collated_df <- cbind(qcf_df, qcf_output_df)
  
  ## Save output
  qCF_collated_df_file <- paste0(output_csv_dir, "qCF_values.csv")
  write.csv(qcf_collated_df, file = qCF_collated_df_file, row.names = FALSE)
}






