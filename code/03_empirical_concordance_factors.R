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
control <- list(run.cf.analyses = FALSE)


#### 2. Prepare functions, variables and packages ####
# Source functions and dataset information
source(paste0(repo_dir, "code/func_empirical_tree_estimation.R"))
source(paste0(repo_dir, "code/data_dataset_info.R"))
alignment_taxa_df <- read.table(paste0(repo_dir, "output/dataset_included_taxa.tsv"), header = T)

# Remove unneeded dataset information
rm(all_taxa, borowiec2015_list, chang2015_list, dunn2008_list, hejnol2009_list, laumer2018_list, laumer2019_list,
   models_list, moroz2014_list, nosenko2013_list, philippe2009_list, philippe2011_list, pick2010_list, ryan2013_list,
   simion2017_list, whelan2015_list, whelan2017_list)

# Open the required dataframe with dataset information
input_df <- read.csv(paste0(output_csv_dir, "cf_analysis_input_paths.csv"), stringsAsFactors = FALSE)



#### 3. Create command lines for calculating gCF and qCF ####
if (control$run.cf.analyses == TRUE){
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
# Extract files from cf_analyses folder
all_files <- list.files(cf_dir, recursive = TRUE)
# Extract only gCF files
gCF_files <- grep("gcf", all_files, value = T)
# Extract 




#### 6. Extract qCF values from key clades  ####
qCF_df_file <- paste0(output_csv_dir, "qCF_tree_files.csv")
if (file.exists(qCF_df_file)){
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
  write.csv(qcf_df, file = qCF_df_file)
} else {
  qcf_df <- read.csv(qCF_df_file, stringsAsFactors = TRUE)
}

# Iterate one row at a time
i = 1

extract.qcf.wrapper <- function(i, qcf_df, 
                                matrix_taxa = matrix_taxa, dataset_info = dataset_info, alignment_taxa_df = alignment_taxa_df){
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
  
  # Extract qCF depending on topology
  if (row$tree_topology == "CTEN_PORI"){
    qcf_extracted <- extract.CTEN_PORI.qcf(dataset = row_dataset, matrix_name = row$matrix_name, tree_file = row$qcf_tree_file,
                                           constraint_clades = constraint_clades)
  }
  
  
}

extract.CTEN_PORI.qcf <- function(dataset, matrix_name, tree_file, constraint_clades){
  # Function to extract qCF for CTEN_PORI tree
  
  
  
  
}
