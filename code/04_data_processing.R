## caitlinch/ancient_ILS/code/04_data_processing.R
# This script estimates maximum likelihood trees under different models of substitution for 14 empirical data sets
# Caitlin Cherryh 2023



#### 1. Input parameters ####
## Specify parameters:
# location                    <- To run locally: location = "Local" 
#                                   or to run on server: location = "Dayhoff"
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# output_dir                  <- Directory for saving results files (csv and text)
# hypothesis_tree_dir         <- Location of constrained ML and ASTRAL trees
# hypothesis_gene_dir         <- Location of constrained ML gene trees

location = "Local"
repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_dir              <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/05_output_files/"
hypothesis_tree_dir     <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_02_empirical_tree_estimation/03_hypothesis_trees/"
hypothesis_gene_dir     <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_03_empirical_genes_constrained_trees/"


#### 2. Prepare libraries and functions ####
# Open packages
library(stringr)

# Open functions
source(paste0(repo_dir, "code/func_data_processing.R"))

# Extract list of all output files
all_output_files <- list.files(output_dir)



#### 3. Prepare tree and constrained tree output for figures ####
# Open results file with constrained tree results
raw_tree_df <- read.csv(paste0(output_dir, grep("ConstrainedTreeResults", all_output_files, value = T)), stringsAsFactors = FALSE)
# Add year column
raw_tree_df$year <- as.numeric(str_extract(raw_tree_df$dataset, "(\\d)+"))
# Add column with model code
raw_tree_df$ML_model_code <- unlist(lapply(raw_tree_df$unconstrained_tree_best_model, function(x){strsplit(x, "\\+")[[1]][1]}))
# Trim unnecessary columns
tree_df <- raw_tree_df[, c("dataset", "matrix_name", "dataset_id", "gene_name", "gene_id",
                           "ML_model_code", "unconstrained_tree_best_model", "unconstrained_tree_alisim_model",
                           "unconstrained_tree_logl", "unconstrained_tree_unconstrained_logl", "unconstrained_tree_numFreeParams",
                           "unconstrained_tree_BIC", "unconstrained_tree_length", "unconstrained_tree_sumInternalBL",
                           "CTEN_LogL", "CTEN_Unconstrained_LogL", "CTEN_NumFreeParams",
                           "CTEN_BIC", "CTEN_TreeLength", "CTEN_SumInternalBranchLengths",
                           "PORI_LogL", "PORI_Unconstrained_LogL", "PORI_NumFreeParams",
                           "PORI_BIC", "PORI_TreeLength", "PORI_SumInternalBranchLengths",
                           "CTEN_PORI_LogL", "CTEN_PORI_Unconstrained_LogL", "CTEN_PORI_NumFreeParams",
                           "CTEN_PORI_BIC", "CTEN_PORI_TreeLength", "CTEN_PORI_SumInternalBranchLengths")]
# Rename columns for plotting
names(tree_df) <- c("dataset", "matrix_name", "dataset_id", "gene_name", "gene_id",
                    "ML_model_code", "ML_best_model", "ML_alisim_model",
                    "ML_logl", "ML_unconstrained_logl", "ML_numFreeParams",
                    "ML_BIC", "ML_tree_length", "ML_sumInternalBL",
                    "CTEN_LogL", "CTEN_Unconstrained_LogL", "CTEN_NumFreeParams",
                    "CTEN_BIC", "CTEN_TreeLength", "CTEN_SumInternalBL",
                    "PORI_LogL", "PORI_Unconstrained_LogL", "PORI_NumFreeParams",
                    "PORI_BIC", "PORI_TreeLength", "PORI_SumInternalBL",
                    "CTEN_PORI_LogL", "CTEN_PORI_Unconstrained_LogL", "CTEN_PORI_NumFreeParams",
                    "CTEN_PORI_BIC", "CTEN_PORI_TreeLength", "CTEN_PORI_SumInternalBL")
# Write the output
tree_output_file <- paste0(output_dir, "results_tree_likelihood.csv")
write.csv(tree_df, file = tree_output_file, row.names = FALSE)



#### 4. Prepare scf output for figures ####




#### 5. Prepare AU test output for figures ####
# Open raw AU test output
raw_au_df <- read.csv(paste0(output_dir, grep("TreeComparisonResults", all_output_files, value = T)), stringsAsFactors = FALSE)
# Extract AU results per ID
au_test_df <- as.data.frame(do.call(rbind, lapply(unique(raw_au_df$ID), summarise.AU.test.results, raw_au_df)))
# Remove columns consisting only of NA
au_test_df <- Filter(function(x)!all(is.na(x)), au_test_df)
# Sort output by year
au_test_df <- au_test_df[order(au_test_df$year, au_test_df$dataset, au_test_df$matrix, au_test_df$gene_id),]
# Reorder columns
au_test_df <- au_test_df[, c("gene_id", "dataset", "matrix", "gene", "year",
                             "topology_test", "tree_1", "tree_2", "tree_3")]
# Write the output
au_test_file <- paste0(output_dir, "results_gene_AU_test.csv")
write.csv(au_test_df, file = au_test_file, row.names = FALSE)



#### 5. Prepare ELW output for figures ####
# Open raw AU test output
raw_au_df <- read.csv(paste0(output_dir, grep("TreeComparisonResults", all_output_files, value = T)), stringsAsFactors = FALSE)
# Extract AU results per ID
elw_df <- as.data.frame(do.call(rbind, lapply(unique(raw_au_df$ID), summarise.eLW, raw_au_df)))
# Remove columns consisting only of NA
elw_df <- Filter(function(x)!all(is.na(x)), elw_df)
# Sort output by year
elw_df <- elw_df[order(elw_df$year, elw_df$dataset, elw_df$matrix, elw_df$gene_id),]
# Reorder columns
elw_df <- elw_df[, c("gene_id", "dataset", "matrix", "gene", "year",
                     "topology_test", "tree_1", "tree_2", "tree_3")]
# Write the output
elw_file <- paste0(output_dir, "results_gene_elw.csv")
write.csv(elw_df, file = elw_file, row.names = FALSE)



#### 6. Collate hypothesis trees ####
# List all files in the directory
all_ht_files <- list.files(hypothesis_tree_dir)
all_ht_files <- grep("00_", all_ht_files, value = TRUE, invert = TRUE)
# Create ID for each set of trees
dataset_id <- unique(unlist(lapply(strsplit(all_ht_files, "\\."), function(x){paste0(x[1], ".", x[2])})))
dataset_id <- grep("Simion|Hejnol", dataset_id, value = TRUE, invert = TRUE)
# Call each ID and collate the tree files
for (id in dataset_id){
  # Extract output files for this dataset
  id_files <- grep(id, all_ht_files, value = T)
  # Collate ML files
  ml_tree_files <- grep("\\.treefile", grep("_ML_", id_files, value = T), value = T)
  ml_tree_files_ordered <- paste0(hypothesis_tree_dir, c(grep("\\.CTEN_ML", ml_tree_files, value = T), 
                                                         grep("\\.PORI_ML", ml_tree_files, value = T), 
                                                         grep("\\.CTEN_PORI_ML", ml_tree_files, value = T)))
  ml_trees <- c(unlist(lapply(ml_tree_files_ordered, readLines)), "")
  collated_ml_tree_file <- paste0(hypothesis_tree_dir, id, ".ModelFinder.ML.collated_hypothesis_trees.treefile")
  write(ml_trees, file = collated_ml_tree_file)
  # Collate ASTRAL files
  astral_tree_files <- grep("\\.tre", grep("_ASTRAL_", id_files, value = T), value = T)
  astral_tree_files_ordered <- paste0(hypothesis_tree_dir, c(grep("\\.CTEN_ASTRAL", astral_tree_files, value = T), 
                                                             grep("\\.PORI_ASTRAL", astral_tree_files, value = T), 
                                                             grep("\\.CTEN_PORI_ASTRAL", astral_tree_files, value = T)))  
  astral_trees <- c(unlist(lapply(astral_tree_files_ordered, readLines)), "")
  collated_astral_tree_file <- paste0(hypothesis_tree_dir, id, ".ModelFinder.ASTRAL.collated_hypothesis_trees.treefile")
  write(astral_trees, file = collated_astral_tree_file)
  
}



#### 7. Collate constrained gene trees ####
# List all files within the hypothesis_gene_dir
all_datasets <- list.files(hypothesis_gene_dir)
all_datasets <- grep("_hypothesis_trees|all_collated_constrained_gene_trees", all_datasets, value = TRUE, invert = TRUE)
all_datasets <- grep("ModelFinder.CTEN.gene_trees|ModelFinder.PORI.gene_trees|ModelFinder.CTEN_PORI.gene_trees", all_datasets, value = TRUE, invert = TRUE)
# Collate the unconstrained gene trees from each dataset
for (d in all_datasets){
  # List all files in the d folder
  d_dir       <- paste0(hypothesis_gene_dir, d, "/")
  d_all_files <- list.files(d_dir)
  # Extract list of unconstrained gene trees
  d_MFP_tree_files      <- grep("\\.MFP\\.treefile", d_all_files, value = TRUE)
  # Open list of trees
  d_MFP_trees      <- c(unlist(lapply(paste0(d_dir, d_MFP_tree_files), readLines)), "")
  # Create output files for collated constrained gene trees
  d_MFP_collated_file      <- paste0(hypothesis_gene_dir, d, ".ModelFinder.MFP.gene_trees.treefile")
  # Save list of trees
  write(d_MFP_trees, file = d_MFP_collated_file)
}
# Collate the constrained gene trees from each dataset
for (d in all_datasets){
  # List all files in the d folder
  d_dir       <- paste0(hypothesis_gene_dir, d, "/")
  d_all_files <- list.files(d_dir)
  # Extract list of trees for each topology
  d_CTEN_tree_files      <- grep("\\.CTEN_tree.treefile", d_all_files, value = TRUE)
  d_PORI_tree_files      <- grep("\\.PORI_tree.treefile", d_all_files, value = TRUE)
  d_CTEN_PORI_tree_files <- grep("\\.CTEN_PORI_tree.treefile", d_all_files, value = TRUE)
  # Open list of trees
  d_CTEN_trees      <- c(unlist(lapply(paste0(d_dir, d_CTEN_tree_files), readLines)), "")
  d_PORI_trees      <- c(unlist(lapply(paste0(d_dir, d_PORI_tree_files), readLines)), "")
  d_CTEN_PORI_trees <- c(unlist(lapply(paste0(d_dir, d_CTEN_PORI_tree_files), readLines)), "")
  # Create output files for collated constrained gene trees
  d_CTEN_collated_file      <- paste0(hypothesis_gene_dir, d, ".ModelFinder.CTEN.gene_trees.treefile")
  d_PORI_collated_file      <- paste0(hypothesis_gene_dir, d, ".ModelFinder.PORI.gene_trees.treefile")
  d_CTEN_PORI_collated_file <- paste0(hypothesis_gene_dir, d, ".ModelFinder.CTEN_PORI.gene_trees.treefile")
  # Save list of trees
  write(d_CTEN_trees, file = d_CTEN_collated_file)
  write(d_PORI_trees, file = d_PORI_collated_file)
  write(d_CTEN_PORI_trees, file = d_CTEN_PORI_collated_file)
}



#### 8. Create dataframe for constrained concordance factors analysis ####
if (location == "Dayhoff"){
  # All Dayhoff paths
  repo_dir              <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  alignment_dir         <- paste0(repo_dir, "data_all/")
  constraint_tree_dir   <- paste0(repo_dir, "output/constraint_trees/")
  hypothesis_tree_dir   <- paste0(repo_dir, "output/hypothesis_trees/")
  collated_genes_dir    <- paste0(repo_dir, "cf_constrained/genes_collated/")
  output_csv_dir        <- paste0(repo_dir, "cf_constrained/")
  gcf_dir               <- paste0(repo_dir, "cf_constrained/gcf/")
  scf_dir               <- paste0(repo_dir, "cf_constrained/scf/")
  qcf_dir               <- paste0(repo_dir, "cf_constrained/qcf/")
  # Add executable program paths
  iqtree2 <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2"
  ASTRAL  <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/astral/Astral/astral.5.7.8.jar"
  # Extract alignment paths
  remove_als                  <- "Simion2017.supermatrix_97sp_401632pos_1719genes|Hejnol2009.Hejnol_etal_2009_FixedNames|Laumer2019.nonbilateria_MARE_BMGE"
  remove_als_partition        <- "Simion2017|Hejnol2009|Laumer2019"
  all_dataset_files           <- list.files(alignment_dir)
  alignment_files             <- paste0(alignment_dir, grep(remove_als, grep("aa.alignment", all_dataset_files, value = T), value = T, invert = T))
  modelfinder_partition_files <- paste0(alignment_dir, grep(remove_als_partition, grep("\\.ModelFinder.partition.nex", all_dataset_files, value = T), value = T, invert = T))
  original_partition_files    <- paste0(alignment_dir, grep(remove_als_partition, grep("_partitions_formatted.nex", all_dataset_files, value = T), value = T, invert = T))
  # Extract constraint trees
  all_constraint_trees        <- grep(remove_als, list.files(constraint_tree_dir), value = T, invert = T)
  CTEN_constraint_trees       <- paste0(constraint_tree_dir, grep("constraint_tree_1.nex", all_constraint_trees, value = T))
  PORI_constraint_trees       <- paste0(constraint_tree_dir, grep("constraint_tree_2.nex", all_constraint_trees, value = T))
  CTEN_PORI_constraint_trees  <- paste0(constraint_tree_dir, grep("constraint_tree_3.nex", all_constraint_trees, value = T))
  # Extract gene trees
  all_gene_collated_trees   <- grep(remove_als, list.files(collated_genes_dir), value = T, invert = T)
  CTEN_gene_trees           <- paste0(collated_genes_dir, grep("\\.ModelFinder\\.CTEN\\.gene_trees\\.treefile", all_gene_collated_trees, value = T))
  PORI_gene_trees           <- paste0(collated_genes_dir, grep("\\.ModelFinder\\.PORI\\.gene_trees\\.treefile", all_gene_collated_trees, value = T))
  CTEN_PORI_gene_trees      <- paste0(collated_genes_dir, grep("\\.ModelFinder\\.CTEN_PORI\\.gene_trees\\.treefile", all_gene_collated_trees, value = T))
  MFP_gene_trees            <- paste0(collated_genes_dir, grep("\\.ModelFinder\\.MFP\\.gene_trees\\.treefile", all_gene_collated_trees, value = T))
  # Extract information about the datasets
  dataset     <- unlist(lapply(strsplit(basename(alignment_files), "\\."), function(x){paste0(x[[1]])}))
  matrix      <- unlist(lapply(strsplit(basename(alignment_files), "\\."), function(x){paste0(x[[2]])}))
  dataset_id  <- unlist(lapply(strsplit(basename(alignment_files), "\\."), function(x){paste0(x[[1]], ".", x[[2]])}))
  # Extract hypothesis tree paths
  all_h_trees                 <- grep(remove_als_partition, list.files(hypothesis_tree_dir), value = T, invert = T)
  CTEN_ML_tree_path           <- paste0(hypothesis_tree_dir, grep("\\.CTEN_ML_tree\\.treefile", all_h_trees, value = T))
  CTEN_ASTRAL_tree_path       <- paste0(hypothesis_tree_dir, grep("\\.CTEN_ASTRAL_tree\\.tre", all_h_trees, value = T))
  PORI_ML_tree_path           <- paste0(hypothesis_tree_dir, grep("\\.PORI_ML_tree\\.treefile", all_h_trees, value = T))
  PORI_ASTRAL_tree_path       <- paste0(hypothesis_tree_dir, grep("\\.PORI_ASTRAL_tree\\.tre", all_h_trees, value = T))
  CTEN_PORI_ML_tree_path      <- paste0(hypothesis_tree_dir, grep("\\.CTEN_PORI_ML_tree\\.treefile", all_h_trees, value = T))
  CTEN_PORI_ASTRAL_tree_path  <- paste0(hypothesis_tree_dir, grep("\\.CTEN_PORI_ASTRAL_tree\\.tre", all_h_trees, value = T))
  # Assemble all this information into a csv file
  constrained_cf_df <- data.frame(dataset = dataset,
                                  matrix = matrix,
                                  dataset_id = dataset_id,
                                  iqtree2 = iqtree2,
                                  astral = ASTRAL, 
                                  alignment_dir = alignment_dir,
                                  alignment_path = alignment_files,
                                  original_partition_path = original_partition_files, 
                                  modelfinder_partition_path = modelfinder_partition_files,
                                  constraint_tree_dir = constraint_tree_dir,
                                  CTEN_constraint_tree_path = CTEN_constraint_trees,
                                  PORI_constraint_tree_path = PORI_constraint_trees,
                                  hypothesis_tree_dir = hypothesis_tree_dir,
                                  CTEN_PORI_constraint_tree_path = CTEN_PORI_constraint_trees,
                                  CTEN_ML_tree_path = CTEN_ML_tree_path,
                                  CTEN_ASTRAL_tree_path = CTEN_ASTRAL_tree_path,
                                  PORI_ML_tree_path = PORI_ML_tree_path,
                                  PORI_ASTRAL_tree_path = PORI_ASTRAL_tree_path,
                                  CTEN_PORI_ML_tree_path = CTEN_PORI_ML_tree_path,
                                  CTEN_PORI_ASTRAL_tree_path = CTEN_PORI_ASTRAL_tree_path,
                                  collated_genes_dir = collated_genes_dir,
                                  CTEN_gene_trees = CTEN_gene_trees,
                                  PORI_gene_trees = PORI_gene_trees,
                                  CTEN_PORI_gene_trees = CTEN_PORI_gene_trees,
                                  MFP_gene_trees = MFP_gene_trees,
                                  output_dir = output_csv_dir,
                                  gcf_dir = gcf_dir,
                                  scf_dir = scf_dir,
                                  qcf_dir = qcf_dir)
  # Save the csv file
  constrained_cf_csv_path <- paste0(output_csv_dir, "constrained_cf_input.csv")
  write.csv(constrained_cf_df, file = constrained_cf_csv_path, row.names = FALSE)
}



#### 9. Nicely format ASTRAL output for plotting ####
# Open constrained_unconstrained_cf raw output
raw_in_df_file <- paste0(output_dir, "constrained_unconstrained_cf_output_raw.csv")
raw_in_df <- read.csv(raw_in_df_file, header = T, stringsAsFactors = F)
raw_in_df <- raw_in_df[ , c(1:30,37)]
# Separate into rows with and without NA for ASTRAL output
raw_na_df <- raw_in_df[which(is.na(raw_in_df$raw_ASTRAL_branch_labels)), ]
raw_astral_df <- raw_in_df[which(is.na(raw_in_df$raw_ASTRAL_branch_labels) == FALSE), ]
# Separate raw ASTRAL output into columns
raw_astral_df$raw_ASTRAL_branch_labels <- gsub("'", "", raw_astral_df$raw_ASTRAL_branch_label)
raw_astral_df$raw_ASTRAL_branch_labels <- gsub("\\[", "", gsub("\\]", "", raw_astral_df$raw_ASTRAL_branch_labels))
split_astral_vals <- strsplit(raw_astral_df$raw_ASTRAL_branch_labels, "\\;")
# Create new columns for the ASTRAL output values
raw_astral_df$q1  <- as.numeric(gsub("q1=", "", unlist(lapply(split_astral_vals, function(x){x[[1]]}))))
raw_astral_df$q2  <- as.numeric(gsub("q2=", "", unlist(lapply(split_astral_vals, function(x){x[[2]]}))))
raw_astral_df$q3  <- as.numeric(gsub("q3=", "", unlist(lapply(split_astral_vals, function(x){x[[3]]}))))
raw_astral_df$f1  <- as.numeric(gsub("f1=", "", unlist(lapply(split_astral_vals, function(x){x[[4]]}))))
raw_astral_df$f2  <- as.numeric(gsub("f2=", "", unlist(lapply(split_astral_vals, function(x){x[[5]]}))))
raw_astral_df$f3  <- as.numeric(gsub("f3=", "", unlist(lapply(split_astral_vals, function(x){x[[6]]}))))
raw_astral_df$pp1 <- as.numeric(gsub("pp1=", "", unlist(lapply(split_astral_vals, function(x){x[[7]]}))))
raw_astral_df$pp2 <- as.numeric(gsub("pp2=", "", unlist(lapply(split_astral_vals, function(x){x[[8]]}))))
raw_astral_df$pp3 <- as.numeric(gsub("pp3=", "", unlist(lapply(split_astral_vals, function(x){x[[9]]}))))
raw_astral_df$QC  <- as.numeric(gsub("QC=", "", unlist(lapply(split_astral_vals, function(x){x[[10]]}))))
raw_astral_df$EN  <- as.numeric(gsub("EN=", "", unlist(lapply(split_astral_vals, function(x){x[[11]]}))))
# Create new columns for the NA rows with all NA for all ASTRAL output values
raw_na_df$q1  <- NA
raw_na_df$q2  <- NA
raw_na_df$q3  <- NA
raw_na_df$f1  <- NA
raw_na_df$f2  <- NA
raw_na_df$f3  <- NA
raw_na_df$pp1 <- NA
raw_na_df$pp2 <- NA
raw_na_df$pp3 <- NA
raw_na_df$QC  <- NA
raw_na_df$EN  <- NA
# Reassemble dataframes
formatted_in_df <- rbind(raw_astral_df, raw_na_df)
# Rearrange row order
formatted_in_df <- formatted_in_df[order(formatted_in_df$dataset, formatted_in_df$matrix,
                                         formatted_in_df$hypothesis_tree, formatted_in_df$branch_num,
                                         formatted_in_df$branch_description, formatted_in_df$gene_type) , ]
rownames(formatted_in_df) <- 1:nrow(formatted_in_df)
# Rearrange column order
formatted_in_df <- formatted_in_df[ , c(1:29, 32:42, 31)]
# Write the nicely formatted csv out
formatted_in_df_path <- paste0(output_dir, "constrained_unconstrained_cf_output_formatted.csv")
write.csv(formatted_in_df, file = formatted_in_df_path, row.names = FALSE)



