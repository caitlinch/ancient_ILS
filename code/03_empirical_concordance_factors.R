## caitlinch/ancient_ILS/code/03_empirical_concordance_factors.R
# This script estimates gene, and quartet concordance factors for empirical gene trees and constrained hypothesis trees.
# Caitlin Cherryh 2024



#### 1. Input parameters ####
## Specify parameters:
# location                    <- Where the script is being run
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# output_csv_dir              <- Location of csv input/output/results files
# cf_dir                      <- Location of concordance factor output files
# iqtree2_server              <- Location of IQ-Tree2 executable
# iqtree2_num_threads         <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
# astral_server               <- Location of ASTRAL executable
# noPlac_dataset_ids          <- datasets to remove Placozoa tips from trees before calculating qCF for Partition models 
#                                   (constrained trees have different Placozoa placement, which throws out quartet comparisons)
# control                     <- Control which sections of the code run

repo_dir              <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_csv_dir        <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/07_output_files/"
cf_dir                <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/06_cf_analyses/"

# Server file paths to run qCF/gCF analyses
iqtree2_server        <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2"
iqtree2_num_threads   <- 30
astral_server         <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/astral/Astral/astral.5.7.8.jar"

# Specify datasets to remove Placozoa before calculating qCF for Partition models 
#     (due to alternative placements of the Placozoa taxa in the different constrained tree topologies)
noPlac_dataset_ids <- c("Chang2015.Chang_AA", "Laumer2018.Tplx_BUSCOeuk", "Nosenko2013.ribosomal_14615_smatrix",
                        "Philippe2009.Philippe_etal_superalignment_FixedNames", "Whelan2015.Dataset10")

# Control commands
control <- list(remove.Plac = FALSE,
                run.cf.analyses = FALSE,
                call.executables = FALSE,
                extract.gcf = FALSE,
                extract.qcf = FALSE,
                reformat.dataframes = FALSE)



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



#### 4. Copy data but remove Placozoa taxa for the problem datasets with inconsistent PLAC placement ####
if (control$remove.Plac == TRUE){
  # Open the required dataframe with dataset information
  input_df <- read.csv(paste0(output_csv_dir, "cf_analysis_input_paths.csv"), stringsAsFactors = FALSE)
  # Duplicate df and add model
  input_df <- rbind(input_df, input_df)
  input_df$model <- c(rep("C60", (nrow(input_df)/2)), rep("Partition", (nrow(input_df)/2)) )
  # Identify which files to update
  noPlac_rows <- which(input_df$dataset_id %in% noPlac_dataset_ids & input_df$model == "Partition")
  # Remove Placozoa taxa from these dataset (ML tree and gene trees) and update file path to "noPlac" path
  input_df$c60_gene_trees[noPlac_rows] <- unlist(lapply(input_df$c60_gene_trees[noPlac_rows], gene.trees.remove.Plac, 
                                                        all_datasets = all_datasets, matrix_taxa = matrix_taxa, alignment_taxa_df = alignment_taxa_df))
  input_df$C60_CTEN_tree[noPlac_rows] <- unlist(lapply(input_df$C60_CTEN_tree[noPlac_rows], tree.remove.Plac, 
                                                       all_datasets = all_datasets, matrix_taxa = matrix_taxa, alignment_taxa_df = alignment_taxa_df))
  input_df$C60_PORI_tree[noPlac_rows] <- unlist(lapply(input_df$C60_PORI_tree[noPlac_rows], tree.remove.Plac, 
                                                       all_datasets = all_datasets, matrix_taxa = matrix_taxa, alignment_taxa_df = alignment_taxa_df))
  input_df$C60_CTEN_PORI_tree[noPlac_rows] <- unlist(lapply(input_df$C60_CTEN_PORI_tree[noPlac_rows], tree.remove.Plac, 
                                                            all_datasets = all_datasets, matrix_taxa = matrix_taxa, alignment_taxa_df = alignment_taxa_df))
  input_df$mfp_gene_trees[noPlac_rows] <- unlist(lapply(input_df$mfp_gene_trees[noPlac_rows], gene.trees.remove.Plac, 
                                                        all_datasets = all_datasets, matrix_taxa = matrix_taxa, alignment_taxa_df = alignment_taxa_df))
  input_df$partition_CTEN_tree[noPlac_rows] <- unlist(lapply(input_df$partition_CTEN_tree[noPlac_rows], tree.remove.Plac, 
                                                             all_datasets = all_datasets, matrix_taxa = matrix_taxa, alignment_taxa_df = alignment_taxa_df))
  input_df$partition_PORI_tree[noPlac_rows] <- unlist(lapply(input_df$partition_PORI_tree[noPlac_rows], tree.remove.Plac, 
                                                             all_datasets = all_datasets, matrix_taxa = matrix_taxa, alignment_taxa_df = alignment_taxa_df))
  input_df$partition_CTEN_PORI_tree[noPlac_rows] <- unlist(lapply(input_df$partition_CTEN_PORI_tree[noPlac_rows], tree.remove.Plac, 
                                                                  all_datasets = all_datasets, matrix_taxa = matrix_taxa, alignment_taxa_df = alignment_taxa_df))
  # Update output paths
  input_df$partition_gcf_output_dir <- paste0(dirname(input_df$partition_gcf_output_dir), "/gcf_", tolower(input_df$model), "/")
  input_df$partition_qcf_output_dir <- paste0(dirname(input_df$partition_qcf_output_dir), "/gcf_", tolower(input_df$model), "/")
  # Save the updated input file
  noPlac_input_df_file <- paste0(output_csv_dir, "cf_analysis_input_paths_noPlac.csv")
  write.csv(input_df, file = noPlac_input_df_file, row.names = F)
}



#### 3. Create command lines for calculating gCF and qCF ####
if (control$run.cf.analyses == TRUE){
  # Open the dataframe to set command lines
  noPlac_input_df_file <- paste0(output_csv_dir, "cf_analysis_input_paths_noPlac.csv")
  input_df <- read.csv(noPlac_input_df_file, stringsAsFactors = FALSE)
  
  # Identify which rows to update with "noPlac" marker
  noPlac_rows <- which(input_df$dataset_id %in% noPlac_dataset_ids & input_df$model == "Partition")
  # Add column to dictate output file note for noPlac runs
  input_df$Plac_id <- ""
  input_df$Plac_id[noPlac_rows] <- "noPlac."
  
  # Add gCF commands
  # $ iqtree2 -te concat.treefile --gcf loci.treefile --prefix concord
  input_df$cten_gcf_prefix  <- paste0(input_df$partition_gcf_output_dir, input_df$dataset_id, ".", input_df$model, ".CTEN.gcf")
  input_df$cten_gcf_command <- paste0(iqtree2_server, " -te ", input_df$partition_CTEN_tree,
                                      " --gcf ", input_df$mfp_gene_trees, 
                                      " -pre ", input_df$cten_gcf_prefix,
                                      " -nt ", iqtree2_num_threads)
  input_df$pori_gcf_prefix  <- paste0(input_df$partition_gcf_output_dir, input_df$dataset_id, ".", input_df$model, ".PORI.gcf")
  input_df$pori_gcf_command <- paste0(iqtree2_server, " -te ", input_df$partition_PORI_tree,
                                      " --gcf ", input_df$mfp_gene_trees, 
                                      " -pre ", input_df$pori_gcf_prefix,
                                      " -nt ", iqtree2_num_threads)
  input_df$cten_pori_gcf_prefix   <- paste0(input_df$partition_gcf_output_dir, input_df$dataset_id, ".", input_df$model, ".CTEN_PORI.gcf")
  input_df$cten_pori_gcf_command  <- paste0(iqtree2_server, " -te ", input_df$partition_CTEN_PORI_tree,
                                            " --gcf ", input_df$mfp_gene_trees, 
                                            " -pre ", input_df$cten_pori_gcf_prefix,
                                            " -nt ", iqtree2_num_threads)
  
  # Add qCF commands
  # $ java -jar astral.5.7.8.jar -q test_data/1kp.tre -i test_data/1KP-genetrees.tre -t 2 -o test_data/1kp-scored-t2.tre
  input_df$CTEN_qcf_tree    <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".", input_df$model, ".CTEN.", input_df$Plac_id, "qcf.tre")
  input_df$CTEN_qcf_log     <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".", input_df$model, ".CTEN.", input_df$Plac_id, "qcf.log")
  input_df$cten_qcf_command <- paste0("java -jar ", astral_server, 
                                      " -q ", input_df$partition_CTEN_tree, 
                                      " -i ", input_df$mfp_gene_trees, " -t 2 ",
                                      " -o ", input_df$CTEN_qcf_tree,
                                      " 2> ", input_df$CTEN_qcf_log)
  input_df$PORI_qcf_tree    <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".", input_df$model, ".PORI.", input_df$Plac_id, "qcf.tre")
  input_df$PORI_qcf_log     <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".", input_df$model, ".PORI.", input_df$Plac_id, "qcf.log")
  input_df$pori_qcf_command <- paste0("java -jar ", astral_server, 
                                      " -q ", input_df$partition_PORI_tree, 
                                      " -i ", input_df$mfp_gene_trees, " -t 2 ",
                                      " -o ", input_df$PORI_qcf_tree,
                                      " 2> ", input_df$PORI_qcf_log)
  input_df$CTEN_PORI_qcf_tree     <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".", input_df$model, ".CTEN_PORI.", input_df$Plac_id, "qcf.tre")
  input_df$CTEN_PORI_qcf_log      <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".", input_df$model, ".CTEN_PORI.", input_df$Plac_id, "qcf.log")
  input_df$cten_pori_qcf_command  <- paste0("java -jar ", astral_server, 
                                            " -q ", input_df$partition_CTEN_PORI_tree, 
                                            " -i ", input_df$mfp_gene_trees, " -t 2 ",
                                            " -o ", input_df$CTEN_PORI_qcf_tree,
                                            " 2> ", input_df$CTEN_PORI_qcf_log)
  
  # Write out command lines to text files
  write(c(input_df$cten_gcf_command, input_df$pori_gcf_command, input_df$cten_pori_gcf_command), 
        file = paste0(output_csv_dir, "gcf_commands_noPlac.txt"))
  write(c(input_df$cten_qcf_command, input_df$pori_qcf_command, input_df$cten_pori_qcf_command), 
        file = paste0(output_csv_dir, "qcf_commands_noPlac.txt"))
  
  # Call executables, if desired
  if (control$call.executables == TRUE){
    # Run gCF calculation, and save screen output as an object to assess for errors
    gcf_run_log <- lapply(c(input_df$cten_gcf_command, input_df$pori_gcf_command, input_df$cten_pori_gcf_command), system, intern = TRUE)
    # Run qCF calculation, and save screen output as an object to assess for errors
    qcf_run_log <- lapply(c(input_df$cten_qcf_command, input_df$pori_qcf_command, input_df$cten_pori_qcf_command), system, intern = TRUE)
  }
  
  # Save commands in the csv file
  write.csv(input_df, file = paste0(output_csv_dir, "cf_analysis_commands_noPlac.csv"), row.names = FALSE)
}



#### 4. Extract gCF values from key clades  ####
if (control$extract.qcf == TRUE){
  ## Get gCF output
  # Specify gCF parameters df
  gcf_df_file <- paste0(output_csv_dir, "gCF_tree_files.csv")
  # Open gCF parameters dataframe
  if (file.exists(gcf_df_file) == FALSE){
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
    # Add note of whether Placozoa is present and remove the rows with inconsistent Placozoa placement
    gcf_df$Plac_present <- "Plac"
    gcf_df$Plac_present[grep("noPlac", gcf_df$gcf_stat_files)] <- "noPlac"
    gcf_df <- gcf_df[which( ! (gcf_df$dataset_id %in% noPlac_dataset_ids & gcf_df$model == "Partition" & gcf_df$Plac_present == "Plac") ), ]
    # Rearrange order of columns
    gcf_df <- gcf_df[, c("id", "dataset", "matrix_name", "dataset_id", "model", "tree_topology", "Plac_present", "gcf_branch_files", "gcf_stat_files")]
    # Write qCF_df
    write.csv(gcf_df, file = gcf_df_file, row.names = FALSE)
  } else {
    gcf_df <- read.csv(gcf_df_file, stringsAsFactors = FALSE)
  }
  
  ## Extract gCF values
  # Extract all gCF values
  gcf_output_list <- lapply(1:nrow(gcf_df), extract.gcf.wrapper, gcf_df = gcf_df, 
                            matrix_taxa = matrix_taxa, all_datasets = all_datasets, 
                            alignment_taxa_df = alignment_taxa_df)
  # Format and output qCF dataframe
  gcf_output_df <- as.data.frame(do.call(rbind, gcf_output_list), stringsAsFactors = FALSE)
  # Cbind to the parameters dataframe
  gcf_collated_df <- cbind(gcf_df, gcf_output_df)
  
  ## Save output
  gcf_collated_df_file <- paste0(output_csv_dir, "gCF_values.csv")
  write.csv(gcf_collated_df, file = gcf_collated_df_file, row.names = FALSE)
}



#### 6. Extract qCF values from key clades  ####
if (control$extract.qcf == TRUE){
  ## Get qCF output
  # Specify qCF parameters df
  qcf_df_file <- paste0(output_csv_dir, "qCF_tree_files.csv")
  # Open qCF parameters dataframe
  if (file.exists(qcf_df_file) == FALSE){
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
    # Add note of whether Placozoa is present
    qcf_df$Plac_present <- "Plac"
    qcf_df$Plac_present[grep("noPlac", qcf_df$qcf_tree_file)] <- "noPlac"
    qcf_df <- qcf_df[which( ! (qcf_df$dataset_id %in% noPlac_dataset_ids & qcf_df$model == "Partition" & qcf_df$Plac_present == "Plac") ), ]
    # Rearrange order of columns
    qcf_df <- qcf_df[, c("id", "dataset", "matrix_name", "dataset_id", "model", "tree_topology", "Plac_present", "qcf_tree_file")]
    # Write qCF_df
    write.csv(qcf_df, file = qcf_df_file, row.names = FALSE)
  } else {
    qcf_df <- read.csv(qcf_df_file, stringsAsFactors = FALSE)
  }
  
  ## Extract qCF values
  # Construct id for iterating over analyses
  qcf_df$analysis_id <- paste0(qcf_df$dataset_id, ".", qcf_df$model)
  qcf_params <- unique(qcf_df$analysis_id)
  # Extract qCF for all datasets
  qcf_output_list <- lapply(qcf_params, extract.qcf.wrapper, qcf_df = qcf_df, 
                            matrix_taxa = matrix_taxa, all_datasets = all_datasets, 
                            alignment_taxa_df = alignment_taxa_df)
  # Format and output qCF dataframe
  qcf_collated_df <- as.data.frame(do.call(rbind, qcf_output_list), stringsAsFactors = FALSE)
  # Add note of whether Placozoa is present
  qcf_collated_df$Plac_present <- "Plac"
  qcf_collated_df$Plac_present[which(qcf_collated_df$dataset_id %in% noPlac_dataset_ids & qcf_collated_df$model == "Partition")] <- "noPlac"
  
  ## Save output
  qcf_collated_df_file <- paste0(output_csv_dir, "qCF_values.csv")
  write.csv(qcf_collated_df, file = qcf_collated_df_file, row.names = FALSE)
}



#### 7. Reformat gCF and qCF dataframes ####
if (control$reformat.dataframes == TRUE){
  ## Open dataframes
  gcf_collated_df_file  <- paste0(output_csv_dir, "gCF_values.csv")
  gcf_collated_df       <- read.csv(gcf_collated_df_file, stringsAsFactors = FALSE)
  qcf_collated_df_file  <- paste0(output_csv_dir, "qCF_values.csv")
  qcf_collated_df       <- read.csv(qcf_collated_df_file, stringsAsFactors = FALSE)
  
  ## Reformat gCF df
  # Call function to reformat dataframe
  gcf_clean_df <- reformat.gCF.df(input_df = gcf_collated_df)
  # Save reformatted dataframe
  gcf_clean_df_file  <- paste0(output_csv_dir, "gCF_values_formatted.csv")
  write.csv(gcf_clean_df, file = gcf_clean_df_file, row.names = FALSE)
  # New dataframe to check values
  check_gcf_df <- data.frame(id = gcf_clean_df$id, 
                             sum_gcf = (as.numeric(gcf_clean_df$CTEN.KEY_gCF) + 
                                          as.numeric(gcf_clean_df$PORI.KEY_gCF) + 
                                          as.numeric(gcf_clean_df$CTENPORI.KEY_gCF) + 
                                          as.numeric(gcf_clean_df$CTEN.KEY_gDFP) ) )
  # Save check dataframe
  check_gcf_df_file  <- paste0(output_csv_dir, "gCF_check.csv")
  write.csv(check_gcf_df, file = check_gcf_df_file, row.names = FALSE)
  
  ## Check qCF values
  # Call function to reformat dataframe
  qcf_clean_df <- reformat.qCF.df(input_df = qcf_collated_df)
  # Save reformatted dataframe
  qcf_clean_df_file  <- paste0(output_csv_dir, "qCF_values_formatted.csv")
  write.csv(qcf_clean_df, file = qcf_clean_df_file, row.names = FALSE)
  # CTEN_q1 + PORI_q1 + CTEN_PORI_q1 = 1
  check_qcf_df <- data.frame(id = qcf_clean_df$id, 
                             sum_q = (as.numeric(qcf_clean_df$CTEN.KEY_q1) + 
                                        as.numeric(qcf_clean_df$PORI.KEY_q1) + 
                                        as.numeric(qcf_clean_df$CTENPORI.KEY_q1)) )
  # Save check dataframe
  check_qcf_df_file  <- paste0(output_csv_dir, "qCF_check.csv")
  write.csv(check_qcf_df, file = check_qcf_df_file, row.names = FALSE)
}



