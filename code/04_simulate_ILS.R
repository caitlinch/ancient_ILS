# ancient_ILS/code/01_simulate_ILS.R
## This script prepares and runs simulations of the animal tree of life under ILS
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
## File paths and computational parameters
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# hypothesis_tree_dir       <- location of constrained ML and ASTRAL trees estimated from empirical dataset Whelan et al. (2017)
# output_dir                <- output directory to save test runs for determining appropriate branch lengths
# iqtree2                   <- path to iqtree2 version 2.2.0
# iqtree2_num_threads       <- Maximum number of parallel threads to allow for IQ-Tree2 runs
# astral                    <- path to ASTRAL executable version 5.7.8
# ms                        <- path to ms executable
# tip_name_csv              <- Taxa reconciliation csv file (to make tip names consistent across datasets)

## Phylogenetic parameters
# substitution_model        <- model(s) to use when:
#                                    - estimating the alignment using Alisim (Alisim is included in IQ-Tree2)
#                                    - to use when estimating the ML trees from simulated alignments in IQ-Tree 2
#                                    - Note: Should be either one model (e.g. "LG") or a vector the same length as the number of genes (e.g. 117 models long)
# num_genes                 <- number of genes to simulate

## Control parameters
# create.simulation.parameters  <- flag to run code to generate parameters for simulations i.e., branch lengths, evolutionary hypothesis (to generate parameters, create.simulation.parameters=TRUE)
# generate.alignments           <- flag to run code to generate aligments from simulation parameters (to generate alignments: generate.alignments=TRUE)
# estimate.trees                <- flag to run code to estimate trees from simulated alignments (to estimate trees: estimate.trees=TRUE)
# conduct.analysis              <- flag to run code to to compare trees and calculate concordance factors (to run: conduct.analysis=TRUE)
# copy.completed.files          <- flag to run code to determine whether to copy output trees to a new file (to copy: copy.completed.files=TRUE)

location = "local"
if (location == "local"){
  # Local runs
  repo_dir              <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  hypothesis_tree_dir   <- paste0(repo_dir, "hypothesis_trees/")
  output_dir            <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/05_simulations/"
  iqtree2               <- "iqtree2"
  iqtree2_num_threads   <- 3
  iqtree2_num_ufb       <- 1000
  astral                <- "/Users/caitlincherryh/Documents/Executables/ASTRAL-5.7.8-master/Astral/astral.5.7.8.jar"
  ms                    <- "/Users/caitlincherryh/Documents/Executables/ms_exec/ms"
  tip_name_csv          <- paste0(repo_dir, "output/Cherryh_MAST_metazoa_taxa_reconciliation.csv")
  simion_parameters_csv <- paste0(repo_dir, "output/monophyletic_gene_tree_branch_models.csv")
} else if (location == "dayhoff"){
  # Remote runs
  repo_dir              <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  hypothesis_tree_dir   <- paste0(repo_dir, "hypothesis_trees/")
  output_dir            <- paste0(repo_dir, "simulation_output/")
  iqtree2               <- paste0(repo_dir, "iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2")
  iqtree2_num_threads   <- 20
  iqtree2_num_ufb       <- 1000
  astral                <- paste0(repo_dir, "astral/Astral/astral.5.7.8.jar")
  ms                    <- paste0(repo_dir, "msdir/ms")
  tip_name_csv          <- paste0(repo_dir, "output/Cherryh_MAST_metazoa_taxa_reconciliation.csv")
  simion_parameters_csv <- paste0(repo_dir, "output/monophyletic_gene_tree_branch_models.csv")
}

## Phylogenetic parameters
substitution_model          <- "LG"
substitution_model   <- "LG"
num_genes <- 200
min_branch_length           <- 10/(200*200) # 10/N, where N = total number of sites. Forces 10 substitutions onto each branch - reasonable over large time period/breadth of diversity
min_coalescent_difference   <- 0.001 # Keep to same magnitude as minimum coalescent interval in the ASTRAL species tree (~0.0027)

## Control parameters
control_parameters <- c("create.simulation.parameters" = FALSE,
                        "generate.alignments" = FALSE,
                        "estimate.trees" = FALSE,
                        "conduct.analysis" = TRUE)



#### 2. Open packages and functions ####
## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))
source(paste0(repo_dir, "code/func_simulations.R"))
source(paste0(repo_dir, "code/func_tree_estimation.R"))
source(paste0(repo_dir, "code/func_analysis.R"))

## Open packages
library(parallel)

# Create output file paths
output_files <- paste0(output_dir, "ancientILS_", c("simulation_parameters.csv", "output_generate_alignments.csv","output_generate_trees.csv", 
                                                    "output_analysis.csv", "branch_qcf_actual.csv", "branch_qcf_estimated.csv"))
names(output_files) <- c("simulations", "alignments", "trees", "analysis", "actual_qcf", "estimated_qcf")



#### 3. Create new shorter tip names ####
tip_name_df <- read.csv(tip_name_csv, stringsAsFactors = FALSE)
tip_name_df <- tip_name_df[tip_name_df$dataset == "Simion2017", ]
rownames(tip_name_df) <- 1:nrow(tip_name_df)
tip_name_df$numbered_name <- 1:nrow(tip_name_df)


#### 4. Create dataframe for simulations ####
if ( (file.exists(output_files[["simulations"]]) == FALSE) | (control_parameters[["create.simulation.parameters"]] == TRUE) ){
  # Create dataframe for ILS simulations
  # ILS simulations: vary length of branch a separating SOM and All Other Metazoans (in Simion2017 ASTRAL tree, branch leading to PORI+PLAC+CNID+BILAT)
  sim_df <- rbind(as.data.frame(expand.grid(replicates = 1:5, 
                                            hypothesis_tree = c(1,2),
                                            branch_a_empirical_length = 0.1544, 
                                            branch_b_empirical_length = 4.191,
                                            simulation_type = "ILS",
                                            simulation_number = "sim1") ) )
  # Add the other columns for the dataframes
  set_params <- c("branch_c_length" = 0.2741, "branch_cnidaria_length" = 0.6078, 
                  "branch_bilateria_length" = 0.7649, "branch_porifera_length" = 0.0864,
                  "branch_all_animals_length" = 0.3705, "branch_outgroup_length" = 0.3705,
                  "ML_tree_depth" = NA, "ASTRAL_tree_depth" = 9.089,
                  "proportion_internal_branches" = NA, "minimum_coalescent_time_difference" = min_coalescent_difference,
                  "num_taxa" = 97, "num_genes" = num_genes, "gene_length" = 200,
                  "dataset" = "Simion2017.supermatrix_97sp_401632pos_1719genes", "dataset_type" = "Protein",
                  "ms" = ms, "ASTRAL" = astral, "iqtree2" = iqtree2,
                  "iqtree2_num_threads" = iqtree2_num_threads, "iqtree2_num_ufb" = 0,
                  "alisim_gene_models" = substitution_model, "ML_tree_estimation_models" = substitution_model)
  set_params_df <- as.data.frame(matrix(set_params, nrow = 1, ncol = length(set_params), byrow = T))
  names(set_params_df) <- names(set_params)
  # Bind the columns from the set_params_df to the sim_df
  sim_df <- cbind(sim_df, set_params_df)
  # Add parameters that are dependent on other parameters
  sim_df$num_sites <- (as.numeric(sim_df$gene_length)*as.numeric(sim_df$num_genes))
  # Add the hypothesis tree name in a new column
  sim_df$hypothesis_tree_file <- as.character(sim_df$hypothesis_tree)
  sim_df$hypothesis_tree_file[which(sim_df$hypothesis_tree == 1)] <- paste0(hypothesis_tree_dir, "Simion2017_ASTRAL_hypothesis_tree_1_CTEN.tre")
  sim_df$hypothesis_tree_file[which(sim_df$hypothesis_tree == 2)] <- paste0(hypothesis_tree_dir, "Simion2017_ASTRAL_hypothesis_tree_2_PORI.tre")
  # Add an ID column
  sim_df$ID <- paste0(sim_df$simulation_number, "_h", sim_df$hypothesis_tree, "_a", sim_df$branch_a_length, "_b", sim_df$branch_b_length, "_rep", sim_df$replicates)
  # Add separate output folder for each rep
  sim_df$output_folder <- paste0(output_dir, sim_df$ID, "/")
  # Reorder columns
  # Note: minimum_coalescent_time_difference column is the time (in coalescent units) to adjust node depth by if multiple nodes with the same taxa have the same coalescence time
  sim_df <- sim_df[,c("dataset", "dataset_type", "ID", "simulation_number", "simulation_type", "hypothesis_tree", "hypothesis_tree_file", "replicates",
                      "branch_a_length", "branch_b_length", "branch_c_length", "branch_all_animals_length", "branch_bilateria_length", "branch_cnidaria_length",
                      "branch_outgroup_length", "branch_porifera_length", "proportion_internal_branches", "minimum_coalescent_time_difference", "ASTRAL_tree_depth",
                      "ML_tree_depth", "num_taxa", "num_genes", "gene_length", "num_sites", "output_folder", "ms", "ASTRAL", "iqtree2", "alisim_gene_models",
                      "iqtree2_num_threads", "iqtree2_num_ufb", "ML_tree_estimation_models")]
  # Save the output file
  write.csv(sim_df, file = output_files[["simulations"]], row.names = FALSE)
}



#### 5. Generate simulated alignments ####
if ((control_parameters[["generate.alignments"]] == TRUE) & (file.exists(output_files[["alignments"]]) == FALSE) ){
  # Open dataframe of simulation parameters
  sim_df <- read.csv(output_files[["simulations"]], stringsAsFactors = FALSE)
  # To generate all simulated alignments
  if (location == "local"){
    generate_alignment_list <- lapply(1:nrow(sim_df), generate.one.alignment.wrapper, sim_df = sim_df, 
                                      converted_taxa_names = simulation_taxa_names,rerun = FALSE)
  } else {
    generate_alignment_list <- mclapply(1:nrow(sim_df), generate.one.alignment.wrapper, sim_df = sim_df, 
                                        converted_taxa_names = simulation_taxa_names, rerun = FALSE, 
                                        mc.cores = num_parallel_threads)
  }
  generate_alignment_df <- as.data.frame(do.call(rbind, generate_alignment_list))
  # Save output dataframe
  write.csv(generate_alignment_df, file = output_files[["alignments"]], row.names = FALSE)
}



#### 6. Estimate trees ####
if ( (control_parameters[["estimate.trees"]] == TRUE) & (file.exists(output_files[["trees"]]) == FALSE) ){
  # Read in the dataframe from the previous step
  generate_alignment_df <- read.csv(output_files[["alignments"]], stringsAsFactors = FALSE)
  # To estimate all trees with a set model for all single simulated alignments
  if (location == "local"){
    tree_list <- lapply(1:nrow(generate_alignment_df), estimate.trees, 
                        df = generate_alignment_df, call.executable.programs = TRUE)
  } else {
    tree_list <- mclapply(1:nrow(generate_alignment_df), estimate.trees, 
                          df = generate_alignment_df, call.executable.programs = FALSE,
                          mc.cores = floor(num_parallel_threads/iqtree2_num_threads) )
  }
  tree_df <- as.data.frame(do.call(rbind, tree_list))
  # Save the dataframe of tree estimation calls
  write.csv(tree_df, file = output_files[["trees"]], row.names = FALSE)
}



#### 7. Conduct analyses ####
if ( (control_parameters[["conduct.analysis"]] == TRUE) & (file.exists(output_files[["analysis"]]) == FALSE) ){
  # Read in the output from the previous step
  tree_df <- read.csv(output_files[["trees"]], stringsAsFactors = FALSE)
  # Call the function to calculate gCF and RF distances
  # To run for one row:
  #       analysis.wrapper(1, df = tree_df, ASTRAL_path = astral, hypothesis_tree_dir = hypothesis_tree_dir, converted_taxa_names = simulation_taxa_names)
  if (location == "local"){
    lapply(1:nrow(tree_df), analysis.wrapper, df = tree_df, ASTRAL_path = astral,
           hypothesis_tree_dir = hypothesis_tree_dir, converted_taxa_names = simulation_taxa_names, rerun.ASTRAL = TRUE)
  } else {
    mclapply(1:nrow(tree_df), analysis.wrapper, df = tree_df, ASTRAL_path = astral,
             hypothesis_tree_dir = hypothesis_tree_dir, converted_taxa_names = simulation_taxa_names, rerun.ASTRAL = TRUE,
             mc.cores = num_parallel_threads)
  }
  
  # Collect and save completed output csvs
  all_output_csvs <- paste0(tree_df$output_folder, tree_df$ID, "_analysis_output.csv")
  completed_csvs <- all_output_csvs[which(file.exists(all_output_csvs) == TRUE)]
  completed_csv_list <- lapply(completed_csvs, read.csv, stringsAsFactors = FALSE)
  completed_csv_df <- as.data.frame(do.call(rbind, completed_csv_list))
  write.csv(completed_csv_df, file = output_files[["analysis"]], row.names = FALSE)
  
  # Collect and save actual qcfs csv files
  a_csv_files <- completed_csv_df$actual_qcf_branch_output_csv
  a_qcf_list <- lapply(a_csv_files, read.csv, stringsAsFactors = FALSE)
  a_qcf_df <- as.data.frame(do.call(rbind, a_qcf_list))
  write.csv(a_qcf_df, file = output_files[["actual_qcf"]], row.names = FALSE)
  
  # Collect and save estimated qcfs csv files
  e_csv_files <- completed_csv_df$estimated_qcf_branch_output_csv
  e_qcf_list <- lapply(e_csv_files, read.csv, stringsAsFactors = FALSE)
  e_qcf_df <- as.data.frame(do.call(rbind, e_qcf_list))
  write.csv(e_qcf_df, file = output_files[["estimated_qcf"]], row.names = FALSE)
}


