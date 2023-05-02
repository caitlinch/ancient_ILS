# ancient_ILS/code/01_simulate_ILS_and_LBA.R
## This script prepares and runs simulations of the animal tree of life under LBA and ILS
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# repo_dir            <- location of caitlinch/ancient_ILS github repository
# hypothesis_tree_dir <- location to directory containing the hypothesis trees (estimated in script 00_prepare_trees.R of this repository)
# output_dir          <- directory to save output
# alignment_path      <- path to the alignment from Whelan et al. 2017 for the alignment Metazoa_Choano_RCFV_strict
# partition_path      <- path to the partition file for the alignment Metazoa_Choano_RCFV_strict (created in script 00_prepare_trees.R)
# iqtree2             <- path to iqtree2 version 2.2.2
# ms                  <- path to ms executable

repo_dir              <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
hypothesis_tree_dir   <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_hypothesis_trees/"
output_dir            <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_simulations/"

alignment_path        <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/Whelan2017.Metazoa_Choano_RCFV_strict.aa.alignment.fa"
partition_path        <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_hypothesis_trees/Whelan2017_models_partitions.nex"
iqtree2               <- "iqtree2"
ms                    <- "/Users/caitlincherryh/Documents/Executables/ms_exec/ms"



#### 2. Open packages and functions ####
source(paste0(repo_dir, "code/func_simulations.R"))



#### 3. Prepare simulation dataframe ####
# Prepare the three sets of simulations separately
sim_df_1 <- as.data.frame(expand.grid(replicates = c(1,2,3,4,5,6,7,8,9,10), branch_a_percent_height = 5.8, 
                                    branch_b_percent_height = c(10,20,30,40,50,60,70), hypothesis_tree = c(1,2,3)))
sim_df_1$simulation_type = "LBA" # vary branch b
sim_df_1$simulation_number = "sim1"
sim_df_2 <- as.data.frame(expand.grid(replicates = c(1,2,3,4,5,6,7,8,9,10), branch_a_percent_height = c(5,20,15,20,25,30,35,40), 
                                      branch_b_percent_height = 38.5, hypothesis_tree = c(1,2,3)))
sim_df_2$simulation_type = "ILS" # vary branch a
sim_df_2$simulation_number = "sim2"
sim_df_3 <- as.data.frame(expand.grid(replicates = c(1,2,3,4,5,6,7,8,9,10), branch_a_percent_height = c(5,20,15,20,25,30,35,40), 
                                      branch_b_percent_height = c(10,20,30,40,50,60,70), hypothesis_tree = c(1,2,3)))
sim_df_3$simulation_type = "LBA+ILS" # vary branch a and branch b
sim_df_3$simulation_number = "sim3"
# Collate the three sets of simulations into a single dataframe
sim_df <- rbind(sim_df_1,sim_df_2, sim_df_3)
# Add the other columns for the dataframes
sim_df$tree_length <- 1.28
sim_df$branch_a_empirical_length <- 0.0746
sim_df$branch_a_simulation_length <- unique(sim_df$tree_length) * (sim_df$branch_a_percent_height/100)
sim_df$branch_b_empirical_length <- 0.4927
sim_df$branch_b_simulation_length <- unique(sim_df$tree_length) * (sim_df$branch_b_percent_height/100)
sim_df$num_taxa <- 76
sim_df$num_genes <- 117
sim_df$num_sites <- 49388
sim_df$gene_length <- "From alignment"
sim_df$dataset <- "Whelan2017.Metazoa_Choano_RCFV_strict"
sim_df$dataset_type <- "Protein"
# Add path for executables
sim_df$ms <- ms
sim_df$iqtree2 <- iqtree2
# Add the hypothesis tree name in a new column
sim_df$hypothesis_tree_file <- as.character(sim_df$hypothesis_tree)
sim_df$hypothesis_tree_file[which(sim_df$hypothesis_tree == 1)] <- paste0(hypothesis_tree_dir, "Whelan2017_hypothesis_tree_1_Cten.treefile")
sim_df$hypothesis_tree_file[which(sim_df$hypothesis_tree == 2)] <- paste0(hypothesis_tree_dir, "Whelan2017_hypothesis_tree_2_Pori.treefile")
sim_df$hypothesis_tree_file[which(sim_df$hypothesis_tree == 3)] <- paste0(hypothesis_tree_dir, "Whelan2017_hypothesis_tree_3_CtenPori.treefile")
# Add an ID column
sim_df$ID <- paste0(sim_df$simulation_number, "_h", sim_df$hypothesis_tree, "_a", sim_df$branch_a_percent_height, "_b", sim_df$branch_b_percent_height, "_rep", sim_df$replicates)
# Add separate output folder for each rep
sim_df$output_folder <- paste0(output_dir, sim_df$ID, "/")
# Reorder columns
sim_df <- sim_df[,c("ID", "dataset", "dataset_type", "num_taxa", "num_genes", "gene_length", "num_sites", "tree_length",
                    "branch_a_empirical_length", "branch_b_empirical_length", "simulation_number", "simulation_type",
                    "hypothesis_tree", "hypothesis_tree_file", "branch_a_percent_height", "branch_a_simulation_length",
                    "branch_b_percent_height", "branch_b_simulation_length", "replicates", "output_folder", "ms", "iqtree2")]
# Save the dataframe
sim_df_op_file <- paste0(output_dir, "ancientILS_simulation_parameters.csv")
write.csv(sim_df, file = sim_df_op_file)





                        