# ancient_ILS/code/01_simulate_ILS_and_LBA.R
## This script prepares and runs simulations of the animal tree of life under LBA and ILS
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# repo_dir            <- Location of caitlinch/ancient_ILS github repository
# output_dir          <- directory to save output
# alignment_path      <- path to the alignment from Whelan et al. 2017 for the alignment Metazoa_Choano_RCFV_strict
# partition_path      <- path to the partition file for the alignment Metazoa_Choano_RCFV_strict (created in script 00_prepare_trees.R)
# iqtree2             <- path to iqtree2 version 2.2.2

repo_dir <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_dir <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_simulations/"
alignment_path <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/Whelan2017.Metazoa_Choano_RCFV_strict.aa.alignment.fa"
partition_path <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_hypothesis_trees/Whelan2017_models_partitions.nex"
iqtree2 <- "iqtree2"



#### 2. Open packages and functions ####
source(paste0(repo_dir, "code/func_prepare_trees.R"))



#### 3. Prepare simulation dataframe ####
sim_df <- as.data.frame(expand.grid(replicates = c(1,2,3,4,5,6,7,8,9,10), branch_a = c(5,20,15,20,25,30,35), branch_b = c(10,20,30,40,50,60,70), hypothesis_tree = c(1,2,3)))
sim_df <- sim_df[,c("hypothesis_tree", "branch_a", "branch_b", "replicates")]
                        