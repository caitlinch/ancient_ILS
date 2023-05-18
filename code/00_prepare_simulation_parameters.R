# ancient_ILS/code/00_prepare_simulation_parameters.R
## This script is used to determine realistic simulation parameters
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# repo_dir                              <- Location of caitlinch/ancient_ILS github repository
# constrained_tree_output_dir           <- directory to save constrained ML trees estimated from the alignment
# tree_output_dir                       <- directory to save other ML trees estimated from the alignment
# alignment_path                        <- path to the alignment from Whelan et al. 2017 for the alignment Metazoa_Choano_RCFV_strict
# models_path                           <- path to the models and genes from Whelan et al. 2017 for the alignment Metazoa_Choano_RCFV_strict

location = "dayhoff"
if (location == "local"){
  repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  constrained_tree_output_dir <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_hypothesis_trees/"
  tree_output_dir             <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_ML_tree_estimation/"
  alignment_path              <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/Whelan2017.Metazoa_Choano_RCFV_strict.aa.alignment.fa"
  models_path                 <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/Metazoa_Choano_RCFV_strict_Models.txt"
} else if (location == "dayhoff"){
  repo_dir                    <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  constrained_tree_output_dir <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/02_empirical_hypothesis_trees/"
  tree_output_dir             <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/02_empirical_tree_estimation/"
  alignment_path              <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/01_empirical_data/Whelan2017.Metazoa_Choano_RCFV_strict.aa.alignment.fa"
  models_path                 <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/01_empirical_data/Metazoa_Choano_RCFV_strict_Models.txt"
}



#### 2. Open packages and functions ####
## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))