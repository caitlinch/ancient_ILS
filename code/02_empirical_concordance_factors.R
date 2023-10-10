## caitlinch/ancient_ILS/code/02_empirical_concordance_factors.R
# This script estimates site, gene, and quartet concordance factors for empirical gene trees and constrained hypothesis trees.
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## Specify parameters:
# alignment_dir               <- Directory containing alignments for all data sets
#                                   Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                                   E.g. Cherryh2022.alignment1.aa.fa
# output_dir                  <- Directory for IQ-Tree output (trees and tree mixtures)
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# iqtree2                     <- Location of IQ-Tree2 stable release
# iqtree_num_threads          <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
# ml_tree_bootstraps          <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree
# hypothesis_tree_bootstraps  <- Number of ultrafast bootstraps (UFB) to perform in IQ-Tree when estimating constrained maximum likelihood trees
# number_parallel_processes   <- The number of simultaneous processes to run at once using mclapply(). 
#                                   If 1, then all processes will run sequentially

location = "local"
if (location == "local"){
  alignment_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/"
  repo_dir <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  iqtree2 <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/03_Software_IQ-Tree/iqtree-2.2.0-MacOSX/bin/iqtree2"
  number_parallel_processes <- 1
  
} else if (location == "dayhoff" | location == "rona" ){
  if (location == "dayhoff"){
    repo_dir <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
    iqtree2 <- "/mnt/data/dayhoff/home/u5348329/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  } else if (location == "rona"){
    repo_dir <- "/home/caitlin/ancient_ILS/"
    iqtree2 <- "/home/caitlin/metazoan-mixtures/iqtree/iqtree-2.2.0-Linux/bin/iqtree2"
  }
  alignment_dir <- paste0(repo_dir, "data_all/")
  output_dir <-  paste0(repo_dir, "output/")
  number_parallel_processes = 1
}

# Set parameters that are identical for all run locations
iqtree_num_threads <- 15
ml_tree_bootstraps <- 1000
hypothesis_tree_bootstraps <- 1000

# Set control parameters
control_parameters <- list(estimate.pmsf.gene.trees = FALSE,
                           estimate.ModelFinder.gene.trees = FALSE,
                           estimate.constrained.pmsf.trees = FALSE,
                           estimate.constrained.ModelFinder.trees = FALSE,
                           estimate.partitioned.ML.trees = FALSE)

