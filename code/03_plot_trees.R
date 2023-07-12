# ancient_ILS/code/03_plot_trees.R
## This script plots and outputs phylogenetic tree figures
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# repo_dir            <- Location of caitlinch/ancient_ILS github repository

repo_dir <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"



#### 2. Open packages and prepare variables ####
library(ape)
library(ggtree)
library(patchwork)



#### 4. Plot figures for introduction and methods ####
# Identify the ASTRAL tree
empirical_trees <- list.files(paste0(repo_dir, "empirical_tree"))
astral_tree <- paste0(repo_dir, "empirical_tree/", grep("Whelan2017_ASTRAL_tree.tre", empirical_trees, value = T))
# Open the ASTRAL tree
a_tree <- read.tree(astral_tree)
# Reroot tree
a_tree <- root(a_tree, outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), resolve.root = TRUE)
# Set missing branch lengths to 2
a_tree$edge.length[which(is.na(a_tree$edge.length))] <- 3


#### 4. Plot phylogenetic hypotheses ####
# Open the possible phylogenetic topologies
hyp_trees <- read.tree(file = paste0(repo_dir, "hypothesis_trees/alternative_phylogenetic_hypotheses.nex"))
c_tree <- hyp_trees[[1]]
p_tree <- hyp_trees[[2]]
cp_tree <- hyp_trees[[3]]




