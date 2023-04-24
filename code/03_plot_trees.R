# ancient_ILS/code/03_plot_trees.R
## This script plots and outputs phylogenetic tree figures
# Caitlin Cherryh, 2022

#### 1. Input parameters ####
# repo_dir            <- Location of caitlinch/ancient_ILS github repository

repo_dir <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"



#### 2. Open packages and prepare variables ####
library(ape)
library(patchwork)


#### 3. Plot phylogenetic hypotheses ####
# Open the trees
hyp_trees <- read.tree(file = paste0(repo_dir, "/trees/alternative_phylogenetic_hypotheses.nex"))

# Plot possible topologies
top_trees <- hyp_trees[1:3]
par(mfrow=c(1,1))
plot(hyp_trees[[1]], cex = 2, edge.width = 3)
plot(hyp_trees[[2]], cex = 2, edge.width = 3)
plot(hyp_trees[[3]], cex = 2, edge.width = 3)
plot(hyp_trees[[4]], cex = 2, edge.width = 3)
plot(hyp_trees[[5]], cex = 2, edge.width = 3)
plot(hyp_trees[[6]], cex = 2, edge.width = 3)
plot(hyp_trees[[7]], cex = 2, edge.width = 3)
plot(hyp_trees[[8]], cex = 2, edge.width = 3)
plot(hyp_trees[[9]], cex = 2, edge.width = 3)

