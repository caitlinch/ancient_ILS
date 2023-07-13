# ancient_ILS/code/03_plot_trees.R
## This script plots and outputs phylogenetic tree figures
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# repo_dir            <- Location of caitlinch/ancient_ILS github repository

repo_dir <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"



#### 2. Open packages and prepare variables ####
# Open packages
library(ape)
library(phytools)
library(ggtree)
library(ggplot2)
library(patchwork)

# Create output directory for plots
plot_dir <- paste0(repo_dir, "figures/")



#### 3. Plot trees estimated from Whelan 2017 dataset ####
## Empirical ASTRAL tree
# Identify the ASTRAL tree
empirical_trees <- list.files(paste0(repo_dir, "empirical_tree"))
astral_tree <- paste0(repo_dir, "empirical_tree/", grep("Whelan2017_ASTRAL_tree.tre", empirical_trees, value = T))
# Open the ASTRAL tree
a_tree <- read.tree(astral_tree)
# Reroot tree
a_tree <- root(a_tree, outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), resolve.root = TRUE)
# Set missing branch lengths to 2
a_tree$edge.length[which(is.na(a_tree$edge.length))] <- 3
# Plot the ASTRAL tree

## Empirical IQ-Tree tree
# Identify the IQ-Tree tree
ml_tree <- paste0(repo_dir, "empirical_tree/", grep("Whelan2017_partitioned_ML_tree.treefile", empirical_trees, value = T))



#### 4. Plot phylogenetic hypotheses ####
# Open the possible phylogenetic topologies
hyp_trees <- read.tree(file = paste0(repo_dir, "hypothesis_trees/alternative_phylogenetic_hypotheses.nex"))
c_tree <- hyp_trees[[1]]
p_tree <- hyp_trees[[2]]
cp_tree <- hyp_trees[[3]]
# Plot Ctenophora-sister tree
c_plot <- ggtree(c_tree, size = 2) + 
  geom_rootedge(0.5, linewidth = 2) +
  geom_tiplab(size = 10, offset = 0.15) +
  xlim(-0.5,6.5)
# Plot Porifera-sister tree
p_plot <- ggtree(p_tree, size = 2) + 
  geom_rootedge(0.5, linewidth = 2) +
  geom_tiplab(size = 10, offset = 0.15) +
  xlim(-0.5,6.5)
# Plot Ctenophora+Porifera-sister tree
cp_plot <- ggtree(cp_tree, size = 2) + 
  geom_rootedge(0.5, linewidth = 2) +
  geom_tiplab(size = 10, offset = 0.15) +
  xlim(-0.5,5)
# Assemble plots using patchwork
patch <- c_plot + p_plot + cp_plot + plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
# Save plot
patch_path <- paste0(plot_dir, "methods_hypothesis_topologies.png")
png(file = patch_path, width = 1400, height = 500, units = "px")
patch
dev.off()



#### 5. Plot simulation hypothesis trees with branch length annotations ####
# Open the possible phylogenetic topologies
hyp_trees <- read.tree(file = paste0(repo_dir, "hypothesis_trees/alternative_phylogenetic_hypotheses.nex"))
c_tree <- hyp_trees[[7]]
p_tree <- hyp_trees[[8]]
cp_tree <- hyp_trees[[9]]
# Extend trees to be ultrametric
c_tree <- force.ultrametric(c_tree, method = "extend")
p_tree <- force.ultrametric(p_tree, method = "extend")
cp_tree <- force.ultrametric(cp_tree, method = "extend")

# Plot and annotate the Ctenophora-sister tree
c_plot <- ggtree(c_tree) +
  geom_rootedge(rootedge = 5) +
  geom_text2(aes(x = 0.33, y = 1.8), label = '0.6278', check_overlap = TRUE, color = 'grey50', size = 6) +
  geom_text2(aes(1.5, 3.8), label = '1.647', check_overlap = TRUE, color = 'grey50', size = 6) +
  geom_text2(aes(1.08, 6.7), label = '0.0853', check_overlap = TRUE, color = 'grey50', size = 6) +
  geom_text2(aes(1.6, 7.75), label = '0.737', check_overlap = TRUE, color = 'grey50', size = 6) +
  geom_text2(aes(1.7, 9.75), label = '0.9214', check_overlap = TRUE, color = 'grey50', size = 6) +
  geom_text2(aes(0.95, 8.75), label = '0.4145', check_overlap = TRUE, color = 'grey50', size = 6) +
  geom_text2(aes(0.5, 7.3), label = '0.1729', check_overlap = TRUE, color = 'grey50', size = 6) + 
  geom_text2(aes(0.33, 5.5), label = '0.6278', check_overlap = TRUE, color = 'grey50', size = 6) +
  geom_text2(aes(1.45, 4.2), label = 'Branch "b"', check_overlap = TRUE, color = 'darkgreen', size = 6) +
  geom_text2(aes(0.38, 7.7), label = 'Branch "a"', check_overlap = TRUE, color = 'darkgreen', size = 6) +
  geom_cladelabel(node = 19, label = "Bilateria", fontsize = 9, align = TRUE, geom = "text", color = c("white", "black")) +
  geom_cladelabel(node = 18, label = "Cnidaria", fontsize = 9, align = TRUE, geom = "text", color = c("white", "black")) +
  geom_cladelabel(node = 16, label = "Porifera", fontsize = 9, align = TRUE, geom = "text", color = c("white", "black")) +
  geom_cladelabel(node = 14, label = "Ctenophora", fontsize = 9, align = TRUE, geom = "text", color = c("white", "black")) +
  geom_cladelabel(node = 12, label = "Outgroup", fontsize = 9, align = TRUE, geom = "text", color = c("white", "black")) +
  geom_segment(aes(x = 0.84, y = 6.5, xend = 0.84, yend = 5.6), arrow = arrow(length = unit(0.15, "cm")), color = "grey70") +
  xlim(0,4.5)

# Plot and annotate the Porifera-first tree
p_plot <- ggtree(p_tree) + xlim(0,7) +
  geom_rootedge(rootedge = 5) + 
  geom_cladelabel(node = 19, label = "Bilateria", fontsize = 9, align = TRUE, geom = "text", color = c("white", "black")) +
  geom_cladelabel(node = 18, label = "Cnidaria", fontsize = 9, align = TRUE, geom = "text", color = c("white", "black")) +
  geom_cladelabel(node = 16, label = "Ctenophora", fontsize = 9, align = TRUE, geom = "text", color = c("white", "black")) +
  geom_cladelabel(node = 14, label = "Porifera", fontsize = 9, align = TRUE, geom = "text", color = c("white", "black")) +
  geom_cladelabel(node = 12, label = "Outgroup", fontsize = 9, align = TRUE, geom = "text", color = c("white", "black")) +
  geom_text2(aes(x = 2, y = 1.2), label = '0.6278', check_overlap = TRUE, color = 'grey50', size = 7) +
  geom_text2(aes(0.5, 4.95), label = '0.6278', check_overlap = TRUE, color = 'grey50', size = 7) +
  geom_text2(aes(1.5, 6.75), label = '0.1729', check_overlap = TRUE, color = 'grey50', size = 7) +
  geom_text2(aes(2.5, 8.25), label = '0.4145', check_overlap = TRUE, color = 'grey50', size = 7) +
  geom_text2(aes(2.4, 3.2), label = '0.0853', check_overlap = TRUE, color = 'grey50', size = 7) +
  geom_text2(aes(3.0, 5.20), label = '1.647', check_overlap = TRUE, color = 'grey50', size = 7) +
  geom_text2(aes(3.5, 7.20), label = '0.737', check_overlap = TRUE, color = 'grey50', size = 7) +
  geom_text2(aes(3.5, 9.20), label = '0.9214', check_overlap = TRUE, color = 'grey50', size = 7) +
  geom_text2(aes(2.5, 3.8), label = 'Branch "b"', check_overlap = TRUE, color = 'darkgreen', size = 7) +
  geom_text2(aes(1.2, 7.3), label = 'Branch "a"', check_overlap = TRUE, color = 'darkgreen', size = 7)

# Assemble c and p plots into a single plot
patch <- c_plot + p_plot + plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 30))

# Save plot
patch_path <- paste0(plot_dir, "methods_Simulation_topologies_BranchLengths.png")
png(file = patch_path, width = 1100, height = 500, units = "px")
patch
dev.off()

