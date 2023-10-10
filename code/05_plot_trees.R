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

# Assemble color palette
cols <- c("clade_labels" = "black", "clade_bars" = "white", "ctenophora_color" = "#E69F00", "porifera_color" = "#56B4E9",
          "branch_labels" = "#009E73", "branch_lengths" = "grey50", "arrows" = "grey70", "empirical_clade_labels" = "grey70")



#### 3. Plot trees estimated from Whelan 2017 dataset ####
## Empirical ASTRAL tree
# Identify the ASTRAL tree
empirical_trees <- list.files(paste0(repo_dir, "empirical_tree"))
astral_tree_path <- paste0(repo_dir, "empirical_tree/", grep("Whelan2017_ASTRAL_tree.tre", empirical_trees, value = T))
# Open the ASTRAL tree
a_tree <- read.tree(astral_tree_path)
# Reroot tree
a_tree <- root(a_tree, outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), resolve.root = TRUE)
# Set missing branch lengths to 2
a_tree$edge.length[which(is.na(a_tree$edge.length))] <- 1
# Rename tips
new_tips <- a_tree$tip.label
new_tips <- gsub("_", " ", new_tips)
new_tips <- gsub(" sp", " sp.", new_tips)
a_tree$tip.label <- new_tips
# Plot the ASTRAL tree
a_plot <- ggtree(a_tree, size = 1.2) +
  geom_rootedge(rootedge = 0.25, linewidth = 1.2) +
  geom_tiplab(size = 8, offset = 0.01) +
  geom_cladelabel(node = 92, label = "Bilateria", fontsize = 12, align = TRUE, geom = "text", color = cols[["empirical_clade_labels"]], offset = 0, angle = 90, hjust  = 0.5) +
  geom_cladelabel(node = 82, label = "Cnidaria", fontsize = 12, align = TRUE, geom = "text", color = cols[["empirical_clade_labels"]], offset = 0.5, angle = 90, hjust  = 0.5) +
  geom_cladelabel(node = 149, label = "Outgroup", fontsize = 12, align = TRUE, geom = "text", color = cols[["empirical_clade_labels"]], offset = 0, angle = 90, hjust  = 0.5) +
  geom_cladelabel(node = 98, label = "Porifera", fontsize = 12, align = TRUE, geom = "text", color = cols[["empirical_clade_labels"]], offset = 0, angle = 90, hjust  = 0.5) +
  geom_cladelabel(node = 120, label = "Ctenophora", fontsize = 12, align = TRUE, geom = "text", color = cols[["empirical_clade_labels"]], offset = 0, angle = 90, hjust  = 0.5)
# Save the ASTRAL tree
# Save plot
plot_path <- paste0(plot_dir, "methods_empirical_Whelan2017_ASTRAL_tree.png")
png(file = plot_path, width = 3500, height = 2100, units = "px")
a_plot
dev.off()

## Empirical IQ-Tree tree
# Identify the IQ-Tree tree
ml_tree_path <- paste0(repo_dir, "empirical_tree/", grep("Whelan2017_partitioned_ML_tree.treefile", empirical_trees, value = T))
# Open the IQ-Tree tree
ml_tree <- read.tree(ml_tree_path)
# Reroot tree
ml_tree <- root(ml_tree, outgroup = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"), resolve.root = TRUE)
# Set missing branch lengths to 2
ml_tree$edge.length[which(is.na(ml_tree$edge.length))] <- 1
# Rename tips
new_tips <- ml_tree$tip.label
new_tips <- gsub("_", " ", new_tips)
new_tips <- gsub(" sp", " sp.", new_tips)
ml_tree$tip.label <- new_tips
# Plot the IQ-Tree tree
ml_plot <- ggtree(ml_tree, size = 1.2) +
  geom_rootedge(rootedge = 0.05, linewidth = 1.2) +
  geom_tiplab(size = 8, offset = 0.01) +
  geom_cladelabel(node = 106, label = "Bilateria", fontsize = 12, align = TRUE, geom = "text", color = cols[["empirical_clade_labels"]], offset = 0.0, angle = 90, hjust  = 0.5) +
  geom_cladelabel(node = 92, label = "Cnidaria", fontsize = 12, align = TRUE, geom = "text", color = cols[["empirical_clade_labels"]], offset = 0.0, angle = 90, hjust  = 0.5) +
  geom_cladelabel(node = 129, label = "Outgroup", fontsize = 12, align = TRUE, geom = "text", color = cols[["empirical_clade_labels"]], offset = 0.0, angle = 90, hjust  = 0.5) +
  geom_cladelabel(node = 111, label = "Porifera", fontsize = 12, align = TRUE, geom = "text", color = cols[["empirical_clade_labels"]], offset = 0.05, angle = 90, hjust  = 0.5) +
  geom_cladelabel(node = 88, label = "Ctenophora", fontsize = 12, align = TRUE, geom = "text", color = cols[["empirical_clade_labels"]], offset = 0.15, angle = 90, hjust  = 0.5)
# Save the ASTRAL tree
# Save plot
plot_path <- paste0(plot_dir, "methods_empirical_Whelan2017_ML_partitioned_tree.png")
png(file = plot_path, width = 3000, height = 1700, units = "px")
ml_plot
dev.off()



#### 4. Plot phylogenetic hypotheses ####
# Open the possible phylogenetic topologies
hyp_trees <- read.tree(file = paste0(repo_dir, "hypothesis_trees/alternative_phylogenetic_hypotheses.nex"))
c_tree <- hyp_trees[[1]]
p_tree <- hyp_trees[[2]]
cp_tree <- hyp_trees[[3]]
# Plot Ctenophora-sister tree
c_plot <- ggtree(c_tree, size = 2) + 
  geom_rootedge(0.5, linewidth = 2) +
  geom_tiplab(aes(color = label), size = 10, offset = 0.15) +
  scale_color_manual(values=c(Bilateria = cols[["clade_labels"]], Cnidaria = cols[["clade_labels"]], 
                              Ctenophora = cols[["ctenophora_color"]], Outgroup = cols[["clade_labels"]],
                              Porifera = cols[["porifera_color"]])) +
  theme(legend.position = "none") +
  xlim(-0.5,6.5)
# Plot Porifera-sister tree
p_plot <- ggtree(p_tree, size = 2) + 
  geom_rootedge(0.5, linewidth = 2) +
  geom_tiplab(aes(color = label), size = 10, offset = 0.15) +
  scale_color_manual(values=c(Bilateria = cols[["clade_labels"]], Cnidaria = cols[["clade_labels"]], 
                              Ctenophora = cols[["ctenophora_color"]], Outgroup = cols[["clade_labels"]],
                              Porifera = cols[["porifera_color"]])) +
  theme(legend.position = "none") +
  xlim(-0.5,6.5)
# Plot Ctenophora+Porifera-sister tree
cp_plot <- ggtree(cp_tree, size = 2) + 
  geom_rootedge(0.5, linewidth = 2) +
  geom_tiplab(aes(color = label), size = 10, offset = 0.15) +
  scale_color_manual(values=c(Bilateria = cols[["clade_labels"]], Cnidaria = cols[["clade_labels"]], 
                              Ctenophora = cols[["ctenophora_color"]], Outgroup = cols[["clade_labels"]],
                              Porifera = cols[["porifera_color"]])) +
  theme(legend.position = "none") +
  xlim(-0.5,5)
# Assemble plots using patchwork
patch <- c_plot + p_plot + cp_plot + plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 45))
# Save plot
patch_path <- paste0(plot_dir, "methods_hypothesis_topologies.png")
png(file = patch_path, width = 1500, height = 500, units = "px")
patch
dev.off()



#### 5. Plot simulation hypothesis trees with branch length annotations, proportional branch lengths ####
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
c_plot <- ggtree(c_tree, size = 1) +
  geom_rootedge(rootedge = 0.1, linewidth = 1.1) +
  geom_cladelabel(node = 19, label = "Bilateria", fontsize = 10, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["clade_labels"]])) +
  geom_cladelabel(node = 18, label = "Cnidaria", fontsize = 10, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["clade_labels"]])) +
  geom_cladelabel(node = 14, label = "Ctenophora", fontsize = 10, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["ctenophora_color"]])) +
  geom_cladelabel(node = 16, label = "Porifera", fontsize = 10, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["porifera_color"]])) +
  geom_cladelabel(node = 12, label = "Outgroup", fontsize = 10, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["clade_labels"]])) +
  geom_text2(aes(x = 0.27, y = 7.65), label = 'Branch "a"', check_overlap = TRUE, color = cols[["branch_labels"]], size = 8) +
  geom_text2(aes(x = 1.45, y = 4.2), label = 'Branch "b"', check_overlap = TRUE, color = cols[["branch_labels"]], size = 8) +
  geom_segment(aes(x = 0.84, y = 6.25, xend = 0.84, yend = 5.6), arrow = arrow(length = unit(0.15, "cm")), color = cols[["arrows"]]) +
  geom_text2(aes(x = 0.265, y = 1.25), label = '0.6278', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 0.265, y = 5.5), label = '0.6278', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 1.45, y = 3.8), label = '1.647', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 1.14, y = 6.5), label = '0.0853', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 1.58, y = 7.75), label = '0.7370', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 1.67, y = 9.75), label = '0.9214', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 0.86, y = 8.75), label = '0.4145', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 0.444, y = 7.25), label = '0.1729', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  xlim(-0.2,4)

# Plot and annotate the Porifera-first tree
p_plot <- ggtree(p_tree, size = 1.1) +
  geom_rootedge(rootedge = 0.1, linewidth = 1.1) + 
  geom_cladelabel(node = 19, label = "Bilateria", fontsize = 10, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["clade_labels"]])) +
  geom_cladelabel(node = 18, label = "Cnidaria", fontsize = 10, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["clade_labels"]])) +
  geom_cladelabel(node = 16, label = "Ctenophora", fontsize = 10, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["ctenophora_color"]])) +
  geom_cladelabel(node = 14, label = "Porifera", fontsize = 10, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["porifera_color"]])) +
  geom_cladelabel(node = 12, label = "Outgroup", fontsize = 10, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["clade_labels"]])) +
  geom_text2(aes(x = 0.27, y = 7.65), label = 'Branch "a"', check_overlap = TRUE, color = cols[["branch_labels"]], size = 8) +
  geom_text2(aes(x = 1.65, y = 6.2), label = 'Branch "b"', check_overlap = TRUE, color = cols[["branch_labels"]], size = 8) +
  geom_segment(aes(x = 0.67, y = 4.2, xend = 0.67, yend = 3.55), arrow = arrow(length = unit(0.15, "cm")), color = cols[["arrows"]]) +
  geom_text2(aes(x = 0.265, y = 1.25), label = '0.6278', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 0.265, y = 5.5), label = '0.6278', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 1.6, y = 5.8), label = '1.647', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 0.97, y = 4.45), label = '0.0853', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 1.58, y = 7.75), label = '0.7370', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 1.67, y = 9.75), label = '0.9214', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 0.86, y = 8.75), label = '0.4145', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  geom_text2(aes(x = 0.444, y = 7.25), label = '0.1729', check_overlap = TRUE, color = cols[["branch_lengths"]], size = 8) +
  xlim(-0.2,4)

# Assemble c and p plots into a single plot
patch <- c_plot + p_plot + plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 40))

# Save plot
patch_path <- paste0(plot_dir, "methods_Simulation_topologies_BranchLengths_proportional.png")
png(file = patch_path, width = 1100, height = 600, units = "px")
patch
dev.off()



#### 5. Plot simulation hypothesis trees with branch length annotations, cladograms ####
# Open the possible phylogenetic topologies
hyp_trees <- read.tree(file = paste0(repo_dir, "hypothesis_trees/alternative_phylogenetic_hypotheses.nex"))
c_tree <- hyp_trees[[4]]
p_tree <- hyp_trees[[5]]
cp_tree <- hyp_trees[[6]]

# Plot and annotate the Ctenophora-sister tree
c_plot <- ggtree(c_tree, size = 1.1) +
  geom_rootedge(rootedge = 0.5) +
  geom_text2(aes(1, 7.4), label = 'Branch "a"', check_overlap = TRUE, color = cols[["branch_labels"]], size = 9) +
  geom_text2(aes(2.55, 3.9), label = 'Branch "b"', check_overlap = TRUE, color = cols[["branch_labels"]], size = 9) +
  geom_cladelabel(node = 19, label = "Bilateria", fontsize = 9, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["clade_labels"]])) +
  geom_cladelabel(node = 18, label = "Cnidaria", fontsize = 9, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["clade_labels"]])) +
  geom_cladelabel(node = 16, label = "Porifera", fontsize = 9, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["porifera_color"]])) +
  geom_cladelabel(node = 14, label = "Ctenophora", fontsize = 9, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["ctenophora_color"]])) +
  geom_cladelabel(node = 12, label = "Outgroup", fontsize = 9, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["clade_labels"]])) +
  xlim(-0.50,7)

# Plot and annotate the Porifera-first tree
p_plot <- ggtree(p_tree, size = 1.1) +
  geom_rootedge(rootedge = 0.5) +
  geom_text2(aes(1.0, 7.35), label = 'Branch "a"', check_overlap = TRUE, color = cols[["branch_labels"]], size = 9) +
  geom_text2(aes(3, 5.85), label = 'Branch "b"', check_overlap = TRUE, color = cols[["branch_labels"]], size = 9) +
  geom_cladelabel(node = 19, label = "Bilateria", fontsize = 9, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["clade_labels"]])) +
  geom_cladelabel(node = 18, label = "Cnidaria", fontsize = 9, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["clade_labels"]])) +
  geom_cladelabel(node = 16, label = "Ctenophora", fontsize = 9, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["ctenophora_color"]])) +
  geom_cladelabel(node = 14, label = "Porifera", fontsize = 9, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["porifera_color"]])) +
  geom_cladelabel(node = 12, label = "Outgroup", fontsize = 9, align = TRUE, geom = "text", color = c(cols[["clade_bars"]], cols[["clade_labels"]])) +
  xlim(-0.5,7)

# Assemble c and p plots into a single plot
patch <- c_plot + p_plot + plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 50))

# Save plot
patch_path <- paste0(plot_dir, "methods_Simulation_topologies_BranchModificationLabels.png")
png(file = patch_path, width = 1100, height = 500, units = "px")
patch
dev.off()
