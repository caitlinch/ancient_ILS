# ancient_ILS/code/05_plot_trees.R
## This script plots and outputs phylogenetic tree figures
# Caitlin Cherryh 2024

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
library(grDevices) # for: hcl.colors(), palette.colors()

## Specify colour palettes used within these plots
# Inbuilt colour palettes
okito_palette           <- palette.colors(palette = "Okabe-Ito")
viridis_palette         <- hcl.colors(n = 8, palette = "Viridis")
mako_palette            <- hcl.colors(n = 5, palette = "Mako")
# Colour palettes for variables
metazoan_clade_palette  <- c("Bilateria" = "#CC79A7", "Cnidaria" = "#009E73", "Ctenophora" = "#56B4E9",
                             "Porifera" = "#E69F00", "Outgroup" = "#999999", "Placozoa" = "#000000")
topology_colours        <- c("Ctenophora" = viridis_palette[1], "Porifera" =  viridis_palette[4], 
                             "Ctenophora+Porifera" = viridis_palette[7], "Paraphyly" = "#999999")
cols <- c("clade_labels" = "black", "clade_bars" = "white", "ctenophora_color" = "#E69F00", "porifera_color" = "#56B4E9",
          "branch_labels" = "#009E73", "branch_lengths" = "grey50", "arrows" = "grey70", "empirical_clade_labels" = "grey70")



#### 3. Plot phylogenetic hypotheses ####
# Open the possible phylogenetic topologies
hyp_trees <- read.tree(file = paste0(repo_dir, "constraint_trees/alternative_phylogenetic_hypotheses.nex"))
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


