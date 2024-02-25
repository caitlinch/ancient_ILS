# ancient_ILS/code/05_plot_figures_genes.R
## This script plots result figures from simulation output
# Caitlin Cherryh, 2023

###### 1. Input parameters ######
## File paths
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# output_dir                <- output directory to save figures

repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/06_figures/"

control_parameters <- list(add.extra.color.palettes = FALSE)


###### 2. Open packages and functions ######
## Open packages
library(reshape2)
library(ggplot2)
library(ggpubr)
library(patchwork)



###### 3. Prepare csvs for plotting  ######




# Save plot
gcf_plot_name <- paste0(repo_dir, "figures/", "empirical_gcf_scores.pdf")
ggsave(filename = gcf_plot_name, plot = gcf_plot, width = 8, height = 10, units = "in")