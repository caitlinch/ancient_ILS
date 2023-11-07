# ancient_ILS/code/05_plot_figures_empirical.R
## This script plots result figures from simulation output
# Caitlin Cherryh, 2023

###### 1. Input parameters ######
## File paths
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# input_dir                 <- location of simulation and analysis result dataframes
# output_dir                <- output directory to save figures

location = "local"
if (location == "local"){
  ## File paths
  repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  input_dir                   <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_simulation_output_files/"
  output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/04_figures/"
}



###### 2. Open packages and functions ######
## Source functions
source(paste0(repo_dir, "code/func_plotting.R"))

## Open packages
library(reshape2)
library(ggplot2)



###### 3. Prepare csvs for plotting  ######
# Identify output csv files
all_csvs <- grep("csv", list.files(input_dir), value = T)
# Identify id.variable columns
id.var_cols <- c("ID", "simulation_type", "hypothesis_tree", "branch_a_length", "branch_b_length", "branch_c_length")

# Set color palettes
qcf_type_palette <- c("#5ab4ac", "#d8b365")
names(qcf_type_palette) = c("Actual", "Estimated")

# Set empirical branch lengths for branches 
emp_bl <- c("a" = 0.1729, "b" = 1.6470)


