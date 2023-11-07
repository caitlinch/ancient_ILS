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
  output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/04_figures/"
}



###### 2. Open packages and functions ######
## Open packages
library(reshape2)
library(ggplot2)



###### 3. Prepare csvs for plotting  ######
# Identify output csv files
empirical_cf_path <- paste0(repo_dir, "output/empirical_concordance_factors.csv")
empirical_df <- read.csv(empirical_cf_path, stringsAsFactors = FALSE)
# Identify id.variable columns
id.var_cols <- c("dataset", "matrix", "hypothesis_tree", "branch_description", "branch_note")

# Specify colour palettes used within these plots
metazoan_clade_palette <- c(Bilateria = "#CC79A7", Cnidaria = "#009E73", Ctenophora = "#56B4E9", Porifera = "#E69F00", Outgroup = "#999999")
cividis5 <- c("#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")

# Extra colour palettes (unused)
if (control_parameters$add.extra.color.palettes == TRUE){
  qcf_type_palette <- c(Actual = "#5ab4ac", Estimated = "#d8b365")
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  tree2_palette <- c("#F0F921FF", "#0D0887FF")
  tree5_palette <- c("#e66101", "#fdb863", "#f7f7f7", "#b2abd2", "#5e3c99")
  tree2_cividis <- c(tree5_cividis[1], tree5_cividis[5])
  tree2_tonal <- c("#bdd7e7", "#2171b5")
  model3_tonal <- c("#980043", "#df65b0", "#d4b9da")
}



###### 4. Plot gCF  ######
# Extract gCF scores into long dataframe
gcf_df <- melt(empirical_df,
               id.vars = id.var_cols,
               measure.vars = c("gCF", "gDF1", "gDF2", "gDFP", "gcf_Label"))
# Extract gCF_N values into long dataframe
gcf_N_df <- melt(empirical_df,
                 id.vars = id.var_cols,
                 measure.vars = c("gCF_N", "gDF1_N", "gDF2_N", "gDFP_N"))



###### 5. Plot sCF  ######
# Extract gCF scores into long dataframe
scf_df <- melt(empirical_df,
               id.vars = id.var_cols,
               measure.vars = c("sCF", "sDF1", "sDF2", "sCF_Label"))
# Extract gCF_N values into long dataframe
scf_N_df <- melt(empirical_df,
                 id.vars = id.var_cols,
                 measure.vars = c("sCF_N", "sDF1_N", "sDF2_N", "sN"))



###### 6. Plot quartet scores  ######
# Extract quartet scores into long dataframe
qs_df <- melt(empirical_df,
              id.vars = id.var_cols,
              measure.vars = c("quartet_score"))


