# ancient_ILS/code/02_plot_figures.R
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
  input_dir                   <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_simulations/"
  output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/04_figures/"
}



###### 2. Open packages and functions ######
## Source functions


## Open packages
library(reshape2)
library(ggplot2)

# Create output file paths




###### 3. Prepare csvs for plotting  ######
# Identify output csv files
all_csvs <- grep("csv", list.files(input_dir), value = T)
# Open simulation results
qcf_df_file <- paste0(input_dir, grep("output_analysis", all_csvs, value = T))
qcf_df <- read.csv(qcf_df_file)
# Add some columns to the dataframe
#     Denote estimated or actual qcf
qcf_df$qcf_type <- as.character(qcf_df$variable)
qcf_df$qcf_type[grep("estimated", as.character(qcf_df$variable))] <- "estimated"
qcf_df$qcf_type[grep("actual", as.character(qcf_df$variable))] <- "actual"
qcf_df$hypothesis_test <- as.character(qcf_df$variable)
#     Denote tree used to calculate qCFs
qcf_df$hypothesis_test[grep("testHyp1_Cten", as.character(qcf_df$hypothesis_test))] <- "Hyp1"
qcf_df$hypothesis_test[grep("testHyp2_Pori", as.character(qcf_df$hypothesis_test))] <- "Hyp2"
qcf_df$hypothesis_test[grep("testHyp3_CtenPori", as.character(qcf_df$hypothesis_test))] <- "Hyp3"
# Identify id.variable columns
id.var_cols <- c("ID", "simulation_type", "hypothesis_tree", "branch_a_length", "branch_b_length", "branch_c_length", "qcf_type", "hypothesis_test")



###### 4. Plot qCF results  ######
#### Plot 1: Branch length vs summary qCF. ####
plot1_df <- melt(data = qcf_df,
                 id.vars = c(id.var_cols),
                 measure.vars = c("actual_num_quartet_trees", "actual_final_quartet_score", "actual_final_normalised_quartet_score",
                                  "actual_qcf_mean", "actual_qcf_median", "actual_qcf_min", "actual_qcf_max", "estimated_num_quartet_trees",
                                  "estimated_final_quartet_score", "estimated_final_normalised_quartet_score", "estimated_qcf_mean",
                                  "estimated_qcf_median", "estimated_qcf_min","estimated_qcf_max"))
plot1a_df <- plot1_df[(plot1_df$variable == "actual_final_normalised_quartet_score" | plot1_df$variable == "estimated_final_normalised_quartet_score"), ]

ggplot(plot1a_df, aes(x = hypothesis_tree, x = value)) + geom_boxplot() + geom_jitter()

#### Plot 2: Branch length vs branch "a" qCF. ####
plot2_df <- melt(data = qcf_df,
                 id.vars = c(id.var_cols),
                 measure.vars = c("actual_testHyp1_Cten_branch_a_qcf_value", "actual_testHyp2_Pori_branch_a_qcf_value", "actual_testHyp3_CtenPori_qcf_value",
                                  "expected_testHyp1_Cten_branch_a_qcf_value", "expected_testHyp2_Pori_branch_a_qcf_value", "expected_testHyp3_CtenPori_qcf_value"))
plot2_df$qcf_type <- as.character(plot2_df$variable)
plot2_df$qcf_type[grep("expected", as.character(plot2_df$variable))] <- "estimated"
plot2_df$qcf_type[grep("actual", as.character(plot2_df$variable))] <- "actual"
plot2_df$hypothesis_test <- as.character(plot2_df$variable)
plot2_df$hypothesis_test[grep("testHyp1_Cten", as.character(plot2_df$hypothesis_test))] <- "Hyp1"
plot2_df$hypothesis_test[grep("testHyp2_Pori", as.character(plot2_df$hypothesis_test))] <- "Hyp2"
plot2_df$hypothesis_test[grep("testHyp3_CtenPori", as.character(plot2_df$hypothesis_test))] <- "Hyp3"

ggplot(plot2_df, aes(x = branch_a_length, y = value, group = hypothesis_tree)) + 
  facet_grid(hypothesis_test~qcf_type) +
  geom_boxplot() 


#### Plot 3: Branch length vs branch "b" qCF. ####
plot3_df <- melt(data = qcf_df,
                 id.vars = c(id.var_cols),
                 measure.vars = c("actual_clade_CTEN_qcf_value", "expected_clade_CTEN_qcf_value"))


#### Plot 4: Branch length vs clade qCF. ####
plot4_df <- melt(data = qcf_df,
                 id.vars = c(id.var_cols),
                 measure.vars = c("actual_clade_BILAT_qcf_value", "actual_clade_CNID_qcf_value", "actual_clade_PORI_qcf_value",
                                  "actual_clade_CTEN_qcf_value", "actual_clade_CHOANO_qcf_value", "actual_clade_BILAT_CTEN_qcf_value",
                                  "actual_clade_ANIMALS_qcf_value", "expected_clade_BILAT_qcf_value", "expected_clade_CNID_qcf_value",
                                  "expected_clade_PORI_qcf_value", "expected_clade_CTEN_qcf_value", "expected_clade_CHOANO_qcf_value",
                                  "expected_clade_BILAT_CTEN_qcf_value", "expected_clade_ANIMALS_qcf_value"))






