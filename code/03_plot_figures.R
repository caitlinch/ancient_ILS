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



###### 4. Plot qCF results from "ancientILS" simulations ######
# Open simulation results
qcf_df_file <- paste0(input_dir, grep("ancientILS", grep("output_analysis", all_csvs, value = T), value = T))
qcf_df <- read.csv(qcf_df_file)
qcf_df$hypothesis_tree <- as.factor(qcf_df$hypothesis_tree)

#### Plot 1: Branch length vs summary qCF. ####
# Create df
plot1_df <- melt(data = qcf_df,
                 id.vars = c(id.var_cols),
                 measure.vars = c("actual_final_normalised_quartet_score", "estimated_final_normalised_quartet_score") )
# Create labs
facet_labs <- c("Ctenophora-sister", "Porifera-sister")
names(facet_labs) = c("1", "2")
# Plots
ggplot(plot1_df, aes(x = branch_a_length, y = value, color = variable)) + 
  facet_grid(hypothesis_tree~., labeller = labeller(hypothesis_tree = facet_labs)) +
  geom_point() +
  geom_vline(xintercept = 0.1729) +
  annotate("text", x = 0.19, y = 0.06, label = "Empirical\nASTRAL\nvalue", color = "Black",
           hjust = 0, vjust = 0.5, size = 4) +
  scale_x_continuous(name = "\nLength of branch a (in coalescent units)", trans='log10') +
  scale_y_continuous(name = "Normalised Final qCF Score\n", limits = c(0,1), breaks = seq(0,1.1,0.1)) +
  scale_color_manual(values = c("#5ab4ac", "#d8b365"),
                     labels = c("Actual", "Estimated"),
                     na.value = "grey50") +
  guides(color = guide_legend(title="qCF score")) +
  theme_bw() +
  theme(strip.text = element_text(size = 15))

ggplot(plot1_df, aes(x = branch_b_length, y = value, color = variable)) + 
  facet_grid(hypothesis_tree~., labeller = labeller(hypothesis_tree = facet_labs)) +
  geom_point() +
  geom_vline(xintercept = 1.6470) +
  annotate("text", x = 1.8, y = 0.06, label = "Empirical\nASTRAL\nvalue", color = "Black",
           hjust = 0, vjust = 0.5, size = 4) +
  scale_x_continuous(name = "\nLength of branch b (in coalescent units)", trans='log10') +
  scale_y_continuous(name = "Normalised Final qCF Score\n", limits = c(0,1), breaks = seq(0,1.1,0.1)) +
  theme_bw() +
  scale_color_manual(values = c("#5ab4ac", "#d8b365"),
                     labels = c("Actual", "Estimated"),
                     na.value = "grey50") +
  guides(color = guide_legend(title="qCF score")) +
  theme(strip.text = element_text(size = 15))


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




###### 4. Plot qCF results from "bothBranchesVary" simulations ######
# Open simulation results
qcf_df_file <- paste0(input_dir, grep("bothBranchesVary", grep("output_analysis", all_csvs, value = T), value = T))
qcf_df <- read.csv(qcf_df_file)

