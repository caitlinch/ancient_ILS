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

# Set color palettes
qcf_type_palette <- c("#5ab4ac", "#d8b365")
names(qcf_type_palette) = c("Actual", "Estimated")

# Set empirical branch lengths for branches 
emp_bl <- c("a" = 0.1729, "b" = 1.6470)




###### 4. Plot qCF results from "ancientILS" simulations ######
# Open simulation results
qcf_df_file <- paste0(input_dir, grep("ancientILS", grep("output_analysis", all_csvs, value = T), value = T))
qcf_df <- read.csv(qcf_df_file)
qcf_df$hypothesis_tree <- as.factor(qcf_df$hypothesis_tree)
# Create labs
tree_labs <- c("Ctenophora-sister", "Porifera-sister")
names(tree_labs) = c("1", "2")
hyp_labs <- c("Branch to\nBILAT+CNID+PORI", "Branch to\nBILAT+CNID+CTEN", "Branch to\nCTEN+PORI")
names(hyp_labs) <- c("Hyp1", "Hyp2", "Hyp3")



#### Plot 1: Branch length vs summary qCF. ####
# Create df
plot1_df <- melt(data = qcf_df,
                 id.vars = c(id.var_cols),
                 measure.vars = c("actual_final_normalised_quartet_score", "estimated_final_normalised_quartet_score") )
# Facet variable for plotting 
plot1_df$qcf_type <- factor(plot1_df$variable,
                            levels = c("actual_final_normalised_quartet_score", "estimated_final_normalised_quartet_score"), 
                            labels = c("Actual", "Estimated") )
plot1_df$qcf_type <- factor(plot1_df$variable,
                            levels = c("actual_final_normalised_quartet_score", "estimated_final_normalised_quartet_score"), 
                            labels = c("Actual", "Estimated") )
# Add new columns for plotting
plot1_df$branch_to_vary <- plot1_df$ID
plot1_df$branch_to_vary[which(plot1_df$branch_a_length == emp_bl["a"])] <- "Branch b"
plot1_df$branch_to_vary[which(plot1_df$branch_b_length == emp_bl["b"])] <- "Branch a"
plot1_df$varied_branch_length <- plot1_df$ID
plot1_df$varied_branch_length[which(plot1_df$branch_a_length == emp_bl["a"])] <- plot1_df$branch_b_length[which(plot1_df$branch_a_length == emp_bl["a"])]
plot1_df$varied_branch_length[which(plot1_df$branch_b_length == emp_bl["b"])] <- plot1_df$branch_a_length[which(plot1_df$branch_b_length == emp_bl["b"])]
plot1_df$varied_branch_length <- factor(plot1_df$varied_branch_length,
                                        levels = c("1e-08", "1e-07", "1e-06", "1e-05", "1e-04", "0.001", "0.01", "0.1", "1", "10"),
                                        labels = c("1e-08", "1e-07", "1e-06", "1e-05", "1e-04", "1e-03", "0.01", "0.1", "1", "10"),
                                        ordered = T)

# Plots
p <- ggplot(plot1_df, aes(x = varied_branch_length, y = value, color = qcf_type)) + 
  facet_grid(hypothesis_tree~branch_to_vary, labeller = labeller(hypothesis_tree = tree_labs)) +
  geom_boxplot() +
  scale_x_discrete(name = "\nLength of branch (in coalescent units)") +
  scale_y_continuous(name = "Normalised Final qCF Score\n", limits = c(0,1), breaks = seq(0,1.1,0.1)) +
  scale_color_manual(values = qcf_type_palette, na.value = "grey50") +
  guides(color = guide_legend(title="qCF type")) +
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_text(size = 30),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 25))
p_file <- paste0(output_dir, "exploratory_normalised-qcf_branch-length.png")
ggsave(filename = p_file, plot = p, device = "png")



#### Plot 2: Branch length vs branch "a" qCF. ####
plot2_df <- melt(data = qcf_df,
                 id.vars = c(id.var_cols),
                 measure.vars = c("actual_testHyp1_Cten_branch_a_qcf_value", "actual_testHyp2_Pori_branch_a_qcf_value", "actual_testHyp3_CtenPori_qcf_value",
                                  "estimated_testHyp1_Cten_branch_a_qcf_value", "estimated_testHyp2_Pori_branch_a_qcf_value", "estimated_testHyp3_CtenPori_qcf_value"))
# Remove any rows with branch a at the empirical branch length (keep only branch b stable for this plot)
plot2_df <- plot2_df[(plot2_df$branch_a_length != emp_bl[["a"]]),]
# Add extract columns for faceting
plot2_df$qcf_type <- factor(plot2_df$variable,
                            levels = c("actual_testHyp1_Cten_branch_a_qcf_value", "estimated_testHyp1_Cten_branch_a_qcf_value", 
                                       "actual_testHyp2_Pori_branch_a_qcf_value", "estimated_testHyp2_Pori_branch_a_qcf_value", 
                                       "actual_testHyp3_CtenPori_qcf_value", "estimated_testHyp3_CtenPori_qcf_value"), 
                            labels = rep(c("Actual", "Estimated"), 3) )
plot2_df$hypothesis_test <- factor(plot2_df$variable,
                                   levels = c("actual_testHyp1_Cten_branch_a_qcf_value", "estimated_testHyp1_Cten_branch_a_qcf_value", 
                                              "actual_testHyp2_Pori_branch_a_qcf_value", "estimated_testHyp2_Pori_branch_a_qcf_value", 
                                              "actual_testHyp3_CtenPori_qcf_value", "estimated_testHyp3_CtenPori_qcf_value"), 
                                   labels = c("Hyp1", "Hyp1", "Hyp2", "Hyp2", "Hyp3", "Hyp3") )

p <- ggplot(plot2_df, aes(x = branch_a_length, y = value, color = qcf_type)) + 
  facet_grid(hypothesis_test~hypothesis_tree, labeller = labeller(hypothesis_tree = tree_labs, hypothesis_test = hyp_labs)) +
  geom_smooth() +
  geom_vline(xintercept = emp_bl[["a"]], linetype = 2, color = "darkgrey") +
  scale_x_continuous(name = "\nLength of branch a (in coalescent units)", trans = "log10") +
  scale_y_continuous(name = "qCF Score\n", limits = c(0,1.05), breaks = seq(0,1.1,0.2)) +
  scale_color_manual(values = qcf_type_palette, na.value = "grey50") +
  guides(color = guide_legend(title="qCF type")) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
p_file <- paste0(output_dir, "results_qcf-score_branch-a-length.png")
ggsave(filename = p_file, plot = p, device = "png")



#### Plot 3: Branch length vs branch "b" qCF. ####
plot3_df <- melt(data = qcf_df,
                 id.vars = c(id.var_cols),
                 measure.vars = c("actual_testHyp1_Cten_branch_a_qcf_value", "actual_testHyp2_Pori_branch_a_qcf_value", "actual_testHyp3_CtenPori_qcf_value",
                                  "estimated_testHyp1_Cten_branch_a_qcf_value", "estimated_testHyp2_Pori_branch_a_qcf_value", "estimated_testHyp3_CtenPori_qcf_value"))
# Remove any rows with branch b at the empirical branch length (keep only branch a stable for this plot)
plot3_df <- plot3_df[(plot3_df$branch_b_length != emp_bl[["b"]]),]
# Add extract columns for faceting
plot3_df$qcf_type <- factor(plot3_df$variable,
                            levels = c("actual_testHyp1_Cten_branch_a_qcf_value", "estimated_testHyp1_Cten_branch_a_qcf_value", 
                                       "actual_testHyp2_Pori_branch_a_qcf_value", "estimated_testHyp2_Pori_branch_a_qcf_value", 
                                       "actual_testHyp3_CtenPori_qcf_value", "estimated_testHyp3_CtenPori_qcf_value"), 
                            labels = rep(c("Actual", "Estimated"), 3) )
plot3_df$hypothesis_test <- factor(plot3_df$variable,
                                   levels = c("actual_testHyp1_Cten_branch_a_qcf_value", "estimated_testHyp1_Cten_branch_a_qcf_value", 
                                              "actual_testHyp2_Pori_branch_a_qcf_value", "estimated_testHyp2_Pori_branch_a_qcf_value", 
                                              "actual_testHyp3_CtenPori_qcf_value", "estimated_testHyp3_CtenPori_qcf_value"), 
                                   labels = c("Hyp1", "Hyp1", "Hyp2", "Hyp2", "Hyp3", "Hyp3") )

p <- ggplot(plot3_df, aes(x = branch_b_length, y = value, color = qcf_type)) + 
  facet_grid(hypothesis_test~hypothesis_tree, labeller = labeller(hypothesis_tree = tree_labs, hypothesis_test = hyp_labs)) +
  geom_smooth() +
  geom_vline(xintercept = emp_bl[["b"]], linetype = 2, color = "darkgrey") +
  scale_x_continuous(name = "\nLength of branch b (in coalescent units)", trans='log10') +
  scale_y_continuous(name = "qCF Score\n", limits = c(0,1.05), breaks = seq(0,1.1,0.2)) +
  scale_color_manual(values = qcf_type_palette, na.value = "grey50") +
  guides(color = guide_legend(title="qCF type")) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
p_file <- paste0(output_dir, "results_qcf-score_branch-b-length.png")
ggsave(filename = p_file, plot = p, device = "png")



#### Plot 4: Branch a length vs clade qCF. ####
plot4_df <- melt(data = qcf_df,
                 id.vars = c(id.var_cols),
                 measure.vars = c("actual_clade_BILAT_qcf_value", "actual_clade_CNID_qcf_value", "actual_clade_PORI_qcf_value",
                                  "actual_clade_CTEN_qcf_value", "actual_clade_CHOANO_qcf_value", "actual_clade_BILAT_CTEN_qcf_value",
                                  "actual_clade_ANIMALS_qcf_value", "estimated_clade_BILAT_qcf_value", "estimated_clade_CNID_qcf_value",
                                  "estimated_clade_PORI_qcf_value", "estimated_clade_CTEN_qcf_value", "estimated_clade_CHOANO_qcf_value",
                                  "estimated_clade_BILAT_CTEN_qcf_value", "estimated_clade_ANIMALS_qcf_value"))
# Remove any rows with branch a at the empirical branch length (keep only branch b stable for this plot)
plot4_df <- plot4_df[(plot4_df$branch_a_length != emp_bl[["a"]]),]
# Format columns as factors for nice plotting
plot4_df$clade <- factor(plot4_df$variable,
                         levels = rev(c("actual_clade_CHOANO_qcf_value", "estimated_clade_CHOANO_qcf_value", 
                                        "actual_clade_ANIMALS_qcf_value", "estimated_clade_ANIMALS_qcf_value", 
                                        "actual_clade_BILAT_CTEN_qcf_value", "estimated_clade_BILAT_CTEN_qcf_value",
                                        "actual_clade_BILAT_qcf_value", "estimated_clade_BILAT_qcf_value", 
                                        "actual_clade_CNID_qcf_value", "estimated_clade_CNID_qcf_value", 
                                        "actual_clade_PORI_qcf_value", "estimated_clade_PORI_qcf_value", 
                                        "actual_clade_CTEN_qcf_value", "estimated_clade_CTEN_qcf_value")), 
                         labels = rev(c("Choanoflagellates", "Choanoflagellates", 
                                        "All Animals", "All Animals", 
                                        "Bilateria+Cnidaria", "Bilateria+Cnidaria",
                                        "Bilateria", "Bilateria", 
                                        "Cnidaria", "Cnidaria",
                                        "Porifera", "Porifera", 
                                        "Ctenophora", "Ctenophora")) )
plot4_df$qcf_type <- factor(plot4_df$variable,
                            levels = c("actual_clade_CHOANO_qcf_value", "estimated_clade_CHOANO_qcf_value", 
                                       "actual_clade_ANIMALS_qcf_value", "estimated_clade_ANIMALS_qcf_value", 
                                       "actual_clade_BILAT_CTEN_qcf_value", "estimated_clade_BILAT_CTEN_qcf_value",
                                       "actual_clade_BILAT_qcf_value", "estimated_clade_BILAT_qcf_value", 
                                       "actual_clade_CNID_qcf_value", "estimated_clade_CNID_qcf_value", 
                                       "actual_clade_PORI_qcf_value", "estimated_clade_PORI_qcf_value", 
                                       "actual_clade_CTEN_qcf_value", "estimated_clade_CTEN_qcf_value"), 
                            labels = rep(c("Actual", "Estimated"), 7) )
plot4_df$branch_a_length <- factor(plot4_df$branch_a_length,
                                   levels = c("1e-08", "1e-07", "1e-06", "1e-05", "1e-04", "0.001", "0.01", "0.1", "1", "10"),
                                   labels = c("1e-08", "1e-07", "1e-06", "1e-05", "1e-04", "1e-03", "0.01", "0.1", "1", "10"),
                                   ordered = T)

p <- ggplot(plot4_df, aes(x = branch_a_length, y = value, color = qcf_type)) + 
  facet_wrap("clade") +
  geom_boxplot() +
  scale_x_discrete(name = "\nLength of branch a (in coalescent units)") +
  scale_y_continuous(name = "qCF Score\n", limits = c(0,1.2), breaks = seq(0,1.6,0.2)) +
  scale_color_manual(values = qcf_type_palette, na.value = "grey50") +
  guides(color = guide_legend(title="qCF type")) +
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 30),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 15)),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 25))
p_file <- paste0(output_dir, "results_clade_branch-a-length.png")
ggsave(filename = p_file, plot = p, device = "png")



#### Plot 5: Branch b length vs clade qCF. ####
plot5_df <- melt(data = qcf_df,
                 id.vars = c(id.var_cols),
                 measure.vars = c("actual_clade_BILAT_qcf_value", "actual_clade_CNID_qcf_value", "actual_clade_PORI_qcf_value",
                                  "actual_clade_CTEN_qcf_value", "actual_clade_CHOANO_qcf_value", "actual_clade_BILAT_CTEN_qcf_value",
                                  "actual_clade_ANIMALS_qcf_value", "estimated_clade_BILAT_qcf_value", "estimated_clade_CNID_qcf_value",
                                  "estimated_clade_PORI_qcf_value", "estimated_clade_CTEN_qcf_value", "estimated_clade_CHOANO_qcf_value",
                                  "estimated_clade_BILAT_CTEN_qcf_value", "estimated_clade_ANIMALS_qcf_value"))
# Remove any rows with branch a at the empirical branch length (keep only branch b stable for this plot)
plot5_df <- plot5_df[(plot5_df$branch_b_length != emp_bl[["b"]]),]
# Format columns as factors for nice plotting
plot5_df$clade <- factor(plot5_df$variable,
                         levels = rev(c("actual_clade_CHOANO_qcf_value", "estimated_clade_CHOANO_qcf_value", 
                                        "actual_clade_ANIMALS_qcf_value", "estimated_clade_ANIMALS_qcf_value", 
                                        "actual_clade_BILAT_CTEN_qcf_value", "estimated_clade_BILAT_CTEN_qcf_value",
                                        "actual_clade_BILAT_qcf_value", "estimated_clade_BILAT_qcf_value", 
                                        "actual_clade_CNID_qcf_value", "estimated_clade_CNID_qcf_value", 
                                        "actual_clade_PORI_qcf_value", "estimated_clade_PORI_qcf_value", 
                                        "actual_clade_CTEN_qcf_value", "estimated_clade_CTEN_qcf_value")), 
                         labels = rev(c("Choanoflagellates", "Choanoflagellates", 
                                        "All Animals", "All Animals", 
                                        "Bilateria+Cnidaria", "Bilateria+Cnidaria",
                                        "Bilateria", "Bilateria", 
                                        "Cnidaria", "Cnidaria",
                                        "Porifera", "Porifera", 
                                        "Ctenophora", "Ctenophora")) )
plot5_df$qcf_type <- factor(plot5_df$variable,
                            levels = c("actual_clade_CHOANO_qcf_value", "estimated_clade_CHOANO_qcf_value", 
                                       "actual_clade_ANIMALS_qcf_value", "estimated_clade_ANIMALS_qcf_value", 
                                       "actual_clade_BILAT_CTEN_qcf_value", "estimated_clade_BILAT_CTEN_qcf_value",
                                       "actual_clade_BILAT_qcf_value", "estimated_clade_BILAT_qcf_value", 
                                       "actual_clade_CNID_qcf_value", "estimated_clade_CNID_qcf_value", 
                                       "actual_clade_PORI_qcf_value", "estimated_clade_PORI_qcf_value", 
                                       "actual_clade_CTEN_qcf_value", "estimated_clade_CTEN_qcf_value"), 
                            labels = rep(c("Actual", "Estimated"), 7) )
plot5_df$branch_a_length <- factor(plot5_df$branch_b_length,
                                   levels = c("1e-08", "1e-07", "1e-06", "1e-05", "1e-04", "0.001", "0.01", "0.1", "1", "10"),
                                   labels = c("1e-08", "1e-07", "1e-06", "1e-05", "1e-04", "1e-03", "0.01", "0.1", "1", "10"),
                                   ordered = T)

p <- ggplot(plot5_df, aes(x = branch_b_length, y = value, color = qcf_type)) + 
  facet_wrap("clade") +
  geom_boxplot() +
  scale_x_discrete(name = "\nLength of branch b (in coalescent units)") +
  scale_y_continuous(name = "qCF Score\n", limits = c(0,1.2), breaks = seq(0,1.6,0.2)) +
  scale_color_manual(values = qcf_type_palette, na.value = "grey50") +
  guides(color = guide_legend(title="qCF type")) +
  theme_bw() +
  theme(strip.text = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 30),
        axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 15, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 15)),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 25))
p_file <- paste0(output_dir, "results_clade_branch-b-length.png")
ggsave(filename = p_file, plot = p, device = "png")






###### 4. Plot qCF results from "bothBranchesVary" simulations ######
# Open simulation results
qcf_df_file <- paste0(input_dir, grep("bothBranchesVary", grep("output_analysis", all_csvs, value = T), value = T))
qcf_df <- read.csv(qcf_df_file)

