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
library(patchwork)



###### 3. Prepare csvs for plotting  ######
# Identify output csv files
empirical_cf_path <- paste0(repo_dir, "output/empirical_concordance_factors.csv")
empirical_df <- read.csv(empirical_cf_path, stringsAsFactors = FALSE)
# Remove unusual branaches
empirical_df <- empirical_df[(empirical_df$branch_description == "To_all_animals" | empirical_df$branch_description == "To_all_other_metazoans" |
                                empirical_df$branch_description == "To_CTEN_clade" | empirical_df$branch_description == "To_PORI_clade"), ]
# Identify id.variable columns
id.var_cols <- c("dataset", "matrix", "hypothesis_tree", "branch_description", "branch_note")

# Specify colour palettes used within these plots
metazoan_clade_palette <- c(Bilateria = "#CC79A7", Cnidaria = "#009E73", Ctenophora = "#56B4E9", Porifera = "#E69F00", Outgroup = "#999999")

# Extra colour palettes (unused)
if (control_parameters$add.extra.color.palettes == TRUE){
  cividis5 <- c("#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
  palette2 <- c("Ctenophora-sister" = "#fdb863", "Porifera-sister" = "#b2abd2")
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
               measure.vars = c("gCF", "gDF1", "gDF2", "gDFP"))
# Add labels for plotting
gcf_df$variable_label <- factor(gcf_df$variable,
                                levels = c("gCF", "gDF1", "gDF2", "gDFP"),
                                labels = c("gCF", "gDF1", "gDF2", "gDFP"),
                                ordered = T)
gcf_df$branch_label <- factor(gcf_df$branch_description,
                              levels = c("To_all_animals", "To_all_other_metazoans", "To_CTEN_clade", "To_PORI_clade"),
                              labels = c("All Animals", "Other Metazoans", "Ctenophora", "Porifera"),
                              ordered = T)
gcf_df$hypothesis_label <- factor(gcf_df$hypothesis_tree,
                                  levels = c("CTEN", "PORI"),
                                  labels = c("Ctenophora-sister", "Porifera-sister"),
                                  ordered = T)
# Plot gCF scores
gcf_plot <- ggplot(gcf_df, aes(x = variable_label, y = value, group = variable_label, fill = hypothesis_label)) +
  geom_boxplot(alpha = 1) +
  facet_grid(branch_label~hypothesis_label) +
  scale_x_discrete(name = "Gene concordance factor (gCF) statistics") +
  scale_y_continuous(name = "Value", breaks = seq(0,100,20), labels = seq(0,100,20), minor_breaks = seq(0,100,10)) +
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"), labels = c("Ctenophora-sister", "Porifera-sister")) +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0, unit = "pt")),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt")),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0),
        strip.text = element_text(size = 18),
        legend.position = "none" )
# Save plot
gcf_plot_name <- paste0(repo_dir, "figures/", "empirical_gcf_scores.pdf")
ggsave(filename = gcf_plot_name, plot = gcf_plot, width = 8, height = 10, units = "in")



###### 5. Plot sCF  ######
# Extract gCF scores into long dataframe
scf_df <- melt(empirical_df,
               id.vars = id.var_cols,
               measure.vars = c("sCF", "sDF1", "sDF2"))
# Add labels for plotting
scf_df$variable_label <- factor(scf_df$variable,
                                levels = c("sCF", "sDF1", "sDF2"),
                                labels = c("sCF", "sDF1", "sDF2"),
                                ordered = T)
scf_df$branch_label <- factor(scf_df$branch_description,
                              levels = c("To_all_animals", "To_all_other_metazoans", "To_CTEN_clade", "To_PORI_clade"),
                              labels = c("All Animals", "Other Metazoans", "Ctenophora", "Porifera"),
                              ordered = T)
scf_df$hypothesis_label <- factor(scf_df$hypothesis_tree,
                                  levels = c("CTEN", "PORI"),
                                  labels = c("Ctenophora-sister", "Porifera-sister"),
                                  ordered = T)
# Plot sCF scores
scf_plot <- ggplot(scf_df, aes(x = variable_label, y = value, group = variable_label, fill = hypothesis_label)) +
  geom_boxplot(alpha = 1) +
  facet_grid(branch_label~hypothesis_label) +
  scale_x_discrete(name = "Site concordance factor (sCF) statistics") +
  scale_y_continuous(name = "Value", breaks = seq(0,100,20), labels = seq(0,100,20), minor_breaks = seq(0,100,10), limits = c(0,100)) +
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"), labels = c("Ctenophora-sister", "Porifera-sister")) +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0, unit = "pt")),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt")),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0),
        strip.text = element_text(size = 18),
        legend.position = "none" )
# Save plot
scf_plot_name <- paste0(repo_dir, "figures/", "empirical_scf_scores.pdf")
ggsave(filename = scf_plot_name, plot = scf_plot, width = 8, height = 10, units = "in")



###### 6. Plot quartet scores  ######
# Extract quartet scores into long dataframe
qs_df <- melt(empirical_df,
              id.vars = id.var_cols,
              measure.vars = c("quartet_score", "gcf_Label"))
# Scale UFB from 0 - 1
qs_df[which(qs_df$variable == "gcf_Label"), ]$value <- qs_df[which(qs_df$variable == "gcf_Label"), ]$value/100
# Add labels for plotting
qs_df$branch_label <- factor(qs_df$branch_description,
                             levels = c("To_all_animals", "To_all_other_metazoans", "To_CTEN_clade", "To_PORI_clade"),
                             labels = c("All Animals", "Other Metazoans", "Ctenophora", "Porifera"),
                             ordered = T)
qs_df$variable_label <- factor(qs_df$variable,
                               levels = c("quartet_score", "gcf_Label"),
                               labels = c("Quartet Score", "Ultra-Fast Bootstrap"),
                               ordered = T)
qs_df$hypothesis_label <- factor(qs_df$hypothesis_tree,
                                 levels = c("CTEN", "PORI"),
                                 labels = c("Ctenophora-sister", "Porifera-sister"),
                                 ordered = T)
# Plot quartet scores
qs_plot <- ggplot(qs_df, aes(x = branch_label, y = value, fill = hypothesis_label)) +
  facet_grid(.~variable_label, scales = "fixed") +
  geom_boxplot(alpha = 1) +
  scale_x_discrete(name = "Clade") +
  scale_y_continuous(name = "Value", breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.05)) +
  scale_fill_manual(name = "Tree topology", values = c("#d8b365", "#5ab4ac"), labels = c("Ctenophora-sister", "Porifera-sister")) +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0, unit = "pt")),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt")),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0),
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.key.size = unit(30, "pt") )
# Save plot
qs_plot_name <- paste0(repo_dir, "figures/", "empirical_quartet_scores.pdf")
ggsave(filename = qs_plot_name, plot = qs_plot, width = 10, height = 8, units = "in")



###### 7. Plot branch lengths  ######
# Extract quartet scores into long dataframe
bl_df <- melt(empirical_df,
              id.vars = id.var_cols,
              measure.vars = c("gCF_length"))
# Add labels for plotting
bl_df$branch_label <- factor(bl_df$branch_description,
                             levels = c("To_all_animals", "To_all_other_metazoans", "To_CTEN_clade", "To_PORI_clade"),
                             labels = c("All Animals", "Other Metazoans", "Ctenophora", "Porifera"),
                             ordered = T)
bl_df$hypothesis_label <- factor(bl_df$hypothesis_tree,
                                 levels = c("CTEN", "PORI"),
                                 labels = c("Ctenophora-sister", "Porifera-sister"),
                                 ordered = T)
# Plot quartet scores
bl_plot <- ggplot(bl_df, aes(x = branch_label, y = value, fill = hypothesis_label)) +
  geom_boxplot(alpha = 1) +
  scale_x_discrete(name = "Clade") +
  scale_y_continuous(name = "Branch length (subs. per site)", breaks = seq(0,0.6,0.1), labels = seq(0,0.6,0.1), minor_breaks = seq(0,0.6,0.05), limits = c(0, 0.6)) +
  scale_fill_manual(name = "Tree topology", values = c("#d8b365", "#5ab4ac"), labels = c("Ctenophora-sister", "Porifera-sister")) +
  theme_bw() +
  theme(axis.title = element_text(size = 18),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0, unit = "pt")),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt")),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0),
        strip.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.key.size = unit(30, "pt") )
# Save plot
bl_plot_name <- paste0(repo_dir, "figures/", "empirical_branch_lengths.pdf")
ggsave(filename = bl_plot_name, plot = bl_plot, width = 7, height = 8, units = "in")



###### 8. Collate branch support and branch length plots  ######
# Plot quartet scores
qs_plot <- ggplot(qs_df, aes(x = branch_label, y = value, fill = hypothesis_label)) +
  facet_grid(.~variable_label, scales = "fixed") +
  geom_boxplot(alpha = 1) +
  scale_x_discrete(name = "Clade") +
  scale_y_continuous(name = "Value", breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.05)) +
  scale_fill_manual(name = "Tree topology", values = c("#d8b365", "#5ab4ac"), labels = c("Ctenophora-sister", "Porifera-sister")) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0, unit = "pt")),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt")),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0),
        strip.text = element_text(size = 16),
        legend.position = "none")
# Plot quartet scores
bl_plot <- ggplot(bl_df, aes(x = branch_label, y = value, fill = hypothesis_label)) +
  geom_boxplot(alpha = 1) +
  scale_x_discrete(name = "Clade") +
  scale_y_continuous(name = "Branch length (subs. per site)", breaks = seq(0,0.6,0.1), labels = seq(0,0.6,0.1), minor_breaks = seq(0,0.6,0.05), limits = c(0, 0.6)) +
  scale_fill_manual(name = "Tree topology", values = c("#d8b365", "#5ab4ac"), labels = c("Ctenophora-sister", "Porifera-sister")) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 10, l = 0, unit = "pt")),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 10, unit = "pt")),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1.0),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.key.size = unit(30, "pt") )
# Create quilt
quilt <- qs_plot + bl_plot + plot_layout(widths = c(2, 1), ncol = 2) + plot_annotation(tag_levels = 'a', tag_suffix = ".") & 
  theme(plot.tag = element_text(size = 30))
# Save quilt
quilt_name <- paste0(repo_dir, "figures/", "empirical_collated_branch.pdf")
ggsave(filename = quilt_name, plot = quilt, width = 12, height = 8, units = "in")

