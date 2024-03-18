# ancient_ILS/code/05_plot_figures_cf.R
## This script plots result figures from simulation output
# Caitlin Cherryh, 2024

## These plots focus on the key branch, i.e., the branch with different rearrangements of the four clades
#     (Ctenophora), (Porifera), (Outgroup), and (Cnidaria+Bilateria)
# Key branch to extract:
#   For CTEN-Tree: "All other animals"
#   For PORI-Tree: "All other animals"
#   For CTEN-PORI-Tree: "CTEN+PORI

###### 1. Input parameters ######
## File paths
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# output_dir                <- location of csv files

repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/05_output_files/"



###### 2. Open packages and functions ######
## Open packages
library(reshape2)
library(ggplot2)
library(ggtern) # for ternary plots of sCF
library(patchwork) # for assembling ternary plots of emperical gCF, sCF and quartet scores

# Specify colour palettes used within these plots
metazoan_clade_palette  <- c(Bilateria = "#CC79A7", Cnidaria = "#009E73", Ctenophora = "#56B4E9", Porifera = "#E69F00", Outgroup = "#999999")
clades_colours          <- c("Ctenophora" = "#2171b5", "Porifera" =  "#E69F00", "Ctenophora+Porifera" = "#009E73")
clades_black            <- c("Ctenophora" = "black", "Porifera" =  "black", "Ctenophora+Porifera" = "black")
passfail_palette        <- c("Pass (p > 0.05)" = "#addd8e", "Fail (p \u2264 0.05)" = "#005a32")



###### 3. Open and prepare csvs for plotting  ######
# List all output files
all_output_files <- list.files(output_dir)

# Open sCF results and remove unused datasets
species_scf_df  <- read.csv(paste0(repo_dir, "output/empirical_dataset_concordance_factors.csv"),  stringsAsFactors = FALSE)
species_scf_df  <- species_scf_df[which(species_scf_df$dataset != "Simion2017"), ]
species_scf_df  <- species_scf_df[which(species_scf_df$dataset != "Hejnol2009"), ]

# Open constrained/unconstrained cf results
cf_df <- read.csv(paste0(output_dir, grep("constrained_unconstrained_cf_output_checked.csv", all_output_files, value = T)), stringsAsFactors = F)
# Remove any Plac-sister branches (these do not include the branch of interest)
cf_df <- cf_df[grep("PlacSister", cf_df$branch_description, invert = T), ]



###### 4. Update and format dataframes for plots ######
# Add any missing columns to species_scf_df
species_scf_df$dataset_id <- paste0(species_scf_df$dataset, ".", species_scf_df$matrix) 
species_scf_df$dataset_id_formatted <- factor(species_scf_df$dataset_id,
                                              levels =  c("Dunn2008.Dunn2008_FixedNames", "Philippe2009.Philippe_etal_superalignment_FixedNames", "Philippe2011.UPDUNN_MB_FixedNames", 
                                                          "Nosenko2013.nonribosomal_9187_smatrix", "Nosenko2013.ribosomal_14615_smatrix", "Ryan2013.REA_EST_includingXenoturbella", 
                                                          "Moroz2014.ED3d", "Borowiec2015.Best108", "Chang2015.Chang_AA", 
                                                          "Whelan2015.Dataset10", "Whelan2017.Metazoa_Choano_RCFV_strict", "Laumer2018.Tplx_BUSCOeuk"),
                                              labels = c("Dunn 2008",  "Philippe 2009", "Philippe 2011", 
                                                         "Nosenko 2013\nnonribosomal", "Nosenko 2013\nribosomal", "Ryan 2013", 
                                                         "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                                         "Whelan 2015", "Whelan 2017",  "Laumer 2018"),
                                              ordered = TRUE)
species_scf_df$tree_topology <- species_scf_df$hypothesis_tree
species_scf_df$tree_topology_formatted <- factor(species_scf_df$tree_topology,
                                                 levels =  c("CTEN", "PORI", "CTEN_PORI"),
                                                 labels = c("Ctenophora", "Porifera", "Ctenophora+Porifera"),
                                                 ordered = TRUE)
species_scf_df$branch_to_clade <- factor(species_scf_df$branch_description,
                                         levels = c("To_all_animals", "To_all_other_metazoans", "To_CTEN_clade", "To_PORI_clade",
                                                    "To_CTEN+PLAC_clade", "To_PORI+PLAC_clade", "To_all_animals_except_PLAC"),
                                         labels = c("ALL_ANIMALS", "ALL_OTHER_ANIMALS", "CTEN", "PORI",
                                                    "CTEN_PLAC", "PORI_PLAC", "ALL_ANIMLS_EXCEPT_PLAC"),
                                         ordered = TRUE)
species_scf_df$ultrafast_bootstrap <- species_scf_df$sCF_Label
# Add any missing columns to cf_df
cf_df$dataset_id <- paste0(cf_df$dataset, ".", cf_df$matrix) 
cf_df$dataset_id_formatted <- factor(cf_df$dataset_id,
                                     levels =  c("Dunn2008.Dunn2008_FixedNames", "Philippe2009.Philippe_etal_superalignment_FixedNames", "Philippe2011.UPDUNN_MB_FixedNames", 
                                                 "Nosenko2013.nonribosomal_9187_smatrix", "Nosenko2013.ribosomal_14615_smatrix", "Ryan2013.REA_EST_includingXenoturbella", 
                                                 "Moroz2014.ED3d", "Borowiec2015.Best108", "Chang2015.Chang_AA", 
                                                 "Whelan2015.Dataset10", "Whelan2017.Metazoa_Choano_RCFV_strict", "Laumer2018.Tplx_BUSCOeuk"),
                                     labels = c("Dunn 2008",  "Philippe 2009", "Philippe 2011", 
                                                "Nosenko 2013\nnonribosomal", "Nosenko 2013\nribosomal", "Ryan 2013", 
                                                "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                                "Whelan 2015", "Whelan 2017",  "Laumer 2018"),
                                     ordered = TRUE)
cf_df$tree_topology <- cf_df$hypothesis_tree
cf_df$tree_topology_formatted <- factor(cf_df$tree_topology,
                                        levels =  c("CTEN", "PORI", "CTEN_PORI"),
                                        labels = c("Ctenophora", "Porifera", "Ctenophora+Porifera"),
                                        ordered = TRUE)
cf_df$branch_to_clade <- factor(cf_df$branch_description,
                                levels = c("Metazoa", "All_other_animals", "CTEN_PORI"),
                                labels = c("ALL_ANIMALS", "ALL_OTHER_ANIMALS", "CTEN_PORI"),
                                ordered = TRUE)
cf_df$clade_formatted <- factor(cf_df$branch_description,
                                levels = c("Metazoa", "All_other_animals", "CTEN_PORI"),
                                labels = c("Metazoa", "All other animals", "Ctenophora+Porifera"),
                                ordered = TRUE)
cf_df$ultrafast_bootstrap <- cf_df$gcf_Label
cf_df$gene_type_formatted <- factor(cf_df$gene_type,
                                    levels = c("Unconstrained", "Constrained"),
                                    labels = c("Unconstrained genes", "Constrained genes"),
                                    ordered = T)



###### 4. Species sCF df: Empirical gCF, sCF and quartet values, All Other Metazoans, 2 Topologies ######
# Extract variables of interest
emp2_df <- species_scf_df
emp2_df <- emp2_df[which(emp2_df$branch_to_clade %in% c("ALL_ANIMALS", "ALL_OTHER_ANIMALS", "CTEN", "PORI")), ]
# Add new columns
emp2_df$clade_formatted <- factor(emp2_df$branch_to_clade,
                                  levels = c("ALL_ANIMALS", "ALL_OTHER_ANIMALS", "CTEN", "PORI"),
                                  labels = c("Metazoa", "All other animals", "Ctenophora", "Porifera"))
# Plot gCF
emp2_gcf <- ggtern(emp2_df, aes(x = gDF1, y = gCF, z = gDF2, color = tree_topology_formatted, shape = tree_topology_formatted)) +
  geom_point(size = 4, alpha = 0.6) +
  facet_wrap(clade_formatted ~., nrow = 1, ncol = 4) +
  labs(title = "a.") +
  scale_color_manual(values = clades_colours, name = "Constrained\ntree topology") +
  scale_shape_manual(values = c("Ctenophora" = 16, "Porifera" = 17), name = "Constrained\ntree topology") +
  theme_bw() +
  theme(plot.title = element_text(size = 30),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.5, "lines"),
        tern.panel.grid.major = element_line(colour = "grey80", linewidth = 8),
        tern.panel.grid.minor = element_line(colour = "grey80", linewidth = 8) ) +
  guides( color = guide_legend(override.aes = list(size = 5)) )  
# Plot sCF
emp2_scf <- ggtern(emp2_df, aes(x = sDF1, y = sCF, z = sDF2, color = tree_topology_formatted, shape = tree_topology_formatted)) +
  geom_point(size = 4, alpha = 0.6) +
  facet_wrap(clade_formatted ~., nrow = 1, ncol = 4) +
  labs(title = "b.") +
  scale_color_manual(values = clades_colours, name = "Constrained\ntree topology") +
  scale_shape_manual(values = c("Ctenophora" = 16, "Porifera" = 17), name = "Constrained\ntree topology") +
  theme_bw() +
  theme(plot.title = element_text(size = 30),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.5, "lines"),
        tern.panel.grid.major = element_line(colour = "grey80", linewidth = 8),
        tern.panel.grid.minor = element_line(colour = "grey80", linewidth = 8) ) +
  guides( color = guide_legend(override.aes = list(size = 5)) ) 
# Plot quartet scores
emp2_qs <- ggtern(emp2_df, aes(x = q2, y = q1, z = q3, color = tree_topology_formatted, shape = tree_topology_formatted)) +
  geom_point(size = 4, alpha = 0.6) +
  facet_wrap(clade_formatted ~., nrow = 1, ncol = 4) +
  labs(title = "c.") +
  scale_color_manual(values = clades_colours, name = "Constrained\ntree topology") +
  scale_shape_manual(values = c("Ctenophora" = 16, "Porifera" = 17), name = "Constrained\ntree topology") +
  theme_bw() +
  theme(plot.title = element_text(size = 30),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.5, "lines"),
        tern.panel.grid.major = element_line(colour = "grey80", linewidth = 8),
        tern.panel.grid.minor = element_line(colour = "grey80", linewidth = 8) ) +
  guides( color = guide_legend(override.aes = list(size = 5)) ) +
  Tlab("qCF") +
  Llab("qDF1") +
  Rlab("qDF2")
# Assemble the three ternary plots using ggtern::grid.arrange 
#     (as ternary plots have three axes, patchwork and ggplot::grid.arrange don't work cleanly here)
quilt2 <- ggtern::grid.arrange(emp2_gcf, emp2_scf, emp2_qs, nrow = 3, ncol = 1)
# Save the quilt
emp_tern_file <- paste0(repo_dir, "figures/", "empirical_dataset_gCF_sCF_qs_ternary_plots_coloured.pdf")
ggsave(filename = emp_tern_file, plot = quilt2, width = 18, height = 12, units = "in")
emp_tern_file_png <- paste0(repo_dir, "figures/", "empirical_dataset_gCF_sCF_qs_ternary_plots_coloured.png")
ggsave(filename = emp_tern_file_png, plot = quilt2, width = 18, height = 12, units = "in")



###### 4. CF df: Empirical gCF, sCF and quartet values, Key branch, 3 Topologies ######
# Extract variables of interest
emp3_df <- cf_df
# Remove any rows that are NOT leading to the "key branch"
emp3_df <- emp3_df[which(emp3_df$hypothesis_tree == "CTEN" & emp3_df$branch_to_clade == "ALL_ANIMALS" |
                           emp3_df$hypothesis_tree == "PORI" & emp3_df$branch_to_clade == "ALL_ANIMALS" |
                           emp3_df$hypothesis_tree == "CTEN_PORI" & emp3_df$branch_to_clade == "CTEN_PORI"), ]
# Plot gCF
emp3_gcf <- ggtern(emp3_df, aes(x = gDF1, y = gCF, z = gDF2, color = tree_topology_formatted, shape = tree_topology_formatted)) +
  geom_mask() +  
  geom_point(size = 3, alpha = 0.6) +
  facet_wrap(gene_type_formatted ~., nrow = 1) +
  labs(title = "a.") +
  scale_color_manual(values = clades_colours, name = "Constrained\ntree topology") +
  scale_shape_manual(values = c("Ctenophora" = 16, "Porifera" = 17, "Ctenophora+Porifera" = 15), name = "Constrained\ntree topology") +
  theme_bw() +
  theme_showgrid() +
  theme(plot.title = element_text(size = 30),
        plot.subtitle = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.5, "lines"),
        tern.panel.grid.major = element_line(colour = "grey80", linewidth = 12),
        tern.panel.grid.minor = element_line(colour = "grey80", linewidth = 12),
        legend.position = "none",
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(-5, -5, -5, -5)) 
# Plot sCF
emp3_unconstrained_df <- emp3_df[which(emp3_df$gene_type == "Unconstrained"), ]
emp3_scf <- ggtern(emp3_unconstrained_df, aes(x = sDF1, y = sCF, z = sDF2, color = tree_topology_formatted, shape = tree_topology_formatted)) +
  geom_mask() +  
  geom_point(size = 3, alpha = 0.6) +
  facet_wrap(gene_type_formatted ~., nrow = 1) +
  labs(title = "b.") +
  scale_color_manual(values = clades_colours, name = "Constrained\ntree topology") +
  scale_shape_manual(values = c("Ctenophora" = 16, "Porifera" = 17, "Ctenophora+Porifera" = 15), name = "Constrained\ntree topology") +
  theme_bw() +
  theme_showgrid() +
  theme(plot.title = element_text(size = 30),
        plot.subtitle = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1.5, "lines"),
        legend.margin = margin(l = 85, unit = "pt"),
        tern.panel.grid.major = element_line(colour = "grey80", linewidth = 12),
        tern.panel.grid.minor = element_line(colour = "grey80", linewidth = 12),
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(-5, -5, -5, -5)) +
  guides( color = guide_legend(override.aes = list(size = 5)) )  
# Plot quartet scores
emp3_qs <- ggtern(emp3_df, aes(x = q2, y = q1, z = q3, color = tree_topology_formatted, shape = tree_topology_formatted)) +
  geom_mask() +  
  geom_point(size = 3, alpha = 0.6) +
  facet_wrap(gene_type_formatted ~., nrow = 1) +
  labs(title = "c.") +
  scale_color_manual(values = clades_colours, name = "Constrained\ntree topology") +
  scale_shape_manual(values = c("Ctenophora" = 16, "Porifera" = 17, "Ctenophora+Porifera" = 15), name = "Constrained\ntree topology") +
  theme_bw() +
  theme_showgrid() +
  theme(plot.title = element_text(size = 30),
        plot.subtitle = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.5, "lines"),
        tern.panel.grid.major = element_line(colour = "grey80", linewidth = 12),
        tern.panel.grid.minor = element_line(colour = "grey80", linewidth = 12),
        legend.position = "none",
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(-5, -5, -5, -5)) +
  Tlab("qCF") +
  Llab("qDF1") +
  Rlab("qDF2")

# Assemble the three ternary plots using ggtern::grid.arrange 
#     (as ternary plots have three axes, patchwork and ggplot::grid.arrange don't work cleanly here)
quilt3 <- ggtern::grid.arrange(emp3_gcf, emp3_scf, emp3_qs, nrow = 3, ncol = 1)
# Save the quilt
emp3_tern_file <- paste0(repo_dir, "figures/", "cf_constrained_unconstrained_ternary_coloured.pdf")
ggsave(filename = emp3_tern_file, plot = quilt3, width = 10, height = 12, units = "in")
emp3_tern_file_png <- paste0(repo_dir, "figures/", "cf_constrained_unconstrained_ternary_coloured.png")
ggsave(filename = emp3_tern_file_png, plot = quilt3, width = 10, height = 12, units = "in")



###### 5. CF + Key branch + 3 topologies + unconstrained gene trees ######
# Key branch to extract:
#   CF = CTEN-Tree, "All other animals"
#   DF1 = PORI-Tree, "All other animals"
#   DF2 = CTEN-PORI-Tree, "CTEN+PORI"

## Prepare unconstrained gene df
# Extract unconstrained gene tree rows with the right branch for each tree topology
filtered_uncon_df <- cf_df[c(which(cf_df$hypothesis_tree == "CTEN" & cf_df$branch_description == "All_other_animals" & cf_df$gene_type == "Unconstrained"),
                             which(cf_df$hypothesis_tree == "PORI" & cf_df$branch_description == "All_other_animals" & cf_df$gene_type == "Unconstrained"),
                             which(cf_df$hypothesis_tree == "CTEN_PORI" & cf_df$branch_description == "CTEN_PORI" & cf_df$gene_type == "Unconstrained")) , 
                           c(1:6, 41:48, 7:9, 16:18, 19:21, 26:28,29, 32, 35, 38:39)]
row.names(filtered_uncon_df) <- 1:nrow(filtered_uncon_df)
# Reformat the dataframe to have the separate CF for each topology around the branch of interest
uncon_df <- extract.correct.CF.values(df = filtered_uncon_df)
uncon_df <- convert.to.long.format(df = uncon_df)

## Prepare constrained gene df
# Extract constrained gene tree rows with the right branch for each tree topology
filtered_con_df <- cf_df[c(which(cf_df$hypothesis_tree == "CTEN" & cf_df$branch_description == "All_other_animals" & cf_df$gene_type == "Constrained"),
                           which(cf_df$hypothesis_tree == "PORI" & cf_df$branch_description == "All_other_animals" & cf_df$gene_type == "Constrained"),
                           which(cf_df$hypothesis_tree == "CTEN_PORI" & cf_df$branch_description == "CTEN_PORI" & cf_df$gene_type == "Constrained")) , 
                         c(1:6, 41:48, 7:9, 16:18, 19:21, 26:28,29, 32, 35, 38:39)]
row.names(filtered_con_df) <- 1:nrow(filtered_con_df)
# Reformat the dataframe to have the separate CF for each topology around the branch of interest
con_df <- extract.correct.CF.values(df = filtered_con_df)
con_df <- convert.to.long.format(df = con_df)

## Assemble the dataframes
plot_df <- rbind(uncon_df, con_df)

## Create the plot
# Plot CFs
cf_tern <- ggtern(plot_df, aes(x = cf_val, y = df1_val, z = df2_val)) +
  geom_mask() +  
  geom_point(size = 5, alpha = 0.6) +
  facet_grid(cf_type ~ gene_type_formatted) +
  theme_bw() +
  theme(strip.clip = "off",
        strip.text = element_text(size = 30),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        tern.panel.grid.major = element_line(colour = "grey80", linewidth = 12),
        tern.panel.grid.minor = element_line(colour = "grey80", linewidth = 12),
        panel.spacing = unit(1.5, "lines")) +
  Tlab("CTEN") +
  Llab("PORI") +
  Rlab("CTEN \u002B\nPORI")
# Save the plot
emp_tern_file <- paste0(repo_dir, "figures/", "cf_constrained_unconstrained_keyBranch_ternary.pdf")
ggsave(filename = emp_tern_file, plot = cf_tern, width = 16, height = 16, units = "in")
emp_tern_file_png <- paste0(repo_dir, "figures/", "cf_constrained_unconstrained_keyBranch_ternary.png")
ggsave(filename = emp_tern_file_png, plot = cf_tern, width = 16, height = 16, units = "in")

## Reformat the csv as a nice table in the manuscript
# Remove unneeded columns
op_plot_df <- rbind(extract.correct.CF.values(df = filtered_uncon_df), extract.correct.CF.values(df = filtered_con_df))
output_df_raw_file <- paste0(output_dir, "constrained_unconstrained_cf_KeyBranch_raw.csv")
write.csv(op_plot_df, file = output_df_raw_file, row.names = FALSE)
# Reorder columns and rows
op_plot_df <- op_plot_df[ , c(3, 8, 10:30, 31, 34)]
op_plot_df$gene_type <- factor(op_plot_df$gene_type, levels = c("Unconstrained", "Constrained"), ordered = T)
op_plot_df <- op_plot_df[order(op_plot_df$dataset, op_plot_df$gene_type), ]
# Rename columns
names(op_plot_df) <- c("Dataset", "Gene tree type", 
                       "gCF", "gCF_N", "gDF1", "gDF1_N", "gDF2", "gDF2_N",
                       "sCF", "sCF_N", "sDF1", "sDF1_N", "sDF2", "sDF2_N",
                       "qCF", "qDF1", "qDF2", "qCF_f", "qDF1_f", "qDF2_f", "qCF_pp", "qDF1_pp", "qDF2_pp", "QC", "EN")
# Save the csv as output
output_df_file <- paste0(output_dir, "constrained_unconstrained_cf_KeyBranch.csv")
write.csv(op_plot_df, file = output_df_file, row.names = FALSE)


