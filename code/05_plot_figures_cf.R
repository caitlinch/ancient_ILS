# ancient_ILS/code/05_plot_figures_cf.R
## This script plots result figures from concordance factors output
# Caitlin Cherryh 2024

## These plots focus on the key branch, i.e., the branch with different rearrangements of the four clades
#     (Ctenophora), (Porifera), (Outgroup), and (Cnidaria+Bilateria)
# Key branch to extract:
#   For CTEN-Tree: "All other animals"
#   For PORI-Tree: "All other animals"
#   For CTEN-PORI-Tree: "CTEN+PORI

# Code snippets for patchwork:
# Add title: 
#   plot_annotation(title = "gCF", theme = theme(plot.title = element_text(size = 25, vjust = 0.5, hjust = 0.5)))
# Add panel labels:
#   plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
# Add both: 
#   plot_annotation(title = "gCF", tag_levels = 'a', tag_suffix = ".", theme = theme(plot.title = element_text(size = 25, vjust = 0.5, hjust = 0.5))) & theme(plot.tag = element_text(size = 30))


##### 1. Input parameters #####
## Specify parameters:
# repo_dir          <- Location of caitlinch/metazoan-mixtures github repository
# output_csv_dir    <- Location of csv input/output/results files
# plot_dir          <- Location to save plots

repo_dir          <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_csv_dir    <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/07_output_files/"
plot_dir          <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/08_figures/"



##### 2. Open packages and functions #####
## Open packages
library(reshape2)
library(ggplot2)
library(ggtern) # for ternary plots of sCF
library(ggpubr) # for adding statistics and regression lines: stat_regline_equation()
library(patchwork) # for assembling ternary plots of emperical gCF, sCF and quartet scores
library(grDevices) # for hcl.colors(), palette.colors()

## Source functions
source(paste0(repo_dir, "code/func_plotting.R"))

## Specify colour palettes used within these plots
# Inbuilt colour palettes
okito_palette           <- palette.colors(palette = "Okabe-Ito")
viridis_palette         <- hcl.colors(n = 8, palette = "Viridis")
mako_palette            <- hcl.colors(n = 5, palette = "Mako")
# Colour palettes for variables
metazoan_clade_palette  <- c("Bilateria" = "#CC79A7", "Cnidaria" = "#009E73", "Ctenophora" = "#56B4E9",
                             "Porifera" = "#E69F00", "Outgroup" = "#999999", "Placozoa" = "#000000")
topology_colours        <- c("Ctenophora" = viridis_palette[1], "Porifera" =  viridis_palette[4], 
                             "Ctenophora+Porifera" = viridis_palette[8], "Paraphyly" = "#999999")

##### 3. Open and prepare csvs for plotting  #####
## List all output files
all_output_files <- paste0(repo_dir, "output/", list.files(output_csv_dir))

## Open empirical cf analysis results
gcf_df    <- read.csv(grep("gCF_values.csv", all_output_files, value = TRUE))
gcf_long  <- read.csv(grep("gCF_values_formatted.csv", all_output_files, value = TRUE))
qcf_df    <- read.csv(grep("qCF_values.csv", all_output_files, value = TRUE))
qcf_long  <- read.csv(grep("qCF_values_formatted.csv", all_output_files, value = TRUE))        



##### 4. Update and format data frames for plots #####
## Add nicely formatted dataset labels to the data frames
dataset_labels                <- c("Dunn2008.Dunn2008_FixedNames", "Philippe2009.Philippe_etal_superalignment_FixedNames", "Philippe2011.UPDUNN_MB_FixedNames", 
                                   "Nosenko2013.nonribosomal_9187_smatrix", "Nosenko2013.ribosomal_14615_smatrix", "Ryan2013.REA_EST_includingXenoturbella", 
                                   "Moroz2014.ED3d", "Borowiec2015.Best108", "Chang2015.Chang_AA", 
                                   "Whelan2015.Dataset10", "Whelan2017.Metazoa_Choano_RCFV_strict", "Laumer2018.Tplx_BUSCOeuk")
dataset_labels_formatted      <- c("Dunn 2008",  "Philippe 2009", "Philippe 2011", 
                                   "Nosenko 2013\nnonribosomal", "Nosenko 2013\nribosomal", "Ryan 2013", 
                                   "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                   "Whelan 2015", "Whelan 2017",  "Laumer 2018")
gcf_df$dataset_id_formatted   <- factor(gcf_df$dataset_id, levels =  dataset_labels, labels = dataset_labels_formatted, ordered = TRUE)
gcf_long$dataset_id_formatted <- factor(gcf_long$dataset_id, levels =  dataset_labels, labels = dataset_labels_formatted, ordered = TRUE)
qcf_df$dataset_id_formatted   <- factor(qcf_df$dataset_id, levels =  dataset_labels, labels = dataset_labels_formatted, ordered = TRUE)
qcf_long$dataset_id_formatted <- factor(qcf_long$dataset_id, levels =  dataset_labels, labels = dataset_labels_formatted, ordered = TRUE)

## Add nicely formatted topology labels to the data frames
topology_labels                   <- c("CTEN", "PORI", "CTEN_PORI")
topology_labels_formatted         <- c("Ctenophora", "Porifera", "Ctenophora+Porifera")
topology_labels_short             <- c("CTEN", "PORI", "CTEN+PORI")
gcf_df$tree_topology_formatted    <- factor(gcf_df$tree_topology, levels =  topology_labels, labels = topology_labels_formatted, ordered = TRUE)
qcf_df$tree_topology_formatted    <- factor(qcf_df$topology, levels =  topology_labels, labels = topology_labels_formatted, ordered = TRUE)
gcf_df$tree_topology_short        <- factor(gcf_df$tree_topology, levels =  topology_labels, labels = topology_labels_short, ordered = TRUE)
qcf_df$tree_topology_short        <- factor(qcf_df$topology, levels =  topology_labels, labels = topology_labels_short, ordered = TRUE)

## Add nicely formatted model labels to the data frames
models                    <- c("Partition", "C60")
models_formatted          <- c("Partition", "C60")
gcf_df$model_formatted    <- factor(gcf_df$model, levels =  models, labels = models_formatted, ordered = TRUE)
gcf_long$model_formatted  <- factor(gcf_long$model, levels =  models, labels = models_formatted, ordered = TRUE)
qcf_df$model_formatted    <- factor(qcf_df$model, levels =  models, labels = models_formatted, ordered = TRUE)
qcf_long$model_formatted  <- factor(qcf_long$model, levels =  models, labels = models_formatted, ordered = TRUE)

## Add nicely formatted PLACOZOA label to the data frames
Plac_labels                     <- c("Plac", "noPlac")
Plac_labels_formatted           <- c("Present", "Removed")
gcf_df$Plac_present_formatted   <- factor(gcf_df$Plac_present, levels =  Plac_labels, labels = Plac_labels_formatted, ordered = TRUE)
gcf_long$Plac_present_formatted <- factor(gcf_long$Plac_present, levels =  Plac_labels, labels = Plac_labels_formatted, ordered = TRUE)
qcf_df$Plac_present_formatted   <- factor(qcf_df$Plac_present, levels =  Plac_labels, labels = Plac_labels_formatted, ordered = TRUE)
qcf_long$Plac_present_formatted <- factor(qcf_long$Plac_present, levels =  Plac_labels, labels = Plac_labels_formatted, ordered = TRUE)

## Add nicely formatted column describing concordance factor type to the data frames
gcf_df$analysis   <- "gCF"
gcf_long$analysis <- "gCF"
qcf_df$analysis   <- "qCF"
qcf_long$analysis <- "qCF"

## Create 2 sets of datasets for plotting into nice grids
# column name: dataset_id
dataset_group_1 <- c("Dunn2008.Dunn2008_FixedNames", "Philippe2009.Philippe_etal_superalignment_FixedNames", "Philippe2011.UPDUNN_MB_FixedNames", 
                     "Nosenko2013.nonribosomal_9187_smatrix", "Nosenko2013.ribosomal_14615_smatrix", "Ryan2013.REA_EST_includingXenoturbella")
dataset_group_2 <- c("Moroz2014.ED3d", "Borowiec2015.Best108", "Chang2015.Chang_AA", 
                     "Whelan2015.Dataset10", "Whelan2017.Metazoa_Choano_RCFV_strict", "Laumer2018.Tplx_BUSCOeuk")



##### 5. gene concordance factors (gCF) ternary plots #####
## PLOT gCF 1: gCF values as bars for each dataset
# Create plot panels
gcf_topology_bars_panel1 <- ggplot(data = gcf_df[which(gcf_df$dataset_id %in% dataset_group_1), ], 
                                   aes(x = tree_topology_short, y = KEY_gCF, fill = tree_topology_short)) +
  geom_bar(stat="identity") +
  facet_grid(model_formatted~dataset_id_formatted) +
  scale_x_discrete(name = "Constrained topology") +
  scale_y_continuous(name = "Value", breaks = seq(0,25,5), labels = seq(0,25,5), minor_breaks = seq(0,25,2.5), limits = c(0,25)) +
  scale_fill_manual(name = "Topology", values = c("CTEN" = topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                  "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]])) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title.x = element_text(size = 16, margin = margin(t=10, r=0, b=20, l=0, unit="pt")),
        axis.text.x = element_text(size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 16, margin = margin(t=0, r=10, b=0, l=0, unit="pt")),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 16, margin = margin(t=0, r=0, b=10, l=0, unit="pt")),
        legend.text = element_text(size = 13))
gcf_topology_bars_panel2 <- ggplot(data = gcf_df[which(gcf_df$dataset_id %in% dataset_group_2), ], 
                                   aes(x = tree_topology_short, y = KEY_gCF, fill = tree_topology_short)) +
  geom_bar(stat="identity") +
  facet_grid(model_formatted~dataset_id_formatted) + 
  scale_x_discrete(name = "Constrained topology") +
  scale_y_continuous(name = "Value", breaks = seq(0,25,5), labels = seq(0,25,5), minor_breaks = seq(0,25,2.5), limits = c(0,25)) +
  scale_fill_manual(name = "Topology", values = c("CTEN" = topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                  "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]])) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title.x = element_text(size = 16, margin = margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.text.x = element_text(size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 16, margin = margin(t=0, r=10, b=0, l=0, unit="pt")),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 16, margin = margin(t=0, r=0, b=10, l=0, unit="pt")),
        legend.text = element_text(size = 13))
# Assemble patchwork and save plot
gcf_topology_bars_quilt <- gcf_topology_bars_panel1 / gcf_topology_bars_panel2 + 
  plot_annotation(title = "Gene Concordance Factors", 
                  theme = theme(plot.title = element_text(size = 25, vjust = 0.5, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0, unit = "pt"))))
ggsave(filename = paste0(plot_dir, "cf_gcf_topology_bar.pdf"), plot = gcf_topology_bars_quilt, width = 12, height = 10, units = "in")
ggsave(filename = paste0(plot_dir, "cf_gcf_topology_bar.png"), plot = gcf_topology_bars_quilt, width = 12, height = 10, units = "in")

## PLOT gCF 2: gCF/gDFP values as bars for each dataset
# Create new dataframe with gDFP values included
gcf_trimmed_df <- gcf_df[ , c("analysis", "id", "dataset", "matrix_name", "dataset_id", "dataset_id_formatted", "model", "model_formatted",
                              "tree_topology", "tree_topology_formatted", "tree_topology_short", "Plac_present", "Plac_present_formatted",
                              "KEY_gCF", "KEY_Length", "CTEN_Length", "CTEN_monophyly", "PORI_Length", "PORI_monophyly")]
gdfp_df <- gcf_df[which(gcf_df$tree_topology == "CTEN"), 
                  c("analysis", "id", "dataset", "matrix_name", "dataset_id", "dataset_id_formatted", "model", "model_formatted",
                    "tree_topology", "tree_topology_formatted", "tree_topology_short", "Plac_present", "Plac_present_formatted",
                    "KEY_gDFP", "KEY_Length", "CTEN_Length", "CTEN_monophyly", "PORI_Length", "PORI_monophyly")]
gdfp_df$tree_topology <- "Paraphyly"
gdfp_df$tree_topology <- "Paraphyly"
names(gdfp_df) <- c("analysis", "id", "dataset", "matrix_name", "dataset_id", "dataset_id_formatted", "model", "model_formatted",
                    "tree_topology", "tree_topology_formatted", "tree_topology_short", "Plac_present", "Plac_present_formatted",
                    "KEY_gCF", "KEY_Length", "CTEN_Length", "CTEN_monophyly", "PORI_Length", "PORI_monophyly")
gdfp_df <- rbind(gcf_trimmed_df, gdfp_df)
# Reformat topology labels
topology_labels_gdfp                <- c("CTEN", "PORI", "CTEN_PORI", "Paraphyly")
topology_labels_formatted_gdfp      <- c("Ctenophora", "Porifera", "Ctenophora+Porifera", "Paraphyly")
topology_labels_short_gdfp          <- c("CTEN", "PORI", "CTEN+PORI", "PARAPHYLY")
gdfp_df$tree_topology_formatted     <- factor(gdfp_df$tree_topology, levels =  topology_labels_gdfp, labels = topology_labels_formatted_gdfp, ordered = TRUE)
gdfp_df$tree_topology_short         <- factor(gdfp_df$tree_topology, levels =  topology_labels_gdfp, labels = topology_labels_short_gdfp, ordered = TRUE)
# Create plot panels
gdfp_topology_bars_panel1 <- ggplot(data = gdfp_df[which(gdfp_df$dataset_id %in% dataset_group_1), ], 
                                    aes(x = tree_topology_short, y = KEY_gCF, fill = tree_topology_short)) +
  geom_bar(stat="identity") +
  facet_grid(model_formatted~dataset_id_formatted) +
  scale_x_discrete(name = "Constrained topology") +
  scale_y_continuous(name = "Value", breaks = seq(0,100,20), labels = seq(0,100,20), minor_breaks = seq(0,100,10), limits = c(0,100)) +
  scale_fill_manual(name = "Topology", values = c("CTEN" = topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                  "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]], "PARAPHYLY" = topology_colours[["Paraphyly"]])) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title.x = element_text(size = 16, margin = margin(t=10, r=0, b=20, l=0, unit="pt")),
        axis.text.x = element_text(size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 16, margin = margin(t=0, r=10, b=0, l=0, unit="pt")),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 16, margin = margin(t=0, r=0, b=10, l=0, unit="pt")),
        legend.text = element_text(size = 13))
gdfp_topology_bars_panel2 <- ggplot(data = gdfp_df[which(gdfp_df$dataset_id %in% dataset_group_2), ], 
                                    aes(x = tree_topology_short, y = KEY_gCF, fill = tree_topology_short)) +
  geom_bar(stat="identity") +
  facet_grid(model_formatted~dataset_id_formatted) + 
  scale_x_discrete(name = "Constrained topology") +
  scale_y_continuous(name = "Value", breaks = seq(0,100,20), labels = seq(0,100,20), minor_breaks = seq(0,100,10), limits = c(0,100)) +
  scale_fill_manual(name = "Topology", values = c("CTEN" = topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                  "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]], "PARAPHYLY" = topology_colours[["Paraphyly"]])) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title.x = element_text(size = 16, margin = margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.text.x = element_text(size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 16, margin = margin(t=0, r=10, b=0, l=0, unit="pt")),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 16, margin = margin(t=0, r=0, b=10, l=0, unit="pt")),
        legend.text = element_text(size = 13))
# Assemble patchwork and save plot
gdfp_topology_bars_quilt <- gdfp_topology_bars_panel1 / gdfp_topology_bars_panel2 + 
  plot_annotation(title = "Gene Concordance Factors",
                  theme = theme(plot.title = element_text(size = 25, vjust = 0.5, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0, unit = "pt"))))
ggsave(filename = paste0(plot_dir, "cf_gdfp_topology_bar.pdf"), plot = gdfp_topology_bars_quilt, width = 12, height = 10, units = "in")
ggsave(filename = paste0(plot_dir, "cf_gdfp_topology_bar.png"), plot = gdfp_topology_bars_quilt, width = 12, height = 10, units = "in")

## PLOT gCF 3: gCF/gDFP boxplots across datasets
# Create plot panels
gcf_boxplot <- ggplot(data = gdfp_df, 
                      aes(x = tree_topology_short, y = KEY_gCF, fill = tree_topology_short)) +
  geom_violin(draw_quantiles = TRUE) + 
  geom_boxplot(width = 0.07) +
  scale_x_discrete(name = "Constrained topology") +
  scale_y_continuous(name = "Value", breaks = seq(0,100,20), labels = seq(0,100,20), minor_breaks = seq(0,100,5), limits = c(0,100)) +
  scale_fill_manual(name = "Topology", values = c("CTEN" = topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                  "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]], "PARAPHYLY" = topology_colours[["Paraphyly"]]), guide = NULL) +
  labs(title = "Gene Concordance Factors") +
  theme_bw() +
  theme(plot.title = element_text(size = 25, hjust = 0.5),
        axis.title.x = element_text(size = 20, margin = margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.text.x = element_text(size = 16, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 20, margin = margin(t=0, r=10, b=0, l=0, unit="pt")),
        axis.text.y = element_text(size = 16))
ggsave(filename = paste0(plot_dir, "cf_gdfp_boxplot.pdf"), plot = gcf_boxplot)
ggsave(filename = paste0(plot_dir, "cf_gdfp_boxplot.png"), plot = gcf_boxplot)



##### 6. quartet concordance factors (qCF) ternary plots #####
## PLOT qCF 1: qCF values as ternary plots across dataset
# aes: x = CTEN, y = PORI, z = CTENPORI
# labs: L = CTEN, T = PORI, R = CTENPORI
qcf_tern <- ggtern(qcf_long, mapping = aes(x = CTEN.KEY_q1, y = PORI.KEY_q1, z = CTENPORI.KEY_q1)) +
  geom_Lline(Lintercept = 0.33, color = "black", linetype = "dashed") +
  geom_Tline(Tintercept = 0.33, color = "black", linetype = "dashed") +
  geom_Rline(Rintercept = 0.33, color = "black", linetype = "dashed") +
  geom_point(size = 5, alpha = 0.6, color = mako_palette[2]) +
  scale_L_continuous(name = "CTEN", breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), limits = c(0,1)) +
  scale_T_continuous(name = "PORI", breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), limits = c(0,1)) +
  scale_R_continuous(name = "CTEN\u002B\nPORI", breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), limits = c(0,1)) +
  labs(title = "qCF") +
  theme_bw() +
  theme(plot.title = element_text(size = 25, hjust = 0.5, face = "bold"),
        tern.panel.grid.major = element_line(colour = "darkgrey", linewidth = 0.6),
        tern.panel.grid.minor = element_line(colour = "darkgrey", linewidth = 0.3),
        tern.axis.title.L = element_text(size = 22, hjust = 0, vjust = 2),
        tern.axis.title.T = element_text(size = 22, hjust = 0.5, vjust = -0.5),
        tern.axis.title.R = element_text(size = 22, hjust = 0.5, vjust = 1),
        tern.axis.text = element_text(size = 18) )
ggsave(filename = paste0(plot_dir, "cf_qcf_ternary.pdf"), plot = qcf_tern, width = 9, height = 9, units = "in")
ggsave(filename = paste0(plot_dir, "cf_qcf_ternary.png"), plot = qcf_tern, width = 9, height = 9, units = "in")

## PLOT qCF 2: qCF values as ternary plots for the three constrained tree topologies
qcf_topology_tern <- ggtern(qcf_df, mapping = aes(x = KEY_q2, y = KEY_q1, z = KEY_q3, color = tree_topology_short)) +
  geom_Lline(Lintercept = 0.33, color = "black", linetype = "dashed") +
  geom_Tline(Tintercept = 0.33, color = "black", linetype = "dashed") +
  geom_Rline(Rintercept = 0.33, color = "black", linetype = "dashed") +
  geom_point(size = 6, alpha = 0.6) +
  scale_L_continuous(name = "qDF1", breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), limits = c(0,1)) +
  scale_T_continuous(name = "qCF", breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), limits = c(0,1)) +
  scale_R_continuous(name = "qDF2", breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), limits = c(0,1)) +
  scale_color_manual(name = "Topology", values = c("CTEN" = topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                   "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]])) +
  theme_bw() +
  theme(tern.panel.grid.major = element_line(colour = "darkgrey", linewidth = 0.6),
        tern.panel.grid.minor = element_line(colour = "darkgrey", linewidth = 0.3),
        tern.axis.title.L = element_text(size = 22, hjust = 0, vjust = 2),
        tern.axis.title.T = element_text(size = 22, vjust = -0.5),
        tern.axis.title.R = element_text(size = 22, hjust = 0.6, vjust = 2),
        tern.axis.text = element_text(size = 18),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 18))
ggsave(filename = paste0(plot_dir, "cf_qcf_topology_ternary.pdf"), plot = qcf_topology_tern)
ggsave(filename = paste0(plot_dir, "cf_qcf_topology_ternary.png"), plot = qcf_topology_tern)

## PLOT qCF 3: qCF values as bars for each dataset
# Create plot panels
qcf_topology_bars_panel1 <- ggplot(data = qcf_df[which(qcf_df$dataset_id %in% dataset_group_1), ], 
                                   aes(x = tree_topology_short, y = KEY_q1, fill = tree_topology_short)) +
  geom_bar(stat="identity") +
  facet_grid(model_formatted~dataset_id_formatted) +
  scale_x_discrete(name = "Constrained topology") +
  scale_y_continuous(name = "Value", breaks = seq(0,0.6,0.1), labels = seq(0,0.6,0.1), minor_breaks = seq(0,0.6,0.05), limits = c(0,0.5)) +
  scale_fill_manual(name = "Topology", values = c("CTEN" = topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                  "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]])) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title.x = element_text(size = 16, margin = margin(t=10, r=0, b=20, l=0, unit="pt")),
        axis.text.x = element_text(size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 16, margin = margin(t=0, r=10, b=0, l=0, unit="pt")),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 16, margin = margin(t=0, r=0, b=10, l=0, unit="pt")),
        legend.text = element_text(size = 13))
qcf_topology_bars_panel2 <- ggplot(data = qcf_df[which(qcf_df$dataset_id %in% dataset_group_2), ], 
                                   aes(x = tree_topology_short, y = KEY_q1, fill = tree_topology_short)) +
  geom_bar(stat="identity") +
  facet_grid(model_formatted~dataset_id_formatted) + 
  scale_x_discrete(name = "Constrained topology") +
  scale_y_continuous(name = "Value", breaks = seq(0,0.6,0.1), labels = seq(0,0.6,0.1), minor_breaks = seq(0,0.6,0.05), limits = c(0,0.5)) +
  scale_fill_manual(name = "Topology", values = c("CTEN" =  topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                  "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]])) +
  theme_bw() +
  theme(strip.text = element_text(size = 15),
        axis.title.x = element_text(size = 16, margin = margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.text.x = element_text(size = 13, angle = 45, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 16, margin = margin(t=0, r=10, b=0, l=0, unit="pt")),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 16, margin = margin(t=0, r=0, b=10, l=0, unit="pt")),
        legend.text = element_text(size = 13))
# Assemble patchwork and save plot
qcf_topology_bars_quilt <- qcf_topology_bars_panel1 / qcf_topology_bars_panel2 + 
  plot_annotation(title = "Quartet Concordance Factors",
                  theme = theme(plot.title = element_text(size = 25, vjust = 0.5, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0, unit = "pt"))))
ggsave(filename = paste0(plot_dir, "cf_qcf_topology_bar.pdf"), plot = qcf_topology_bars_quilt, width = 12, height = 10, units = "in")
ggsave(filename = paste0(plot_dir, "cf_qcf_topology_bar.png"), plot = qcf_topology_bars_quilt, width = 12, height = 10, units = "in")



##### 7. Concordance factors and branch lengths #####
## Construct data frame of qCF/gCF and branch lengths
gcfdfp_bl_df          <- gdfp_df[, c("analysis", "id", "dataset", "matrix_name", "dataset_id", "dataset_id_formatted", "model", "model_formatted",
                                     "tree_topology", "tree_topology_formatted", "tree_topology_short", "Plac_present", "Plac_present_formatted",
                                     "KEY_gCF","KEY_Length", "CTEN_Length", "CTEN_monophyly", "PORI_Length", "PORI_monophyly")]
names(gcfdfp_bl_df)   <- c("analysis", "id", "dataset", "matrix_name", "dataset_id", "dataset_id_formatted", "model", "model_formatted",
                           "tree_topology", "tree_topology_formatted", "tree_topology_short", "Plac_present", "Plac_present_formatted",
                           "KEY_CF","KEY_Length", "CTEN_Length", "CTEN_monophyly", "PORI_Length", "PORI_monophyly")
gcf_bl_df <- gcfdfp_bl_df[which(gcfdfp_bl_df$tree_topology != "Paraphyly"), ]
qcf_bl_df             <- qcf_df[ , c("analysis", "id", "dataset", "matrix_name", "dataset_id", "dataset_id_formatted", "model", "model_formatted",
                                     "topology", "tree_topology_formatted", "tree_topology_short", "Plac_present", "Plac_present_formatted",
                                     "KEY_q1", "KEY_branch_length", "CTEN_branch_length", "CTEN_monophyly", "PORI_branch_length", "PORI_monophyly")]
names(qcf_bl_df)      <- c("analysis", "id", "dataset", "matrix_name", "dataset_id", "dataset_id_formatted", "model", "model_formatted",
                           "tree_topology", "tree_topology_formatted", "tree_topology_short", "Plac_present", "Plac_present_formatted",
                           "KEY_CF","KEY_Length", "CTEN_Length", "CTEN_monophyly", "PORI_Length", "PORI_monophyly")
bl_df                 <- rbind(gcf_bl_df, qcf_bl_df)

# Set theming
branch_length_theming <- theme_bw() +
  theme(strip.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, margin = margin(t=10, r=0, b=0, l=0, unit="pt")),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_text(size = 16, margin = margin(t=0, r=10, b=0, l=0, unit="pt")),
        axis.text.y = element_text(size = 13),
        legend.title = element_text(size = 16, margin = margin(t=0, r=0, b=10, l=0, unit="pt")),
        legend.text = element_text(size = 13),
        legend.position = "none",
        panel.spacing = unit(8, "pt") )

## PLOT 1: CTEN branch lengths against CF value
# Create gcf panel
gcf_cten_panel <- ggplot(gcf_bl_df, aes(x = CTEN_Length, y = KEY_CF, color = tree_topology_short)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y~x) +
  facet_grid(model_formatted~tree_topology_short) +
  scale_x_continuous(name = "Ctenophora branch length", breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25), minor_breaks =  seq(0, 1, 0.125), limits = c(0, 1)) +
  scale_y_continuous(name = "gCF value", breaks = seq(0, 25, 5), labels = seq(0, 25, 5), minor_breaks =  seq(0, 25, 2.5), limits = c(0, 25)) +
  scale_color_manual(name = "Topology", values = c("CTEN" =  topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                   "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]])) +
  branch_length_theming
gcf_cten_stats <- gcf_cten_panel + stat_cor(label.y = 24, color = "black")
# Create qcf panel
qcf_cten_panel <- ggplot(qcf_bl_df, aes(x = CTEN_Length, y = KEY_CF, color = tree_topology_short)) +
  geom_smooth(method = "lm", formula = y~x) +
  geom_point(size = 2, alpha = 0.7) +
  facet_grid(model_formatted~tree_topology_short) +
  scale_x_continuous(name = "Ctenophora branch length", breaks = seq(0, 4.5, 1), labels = seq(0, 4.5, 1), minor_breaks =  seq(0, 4.5, 0.5), limits = c(0, 4.5)) +
  scale_y_continuous(name = "qCF value", breaks = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1), minor_breaks =  seq(0, 0.6, 0.1), limits = c(0, 0.6)) +
  scale_color_manual(name = "Topology", values = c("CTEN" =  topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                   "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]])) +
  branch_length_theming
qcf_cten_stats <- qcf_cten_panel + stat_cor(label.y = 0.57, color = "black")
# Assemble plots patchwork and save plot
bl_cten_quilt <- gcf_cten_panel / qcf_cten_panel + 
  plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
ggsave(filename = paste0(plot_dir, "cf_bl_cten.pdf"), plot = bl_cten_quilt, width = 10, units = "in")
ggsave(filename = paste0(plot_dir, "cf_bl_cten.png"), plot = bl_cten_quilt, width = 10, units = "in")
# Assemble stats patchwork and save plot
bl_cten_stat <- gcf_cten_stats / qcf_cten_stats + 
  plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
ggsave(filename = paste0(plot_dir, "cf_bl_cten_stats.pdf"), plot = bl_cten_stat, width = 10, units = "in")
ggsave(filename = paste0(plot_dir, "cf_bl_cten_stats.png"), plot = bl_cten_stat, width = 10, units = "in")

## PLOT 2: KEY branch lengths against CF value
# Create gcf panel
gcf_key_panel <- ggplot(gcf_bl_df, aes(x = KEY_Length, y = KEY_CF, color = tree_topology_short)) +
  geom_smooth(method = "lm", formula = y~x) +
  geom_point(size = 2, alpha = 0.7) +
  facet_grid(model_formatted~tree_topology_short) +
  scale_x_continuous(name = "Key branch length", seq(0, 0.08, 0.02), labels = seq(0, 0.08, 0.02), minor_breaks =  seq(0, 0.08, 0.01), limits = c(0, 0.08, 0.01)) +
  scale_y_continuous(name = "gCF value", breaks = seq(0, 25, 5), labels = seq(0, 25, 5), minor_breaks =  seq(0, 25, 2.5), limits = c(0, 25)) +
  scale_color_manual(name = "Topology", values = c("CTEN" =  topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                   "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]])) +
  branch_length_theming
gcf_key_stats <- gcf_key_panel + stat_cor(label.y = 24, color = "black")
# Create qcf panel
qcf_key_panel <- ggplot(qcf_bl_df, aes(x = KEY_Length, y = KEY_CF, color = tree_topology_short)) +
  geom_smooth(method = "lm", formula = y~x) +
  geom_point(size = 2, alpha = 0.7) +
  facet_grid(model_formatted~tree_topology_short) +
  scale_x_continuous(name = "Key branch length", breaks = seq(0, 0.25, 0.1), labels = seq(0, 0.25, 0.1), minor_breaks =  seq(0, 0.25, 0.025), limits = c(0, 0.25)) +
  scale_y_continuous(name = "qCF value", breaks = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1), minor_breaks =  seq(0, 0.6, 0.1), limits = c(0, 0.6)) +
  scale_color_manual(name = "Topology", values = c("CTEN" =  topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                   "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]])) +
  branch_length_theming
qcf_key_stats <- qcf_key_panel + stat_cor(label.y = 0.57, color = "black")
# Assemble patchwork and save plot
bl_key_quilt <- gcf_key_panel / qcf_key_panel + 
  plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
ggsave(filename = paste0(plot_dir, "cf_bl_key.pdf"), plot = bl_key_quilt, width = 10, units = "in")
ggsave(filename = paste0(plot_dir, "cf_bl_key.png"), plot = bl_key_quilt, width = 10, units = "in")
# Assemble stats patchwork and save plot
bl_key_stat <- gcf_key_stats / qcf_key_stats + 
  plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
ggsave(filename = paste0(plot_dir, "cf_bl_key_stats.pdf"), plot = bl_key_stat, width = 10, units = "in")
ggsave(filename = paste0(plot_dir, "cf_bl_key_stats.png"), plot = bl_key_stat, width = 10, units = "in")

## PLOT 3: PORI branch lengths against CF value
# Create gcf panel
gcf_pori_panel <- ggplot(gcf_bl_df, aes(x = PORI_Length, y = KEY_CF, color = tree_topology_short)) +
  geom_smooth(method = "lm", formula = y~x) +
  geom_point(size = 2, alpha = 0.7) +
  facet_grid(model_formatted~tree_topology_short) +
  scale_x_continuous(name = "Porifera branch length", seq(0, 0.08, 0.04), labels = seq(0, 0.08, 0.04), minor_breaks =  seq(0, 0.08, 0.01), limits = c(0, 0.08)) +
  scale_y_continuous(name = "gCF value", breaks = seq(0, 25, 5), labels = seq(0, 25, 5), minor_breaks =  seq(0, 25, 2.5), limits = c(0, 25)) +
  scale_color_manual(name = "Topology", values = c("CTEN" =  topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                   "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]])) +
  branch_length_theming
gcf_pori_stats <- gcf_pori_panel + stat_cor(label.y = 24, color = "black")
# Create qcf panel
qcf_pori_panel <- ggplot(qcf_bl_df, aes(x = PORI_Length, y = KEY_CF, color = tree_topology_short)) +
  geom_smooth(method = "lm", formula = y~x) +
  geom_point(size = 2, alpha = 0.7) +
  facet_grid(model_formatted~tree_topology_short) +
  scale_x_continuous(name = "Porifera branch length", breaks = seq(0, 0.25, 0.1), labels = seq(0, 0.25, 0.1), minor_breaks =  seq(0, 0.25, 0.025), limits = c(0, 0.25)) +
  scale_y_continuous(name = "qCF value", breaks = seq(0, 0.6, 0.1), labels = seq(0, 0.6, 0.1), minor_breaks =  seq(0, 0.6, 0.1), limits = c(0, 0.6)) +
  scale_color_manual(name = "Topology", values = c("CTEN" =  topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                   "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]])) +
  branch_length_theming
qcf_pori_stats <- qcf_pori_panel + stat_cor(label.y = 0.57, color = "black")
# Assemble patchwork and save plot
bl_pori_quilt <- gcf_pori_panel / qcf_pori_panel + 
  plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
ggsave(filename = paste0(plot_dir, "cf_bl_pori.pdf"), plot = bl_pori_quilt, width = 10, units = "in")
ggsave(filename = paste0(plot_dir, "cf_bl_pori.png"), plot = bl_pori_quilt, width = 10, units = "in")
# Assemble stats patchwork and save plot
bl_pori_stat <- gcf_pori_stats / qcf_pori_stats + 
  plot_annotation(tag_levels = 'a', tag_suffix = ".") & theme(plot.tag = element_text(size = 30))
ggsave(filename = paste0(plot_dir, "cf_bl_pori_stats.pdf"), plot = bl_pori_stat, width = 10, units = "in")
ggsave(filename = paste0(plot_dir, "cf_bl_pori_stats.png"), plot = bl_pori_stat, width = 10, units = "in")





#################### OLD CODE
##### 4. CF df: Empirical gCF, sCF and quartet values, Key branch, 3 Topologies #####
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



##### 5. CF + Key branch + 3 topologies + unconstrained gene trees #####
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


