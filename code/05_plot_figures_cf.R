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
library(ggpubr) # for adding statistics and regression lines: stat_regline_equation(), stat_cor()
library(patchwork) # for assembling ternary plots of emperical gCF, sCF and quartet scores
library(grDevices) # for: hcl.colors(), palette.colors()
library(stringr) # for: str_extract_all()

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

##### 3. Open and prepare csvs for plotting  #####
## List all output files
all_output_files <- paste0(repo_dir, "output/", list.files(output_csv_dir))

## Open empirical cf analysis results
gcf_df    <- read.csv(grep("gCF_values.csv", all_output_files, value = TRUE))
gcf_long  <- read.csv(grep("gCF_values_formatted.csv", all_output_files, value = TRUE))
qcf_df    <- read.csv(grep("qCF_values.csv", all_output_files, value = TRUE))
qcf_long  <- read.csv(grep("qCF_values_formatted.csv", all_output_files, value = TRUE))
scf_df    <- read.csv(grep("gene_scf.csv", all_output_files, value = TRUE))
scf_df    <- scf_df[which(scf_df$dataset != "Hejnol2009"), ]
logl_df   <- read.csv(grep("gene_tree_likelihood.csv", all_output_files, value = TRUE))
logl_df   <- logl_df[which(logl_df$dataset_id != "Hejnol2009.Hejnol_etal_2009_FixedNames" & 
                             logl_df$dataset_id != "Simion2017.supermatrix_97sp_401632pos_1719genes"), ]



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
dataset_labels_oneline        <- c("Dunn 2008",  "Philippe 2009", "Philippe 2011", 
                                   "Nosenko 2013 nonribo.", "Nosenko 2013 ribo.", "Ryan 2013", 
                                   "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                   "Whelan 2015", "Whelan 2017",  "Laumer 2018")
gcf_df$dataset_id_formatted   <- factor(gcf_df$dataset_id, levels =  dataset_labels, labels = dataset_labels_formatted, ordered = TRUE)
gcf_long$dataset_id_formatted <- factor(gcf_long$dataset_id, levels =  dataset_labels, labels = dataset_labels_formatted, ordered = TRUE)
qcf_df$dataset_id_formatted   <- factor(qcf_df$dataset_id, levels =  dataset_labels, labels = dataset_labels_formatted, ordered = TRUE)
qcf_long$dataset_id_formatted <- factor(qcf_long$dataset_id, levels =  dataset_labels, labels = dataset_labels_formatted, ordered = TRUE)
scf_df$dataset_id_formatted   <- factor(scf_df$dataset_id, levels =  dataset_labels, labels = dataset_labels_formatted, ordered = TRUE)
logl_df$dataset_id_formatted  <- factor(logl_df$dataset_id, levels =  dataset_labels, labels = dataset_labels_oneline, ordered = TRUE)

## Add nicely formatted topology labels to the data frames
topology_labels                   <- c("CTEN", "PORI", "CTEN_PORI")
topology_labels_formatted         <- c("Ctenophora", "Porifera", "Ctenophora+Porifera")
topology_labels_short             <- c("CTEN", "PORI", "CTEN+PORI")
gcf_df$tree_topology_formatted    <- factor(gcf_df$tree_topology, levels =  topology_labels, labels = topology_labels_formatted, ordered = TRUE)
qcf_df$tree_topology_formatted    <- factor(qcf_df$topology, levels =  topology_labels, labels = topology_labels_formatted, ordered = TRUE)
scf_df$tree_topology_formatted    <- factor(scf_df$tree_topology, levels =  topology_labels, labels = topology_labels_formatted, ordered = TRUE)
gcf_df$tree_topology_short        <- factor(gcf_df$tree_topology, levels =  topology_labels, labels = topology_labels_short, ordered = TRUE)
qcf_df$tree_topology_short        <- factor(qcf_df$topology, levels =  topology_labels, labels = topology_labels_short, ordered = TRUE)
scf_df$tree_topology_short        <- factor(scf_df$tree_topology, levels =  topology_labels, labels = topology_labels_short, ordered = TRUE)

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
scf_df$analysis   <- "sCF"

## Create 2 sets of datasets for plotting into nice grids
# Column: dataset_id
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
  geom_smooth(method = "lm", formula = y~x, alpha = 0.25, linewidth = 0,
              aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 25, 25, after_stat(ymax)))) +
  stat_smooth(method = "lm" , formula = y~x, geom = "line", linewidth = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  facet_grid(model_formatted~tree_topology_short) +
  scale_x_continuous(name = "Ctenophora branch length", breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25), minor_breaks =  seq(0, 1, 0.125), limits = c(0, 1)) +
  scale_y_continuous(name = "gCF value", breaks = seq(0, 25, 5), labels = seq(0, 25, 5), minor_breaks =  seq(0, 25, 2.5), limits = c(0, 25)) +
  scale_color_manual(name = "Topology", values = c("CTEN" =  topology_colours[["Ctenophora"]], "PORI" =  topology_colours[["Porifera"]], 
                                                   "CTEN+PORI" = topology_colours[["Ctenophora+Porifera"]])) +
  branch_length_theming
gcf_cten_stats <- gcf_cten_panel + stat_cor(label.y = 24, color = "black")
# Create qcf panel
qcf_cten_panel <- ggplot(qcf_bl_df, aes(x = CTEN_Length, y = KEY_CF, color = tree_topology_short)) +
  geom_smooth(method = "lm", formula = y~x, alpha = 0.25, linewidth = 0,
              aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 0.6, 0.6, after_stat(ymax)))) +
  stat_smooth(method = "lm" , formula = y~x, geom = "line", linewidth = 1, alpha = 0.7) +
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
  geom_smooth(method = "lm", formula = y~x, alpha = 0.25, linewidth = 0,
              aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 25, 25, after_stat(ymax)))) +
  stat_smooth(method = "lm" , formula = y~x, geom = "line", linewidth = 1, alpha = 0.7) +
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
  geom_smooth(method = "lm", formula = y~x, alpha = 0.25, linewidth = 0,
              aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 0.6, 0.6, after_stat(ymax)))) +
  stat_smooth(method = "lm" , formula = y~x, geom = "line", linewidth = 1, alpha = 0.7) +
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
  geom_smooth(method = "lm", formula = y~x, alpha = 0.25, linewidth = 0,
              aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 25, 25, after_stat(ymax)))) +
  stat_smooth(method = "lm" , formula = y~x, geom = "line", linewidth = 1, alpha = 0.7) +
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
  geom_smooth(method = "lm", formula = y~x, alpha = 0.25, linewidth = 0,
              aes(ymin = ifelse(after_stat(ymin) < 0, 0, after_stat(ymin)), ymax = ifelse(after_stat(ymax) > 0.6, 0.6, after_stat(ymax)))) +
  stat_smooth(method = "lm" , formula = y~x, geom = "line", linewidth = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 0.7) +
  facet_grid(model_formatted~tree_topology_short) +
  scale_x_continuous(name = "Porifera branch length", breaks = seq(0, 0.35, 0.1), labels = seq(0, 0.35, 0.1), minor_breaks =  seq(0, 0.35, 0.025), limits = c(0, 0.35)) +
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



##### 8. Creating a nice output table #####
## Reformat the gCF results for a nice table in the manuscript
pretty_gcf <- gcf_long[ c("dataset", "matrix_name", "model",
                          "CTEN.KEY_gCF",  "PORI.KEY_gCF", "CTENPORI.KEY_gCF", "CTEN.KEY_gDFP", 
                          "CTEN.KEY_gCF_N",  "PORI.KEY_gCF_N", "CTENPORI.KEY_gCF_N", "CTEN.KEY_gDFP_N",
                          "CTEN.CTEN_gCF", "PORI.CTEN_gCF", "CTENPORI.CTEN_gCF", 
                          "CTEN.CTEN_gCF_N", "PORI.CTEN_gCF_N", "CTENPORI.CTEN_gCF_N",
                          "CTEN.PORI_gCF", "PORI.PORI_gCF", "CTENPORI.PORI_gCF",
                          "CTEN.PORI_gCF_N", "PORI.PORI_gCF_N", "CTENPORI.PORI_gCF_N",
                          "Plac_present", "analysis")]
names(pretty_gcf) <- c("dataset", "matrix_name", "model",
                       "KEY_gCF_CTEN", "KEY_gCF_PORI", "KEY_gCF_CTENPORI", "KEY_gDFP",
                       "KEY_gCF_N_CTEN", "KEY_gCF_N_PORI", "KEY_gCF_N_CTENPORI", "KEY_gDFP_N",
                       "CTEN_gCF_CTEN", "CTEN_gCF_PORI", "CTEN_gCF_CTENPORI",
                       "CTEN_gCF_N_CTEN", "CTEN_gCF_N_PORI", "CTEN_gCF_N_CTENPORI",
                       "PORI_gCF_CTEN", "PORI_gCF_PORI", "PORI_gCF_CTENPORI",
                       "PORI_gCF_N_CTEN", "PORI_gCF_N_PORI",  "PORI_gCF_N_CTENPORI",
                       "Plac_present", "analysis")
pretty_gcf$year <- as.numeric(unlist(stringr::str_extract_all(pretty_gcf$dataset, "\\d{4}")))
pretty_gcf$model_code <- factor(pretty_gcf$model, levels = c("Partition", "C60"), labels = c(1,2), ordered = T)
pretty_gcf <- pretty_gcf[order(pretty_gcf$year, pretty_gcf$dataset, pretty_gcf$matrix_name, pretty_gcf$model_code, decreasing = F), ]
pretty_gcf_file <- paste0(output_csv_dir, "ms_formatted_gCF_table.csv")
write.csv(pretty_gcf, file = pretty_gcf_file, row.names = FALSE)

## Reformat the qCF results for a nice table in the manuscript
pretty_qcf <- qcf_long[ c("dataset", "matrix_name", "model",
                          "CTEN.KEY_q1", "PORI.KEY_q1", "CTENPORI.KEY_q1",
                          "CTEN.KEY_QC", "PORI.KEY_QC", "CTENPORI.KEY_QC",
                          "CTEN.KEY_EN", "PORI.KEY_EN", "CTENPORI.KEY_EN",
                          "CTEN.CTEN_q1", "PORI.CTEN_q1", "CTENPORI.CTEN_q1", 
                          "CTEN.CTEN_QC", "PORI.CTEN_QC", "CTENPORI.CTEN_QC", 
                          "CTEN.CTEN_EN", "PORI.CTEN_EN", "CTENPORI.CTEN_EN", 
                          "CTEN.PORI_q1", "PORI.PORI_q1", "CTENPORI.PORI_q1",
                          "CTEN.PORI_QC", "PORI.PORI_QC", "CTENPORI.PORI_QC",
                          "CTEN.PORI_EN", "PORI.PORI_EN", "CTENPORI.PORI_EN",
                          "Plac_present", "analysis")]
names(pretty_qcf) <- c("dataset", "matrix_name", "model",
                       "KEY_qCF_CTEN", "KEY_qCF_PORI", "KEY_qCF_CTENPORI",
                       "KEY_nQuartets_CTEN", "KEY_nQuartets_PORI", "KEY_nQuartets_CTENPORI",
                       "KEY_enGenes_CTEN", "KEY_enGenes_PORI", "KEY_enGenes_CTENPORI",
                       "CTEN_qCF_CTEN", "CTEN_qCF_PORI", "CTEN_qCF_CTENPORI", 
                       "CTEN_nQuartets_CTEN", "CTEN_nQuartets_PORI", "CTEN_nQuartets_CTENPORI",
                       "CTEN_enGenes_CTEN", "CTEN_enGenes_PORI", "CTEN_enGenes_CTENPORI",
                       "PORI_qCF_CTEN", "PORI_qCF_PORI", "PORI_qCF_CTENPORI",
                       "PORI_nQuartets_CTEN", "PORI_nQuartets_PORI", "PORI_nQuartets_CTENPORI",
                       "PORI_enGenes_CTEN", "PORI_enGenes_PORI", "PORI_enGenes_CTENPORI",
                       "Plac_present", "analysis")
pretty_qcf$year <- as.numeric(unlist(stringr::str_extract_all(pretty_qcf$dataset, "\\d{4}")))
pretty_qcf$model_code <- factor(pretty_qcf$model, levels = c("Partition", "C60"), labels = c(1,2), ordered = T)
pretty_qcf <- pretty_qcf[order(pretty_qcf$year, pretty_qcf$dataset, pretty_qcf$matrix_name, pretty_qcf$model_code, decreasing = F), ]
pretty_qcf_file <- paste0(output_csv_dir, "ms_formatted_qCF_table.csv")
write.csv(pretty_qcf, file = pretty_qcf_file, row.names = FALSE)



###### 9. Plot gene sCF  ######
# Extract gCF scores into long dataframe
key_scf       <- scf_df[which(scf_df$branch_to_clade == "ALL_OTHER_ANIMALS"), ]
key_scf       <- key_scf[which(key_scf$tree_topology != "CTEN_PORI"), ]
# Plot sCF ternary plots
scf_tern_1 <- ggtern(key_scf[which(key_scf$dataset_id %in% unique(key_scf$dataset_id)[1:4]), ], 
                     aes(x = sDF1, y = sCF, z = sDF2, color = tree_topology_formatted, shape = tree_topology_formatted)) +
  geom_point(size = 4, alpha = 0.6) +
  facet_grid(tree_topology_formatted~dataset_id_formatted) +
  scale_color_manual(name = "Constrained\ntree topology", values = c("Ctenophora" = topology_colours[["Ctenophora"]], "Porifera" = topology_colours[["Porifera"]])) +
  scale_shape_manual(name = "Constrained\ntree topology", values = c("Ctenophora" = 16, "Porifera" = 17)) +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        tern.panel.grid.major = element_line(colour = "darkgrey", linewidth = 0.4),
        tern.axis.title.L = element_text(size = 8),
        tern.axis.title.T = element_text(size = 8),
        tern.axis.title.R = element_text(size = 8),
        tern.axis.text = element_text(size = 8)) +
  guides(color = "none", shape = "none") +
  Tlab("CF") +
  Llab("DF1") +
  Rlab("DF2")
scf_tern_2 <- ggtern(key_scf[which(key_scf$dataset_id %in% unique(key_scf$dataset_id)[5:8]), ], 
                     aes(x = sDF1, y = sCF, z = sDF2, color = tree_topology_formatted, shape = tree_topology_formatted)) +
  geom_point(size = 4, alpha = 0.6) +
  facet_grid(tree_topology_formatted~dataset_id_formatted) +
  scale_color_manual(name = "Constrained\ntree topology", values = c("Ctenophora" = topology_colours[["Ctenophora"]], "Porifera" = topology_colours[["Porifera"]])) +
  scale_shape_manual(name = "Constrained\ntree topology", values = c("Ctenophora" = 16, "Porifera" = 17)) +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        tern.panel.grid.major = element_line(colour = "darkgrey", linewidth = 0.4),
        tern.axis.title.L = element_text(size = 8),
        tern.axis.title.T = element_text(size = 8),
        tern.axis.title.R = element_text(size = 8),
        tern.axis.text = element_text(size = 8)) +
  guides(color = "none", shape = "none") +
  Tlab("CF") +
  Llab("DF1") +
  Rlab("DF2")
scf_tern_3 <- ggtern(key_scf[which(key_scf$dataset_id %in% unique(key_scf$dataset_id)[9:12]), ], 
                     aes(x = sDF1, y = sCF, z = sDF2, color = tree_topology_formatted, shape = tree_topology_formatted)) +
  geom_point(size = 4, alpha = 0.6) +
  facet_grid(tree_topology_formatted~dataset_id_formatted) +
  scale_color_manual(name = "Constrained\ntree topology", values = c("Ctenophora" = topology_colours[["Ctenophora"]], "Porifera" = topology_colours[["Porifera"]])) +
  scale_shape_manual(name = "Constrained\ntree topology", values = c("Ctenophora" = 16, "Porifera" = 17)) +
  theme_bw() +
  theme(strip.text = element_text(size = 12),
        tern.panel.grid.major = element_line(colour = "darkgrey", linewidth = 0.4),
        tern.axis.title.L = element_text(size = 8),
        tern.axis.title.T = element_text(size = 8),
        tern.axis.title.R = element_text(size = 8),
        tern.axis.text = element_text(size = 8)) +
  guides(color = "none", shape = "none") +
  Tlab("CF") +
  Llab("DF1") +
  Rlab("DF2")
# Assemble the three ternary plots using ggtern::grid.arrange 
#     (as ternary plots have three axes, patchwork and ggplot::grid.arrange don't work cleanly here)
scf_quilt <- ggtern::grid.arrange(scf_tern_1, scf_tern_2, scf_tern_3, nrow = 3, ncol = 1)
ggsave(filename = paste0(plot_dir, "cf_gene_scf.pdf"), plot = scf_quilt, width = 10, height = 12, units = "in")
ggsave(filename = paste0(plot_dir, "cf_gene_scf.png"), plot = scf_quilt, width = 10, height = 12, units = "in")



###### 9. Plot gene log-likelihood  ######
# Extract the best logl for each row
extract.best.BIC.topology <- function(i, logl_df){
  i_row <- logl_df[i, ]
  if (is.na(i_row$CTEN_BIC) & is.na(i_row$PORI_BIC) & is.na(i_row$CTEN_PORI_BIC)){
    # Missing all BIC
    new_i_row <- c(as.character(i_row[1:14]), as.character(i_row$dataset_id_formatted), "NONE",
                   NA, NA, NA, NA, NA, NA)
  } else if ((i_row$CTEN_BIC < i_row$PORI_BIC) & (i_row$CTEN_BIC < i_row$CTEN_PORI_BIC)){
    # CTEN = best (lowest) BIC
    new_i_row <- as.character(c(i_row[1:14], as.character(i_row$dataset_id_formatted), "PORI",
                                i_row$CTEN_LogL, i_row$CTEN_Unconstrained_LogL, i_row$CTEN_NumFreeParams,
                                i_row$CTEN_BIC, i_row$CTEN_TreeLength, i_row$CTEN_SumInternalBL))
  } else if ((i_row$PORI_BIC < i_row$CTEN_BIC) & (i_row$PORI_BIC < i_row$CTEN_PORI_BIC)){
    # PORI = best (lowest) BIC
    new_i_row <- as.character(c(i_row[1:14], as.character(i_row$dataset_id_formatted), "PORI",
                                i_row$PORI_LogL, i_row$PORI_Unconstrained_LogL, i_row$PORI_NumFreeParams,
                                i_row$PORI_BIC, i_row$PORI_TreeLength, i_row$PORI_SumInternalBL))
  } else if ((i_row$CTEN_PORI_BIC < i_row$PORI_BIC) & (i_row$CTEN_PORI_BIC < i_row$CTEN_BIC)){
    # CTEN_PORI = best (lowest) BIC
    new_i_row <- as.character(c(i_row[1:14], as.character(i_row$dataset_id_formatted), "PORI",
                                i_row$CTEN_PORI_LogL, i_row$CTEN_PORI_Unconstrained_LogL, i_row$CTEN_PORI_NumFreeParams,
                                i_row$CTEN_PORI_BIC, i_row$CTEN_PORI_TreeLength, i_row$CTEN_PORI_SumInternalBL))
    
  } else if ((i_row$CTEN_BIC == i_row$PORI_BIC) & (i_row$CTEN_BIC < i_row$CTEN_PORI_BIC)){
    # Tie: CTEN and PORI
    new_i_row <- as.character(c(i_row[1:14], as.character(i_row$dataset_id_formatted), "CTEN_and_PORI",
                                i_row$CTEN_PORI_LogL, i_row$CTEN_PORI_Unconstrained_LogL, i_row$CTEN_PORI_NumFreeParams,
                                i_row$CTEN_PORI_BIC, i_row$CTEN_PORI_TreeLength, i_row$CTEN_PORI_SumInternalBL))
  } else if ((i_row$CTEN_BIC == i_row$CTEN_PORI_BIC) & (i_row$CTEN_BIC < i_row$PORI_BIC)){
    # Tie: CTEN and CTEN_PORI
    new_i_row <- as.character(c(i_row[1:14], as.character(i_row$dataset_id_formatted), "CTEN_and_CTENPORI",
                                i_row$CTEN_PORI_LogL, i_row$CTEN_PORI_Unconstrained_LogL, i_row$CTEN_PORI_NumFreeParams,
                                i_row$CTEN_PORI_BIC, i_row$CTEN_PORI_TreeLength, i_row$CTEN_PORI_SumInternalBL))
  } else if ((i_row$PORI_BIC == i_row$CTEN_PORI_BIC) & (i_row$PORI_BIC < i_row$CTEN_BIC)) {
    # Tie: PORI and CTEN_PORI
    new_i_row <- as.character(c(i_row[1:14], as.character(i_row$dataset_id_formatted), "PORI_and_CTENPORI",
                                i_row$CTEN_PORI_LogL, i_row$CTEN_PORI_Unconstrained_LogL, i_row$CTEN_PORI_NumFreeParams,
                                i_row$CTEN_PORI_BIC, i_row$CTEN_PORI_TreeLength, i_row$CTEN_PORI_SumInternalBL))
  } else if ((i_row$CTEN_BIC == i_row$PORI_BIC) & (i_row$CTEN_BIC == i_row$CTEN_PORI_BIC) & (i_row$PORI_BIC == i_row$CTEN_PORI_BIC)) {
    # Tie: CTEN and PORI and CTEN_PORI
    new_i_row <- as.character(c(i_row[1:14], as.character(i_row$dataset_id_formatted), "CTEN_and_PORI_and_CTENPORI",
                                i_row$CTEN_PORI_LogL, i_row$CTEN_PORI_Unconstrained_LogL, i_row$CTEN_PORI_NumFreeParams,
                                i_row$CTEN_PORI_BIC, i_row$CTEN_PORI_TreeLength, i_row$CTEN_PORI_SumInternalBL))
  } else {
    # No single best BIC
    new_i_row <- c(as.character(i_row[1:14]), as.character(i_row$dataset_id_formatted), "TIE",
                   NA, NA, NA, NA, NA, NA)
  }
  names(new_i_row) <- c(names(i_row[1:14]), "dataset_id_formatted", "best_topology", 
                        "BEST_LogL", "BEST_Unconstrained_LogL", "BEST_NumFreeParams", 
                        "BEST_BIC", "BEST_TreeLength", "BEST_SumInternalBL")
  return(new_i_row)
}
# Call function
best_list <- lapply(1:nrow(logl_df), extract.best.BIC.topology, logl_df = logl_df)
best_df <- as.data.frame(do.call(rbind, best_list))
# Format best_df
best_df$best_topology_formatted <- factor(best_df$best_topology,
                                          levels = c("PORI", "CTEN_and_PORI", "CTEN_and_CTENPORI", "PORI_and_CTENPORI", "CTEN_and_PORI_and_CTENPORI", "NONE"),
                                          label = c("PORI", "CTEN & PORI", "CTEN & CTEN+PORI", "PORI & CTEN+PORI", "All equal", "None"),
                                          ordered = TRUE)
best_df$BEST_BIC    <- as.numeric(best_df$BEST_BIC)
best_df$BEST_LogL   <- as.numeric(best_df$BEST_LogL)
best_df$ML_logl     <- as.numeric(best_df$ML_logl)
best_df$LogL_diff   <- best_df$ML_logl - best_df$BEST_LogL
# Remove "None" rows
best_df <- best_df[which(best_df$best_topology != "NONE"), ]
# Plot best_df as bar charts
best_topology_plot <- ggplot(data = best_df, aes(x = dataset_id_formatted, fill = best_topology_formatted)) +
  geom_bar(position = "fill") + 
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "Percent of genes (%)", breaks = seq(0,1,0.2), labels = seq(0,1,0.2), minor_breaks = seq(0,1,0.1)) +
  scale_fill_viridis_d(name = "Best topology\n(by BIC)", option = "H") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 16, margin = margin(r = 10)),
        axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1.0, margin = margin(b = 5)),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16, hjust = 0.5),
        legend.text = element_text(size = 12) )
ggsave(filename = paste0(plot_dir, "cf_gene_BIC_topology.pdf"), plot = best_topology_plot)
ggsave(filename = paste0(plot_dir, "cf_gene_BIC_topology.png"), plot = best_topology_plot)

# Plot difference in BIC
bic_long <- melt(logl_df, 
                 id.vars = names(logl_df)[1:14],
                 measure.vars = c("CTEN_BIC", "PORI_BIC", "CTEN_PORI_BIC") )
bic_long$topology <- factor(bic_long$variable,
                            levels = c("CTEN_BIC", "PORI_BIC", "CTEN_PORI_BIC"),
                            labels = c("CTEN", "PORI", "CTEN+PORI"),
                            ordered = T)
# Plot histogram of BIC values
bic_histogram <- ggplot(data = bic_long, aes(x = value, fill = topology)) +
  geom_histogram(bins = 50) +
  facet_grid(topology~.) +
  scale_x_continuous(name = "BIC value", breaks = c(0, 50000, 100000, 150000), labels = c("0", "50K", "100K", "150K"), limits = c(0,150000)) +
  scale_y_continuous(name = "Count") +
  scale_fill_manual(name = "Constrained\ngene topology", values = c("CTEN" = topology_colours[["Ctenophora"]], 
                                                                    "PORI" = topology_colours[["Porifera"]], 
                                                                    "CTEN+PORI"= topology_colours[["Ctenophora+Porifera"]]) ) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14, margin = margin(t=10)),
        axis.title.y = element_text(size = 14, margin = margin(r=10)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14, hjust = 0.5, margin = margin(b=10)),
        legend.text = element_text(size = 12) )
ggsave(filename = paste0(plot_dir, "cf_gene_BIC_histogram.pdf"), plot = bic_histogram)
ggsave(filename = paste0(plot_dir, "cf_gene_BIC_histogram.png"), plot = bic_histogram)



