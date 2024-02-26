# ancient_ILS/code/05_plot_figures_genes.R
## This script plots result figures from simulation output
# Caitlin Cherryh, 2023

###### 1. Input parameters ######
## File paths
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# output_dir                <- output directory to save figures

repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/06_figures/"

control_parameters <- list(add.extra.color.palettes = FALSE,
                           plot.ternary = TRUE,
                           plot.boxplots = TRUE,
                           plot.branch.lengths = TRUE)



###### 2. Open packages and functions ######
## Open packages
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggtern)
library(patchwork)

# Specify colour palettes used within these plots
metazoan_clade_palette  <- c(Bilateria = "#CC79A7", Cnidaria = "#009E73", Ctenophora = "#56B4E9", Porifera = "#E69F00", Outgroup = "#999999")
tonal_palette           <- c("Ctenophora" = "#bdd7e7", "Porifera" =  "#2171b5")
boxplot_bars            <- c("Ctenophora" = "#e7d4e8", "Porifera" =  "#d9f0d3")
boxplot_points          <- c("Ctenophora" = "#762a83", "Porifera" =  "#1b7837")
bl_bars                 <- c("Ctenophora" = "#2171b5", "Porifera" =  "#E69F00")
bl_points               <- c("Ctenophora" = "black", "Porifera" =  "black")

# Extra colour palettes (unused)
if (control_parameters$add.extra.color.palettes == TRUE){
  cividis5          <- c("#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
  palette2          <- c("Ctenophora-sister" = "#fdb863", "Porifera-sister" = "#b2abd2")
  qcf_type_palette  <- c(Actual = "#5ab4ac", Estimated = "#d8b365")
  cbPalette         <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  tree2_palette     <- c("#F0F921FF", "#0D0887FF")
  tree5_palette     <- c("#e66101", "#fdb863", "#f7f7f7", "#b2abd2", "#5e3c99")
  tree2_cividis     <- c(tree5_cividis[1], tree5_cividis[5])
  tree2_tonal       <- c("#bdd7e7", "#2171b5")
  model3_tonal      <- c("#980043", "#df65b0", "#d4b9da")
}



###### 3. Open and prepare csvs for plotting  ######
# Open gene csv results
ll_df           <- read.csv(paste0(repo_dir, "output/results_gene_tree_likelihood.csv"), stringsAsFactors = FALSE)
elw_df          <- read.csv(paste0(repo_dir, "output/results_gene_AU_test.csv"), stringsAsFactors = FALSE)
au_df           <- read.csv(paste0(repo_dir, "output/results_gene_elw.csv"), stringsAsFactors = FALSE)
scf_df          <- read.csv(paste0(repo_dir, "output/results_gene_scf.csv"), stringsAsFactors = FALSE)
species_scf_df  <- read.csv(paste0(repo_dir, "output/empirical_dataset_concordance_factors.csv"),  stringsAsFactors = FALSE)

# Remove Simion 2017 from all dfs (haven't successfully extracted sCF per gene yet OR estimated species trees)
ll_df           <- ll_df[which(ll_df$dataset != "Simion2017"), ]
elw_df          <- elw_df[which(elw_df$dataset != "Simion2017"), ]
au_df           <- au_df[which(au_df$dataset != "Simion2017"), ]
scf_df          <- scf_df[which(scf_df$dataset != "Simion2017"), ]
species_scf_df  <- species_scf_df[which(species_scf_df$dataset != "Simion2017"), ]

# Remove Hejnol 2009 from all dfs (haven't successfully estimated species trees yet)
ll_df           <- ll_df[which(ll_df$dataset != "Hejnol2009"), ]
elw_df          <- elw_df[which(elw_df$dataset != "Hejnol2009"), ]
au_df           <- au_df[which(au_df$dataset != "Hejnol2009"), ]
scf_df          <- scf_df[which(scf_df$dataset != "Hejnol2009"), ]
species_scf_df  <- species_scf_df[which(species_scf_df$dataset != "Hejnol2009"), ]



###### 4. Update and format sCF files for plots ######
if (control_parameters$plot.ternary == TRUE | control_parameters$plot.boxplots == TRUE | control_parameters$plot.branch.lengths == TRUE){
  # Add any missing columns
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
  species_scf_df$dataset_type <- "species"
  scf_df$dataset_type         <- "gene"
  scf_df$ultrafast_bootstrap <- scf_df$ultafast_bootstrap # typo in column name
  
  # Add new column for plot output dataset id
  scf_df$dataset_id_formatted <- factor(scf_df$dataset_id,
                                        levels =  c("Dunn2008.Dunn2008_FixedNames", "Philippe2009.Philippe_etal_superalignment_FixedNames", "Philippe2011.UPDUNN_MB_FixedNames", 
                                                    "Nosenko2013.nonribosomal_9187_smatrix", "Nosenko2013.ribosomal_14615_smatrix", "Ryan2013.REA_EST_includingXenoturbella",
                                                    "Moroz2014.ED3d", "Borowiec2015.Best108", "Chang2015.Chang_AA", 
                                                    "Whelan2015.Dataset10", "Whelan2017.Metazoa_Choano_RCFV_strict", "Laumer2018.Tplx_BUSCOeuk"),
                                        labels = c("Dunn 2008",  "Philippe 2009", "Philippe 2011", 
                                                   "Nosenko 2013\nnonribosomal", "Nosenko 2013\nribosomal", "Ryan 2013",
                                                   "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                                   "Whelan 2015", "Whelan 2017",  "Laumer 2018"),
                                        ordered = TRUE)
  scf_df$tree_topology_formatted <- factor(scf_df$tree_topology,
                                           levels =  c("CTEN", "PORI", "CTEN_PORI"),
                                           labels = c("Ctenophora", "Porifera", "Ctenophora+Porifera"),
                                           ordered = TRUE)
  
  # Remove any branches without sCF, sDF1 and sDF2 values
  scf_trimmed_df <- scf_df[which(is.na(scf_df$sCF) == FALSE & 
                                   is.na(scf_df$sDF1) == FALSE & 
                                   is.na(scf_df$sDF2) == FALSE), ]
}



###### 5. Species and gene sCFs: ternary plots for ALL ANIMALS branch, CTEN/PORI topologies for 12 datasets ######
if (control_parameters$plot.ternary == TRUE){
  # Split into two dataframes (two plots)
  plot_scf_1 <- rbind(scf_trimmed_df[which(scf_trimmed_df$dataset_id %in% c("Dunn2008.Dunn2008_FixedNames", "Philippe2009.Philippe_etal_superalignment_FixedNames", 
                                                                            "Philippe2011.UPDUNN_MB_FixedNames", "Ryan2013.REA_EST_includingXenoturbella") &
                                             scf_trimmed_df$branch_to_clade == "ALL_ANIMALS" &
                                             scf_trimmed_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2")],
                      species_scf_df[which(species_scf_df$dataset_id %in% c("Dunn2008.Dunn2008_FixedNames", "Philippe2009.Philippe_etal_superalignment_FixedNames", 
                                                                            "Philippe2011.UPDUNN_MB_FixedNames", "Ryan2013.REA_EST_includingXenoturbella") &
                                             species_scf_df$branch_to_clade == "ALL_ANIMALS" &
                                             species_scf_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2")] )
  plot_scf_2 <- rbind(scf_trimmed_df[which(scf_trimmed_df$dataset_id %in% c("Nosenko2013.nonribosomal_9187_smatrix", "Nosenko2013.ribosomal_14615_smatrix",
                                                                            "Moroz2014.ED3d", "Borowiec2015.Best108") &
                                             scf_trimmed_df$branch_to_clade == "ALL_ANIMALS" &
                                             scf_trimmed_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2")],
                      species_scf_df[which(species_scf_df$dataset_id %in% c("Nosenko2013.nonribosomal_9187_smatrix", "Nosenko2013.ribosomal_14615_smatrix", 
                                                                            "Moroz2014.ED3d", "Borowiec2015.Best108") &
                                             species_scf_df$branch_to_clade == "ALL_ANIMALS" &
                                             species_scf_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2")] )
  plot_scf_3 <- rbind(scf_trimmed_df[which(scf_trimmed_df$dataset_id %in% c("Chang2015.Chang_AA", "Whelan2015.Dataset10", 
                                                                            "Whelan2017.Metazoa_Choano_RCFV_strict", "Laumer2018.Tplx_BUSCOeuk") &
                                             scf_trimmed_df$branch_to_clade == "ALL_ANIMALS" &
                                             scf_trimmed_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2")],
                      species_scf_df[which(species_scf_df$dataset_id %in% c("Chang2015.Chang_AA", "Whelan2015.Dataset10",
                                                                            "Whelan2017.Metazoa_Choano_RCFV_strict", "Laumer2018.Tplx_BUSCOeuk") &
                                             species_scf_df$branch_to_clade == "ALL_ANIMALS" &
                                             species_scf_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2")] )
  # Create a ternary plot of the sCFs - three separate plots of 4 datasets each
  scf1_plot <- ggtern(data = plot_scf_1, aes(x = sDF1, y = sCF, z = sDF2, shape = dataset_type, color = dataset_type)) + 
    geom_point(alpha = 0.7, size = 3) +
    facet_grid(dataset_id_formatted ~ tree_topology_formatted) +
    labs(title = "Clade: Metazoa") +
    scale_color_manual(values = c("species" = "black", "gene" = "grey80"), name = "sCF Type") +
    scale_shape_manual(values = c("species" = 17, "gene" = 19), name = "sCF Type") + 
    theme_bw() +
    theme(plot.title = element_text(size = 26, hjust = 0.5),
          strip.background = element_blank(), 
          strip.text.x = element_text(size = 22, color = "grey40", margin = margin(t = 5, r = 5, b = 15, l = 5, "pt")),
          strip.text.y = element_text(size = 22, color = "grey40", margin = margin(t = 5, r = 15, b = 5, l = 20, "pt")),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 26), 
          legend.text = element_text(size = 22), 
          legend.key.size = unit(2, "lines"),
          panel.spacing = unit(1, "lines")) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  
  scf2_plot <- ggtern(data = plot_scf_2, aes(x = sDF1, y = sCF, z = sDF2, shape = dataset_type, color = dataset_type)) + 
    geom_point(alpha = 0.7, size = 3) +
    facet_grid(dataset_id_formatted ~ tree_topology_formatted) +
    labs(title = "Clade: Metazoa") +
    scale_color_manual(values = c("species" = "black", "gene" = "grey80"), name = "sCF Type") +
    scale_shape_manual(values = c("species" = 17, "gene" = 19), name = "sCF Type") + 
    theme_bw() +
    theme(plot.title = element_text(size = 26, hjust = 0.5),
          strip.background = element_blank(), 
          strip.text.x = element_text(size = 22, color = "grey40", margin = margin(t = 5, r = 5, b = 15, l = 5, "pt")),
          strip.text.y = element_text(size = 22, color = "grey40", margin = margin(t = 5, r = 15, b = 5, l = 20, "pt")),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 26), 
          legend.text = element_text(size = 22), 
          legend.key.size = unit(2, "lines"),
          panel.spacing = unit(1, "lines")) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  scf3_plot <- ggtern(data = plot_scf_3, aes(x = sDF1, y = sCF, z = sDF2, shape = dataset_type, color = dataset_type)) + 
    geom_point(alpha = 0.7, size = 3) +
    facet_grid(dataset_id_formatted ~ tree_topology_formatted) +
    labs(title = "Clade: Metazoa") +
    scale_color_manual(values = c("species" = "black", "gene" = "grey80"), name = "sCF Type") +
    scale_shape_manual(values = c("species" = 17, "gene" = 19), name = "sCF Type") + 
    theme_bw() +
    theme(plot.title = element_text(size = 26, hjust = 0.5),
          strip.background = element_blank(), 
          strip.text.x = element_text(size = 22, color = "grey40", margin = margin(t = 5, r = 5, b = 15, l = 5, "pt")),
          strip.text.y = element_text(size = 22, color = "grey40", margin = margin(t = 5, r = 15, b = 5, l = 20, "pt")),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 26), 
          legend.text = element_text(size = 22), 
          legend.key.size = unit(2, "lines"),
          panel.spacing = unit(1, "lines")) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  # Save plots
  scf1_plot_name <- paste0(repo_dir, "figures/", "gene_scf_AllAnimals_ternary_plot1.pdf")
  ggsave(filename = scf1_plot_name, plot = scf1_plot, width = 12, height = 12, units = "in")
  scf2_plot_name <- paste0(repo_dir, "figures/", "gene_scf_AllAnimals_ternary_plot2.pdf")
  ggsave(filename = scf2_plot_name, plot = scf2_plot, width = 12, height = 12, units = "in")
  scf3_plot_name <- paste0(repo_dir, "figures/", "gene_scf_AllAnimals_ternary_plot3.pdf")
  ggsave(filename = scf3_plot_name, plot = scf3_plot, width = 12, height = 12, units = "in")
}



###### 6. Species and gene sCFs: ternary plots for ALL ANIMALS branch, CTEN/PORI topologies for 12 datasets ######
if (control_parameters$plot.ternary == TRUE){
  # Split into two dataframes (two plots)
  plot_scf_4 <- rbind(scf_trimmed_df[which(scf_trimmed_df$dataset_id %in% c("Dunn2008.Dunn2008_FixedNames", "Philippe2009.Philippe_etal_superalignment_FixedNames", 
                                                                            "Philippe2011.UPDUNN_MB_FixedNames", "Ryan2013.REA_EST_includingXenoturbella") &
                                             scf_trimmed_df$branch_to_clade == "ALL_OTHER_ANIMALS" &
                                             scf_trimmed_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2")],
                      species_scf_df[which(species_scf_df$dataset_id %in% c("Dunn2008.Dunn2008_FixedNames", "Philippe2009.Philippe_etal_superalignment_FixedNames", 
                                                                            "Philippe2011.UPDUNN_MB_FixedNames", "Ryan2013.REA_EST_includingXenoturbella") &
                                             species_scf_df$branch_to_clade == "ALL_OTHER_ANIMALS" &
                                             species_scf_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2")] )
  plot_scf_5 <- rbind(scf_trimmed_df[which(scf_trimmed_df$dataset_id %in% c("Nosenko2013.nonribosomal_9187_smatrix", "Nosenko2013.ribosomal_14615_smatrix",
                                                                            "Moroz2014.ED3d", "Borowiec2015.Best108") &
                                             scf_trimmed_df$branch_to_clade == "ALL_OTHER_ANIMALS" &
                                             scf_trimmed_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2")],
                      species_scf_df[which(species_scf_df$dataset_id %in% c("Nosenko2013.nonribosomal_9187_smatrix", "Nosenko2013.ribosomal_14615_smatrix", 
                                                                            "Moroz2014.ED3d", "Borowiec2015.Best108") &
                                             species_scf_df$branch_to_clade == "ALL_OTHER_ANIMALS" &
                                             species_scf_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2")] )
  plot_scf_6 <- rbind(scf_trimmed_df[which(scf_trimmed_df$dataset_id %in% c("Chang2015.Chang_AA", "Whelan2015.Dataset10", 
                                                                            "Whelan2017.Metazoa_Choano_RCFV_strict", "Laumer2018.Tplx_BUSCOeuk") &
                                             scf_trimmed_df$branch_to_clade == "ALL_OTHER_ANIMALS" &
                                             scf_trimmed_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2")],
                      species_scf_df[which(species_scf_df$dataset_id %in% c("Chang2015.Chang_AA", "Whelan2015.Dataset10",
                                                                            "Whelan2017.Metazoa_Choano_RCFV_strict", "Laumer2018.Tplx_BUSCOeuk") &
                                             species_scf_df$branch_to_clade == "ALL_OTHER_ANIMALS" &
                                             species_scf_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2")] )
  # Create a ternary plot of the sCFs - three separate plots of 4 datasets each
  scf4_plot <- ggtern(data = plot_scf_4, aes(x = sDF1, y = sCF, z = sDF2, shape = dataset_type, color = dataset_type)) + 
    geom_point(alpha = 0.7, size = 3) +
    facet_grid(dataset_id_formatted ~ tree_topology_formatted) +
    labs(title = "Clade: All Other Animals",
         subtitle = "Branch separating first clade to diverge from all other Metazoan species") +
    scale_color_manual(values = c("species" = "black", "gene" = "grey80"), name = "sCF Type") +
    scale_shape_manual(values = c("species" = 17, "gene" = 19), name = "sCF Type") + 
    theme_bw() +
    theme(plot.title = element_text(size = 26, hjust = 0.5),
          plot.subtitle = element_text(size = 20, hjust = 0.5),
          strip.background = element_blank(), 
          strip.text.x = element_text(size = 22, color = "grey40", margin = margin(t = 5, r = 5, b = 15, l = 5, "pt")),
          strip.text.y = element_text(size = 22, color = "grey40", margin = margin(t = 5, r = 15, b = 5, l = 20, "pt")),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 26), 
          legend.text = element_text(size = 22), 
          legend.key.size = unit(2, "lines"),
          panel.spacing = unit(1, "lines")) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  scf5_plot <- ggtern(data = plot_scf_5, aes(x = sDF1, y = sCF, z = sDF2, shape = dataset_type, color = dataset_type)) + 
    geom_point(alpha = 0.7, size = 3) +
    facet_grid(dataset_id_formatted ~ tree_topology_formatted) +
    labs(title = "Clade: All Other Animals",
         subtitle = "Branch separating first clade to diverge from all other Metazoan species") +
    scale_color_manual(values = c("species" = "black", "gene" = "grey80"), name = "sCF Type") +
    scale_shape_manual(values = c("species" = 17, "gene" = 19), name = "sCF Type") + 
    theme_bw() +
    theme(plot.title = element_text(size = 26, hjust = 0.5),
          plot.subtitle = element_text(size = 20, hjust = 0.5),
          strip.background = element_blank(), 
          strip.text.x = element_text(size = 22, color = "grey40", margin = margin(t = 5, r = 5, b = 15, l = 5, "pt")),
          strip.text.y = element_text(size = 22, color = "grey40", margin = margin(t = 5, r = 15, b = 5, l = 20, "pt")),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 26), 
          legend.text = element_text(size = 22), 
          legend.key.size = unit(2, "lines"),
          panel.spacing = unit(1, "lines")) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  scf6_plot <- ggtern(data = plot_scf_6, aes(x = sDF1, y = sCF, z = sDF2, shape = dataset_type, color = dataset_type)) + 
    geom_point(alpha = 0.7, size = 3) +
    facet_grid(dataset_id_formatted ~ tree_topology_formatted) +
    labs(title = "Clade: All Other Animals",
         subtitle = "Branch separating first clade to diverge from all other Metazoan species") +
    scale_color_manual(values = c("species" = "black", "gene" = "grey80"), name = "sCF Type") +
    scale_shape_manual(values = c("species" = 17, "gene" = 19), name = "sCF Type") + 
    theme_bw() +
    theme(plot.title = element_text(size = 26, hjust = 0.5),
          plot.subtitle = element_text(size = 20, hjust = 0.5),
          strip.background = element_blank(), 
          strip.text.x = element_text(size = 22, color = "grey40", margin = margin(t = 5, r = 5, b = 15, l = 5, "pt")),
          strip.text.y = element_text(size = 22, color = "grey40", margin = margin(t = 5, r = 15, b = 5, l = 20, "pt")),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 26), 
          legend.text = element_text(size = 22), 
          legend.key.size = unit(2, "lines"),
          panel.spacing = unit(1, "lines")) +
    guides(shape = guide_legend(override.aes = list(size = 5)))
  # Save plots
  scf4_plot_name <- paste0(repo_dir, "figures/", "gene_scf_AllOtherAnimals_ternary_plot1.pdf")
  ggsave(filename = scf4_plot_name, plot = scf4_plot, width = 12, height = 12, units = "in")
  scf5_plot_name <- paste0(repo_dir, "figures/", "gene_scf_AllOtherAnimals_ternary_plot2.pdf")
  ggsave(filename = scf5_plot_name, plot = scf5_plot, width = 12, height = 12, units = "in")
  scf6_plot_name <- paste0(repo_dir, "figures/", "gene_scf_AllOtherAnimals_ternary_plot3.pdf")
  ggsave(filename = scf6_plot_name, plot = scf6_plot, width = 12, height = 12, units = "in")
}



###### 7. Species and gene sCFs: box plot for branch to ALL ANIMALS ######
if (control_parameters$plot.boxplots == TRUE){
  # Assemble boxplot dataframe
  boxplot_scf_df_1 <- scf_trimmed_df[which(scf_trimmed_df$branch_to_clade == "ALL_ANIMALS" &
                                             scf_trimmed_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", 
                                       "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2", "sN", "ultrafast_bootstrap") ]
  boxplot_scf_df_1b <- species_scf_df[which(species_scf_df$branch_to_clade == "ALL_ANIMALS" &
                                              species_scf_df$tree_topology %in% c("CTEN", "PORI")), 
                                      c("dataset", "matrix", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", 
                                        "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2", "sN", "ultrafast_bootstrap")]
  # Convert boxplot_scf_df to long format
  long_boxplot_df_1 <- melt(boxplot_scf_df_1,
                            id.vars = c("dataset", "matrix", "dataset_id_formatted", "tree_topology", 
                                        "tree_topology_formatted", "branch_to_clade", "dataset_type"),
                            measure.vars = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap") )
  long_boxplot_df_1$variable_labels <- factor(long_boxplot_df_1$variable,
                                              levels = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap"),
                                              labels = c("sCF", "sDF1", "sDF2","UFB"),
                                              ordered = TRUE)
  long_boxplot_df_1b <- melt(boxplot_scf_df_1b,
                             id.vars = c("dataset", "matrix", "dataset_id_formatted", "tree_topology", 
                                         "tree_topology_formatted", "branch_to_clade", "dataset_type"),
                             measure.vars = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap") )
  long_boxplot_df_1b$variable_labels <- factor(long_boxplot_df_1b$variable,
                                               levels = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap"),
                                               labels = c("sCF", "sDF1", "sDF2","UFB"),
                                               ordered = TRUE)
  long_boxplot_df_1b$color <- factor(as.character(long_boxplot_df_1b$tree_topology_formatted),
                                     levels = c("Ctenophora", "Porifera"),
                                     labels = c(boxplot_points[["Ctenophora"]], boxplot_points[["Porifera"]]),
                                     ordered = FALSE)
  
  # Create faceted boxplot - each facet is a different branch, and each boxplot is a different variable, with colors to denote topology
  boxplot1 <- ggplot(long_boxplot_df_1) +
    geom_boxplot(aes(x = variable_labels, y = value, fill = tree_topology_formatted)) +
    geom_point(data = long_boxplot_df_1b, aes(x = variable_labels, y = value, shape = tree_topology_formatted, color = tree_topology_formatted), size = 3) +
    facet_wrap(dataset_id_formatted ~ .) +
    scale_x_discrete(name = "Site Concordance Factors") +
    scale_y_continuous(name = "Value") +
    scale_fill_manual(values = boxplot_bars, name = "Gene sCF") +
    scale_color_manual(values = boxplot_points, name = "Concatenated sCF") +
    scale_shape_discrete(name = "Concatenated sCF") +
    labs(title = "Clade: Metazoan") +
    theme_bw() +
    theme(plot.title = element_text(size = 22, hjust = 0.5),
          plot.subtitle = element_text(size = 18, hjust = 0.5),
          strip.text = element_text(size = 14),
          axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12, vjust = 1.0, hjust = 1.0, angle = 45),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 12),
          legend.key.size = unit(2, "lines"),) +
    guides(fill = guide_legend(override.aes = list(size = 8)),
           color = guide_legend(override.aes = list(size = 8)))
  # Save plot
  boxplot1_plot_name <- paste0(repo_dir, "figures/", "gene_scf_AllAnimals_boxplot.pdf")
  ggsave(filename = boxplot1_plot_name, plot = boxplot1, width = 10, height = 10, units = "in")
}



###### 8. Species and gene sCFs: box plot for branch to ALL OTHER ANIMALS ######
if (control_parameters$plot.boxplots == TRUE){
  # Assemble boxplot dataframe
  boxplot_scf_df_2 <- scf_trimmed_df[which(scf_trimmed_df$branch_to_clade == "ALL_OTHER_ANIMALS" &
                                             scf_trimmed_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", 
                                       "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2", "sN", "ultrafast_bootstrap") ]
  boxplot_scf_df_2b <- species_scf_df[which(species_scf_df$branch_to_clade == "ALL_OTHER_ANIMALS" &
                                              species_scf_df$tree_topology %in% c("CTEN", "PORI")), 
                                      c("dataset", "matrix", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", 
                                        "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2", "sN", "ultrafast_bootstrap")]
  # Convert boxplot_scf_df to long format
  long_boxplot_df_2 <- melt(boxplot_scf_df_2,
                            id.vars = c("dataset", "matrix", "dataset_id_formatted", "tree_topology", 
                                        "tree_topology_formatted", "branch_to_clade", "dataset_type"),
                            measure.vars = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap") )
  long_boxplot_df_2$variable_labels <- factor(long_boxplot_df_2$variable,
                                              levels = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap"),
                                              labels = c("sCF", "sDF1", "sDF2","UFB"),
                                              ordered = TRUE)
  long_boxplot_df_2b <- melt(boxplot_scf_df_2b,
                             id.vars = c("dataset", "matrix", "dataset_id_formatted", "tree_topology", 
                                         "tree_topology_formatted", "branch_to_clade", "dataset_type"),
                             measure.vars = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap") )
  long_boxplot_df_2b$variable_labels <- factor(long_boxplot_df_2b$variable,
                                               levels = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap"),
                                               labels = c("sCF", "sDF1", "sDF2","UFB"),
                                               ordered = TRUE)
  long_boxplot_df_2b$color <- factor(as.character(long_boxplot_df_2b$tree_topology_formatted),
                                     levels = c("Ctenophora", "Porifera"),
                                     labels = c(boxplot_points[["Ctenophora"]], boxplot_points[["Porifera"]]),
                                     ordered = FALSE)
  # Create faceted boxplot - each facet is a different branch, and each boxplot is a different variable, with colors to denote topology
  boxplot2 <- ggplot(long_boxplot_df_2) +
    geom_boxplot(aes(x = variable_labels, y = value, fill = tree_topology_formatted)) +
    geom_point(data = long_boxplot_df_2b, aes(x = variable_labels, y = value, shape = tree_topology_formatted, color = tree_topology_formatted), size = 3) +
    facet_wrap(dataset_id_formatted ~ .) +
    scale_x_discrete(name = "Site Concordance Factors") +
    scale_y_continuous(name = "Value") +
    scale_fill_manual(values = boxplot_bars, name = "Gene sCF") +
    scale_color_manual(values = boxplot_points, name = "Concatenated sCF") +
    scale_shape_discrete(name = "Concatenated sCF") +
    labs(title = "Clade: All Other Animals",
         subtitle = "Branch separating first clade to diverge from all other Metazoan species") +
    theme_bw() +
    theme(plot.title = element_text(size = 22, hjust = 0.5),
          plot.subtitle = element_text(size = 18, hjust = 0.5),
          strip.text = element_text(size = 14),
          axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12, vjust = 1.0, hjust = 1.0, angle = 45),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 12),
          legend.key.size = unit(2, "lines"),) +
    guides(fill = guide_legend(override.aes = list(size = 8)),
           color = guide_legend(override.aes = list(size = 8)))
  # Save plot
  boxplot2_plot_name <- paste0(repo_dir, "figures/", "gene_scf_AllOtherAnimals_boxplot.pdf")
  ggsave(filename = boxplot2_plot_name, plot = boxplot2, width = 10, height = 10, units = "in")
}



###### 9. Species and gene sCFs: box plot for branch to CTEN ######
if (control_parameters$plot.boxplots == TRUE){
  # Assemble boxplot dataframe
  boxplot_scf_df_3 <- scf_trimmed_df[which(scf_trimmed_df$branch_to_clade == "ALL_OTHER_ANIMALS" &
                                             scf_trimmed_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", 
                                       "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2", "sN", "ultrafast_bootstrap") ]
  boxplot_scf_df_3b <- species_scf_df[which(species_scf_df$branch_to_clade == "ALL_OTHER_ANIMALS" &
                                              species_scf_df$tree_topology %in% c("CTEN", "PORI")), 
                                      c("dataset", "matrix", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", 
                                        "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2", "sN", "ultrafast_bootstrap")]
  # Convert boxplot_scf_df to long format
  long_boxplot_df_3 <- melt(boxplot_scf_df_3,
                            id.vars = c("dataset", "matrix", "dataset_id_formatted", "tree_topology", 
                                        "tree_topology_formatted", "branch_to_clade", "dataset_type"),
                            measure.vars = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap") )
  long_boxplot_df_3$variable_labels <- factor(long_boxplot_df_3$variable,
                                              levels = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap"),
                                              labels = c("sCF", "sDF1", "sDF2","UFB"),
                                              ordered = TRUE)
  long_boxplot_df_3b <- melt(boxplot_scf_df_3b,
                             id.vars = c("dataset", "matrix", "dataset_id_formatted", "tree_topology", 
                                         "tree_topology_formatted", "branch_to_clade", "dataset_type"),
                             measure.vars = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap") )
  long_boxplot_df_3b$variable_labels <- factor(long_boxplot_df_3b$variable,
                                               levels = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap"),
                                               labels = c("sCF", "sDF1", "sDF2","UFB"),
                                               ordered = TRUE)
  long_boxplot_df_3b$color <- factor(as.character(long_boxplot_df_3b$tree_topology_formatted),
                                     levels = c("Ctenophora", "Porifera"),
                                     labels = c(boxplot_points[["Ctenophora"]], boxplot_points[["Porifera"]]),
                                     ordered = FALSE)
  # Create faceted boxplot - each facet is a different branch, and each boxplot is a different variable, with colors to denote topology
  boxplot3 <- ggplot(long_boxplot_df_3) +
    geom_boxplot(aes(x = variable_labels, y = value, fill = tree_topology_formatted)) +
    geom_point(data = long_boxplot_df_3b, aes(x = variable_labels, y = value, shape = tree_topology_formatted, color = tree_topology_formatted), size = 3) +
    facet_wrap(dataset_id_formatted ~ .) +
    scale_x_discrete(name = "Site Concordance Factors") +
    scale_y_continuous(name = "Value") +
    scale_fill_manual(values = boxplot_bars, name = "Gene sCF") +
    scale_color_manual(values = boxplot_points, name = "Concatenated sCF") +
    scale_shape_discrete(name = "Concatenated sCF") +
    labs(title = "Clade: Ctenophora") +
    theme_bw() +
    theme(plot.title = element_text(size = 22, hjust = 0.5),
          plot.subtitle = element_text(size = 18, hjust = 0.5),
          strip.text = element_text(size = 14),
          axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12, vjust = 1.0, hjust = 1.0, angle = 45),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 12),
          legend.key.size = unit(2, "lines"),) +
    guides(fill = guide_legend(override.aes = list(size = 8)),
           color = guide_legend(override.aes = list(size = 8)))
  # Save plot
  boxplot3_plot_name <- paste0(repo_dir, "figures/", "gene_scf_Ctenophora_boxplot.pdf")
  ggsave(filename = boxplot3_plot_name, plot = boxplot3, width = 10, height = 10, units = "in")
}



###### 10. Species and gene sCFs: box plot for branch to PORI ######
if (control_parameters$plot.boxplots == TRUE){
  # Assemble boxplot dataframe
  boxplot_scf_df_4 <- scf_trimmed_df[which(scf_trimmed_df$branch_to_clade == "PORI" &
                                             scf_trimmed_df$tree_topology %in% c("CTEN", "PORI")), 
                                     c("dataset", "matrix", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", 
                                       "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2", "sN", "ultrafast_bootstrap") ]
  boxplot_scf_df_4b <- species_scf_df[which(species_scf_df$branch_to_clade == "PORI" &
                                              species_scf_df$tree_topology %in% c("CTEN", "PORI")), 
                                      c("dataset", "matrix", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", 
                                        "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2", "sN", "ultrafast_bootstrap")]
  # Convert boxplot_scf_df to long format
  long_boxplot_df_4 <- melt(boxplot_scf_df_4,
                            id.vars = c("dataset", "matrix", "dataset_id_formatted", "tree_topology", 
                                        "tree_topology_formatted", "branch_to_clade", "dataset_type"),
                            measure.vars = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap") )
  long_boxplot_df_4$variable_labels <- factor(long_boxplot_df_4$variable,
                                              levels = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap"),
                                              labels = c("sCF", "sDF1", "sDF2","UFB"),
                                              ordered = TRUE)
  long_boxplot_df_4b <- melt(boxplot_scf_df_4b,
                             id.vars = c("dataset", "matrix", "dataset_id_formatted", "tree_topology", 
                                         "tree_topology_formatted", "branch_to_clade", "dataset_type"),
                             measure.vars = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap") )
  long_boxplot_df_4b$variable_labels <- factor(long_boxplot_df_4b$variable,
                                               levels = c("sCF", "sDF1", "sDF2","ultrafast_bootstrap"),
                                               labels = c("sCF", "sDF1", "sDF2","UFB"),
                                               ordered = TRUE)
  long_boxplot_df_4b$color <- factor(as.character(long_boxplot_df_4b$tree_topology_formatted),
                                     levels = c("Ctenophora", "Porifera"),
                                     labels = c(boxplot_points[["Ctenophora"]], boxplot_points[["Porifera"]]),
                                     ordered = FALSE)
  # Create faceted boxplot - each facet is a different branch, and each boxplot is a different variable, with colors to denote topology
  boxplot4 <- ggplot(long_boxplot_df_4) +
    geom_boxplot(aes(x = variable_labels, y = value, fill = tree_topology_formatted)) +
    geom_point(data = long_boxplot_df_4b, aes(x = variable_labels, y = value, shape = tree_topology_formatted, color = tree_topology_formatted), size = 3) +
    facet_wrap(dataset_id_formatted ~ .) +
    scale_x_discrete(name = "Site Concordance Factors") +
    scale_y_continuous(name = "Value") +
    scale_fill_manual(values = boxplot_bars, name = "Gene sCF") +
    scale_color_manual(values = boxplot_points, name = "Concatenated sCF") +
    scale_shape_discrete(name = "Concatenated sCF") +
    labs(title = "Clade: Porifera") +
    theme_bw() +
    theme(plot.title = element_text(size = 22, hjust = 0.5),
          plot.subtitle = element_text(size = 18, hjust = 0.5),
          strip.text = element_text(size = 14),
          axis.title.x = element_text(size = 14, margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt")),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12, vjust = 1.0, hjust = 1.0, angle = 45),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 14), 
          legend.text = element_text(size = 12),
          legend.key.size = unit(2, "lines"),) +
    guides(fill = guide_legend(override.aes = list(size = 8)),
           color = guide_legend(override.aes = list(size = 8)))
  # Save plot
  boxplot4_plot_name <- paste0(repo_dir, "figures/", "gene_scf_Porifera_boxplot.pdf")
  ggsave(filename = boxplot4_plot_name, plot = boxplot4, width = 10, height = 10, units = "in")
}



###### 11. Species and gene branch lengths ######
# Collate all branch lengths into a single dataframe
species_scf_df$branch_length <- species_scf_df$sCF_length
bl_df <- scf_df[ , c("dataset", "matrix", "dataset_id", "dataset_id_formatted", 
                     "tree_topology", "tree_topology_formatted", "branch_to_clade",
                     "dataset_type", "branch_length")]
species_bl_df <- species_scf_df[ , c("dataset", "matrix", "dataset_id", "dataset_id_formatted", 
                                     "tree_topology", "tree_topology_formatted", "branch_to_clade",
                                     "dataset_type", "branch_length")]
# Remove branch lengths with NA values
bl_df <- bl_df[which(is.na(bl_df$branch_length) == FALSE), ]
bl_df <- bl_df[which(bl_df$tree_topology_formatted %in% c("Ctenophora", "Porifera")), ]
bl_df <- bl_df[which(bl_df$branch_to_clade %in% c("ALL_ANIMALS", "ALL_OTHER_ANIMALS", "CTEN", "PORI")), ]
species_bl_df <- species_bl_df[which(is.na(species_bl_df$branch_length) == FALSE), ]
species_bl_df <- species_bl_df[which(species_bl_df$tree_topology_formatted %in% c("Ctenophora", "Porifera")), ]
species_bl_df <- species_bl_df[which(species_bl_df$branch_to_clade %in% c("ALL_ANIMALS", "ALL_OTHER_ANIMALS", "CTEN", "PORI")), ]
# Melt datasets to long format
long_bl_df <- melt(bl_df,
                   id.vars = c("dataset", "matrix", "dataset_id", "dataset_id_formatted", 
                               "tree_topology", "tree_topology_formatted", "branch_to_clade",
                               "dataset_type"),
                   measure.vars = c("branch_length"))
long_species_bl_df <- melt(species_bl_df,
                           id.vars = c("dataset", "matrix", "dataset_id", "dataset_id_formatted", 
                                       "tree_topology", "tree_topology_formatted", "branch_to_clade",
                                       "dataset_type"),
                           measure.vars = c("branch_length"))
# Add new nicely formatted branch_to_clade
long_bl_df$clade_formatted <- factor(long_bl_df$branch_to_clade,
                                     levels = c("ALL_ANIMALS", "ALL_OTHER_ANIMALS", "CTEN", "PORI"),
                                     labels = c("Metazoa", "Other animals", "Ctenophora", "Porifera"))
long_species_bl_df$clade_formatted <- factor(long_species_bl_df$branch_to_clade,
                                             levels = c("ALL_ANIMALS", "ALL_OTHER_ANIMALS", "CTEN", "PORI"),
                                             labels = c("Metazoa", "Other animals", "Ctenophora", "Porifera"))
# Plot branch lengths
bl_plot <- ggplot(long_bl_df) +
  geom_boxplot(aes(x = tree_topology_formatted, y = value, fill = tree_topology_formatted), outlier.colour = "grey40", color = "grey40") +
  geom_point(data = long_species_bl_df, aes(x = tree_topology_formatted, y = value, color = tree_topology_formatted, shape = tree_topology_formatted), size = 4) +
  facet_grid(clade_formatted ~ dataset_id_formatted, scale = "free_y") +
  scale_x_discrete(name = "Constrained tree topology") +
  scale_y_continuous(name = "Branch length (substitutions per site)") +
  scale_fill_manual(values = bl_bars, name = "Genes") +
  scale_color_manual(values = bl_points, name = "Concatenated") +
  scale_shape_discrete(name = "Concatenated") +
  theme_bw() +
  theme(plot.title = element_text(size = 22, hjust = 0.5),
        plot.subtitle = element_text(size = 18, hjust = 0.5),
        strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 14, margin = margin(t = 15, r = 0, b = 0, l = 0, unit = "pt")),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 15, b = 0, l = 0, unit = "pt")),
        axis.text.x = element_text(size = 12, vjust = 1.0, hjust = 1.0, angle = 45),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12),
        legend.key.size = unit(2, "lines")) +
  guides(fill = guide_legend(override.aes = list(size = 8)),
         color = guide_legend(override.aes = list(size = 8)))
# Save plot
bl_plot_name <- paste0(repo_dir, "figures/", "GeneSpecies_branch_lengths_boxplot.pdf")
ggsave(filename = bl_plot_name, plot = bl_plot, width = 18, height = 12, units = "in")






###### 12. Species and gene sCFs: number of decisive sites c.f. total number of sites ######




###### 13. Log likelihoods ######




###### 14. AU tests ######





