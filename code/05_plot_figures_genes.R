# ancient_ILS/code/05_plot_figures_genes.R
## This script plots result figures from simulation output
# Caitlin Cherryh, 2023

###### 1. Input parameters ######
## File paths
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# output_dir                <- output directory to save figures

repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/06_figures/"

control_parameters <- list(add.extra.color.palettes = FALSE)



###### 2. Open packages and functions ######
## Open packages
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggtern)
library(patchwork)

# Specify colour palettes used within these plots
metazoan_clade_palette  <- c(Bilateria = "#CC79A7", Cnidaria = "#009E73", Ctenophora = "#56B4E9", Porifera = "#E69F00", Outgroup = "#999999")
palette3                <- c("CTEN" = "#deebf7", "PORI" = "#9ecae1",  "CTEN_PORI" = "#3182bd")

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



###### 4. Update and integrate csv files ######
# Add any missing columns
species_scf_df$dataset_id <- paste0(species_scf_df$dataset, ".", species_scf_df$matrix) 
species_scf_df$dataset_id_formatted <- factor(species_scf_df$dataset_id,
                                              levels =  c("Dunn2008.Dunn2008_FixedNames", "Philippe2009.Philippe_etal_superalignment_FixedNames", "Philippe2011.UPDUNN_MB_FixedNames", 
                                                          "Ryan2013.REA_EST_includingXenoturbella", "Nosenko2013.nonribosomal_9187_smatrix", "Nosenko2013.ribosomal_14615_smatrix",
                                                          "Moroz2014.ED3d", "Borowiec2015.Best108", "Chang2015.Chang_AA", 
                                                          "Whelan2015.Dataset10", "Whelan2017.Metazoa_Choano_RCFV_strict", "Laumer2018.Tplx_BUSCOeuk"),
                                              labels = c("Dunn 2008",  "Philippe 2009", "Philippe 2011", 
                                                         "Ryan 2013", "Nosenko 2013\nnonribosomal", "Nosenko 2013\nribosomal", 
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
                                                  "Ryan2013.REA_EST_includingXenoturbella", "Nosenko2013.nonribosomal_9187_smatrix", "Nosenko2013.ribosomal_14615_smatrix",
                                                  "Moroz2014.ED3d", "Borowiec2015.Best108", "Chang2015.Chang_AA", 
                                                  "Whelan2015.Dataset10", "Whelan2017.Metazoa_Choano_RCFV_strict", "Laumer2018.Tplx_BUSCOeuk"),
                                      labels = c("Dunn 2008",  "Philippe 2009", "Philippe 2011", 
                                                 "Ryan 2013", "Nosenko 2013\nnonribosomal", "Nosenko 2013\nribosomal", 
                                                 "Moroz 2014", "Borowiec 2015", "Chang 2015", 
                                                 "Whelan 2015", "Whelan 2017",  "Laumer 2018"),
                                      ordered = TRUE)
scf_df$tree_topology_formatted <- factor(scf_df$tree_topology,
                                         levels =  c("CTEN", "PORI", "CTEN_PORI"),
                                         labels = c("Ctenophora", "Porifera", "Ctenophora+Porifera"),
                                         ordered = TRUE)



###### 5. Species and gene sCFs: ternary plots for ALL ANIMALS branch, CTEN/PORI topologies for 12 datasets ######
# Remove any branches without sCF, sDF1 and sDF2 values
scf_trimmed_df <- scf_df[which(is.na(scf_df$sCF) == FALSE & 
                                 is.na(scf_df$sDF1) == FALSE & 
                                 is.na(scf_df$sDF2) == FALSE), ]
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



###### 6. Species and gene sCFs: ternary plots for ALL ANIMALS branch, CTEN/PORI topologies for 12 datasets ######
# Remove any branches without sCF, sDF1 and sDF2 values
scf_trimmed_df <- scf_df[which(is.na(scf_df$sCF) == FALSE & 
                                 is.na(scf_df$sDF1) == FALSE & 
                                 is.na(scf_df$sDF2) == FALSE), ]
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



###### 7. Species and gene sCFs: box plot ######
# Remove any branches without sCF, sDF1 and sDF2 values
scf_trimmed_df <- scf_df[which(is.na(scf_df$sCF) == FALSE & 
                                 is.na(scf_df$sDF1) == FALSE & 
                                 is.na(scf_df$sDF2) == FALSE), ]
# Assemble boxplot dataframe
boxplot_scf_df <- rbind(scf_trimmed_df[which(scf_trimmed_df$branch_to_clade == "ALL_OTHER_ANIMALS" &
                                               scf_trimmed_df$tree_topology %in% c("CTEN", "PORI")), 
                                       c("dataset", "matrix", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", 
                                         "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2", "sN", "ultrafast_bootstrap") ],
                        species_scf_df[which(species_scf_df$branch_to_clade == "ALL_OTHER_ANIMALS" &
                                               species_scf_df$tree_topology %in% c("CTEN", "PORI")), 
                                       c("dataset", "matrix", "dataset_id_formatted", "tree_topology", "tree_topology_formatted", 
                                         "branch_to_clade", "dataset_type", "sCF", "sDF1", "sDF2", "sN", "ultrafast_bootstrap")] )


