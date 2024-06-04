# ancient_ILS
#### Evaluating the impact of conflicting phylogenetic signal on hypotheses of early animal evolution

Caitlin Cherryh

June 2024

***
### Summary
The evolutionary relationships between the Metazoa have been particularly difficult to resolve, due to the deep evolutionary timescales and due to a rapid radiation at the root of the animal tree. Both incomplete lineage sorting and long branch attraction have been hypothesised to impact metazoan tree inference. By applying concordance factors, the topological variation at key branches within the metazoan tree can be quantified to determine the evolutionary processes involved. This repository contains the scripts necessary to perform this analyis by calculating and analysing the gene and quartet concordance factors of 12 empirical metazoan phylogenetic datasets.

#### Contents
+ Scripts
    + All scripts necessary to completely replicate this analysis are included in the `code/` folder
    + Each script includes an overview, a list of necessary parameters or file paths,  and a list of software necessary to run that script
+ Output
    + All `.csv` files containing output from each stage of the analysis can be found in the `output/` folder
    + 'dataset_included_taxa.tsv`: Table with one column per empirical dataset. Column names are dataset names, and all taxon names within that dataset are listed one per row
    + `empirical_alignment_files`: Table listing dataset name, dataset ID, alignment file and partition file for each empirical dataset
    + `input_partition_gene_tree_estimation.csv`: Input file for estimating gene trees
    + `input_c60_gene_tree_estimation.csv`: Input file for estimating gene trees
    + `results_gene_tree_likelihood.csv`: Table with log likelihood and BIC score for all gene trees estimated under partition models
    + `input_paths_cf_analysis.csv`: Input file for concordance factor analysis
    + `results_cf_gCF_values.csv`: Gene concordance factor results for all datasets
    + `results_cf_gCF_values_formatted.csv`: Gene concordance factor results for all datasets, formatted differently for plotting
    + `results_cf_qCF_values.csv`: Quartet concordance factor results for all datasets
    + `results_cf_qCF_values_formatted.csv`: Quartet concordance factor results for all datasets, formatted differently for plotting
+ Constraint trees
    + The `constraint_trees` directory contains a file `alternative_phylogenetic_hypotheses.nxe`, which contains Newick trees for the three evolutionary hypotheses considered in this study
+ Instructions for replication
    + Instructions for replicating these analyses are in this `README.md` file.

***
### Instructions to reproduce the analyses:
To fully replicate the analyses, follow these steps:

1. Download and install the software programs and R packages required to run each of the treelikeness metrics
    + IQ-Tree2
    + SplitsTree4 (v4.17.1 or higher)
    + Phylogemetric
    + fast_TIGER
2. Create the conda environment `gene_filtering` using the `environment.yaml` file
3. Prepare simulated alignments
    + The script `code/01_simulations.R` will generate both sets of simulated alignments (random tree simulations and introgression simulations)
    + The random tree simulations mimic decreased treelikeness by increasing the number of evolutionary histories present within an alignment
    + The introgression simulations mimic decreased treelikeness by increasing the proportion of sites within the alignment impacted by a single introgression event
4. Apply treelikeness metrics
    + The script `code/02_apply_treelikeness_metrics.R` will fully replicate our analysis by applying each treelikeness metric to each simulated alignment
    + To apply specific treelikeness metrics to a single alignment, use the functions in the script `code/func_metrics.R`
    + To apply the tree proportion test to a single alignment, use the functions in the script `code/func_tree_proportion.R`
5. Process and analyse results
    + The script `code/03_data_analysis.R` will plot results from the treelikeness metrics for both sets of simulations

***
