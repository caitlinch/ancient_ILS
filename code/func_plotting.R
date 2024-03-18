# ancient_ILS/code/func_plotting.R
## This script contains functions to make plotting easier or prettier
# Caitlin Cherryh, 2024



extract.correct.CF.values <- function(df){
  # Function that given a dataframe, turns the CF values for 3 tree topologies into the
  #   comparable CF/DF1/DF2 values for 1 tree, assuming that the 3 trees are each unique 
  #   arrangements of the same 4 clades around the same branch
  
  # Reformat dataframe into new one with CF/DF1/DF2 values
  new_df <- data.frame(dataset = df$dataset[1:12], 
                       matrix = df$matrix[1:12], 
                       dataset_id_formatted = df$dataset_id_formatted[1:12], 
                       hypothesis_tree = df$hypothesis_tree[1:12],
                       tree_topology_formatted = df$tree_topology_formatted[1:12],
                       branch_to_clade = df$branch_to_clade[1:12],
                       clade_formatted = df$clade_formatted[1:12],
                       gene_type = df$gene_type[1:12],
                       gene_type_formatted = df$gene_type_formatted[1:12],
                       gCF = df$gCF[1:12],
                       gCF_N = df$gCF_N[1:12],
                       gDF1 = df$gCF[13:24],
                       gDF1_N = df$gCF_N[13:24],
                       gDF2 = df$gCF[25:36],
                       gDF2_N = df$gCF_N[25:36],
                       sCF = df$sCF[1:12],
                       sCF_N = df$sCF_N[1:12],
                       sDF1 = df$sCF[13:24],
                       sDF1_N = df$sCF_N[13:24],
                       sDF2 = df$sCF[25:36],
                       sDF2_N = df$sCF_N[25:36],
                       q1 = df$q1[1:12],
                       q2 = df$q1[13:24],
                       q3 = df$q1[25:36],
                       f1 = df$f1[1:12],
                       f2 = df$f1[13:24],
                       f3 = df$f1[25:36],
                       pp1 = df$pp1[1:12],
                       pp2 = df$pp1[13:24],
                       pp3 = df$pp1[25:36],
                       QC1 = df$QC[1:12],
                       QC2 = df$QC[13:24],
                       QC3 = df$QC[25:36],
                       EN1 = df$EN[1:12],
                       EN2 = df$EN[13:24],
                       EN3 = df$EN[25:36])
  # Return the new dataframe
  return(new_df)
}