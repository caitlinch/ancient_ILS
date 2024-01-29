# ancient_ILS/code/func_data_analysis.R
## This script includes functions to analyse concordance factors, quartet scores, and phylogenetic trees
# Caitlin Cherryh, 2023

#### Extract details from IQ-Tree output files ####
extract.unconstrained.tree.details <- function(row_id, dataframe){
  # Extract best model, rates, and results of IQ-Tree run
  
  # Identify run and extract row
  temp_row <- dataframe[1,]
  # Identify IQ-Tree file from intial run
  initial_iqtree_file <- paste0(temp_row$gene_directory, temp_row$unconstrained_tree_prefix, ".iqtree")
  # Run the functions that extract details from the iqtree file
  best_model_op   <- extract.best.model(initial_iqtree_file)
  rates_op        <- extract.rates(initial_iqtree_file)
  tree_ll_op      <- extract.tree.log.likelihood(initial_iqtree_file, var = "All")
  gamma1_op       <- extract.gamma.values(initial_iqtree_file, gamma.parameter = "List")
  gamma2_op       <- extract.gamma.values(initial_iqtree_file, gamma.parameter = "Shape")
  statfreq_op     <- extract.state.frequencies(initial_iqtree_file)
  alisim_model_op <- extract.alisim.model(iqtree_file)
  # Assemble output into the nice temp row
  op_row <- temp_row
  op_row$unconstrained_tree_best_model            <- best_model_op
  op_row$unconstrained_tree_alisim_model          <- alisim_model_op
  op_row$unconstrained_tree_gamma_shape           <- gamma2_op
  op_row$unconstrained_tree_gamma_categories      <- gamma1_op
  op_row$unconstrained_tree_rates                 <- rates_op
  op_row$unconstrained_tree_state_freqs           <- statfreq_op
  op_row$unconstrained_tree_logl                  <- tree_ll_op[1]
  op_row$unconstrained_tree_unconstrained_logl    <- tree_ll_op[2]
  op_row$unconstrained_tree_numFreeParams         <- tree_ll_op[3]
  op_row$unconstrained_tree_BIC                   <- tree_ll_op[4]
  op_row$unconstrained_tree_length                <- tree_ll_op[5]
  op_row$unconstrained_tree_sumInternalBL         <- tree_ll_op[6]
  # Return the output
  return(op_row)
}



extract.best.model <- function(iqtree_file){
  # Function that will extract the best model of sequence evolution or the model of sequence evolution used,
  #   given a .iqtree file
  
  # Check if the file is a PMSF file
  # The PMSF models do not include a model search or a BIC score for the model
  pmsf_check <- grepl("PMSF", iqtree_file)
  if (pmsf_check == FALSE){
    # The model is not a PMSF model. Continue to extract best model
    # Check if the iqtree file exists
    if (file.exists(iqtree_file) == TRUE){
      # If the iqtree_file does exist:
      ## Open the .iqtree file:
      iq_lines <- readLines(iqtree_file)
      
      ## Check for a ModelFinder section:
      # Determine whether there is a ModelFinder section
      mf_ind <- grep("ModelFinder", iq_lines)
      # Determine whether there is a line detailing the best model
      bm_ind <- grep("Best-fit model according to", iq_lines)
      
      ## Check for a Substitution Process section:
      # Determine the starting line of this section
      sp_ind <- grep("SUBSTITUTION PROCESS", iq_lines)
      # Determine the line detailing the model used
      mos_ind <- grep("Model of substitution", iq_lines)
      
      ## Extract the best fit model from the .iqtree file:
      if ((identical(mf_ind, integer(0)) == FALSE) & (identical(bm_ind, integer(0)) == FALSE)){
        # If ModelFinder was run, extract the best model from the ModelFinder section of the .iqtree file
        # Extract the line containing the best fit model
        m_line <- iq_lines[bm_ind]
      } else if ((identical(sp_ind, integer(0)) == FALSE) & (identical(mos_ind, integer(0)) == FALSE)) {
        # If there is no ModelFinder section, extract the model used from the substitution process section
        m_line <- iq_lines[mos_ind]
      } else {
        m_line <- "NA:NA"
      }
      
      ## Format the model nicely for output: 
      # Split the line at the colon into two parts
      m_line_split <- strsplit(m_line, ":")[[1]]
      # If the best model is a single model, the length of m_line_split will be 2
      #     One section for the explanatory text and one for the model
      # If the best model is a partition model, it will have more than two sections when split by colons
      # Extract the second part of the line onwards (contains the best fit model)
      best_model <- m_line_split[2:length(m_line_split)]
      # If best_model is longer than 1, paste it together again using colons
      if (length(best_model) >1){
        best_model <- paste(best_model, collapse = ":")
      }
      # Remove any white space from the best model
      best_model <- gsub(" ", "", best_model)
    } else if (file.exists(iqtree_file) == FALSE){
      # If the iqtree_file doesn't exist, return NA
      best_model = NA
    } # end if (file.exists(iqtree_file) == TRUE){
  } else if (pmsf_check == TRUE){
    # If the PMSF model was used, there is no modelfinder output to search and no comparison BIC
    # Return NA
    best_model = NA
  } # end if (pmsf_check == FALSE){
  # Return the best model from the iqtree_file (if the file exists)
  return(best_model)
}




extract.rates <- function(iqtree_file){
  # Function to extract the rate parameters of a model from in IQ-Tree
  
  # Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    # Check for +R parameters 
    rate_ind <- grep("Site proportion and rates\\:", iq_lines)
    # Check if the rate_ind is present (if yes, that means there's a line containing the rates and weights)
    if (identical(rate_ind,integer(0)) == FALSE & class(rate_ind) == "integer"){
      ## Rates (+R)
      # The site proportion and weights values are present
      # Extract the rate weights and parameters from the iqtree file
      rates_line <- iq_lines[rate_ind]
      # Remove text from the beginning of the line
      rates_raw <- strsplit(rates_line, "\\:")[[1]][2]
      # Split up by the spaces
      split_rates <- strsplit(rates_raw, " ")[[1]]
      # Remove any entries that are empty characters (i.e. "")
      split_rates <- split_rates[split_rates != ""]
      # Split rates again by the commas (",")
      split2_rates <- unlist(strsplit(split_rates, ","))
      # Remove brackets from values
      rate_vals <- unlist(strsplit(unlist(strsplit(split2_rates, "\\(")), "\\)"))
      # Make output (for attaching into IQ-Tree2 command line)
      # Order is: proportion_1, rate_1, proportion_2, rate_2, ....., proportion_n, rate_n
      rates_op <- paste(rate_vals, collapse = ",")
    } else {
      # The site proportion and weights values are missing
      # Return NA for this file
      rates_op <- NA
    } # end if (identical(rate_ind,logical(0)) == FALSE & class(rate_ind) == "integer")
  } # end if (file.exists(iqtree_file) == TRUE)
  
  # Return the output
  return(rates_op)
}




extract.tree.log.likelihood <- function(iqtree_file, var = "LogL"){
  # Function to extract the log likelihood of a model from in IQ-Tree
  # Can extract either the log likelihood (var = "LogL), the unconstrained log likelihood (without tree) (var = "ULL"),
  #   the number of free parameters (var = "NFP), the BIC (var = "BIC"), the total tree length (var = "TrLen"), 
  #   the sum of internal branch lengths (var = "SIBL"),  or all of the above (var = "All")
  
  # Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    ## Extract the variables from the iqtree file
    if (var == "LogL"){
      # Log Likelihood
      output <- gsub(" ", "", strsplit(strsplit(iq_lines[grep("Log-likelihood of the tree\\:", iq_lines)], "\\:")[[1]][2], "\\(")[[1]][1])
    } else if (var == "ULL"){
      # Unconstrained log likelihood (without tree)
      output <- gsub(" ", "", strsplit(iq_lines[grep("Unconstrained log-likelihood \\(without tree\\)\\:", iq_lines)], "\\:")[[1]][2])
    } else if (var == "NFP"){
      # Number of free parameters
      output <- gsub(" ", "", strsplit(iq_lines[grep("Number of free parameters \\(\\#branches \\+ \\#model parameters\\)\\:", iq_lines)], "\\:")[[1]][2])
    } else if (var == "BIC"){
      # BIC
      output <- gsub(" ", "", strsplit(iq_lines[grep("Bayesian information criterion \\(BIC\\) score\\:", iq_lines)], "\\:")[[1]][2])
    } else if (var == "TrLen"){
      # Tree length
      output <- gsub(" ", "", strsplit(iq_lines[grep("Total tree length \\(sum of branch lengths\\)\\:", iq_lines)], "\\:")[[1]][2])
    } else if (var == "SIBL"){
      # Sum of internal branch lengths
      output <- gsub(" ", "", strsplit(strsplit(iq_lines[grep("Sum of internal branch lengths\\:", iq_lines)], "\\:")[[1]][2], "\\(")[[1]][1])
    } else if (var == "All"){
      # Extract all values
      logl <- gsub(" ", "", strsplit(strsplit(iq_lines[grep("Log-likelihood of the tree\\:", iq_lines)], "\\:")[[1]][2], "\\(")[[1]][1])
      ull <- gsub(" ", "", strsplit(iq_lines[grep("Unconstrained log-likelihood \\(without tree\\)\\:", iq_lines)], "\\:")[[1]][2])
      nfp <- gsub(" ", "", strsplit(iq_lines[grep("Number of free parameters \\(\\#branches \\+ \\#model parameters\\)\\:", iq_lines)], "\\:")[[1]][2])
      bic <- gsub(" ", "", strsplit(iq_lines[grep("Bayesian information criterion \\(BIC\\) score\\:", iq_lines)], "\\:")[[1]][2])
      treelen <- gsub(" ", "", strsplit(iq_lines[grep("Total tree length \\(sum of branch lengths\\)\\:", iq_lines)], "\\:")[[1]][2])
      sibl <- gsub(" ", "", strsplit(strsplit(iq_lines[grep("Sum of internal branch lengths\\:", iq_lines)], "\\:")[[1]][2], "\\(")[[1]][1])
      # Assemble all variables into a vector (for outputting all informating at once)
      output <- c("LogL" = logl, "Unconstrained_LogL" = ull, "NumFreeParams" = nfp, 
                  "BIC" = bic, "TreeLength" = treelen, "SumInternalBranchLengths" = sibl)
    } # end if (var == "LogL")
  } # end if (file.exists(iqtree_file) == TRUE)
  
  # Return the output
  return(output)
}




extract.gamma.values <- function(iqtree_file, gamma.parameter = "List"){
  # Function to extract the gamma parameters of a model from in IQ-Tree
  #   gamma.parameter = "List" - if model has gamma rate, returns the site relative weights and proportions
  #   gamma.parameter = "Shape" - if model has gamma rate, returns the gamma shape alpha 
  
  # Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    # Check for +G parameters 
    gamma_shape_ind <- grep("Gamma shape alpha\\:", iq_lines)
    gamma_list_ind <- grep("Model of rate heterogeneity\\: Gamma", iq_lines)
    # Check if the gamma parameters are present (if yes, that means there's a line containing the rates and weights)
    if ( (identical(gamma_shape_ind,integer(0)) == FALSE & class(gamma_shape_ind) == "integer") |
         (identical(gamma_list_ind,integer(0)) == FALSE & class(gamma_list_ind) == "integer") ){
      ## Gamma (+G)
      if (gamma.parameter == "Shape"){
        ## Extract the value of the gamma shape alpha
        # Extract the line containing the gamma alpha value
        gamma_line <- iq_lines[gamma_shape_ind]
        # Process the line to extract just the gamma alpha value
        raw_gamma_alpha <- strsplit(gamma_line, "\\:")[[1]][2]
        gamma_alpha <- as.numeric(gsub(" ", "", raw_gamma_alpha))
        # Set the output 
        gamma_op <- gamma_alpha
      } else if (gamma.parameter == "List"){
        ## Extract the gamma weights and proportions
        #Create the matrix for discrete gamma categories
        g_start <- grep(" Category", iq_lines) + 1 # get the index for the first line of the gamma categories matrix
        empty   <- which(iq_lines == "") # get indexes of all empty lines
        empty   <- empty[empty > g_start] # get empty lines above gamma categories matrix
        g_end   <- empty[1] - 1 # get end index for gamma categories matrix (one less than next empty line)
        end_line <- iq_lines[g_end]
        # if the end isn't an empty line, subtract one from the end count 
        # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
        # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
        check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1], " ")[[1]])
        if (check_line > 3){
          # If the check_line is longer than 3 characters, it won't be a group for the gamma categories but an instruction
          # Instructions can be excluded from the gamma matrix (but categories can't)
          g_end = g_end - 1
        }
        # Extract the lines of interest
        g_lines <- iq_lines[g_start:g_end]
        # Split the g_lines at the spaces
        g_lines_split <- strsplit(g_lines, " ")
        # Remove empty values
        g_lines_neat <- lapply(1:length(g_lines_split), function(i){g_lines_split[[i]] <- g_lines_split[[i]][which(g_lines_split[[i]] != "")]})
        # Collect the values in the right order for the output
        raw_gamma_vals <- unlist(lapply(1:length(g_lines_neat), function(i){c(g_lines_neat[[i]][3], g_lines_neat[[i]][2])}))
        # Create a nice output format
        gamma_vals <- paste(gsub(" ", "", raw_gamma_vals), collapse = ",")
        # Set the output
        gamma_op <- gamma_vals
      } # end if (gamma.parameter == "Shape"){
    }else {
      # The site proportion and weights values are missing
      # Return NA for this file
      gamma_op <- NA
    } # end if ( (identical(gamma_shape_ind,integer(0)) == FALSE & class(gamma_shape_ind) == "integer") | ( ...) ){
  } # end if (file.exists(iqtree_file) == TRUE)
  
  # Return the output
  return(gamma_op)
}




extract.state.frequencies <- function(iqtree_file){
  # Given an iqtree file, this function will extract the state frequencies for the alignment
  
  ## Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    ## Extract state frequency details:
    # Check for presence of state frequencies details
    ind <- grep("State frequencies:", iq_lines)
    if (identical(ind, integer(0)) == FALSE){
      ## Determine whether state frequencies were empirically determined (or not!)
      # If there is a section for state frequencies, extract and output details 
      sf1 <- strsplit(iq_lines[[ind]], ":")[[1]][2]
      # Check whether state frequencies are needed
      sf1_squashed <- gsub(" ", "", sf1)
      if (sf1_squashed == "(empiricalcountsfromalignment)"){
        ## If state frequencies were determined from the alignment, 
        #    extract the lines with state frequencies from the .iqtree file
        # Get starting line for frequencies
        start_ind <- grep("State frequencies:", iq_lines) + 2
        # Take the 20 lines containing AA frequencies
        freq_lines <- iq_lines[start_ind:(start_ind+19)]
        # Split up the frequency lines into the label and the frequency
        freq_split <- unlist(strsplit(freq_lines, "="))
        
        ## Process the frequencies 
        # Get the frequency
        freq_nums <- freq_split[c(FALSE, TRUE)]
        # Remove any spaces (from IQTree formatting)
        freq_nums <- gsub(" ","",freq_nums)
        
        ## Process the labels
        # Get corresponding AA letter
        freq_names <- freq_split[c(TRUE, FALSE)]
        # Remove IQTree formatting
        freq_names <- gsub("pi\\(", "", freq_names)
        freq_names <- gsub("\\)", "", freq_names)
        freq_names <- gsub(" ", "", freq_names)
        
        ## Make the output pretty
        # Create a nice output by pasting together the frequencies in order
        f_op = paste(freq_nums, collapse = ",")
        
      } else if (sf1_squashed == "(equalfrequencies)"){
        ## If state frequencies were equal, generate equal state frequencies
        f_op <- paste(as.character(rep(1/20, 20)), collapse = ",")
      } else {
        f_op <- "State frequencies from model"
      }
    } else if (identical(ind, integer(0)) == TRUE){
      # If no details on state frequencies are needed, return an empty dataframe
      f_op <- NA
    }
  }
  
  ## Output the frequencies
  return(f_op)
}




extract.alisim.model <- function(iqtree_file){
  # Given an iqtree file, this function will extract the alisim model (if present)
  
  ## Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    ## Check for the phrase "ALISIM COMMAND"
    line_check <- grep("ALISIM COMMAND", iq_lines)
    if (identical(line_check, integer(0)) == FALSE){
      ## Extract the alisim model
      # Extract the line with the alisim command
      alisim_check <- grep("--alisim", iq_lines)
      # Extract the first line that contains the alisim_check - this will include the best model
      alisim_command_line <- iq_lines[ alisim_check[1] ]
      # Remove double quotes and back slashes from alisim command line
      alisim_command_line_formatted  <- gsub("\"", "'", alisim_command_line)
      # Split the string up to extract the model
      split1 <- strsplit(alisim_command_line_formatted, "-m")
      split2 <- strsplit(split1[[1]][2], "--length")
      model_chunk <- split2[[1]][1]
      # Format the model chunk
      alisim_model <- gsub("'| ", "", model_chunk)
    } else if (identical(line_check, integer(0)) == TRUE){
      ## Return NA
      # If no alisim command, return empty output
      alisim_model <- NA
    }
  } else {
    ## Return NA
    # If no IQ-Tree file, return empty output
    alisim_model <- NA
  }
  
  ## Return the alisim model
  return(alisim_model)
}

