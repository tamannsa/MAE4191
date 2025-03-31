#######################################
#         MASTER THESIS               #
#######################################


#######################################
#             PACKAGES                #
#######################################

library(haven)   # Import SPSS, Stata, and SAS files
library(dplyr)   # Data manipulation 
library(naniar)  # Handling missing data visualization
library(survey)  # Complex survey analysis 
library(VIM)     # Visualization of missing values
library(mice)    # Multiple Imputation by Chained Equations
library(gridExtra) # Arrange multiple ggplot2 plots
library(ggplot2) # Data visualization
library(mitools) # Tools for analyzing multiply imputed datasets
library(tidyverse)
library(corrplot)
library(patchwork) # For combining plots
library(car)      # For VIF calculations

# =============================================================================
# Citations: Retrieve package citations & R version 
# =============================================================================
citation("survey") 
citation("mice")
citation("mitools")

version$version.string # Displays the current R version
# R version 4.4.2

#######################################
#             LOAD DATA               #
#######################################

# Source: PISA 2022 Student Questionnaire Data (OECD)
# Link: https://www.oecd.org/en/data/datasets/pisa-2022-database.html#data
# Data dimensions: 613,744 observations x 1,278 variables

setwd("~/UiO") # Set working directory 

# Load PISA 2022 Student Questionnaire data from SPSS file
CY08MSP_STU_QQQ <- read_sav("MAE4191/PISA/2022data/CY08MSP_STU_QQQ.SAV", user_na=TRUE)

# Store original dataset in a master copy for further processing
Master_data <- CY08MSP_STU_QQQ 


#######################################
#         VARIABLE SELECTION          #
#######################################

# =============================================================================
# Data preparation
# =============================================================================

# Function to prepare PISA 2022 student data
prepare_pisa_data <- function(data, countries, predictors, convert_to_na = TRUE) {
  # Define PISA missing value codes and their meanings
  invalid_codes <- list(
    missing = c(99, 999, 9999, 99999, 999999, 9999999),         # Missing
    invalid = c(97, 997, 9997, 99997, 999997, 9999997),         # N/A (Not Applicable)
    not_reached = c(98, 998, 9998, 99998, 999998, 9999998),     # Not Reached
    not_admin = c(95, 995, 9995, 99995, 999995, 9999995),       # Quest Not Admin
    logical_skip = c(96, 996, 9996, 99996, 999996, 9999996)     # Logical Skip
  )
  
  # Identify columns of interest
  pv_cols <- paste0("PV", 1:10, "MATH") # Mathematics plausible values
  weight_cols <- grep("^W_FSTR", names(data), value = TRUE) # Replicate weight columns
  
  # Subset data: Retain only selected countries and relevant variables
  data_subset <- data %>%
    filter(CNT %in% countries) %>% # Keep only specified countries
    select(CNT, W_FSTUWT, all_of(weight_cols), 
           all_of(predictors), all_of(pv_cols)) # Retain key variables
  
  # Convert labelled variables (from SPSS) to standard numeric values
  data_subset <- data_subset %>%
    mutate(across(everything(), haven::zap_labels)) # Remove SPSS-style labels
  
  # Function to handle different types of missing values
  recode_missing <- function(x) {
    if(convert_to_na) {
      case_when(
        x %in% invalid_codes$missing ~ NA_real_,      
        x %in% invalid_codes$invalid ~ NA_real_,      
        x %in% invalid_codes$not_reached ~ NA_real_,  
        x %in% invalid_codes$not_admin ~ NA_real_,    
        x %in% invalid_codes$logical_skip ~ NA_real_, 
        TRUE ~ as.numeric(x) # Convert remaining values to numeric
      )
    } else {
      as.numeric(x)  # Convert without altering missing values
    }
  }
  
  # Apply missing value recoding to predictors and plausible values
  data_subset <- data_subset %>%
    mutate(across(c(all_of(predictors), all_of(pv_cols)), recode_missing))
  
  return(data_subset) # Return cleaned and subsetted dataset
}

# Subset and clean the PISA dataset for selected countries and variables
data_subset <- prepare_pisa_data(
  data = Master_data,
  countries = c("NOR", "CHE", "FIN", "SWE", "DNK"), # Selected countries
  predictors = c("IMMIG", "RELATST", "BELONG", "BULLIED", "FEELSAFE", "SCHRISK", 
                 "PERSEVAGR", "CURIOAGR", "COOPAGR", "GROSAGR", "DISCLIM", "TEACHSUP", 
                 "COGACRCO", "COGACMCO", "MATHEFF", "MATHEF21", "MATHPERS", "ANXMAT", 
                 "CREATEFF", "CREATSCH") # Selected student background and attitude variables
)

# =============================================================================
# Weighted predictor analysis
# =============================================================================

# Function to analyze weighted predictor characteristics
analyze_weighted_predictors <- function(data, predictors) {
  results <- list() # Initialize result storage
  
  # Identify replicate weight columns used for variance estimation
  rep_weights <- grep("^W_FSTR", names(data), value = TRUE)
  
  ## Step 1: MISSING DATA ANALYSIS (WEIGHTED)
  missing_analysis <- lapply(predictors, function(pred) {
    # Define survey design for PISA's complex weighting
    design <- svydesign(
      id = ~1, # No clustering assumed at this stage
      weights = ~W_FSTUWT, # Final student weight
      data = data,
      repweights = data[rep_weights], # Replicate weights
      type = "Fay", # Fay's method for variance estimation
      rho = 0.5 # Common choice for Fay's adjustment
    )
    
    # Compute the weighted proportion of missing values for each predictor
    missing_prop <- svymean(~is.na(get(pred)), design, na.rm = TRUE)
    
    # Store results
    data.frame(
      variable = pred,
      missing_percent = as.numeric(missing_prop) * 100, # Convert to percentage
      se = sqrt(vcov(missing_prop)) * 100 # Standard error
    )
  }) %>% bind_rows()
  
  # Store missing data analysis results
  results$missing <- missing_analysis %>%
    arrange(desc(missing_percent)) # Order variables by missingness
  
  ## Step 2: WEIGHTED CORRELATIONS WITH MATH ACHIEVEMENT
  pv_correlations <- lapply(predictors, function(pred) {
    # Compute weighted correlation between each predictor and math achievement
    pv_cors <- sapply(1:10, function(i) {
      pv_col <- paste0("PV", i, "MATH") # Math plausible value column
      # Compute weighted correlation using a custom function
      cor_est <- weighted_cor(
        x = data[[pred]],
        y = data[[pv_col]],
        weights = data$W_FSTUWT
      )
      return(cor_est)
    })
    
    # Store results: average correlation across plausible values
    data.frame(
      predictor = pred,
      avg_correlation = mean(pv_cors, na.rm = TRUE),
      min_correlation = min(pv_cors, na.rm = TRUE),
      max_correlation = max(pv_cors, na.rm = TRUE)
    )
  }) %>% bind_rows()
  
  # Store and sort correlation results
  results$pv_correlations <- pv_correlations %>%
    arrange(desc(abs(avg_correlation))) # Order by strength of relationship
  
  ## Step 3: WEIGHTED CORRELATIONS BETWEEN PREDICTORS
  # Create an empty matrix to store correlations
  pred_cors <- matrix(NA, length(predictors), length(predictors))
  rownames(pred_cors) <- colnames(pred_cors) <- predictors # Label matrix
  
  # Compute pairwise weighted correlations
  for(i in 1:length(predictors)) {
    for(j in i:length(predictors)) {
      if(i != j) {
        cor_est <- weighted_cor(
          x = data[[predictors[i]]],
          y = data[[predictors[j]]],
          weights = data$W_FSTUWT
        )
        pred_cors[i,j] <- pred_cors[j,i] <- cor_est # Fill symmetric matrix
      } else {
        pred_cors[i,i] <- 1 # Perfect correlation with itself
      }
    }
  }
  
  # Convert matrix to dataframe format for easier handling
  results$predictor_correlations <- as.data.frame(as.table(pred_cors)) %>%
    rename(Parameter1 = Var1, Parameter2 = Var2, r = Freq) %>%
    filter(Parameter1 != Parameter2) # Remove self-correlations
  
  return(results) # Return all analysis results
}

# Function to calculate weighted correlation between two variables
weighted_cor <- function(x, y, weights) {
  # Compute weighted deviations from mean
  wx <- weights * (x - weighted.mean(x, weights, na.rm = TRUE))
  wy <- weights * (y - weighted.mean(y, weights, na.rm = TRUE))
  # Compute weighted Pearson correlation
  sum(wx * wy, na.rm = TRUE) / sqrt(sum(wx * wx, na.rm = TRUE) * sum(wy * wy, na.rm = TRUE))
}

# Run weighted analyses
results <- analyze_weighted_predictors(
  data = data_subset,
  predictors = c("IMMIG", "RELATST", "BELONG", "BULLIED", "FEELSAFE", "SCHRISK", 
                 "PERSEVAGR", "CURIOAGR", "COOPAGR", "GROSAGR", "DISCLIM", "TEACHSUP", 
                 "COGACRCO", "COGACMCO", "MATHEFF", "MATHEF21", "MATHPERS", "ANXMAT", 
                 "CREATEFF", "CREATSCH")
)

# Function to filter through variables based on already decided criteria
suggest_predictors <- function(analysis_results, 
                               missing_threshold = 30,
                               correlation_threshold = 0.7,
                               min_pv_correlation = 0.1) {
  
  ## Step 1: Remove predictors with too many missing values
  keep_predictors <- analysis_results$missing %>%
    filter(missing_percent < missing_threshold) %>%
    pull(variable)
  
  ## Step 2: Identify highly correlated predictor pairs
  high_cors <- analysis_results$predictor_correlations %>%
    filter(abs(r) > correlation_threshold & Parameter1 != Parameter2)
  
  ## Step 3: Remove one predictor from each highly correlated pair
  if(nrow(high_cors) > 0) {
    for(i in 1:nrow(high_cors)) {
      pred1 <- high_cors$Parameter1[i]
      pred2 <- high_cors$Parameter2[i]
      # Step 3a: Retrieve average correlations with the outcome
      cor1 <- analysis_results$pv_correlations %>%
        filter(predictor == pred1) %>%
        pull(avg_correlation)
      cor2 <- analysis_results$pv_correlations %>%
        filter(predictor == pred2) %>%
        pull(avg_correlation)
      # Step 3b: Keep the predictor with the stronger relationship to the outcome
      if(abs(cor1) < abs(cor2)) {
        keep_predictors <- keep_predictors[keep_predictors != pred1]
      } else {
        keep_predictors <- keep_predictors[keep_predictors != pred2]
      }
    }
  }
  
  # Step 4: Retain only predictors with a meaningful correlation to the outcome
  final_predictors <- analysis_results$pv_correlations %>%
    filter(predictor %in% keep_predictors,
           abs(avg_correlation) >= min_pv_correlation) %>%
    pull(predictor)
  return(final_predictors)
}

# Get suggested predictors based on missingness, multicollinearity, and correlation thresholds
final_predictors <- suggest_predictors(
  results,
  missing_threshold = 30, # Remove predictors with more than 30% missing values
  correlation_threshold = 0.7, # Remove one predictor from pairs with correlation > 0.7
  min_pv_correlation = 0.1 # Keep predictors with at least 0.1 correlation with the outcome
)
# Display the final list of selected predictors
final_predictors

# Identify variables that are completely missing for each country
missing_items_by_country <- data_subset %>%
  group_by(CNT) %>%
  summarise(across(everything(), ~ all(is.na(.)), .names = "missing_{.col}"))
# View the results to assess patterns of complete missingness
print(missing_items_by_country)


#######################################
#       MISSINGNESS DIAGNOSTICS       #
#######################################

# Subset dataset for relevant countries
countries_of_interest <- c("NOR", "SWE", "DNK", "FIN", "CHE")
Master_data <- CY08MSP_STU_QQQ %>%
  filter(CNT %in% countries_of_interest)

# Remove CY08MSP_STU_QQQ to open up memory
rm(CY08MSP_STU_QQQ)

# Summary of missing values for key variables
missing_summary <- colSums(
  is.na(Master_data[, c(
      "RELATST", "GROSAGR", 
      "DISCLIM", "MATHEFF", 
      "MATHEF21", "MATHPERS", 
      "ANXMAT")
      ]
    )
  )
print(missing_summary)

# Function to check missingness patterns, including special missing codes
check_specific_codes <- function(variable) {
  # Ensure proper handling of labelled data
  if (inherits(variable, "haven_labelled")) {
    variable_numeric <- as.numeric(variable) # Convert labels to numeric while preserving NA
  } else {
    variable_numeric <- variable
  }
  # Define special missing codes of interest
  special_codes <- c(95, 97, 98, 99)
  # Count occurrences of each special missing code
  code_results <- sapply(special_codes, function(code) {
    count <- sum(variable_numeric == code, na.rm = TRUE) # Count occurrences
    pct <- (count / length(variable_numeric)) * 100 # Calculate percentage
    return(c(count = count, percentage = pct))
  })
  # Display results for special codes
  cat("Special Code Analysis:\n")
  for (i in seq_along(special_codes)) {
    cat(sprintf("Code %d: Count = %d (%.2f%%)\n", 
                special_codes[i], 
                code_results[1, i], 
                code_results[2, i]))
  }
  # Count and display regular NA occurrences
  na_count <- sum(is.na(variable_numeric))
  na_pct <- (na_count / length(variable_numeric)) * 100
  cat(sprintf("Regular NA: Count = %d (%.2f%%)\n", na_count, na_pct))
}

# =============================================================================
# Investigating missingness in scales
# =============================================================================

# Apply missingness check function to each variable and print results
variables_to_check <- c("RELATST", "GROSAGR", "DISCLIM", 
                        "MATHEFF", "MATHEF21", "MATHPERS", "ANXMAT")

# Print results of missingness by type 
for (var in variables_to_check) {
  cat(sprintf("\n%s missingness pattern:\n", var))
  check_specific_codes(Master_data[[var]])
}

# =============================================================================
# Investigating missingness in part-variables
# =============================================================================


# Create item values
ST267_vars <- paste0("ST267Q", sprintf("%02d", 1:08), "JA") # RELATST
ST263_vars <- paste0("ST263Q", sprintf("%02d", seq(2, 8, by = 2)), "JA") # GROSAGR
ST273_vars <- paste0("ST273Q", sprintf("%02d", 1:07), "JA") # DISCLIM
ST290_vars <- paste0("ST290Q", sprintf("%02d", 1:09), "WA") # MATHEFF
ST291_vars <- paste0("ST291Q", sprintf("%02d", 1:10), "JA") # MATHEFF21
ST293_vars <- paste0("ST293Q", sprintf("%02d", 1:09), "JA") # MATHPERS
ST292_vars <- paste0("ST292Q", sprintf("%02d", 1:06), "JA") # ANXMAT

# Loop through and analyze all items in RELATST
for (var in ST267_vars) {
  cat("\n--- Analysis for", var, "---\n")
  check_specific_codes(Master_data[[var]])
}
# Loop through and analyze all items in GROSAGR
for (var in ST263_vars) {
  cat("\n--- Analysis for", var, "---\n")
  check_specific_codes(Master_data[[var]])
}
# Loop through and analyze all items in DISCLIM
for (var in ST273_vars) {
  cat("\n--- Analysis for", var, "---\n")
  check_specific_codes(Master_data[[var]])
}
# Loop through and analyze all items in MATHEFF
for (var in ST290_vars) {
  cat("\n--- Analysis for", var, "---\n")
  check_specific_codes(Master_data[[var]])
}
# Loop through and analyze all items in MATHEFF21
for (var in ST291_vars) {
  cat("\n--- Analysis for", var, "---\n")
  check_specific_codes(Master_data[[var]])
}
# Loop through and analyze all items in MATHPERS
for (var in ST293_vars) {
  cat("\n--- Analysis for", var, "---\n")
  check_specific_codes(Master_data[[var]])
}
# Loop through and analyze all items in ANXMAT
for (var in ST292_vars) {
  cat("\n--- Analysis for", var, "---\n")
  check_specific_codes(Master_data[[var]])
}

# =============================================================================
# Convert to NA
# =============================================================================

# Replace specific values (95, 97, 98, 99) with NA in the entire dataset
Master_data <- Master_data %>%
  mutate(across(everything(), ~ replace(., . %in% c(95, 97, 98, 99, 
                                                    995, 997, 998, 999, 
                                                    9995, 9997, 9998, 9999, 
                                                    99995, 99997, 99998, 99999, 
                                                    999995, 999997, 999998, 999999, 
                                                    9999995, 9999997, 9999998, 9999999), NA)))

# =============================================================================
# Complete cases
# =============================================================================

# Subset for variables of interest
df_select <- Master_data[,c("RELATST", "GROSAGR", "DISCLIM", "MATHEFF", "MATHEF21", "MATHPERS", "ANXMAT")]
# Count the number of rows with complete data
num_complete_cases <- sum(complete.cases(df_select))
# Calculate the percentage of complete cases
total_students <- nrow(df_select)
percentage_complete <- (num_complete_cases / total_students) * 100
# Print the result
cat(sprintf("Number of students with complete data: %d (%.2f%%)\n", 
            num_complete_cases, percentage_complete))

# =============================================================================
# Missing Data Visualization and Tests
# =============================================================================

# Aggregate missing data pattern
aggr(Master_data[, c("RELATST", "GROSAGR", "DISCLIM", "MATHEFF", "MATHEF21", "MATHPERS", "ANXMAT")], 
     numbers = TRUE, prop = FALSE)

# Little`s MCAR test on selected variables
mcar_test(data=df_select)

# =============================================================================
# Running the Survey-Weighted Analysis of Missing Data
# =============================================================================

# Cleaning the data by removing haven labels
Master_data <- as.data.frame(lapply(Master_data, haven::zap_labels))

# Missing data analysis
analyze_missing_data_survey <- function(data, vars, weights, strata_var, psu_var) {
  ## Step 1: Create working dataset by selecting relevant variables
  df_select <- data[, c(vars, weights, strata_var, psu_var)]
  ## Step 2: Count the number of PSUs (Primary Sampling Units) per stratum
  psu_counts <- table(df_select[[strata_var]])
  small_strata <- names(psu_counts[psu_counts < 3]) # Identify strata with fewer than 3 PSUs
  ## Step 3: Combine small strata into one category for better analysis
  df_select$stratum_modified <- df_select[[strata_var]]
  df_select$stratum_modified[df_select[[strata_var]] %in% small_strata] <- "combined_small_strata"
  ## Step 4: Create missingness indicators for each variable
  for(var in vars) {
    df_select[[paste0(var, "_missing")]] <- as.numeric(is.na(df_select[[var]])) # 1 if missing, 0 if not
  }
  ## Step 5: Define survey design with modified strata
  options(survey.lonely.psu = "certainty") # Change handling of single PSUs to "certainty"
  survey_design <- svydesign(
    id = as.formula(paste0("~", psu_var)), # Primary sampling unit (PSU)
    strata = ~stratum_modified, # Stratification variable (modified for small strata)
    weights = as.formula(paste0("~", weights)), # Survey weights (Final trimmed nonresponse adjusted student weight)
    data = df_select # The data to use
  )
  ## Step 6: Loop through all variables to perform analysis
  for(var in vars) {
    explanatory <- vars[vars != var] # Explanatory variables (all but the current variable)
    results <- list() # Store results for each explanatory variable
    # Step 6a: Loop through explanatory variables
    for(exp_var in explanatory) {
      # Define the formula for the analysis (mean calculation)
      form <- as.formula(paste0("~", exp_var))
      # Try-catch for safe calculation of means and tests
      tryCatch({
        # Calculate the mean for non-missing values
        resp_mean <- svymean(form, 
                             subset(survey_design, get(paste0(var, "_missing")) == 0), # Only non-missing values for 'var'
                             na.rm = TRUE)
        # Calculate the mean for missing values
        nonresp_mean <- svymean(form, 
                                subset(survey_design, get(paste0(var, "_missing")) == 1), # Only missing values for 'var'
                                na.rm = TRUE)
        # Perform Wilcoxon rank-sum test (Robust variance estimation)
        test <- svyranktest(
          as.formula(paste0(exp_var, " ~ factor(", var, "_missing)")),
          design = survey_design,
          test.statistic = "Wilcoxon",
          variance = "bootstrapped"
        )
        # Store results in a data frame
        results[[exp_var]] <- data.frame(
          Variable = exp_var,
          Not_missing = paste0(round(coef(resp_mean), 2), " (", 
                               round(sqrt(vcov(resp_mean)), 2), ")"),
          Missing = paste0(round(coef(nonresp_mean), 2), " (", 
                           round(sqrt(vcov(nonresp_mean)), 2), ")"),
          p = round(test$p.value, 4) # P-value for the statistical test
        )
      }, error = function(e) {
        # Handle any errors in calculations
        results[[exp_var]] <- data.frame(
          Variable = exp_var,
          Not_missing = "Error in calculation",
          Missing = "Error in calculation",
          p = NA
        )
      })
    }
    # Step 7: Combine results for all explanatory variables into one table
    final_table <- do.call(rbind, results)
    print(knitr::kable(final_table, 
                       caption = paste("Survey-weighted missing data analysis for", var),
                       row.names = FALSE))
  }
}

# Survey-weighted analysis of missing data 
analyze_missing_data_survey(
  data = Master_data,
  vars = c("RELATST", "GROSAGR", "DISCLIM", "MATHEFF", "MATHEF21", "MATHPERS", "ANXMAT"),
  weights = "W_FSTUWT",      # PISA final student weight
  strata_var = "STRATUM",    # Stratification variable
  psu_var = "CNTSCHID"       # Primary sampling unit (school)
)

#######################################
#            Save dataset             #
#######################################

# Save dataset as SPSS .sav file
write_sav(Master_data, "~/UiO/MAE4191/Coding/Master_data.sav")

#######################################
#         MULTIPLE IMPUTATION         #
#######################################

# Multiple Imputation with Plausible Values
# =============================================================================
# This script performs the following steps:
# 1. Identifies auxiliary variables based on correlation with target variables
# 2. Processes PISA 2022 data by creating separate datasets for each plausible value
# 3. Creates test subsets by country (NOR, SWE, DNK, FIN, CHE)
# 4. Performs multiple imputation on test and full datasets
# =============================================================================

# Load prepared dataset
Master_data <- read_sav("Master_data.sav")

# Define countries of interest
countries <- c("NOR", "SWE", "DNK", "FIN", "CHE")

# Define target variables
target_vars <- c("ANXMAT", "MATHPERS", "MATHEFF", 
                 "DISCLIM", "GROSAGR", "RELATST",
                 "ESCS")

# Define survey variables (IDs, weights, etc.)
survey_vars <- c("CNT", "CNTSCHID", "CNTSTUID", 
                 "W_FSTUWT", paste0("W_FSTURWT", 1:80))

# =============================================================================
# Functions
# =============================================================================

# Function 1: Identify auxiliary variables
# =============================================================================
identify_auxiliary_vars <- function(data, target_vars, threshold_cor = 0.5, 
                                    threshold_miss = 0.05, additional_vars = NULL) {
  result <- list()
  # Calculate missing data proportions for all variables
  missing_prop <- colMeans(is.na(data))
  # Identify all PV variables to exclude them
  pv_vars <- grep("^PV[0-9]+", names(data), value = TRUE)
  # Only look at variables with acceptable missing data levels
  potential_auxiliaries <- names(missing_prop)[missing_prop <= threshold_miss & 
                                                 !(names(missing_prop) %in% target_vars) &
                                                 !(names(missing_prop) %in% pv_vars)]  # Add this line to exclude PVs
  # Include manual auxiliary variables regardless of missingness
  if (!is.null(additional_vars)) {
    additional_vars <- additional_vars[additional_vars %in% names(data)]
    potential_auxiliaries <- unique(c(potential_auxiliaries, additional_vars))
  }
  # Calculate correlations with target variables
  correlations <- data.frame(variable = character(),
                             target = character(),
                             correlation = numeric())
  for (var in potential_auxiliaries) {
    for (target in target_vars) {
      if (var != target) {
        # Use only complete cases for correlation
        complete_data <- data[complete.cases(data[, c(var, target)]), c(var, target)]
        # Only calculate if we have enough data and variable has variation
        if (nrow(complete_data) > 30 && 
            length(unique(complete_data[[var]])) > 1 && 
            length(unique(complete_data[[target]])) > 1) {
          tryCatch({
            cor_value <- cor(complete_data[[var]], complete_data[[target]])
            if (!is.na(cor_value) && abs(cor_value) >= threshold_cor) {
              correlations <- rbind(correlations, 
                                    data.frame(variable = var, 
                                               target = target, 
                                               correlation = cor_value))
            }
          }, error = function(e) {
            # Skip if correlation cannot be calculated
          })
        }
      }
    }
  }
  # Sort by absolute correlation
  correlations$abs_cor <- abs(correlations$correlation)
  correlations <- correlations[order(-correlations$abs_cor), ]
  # Identify recommended auxiliary variables
  selected_auxiliaries <- unique(correlations$variable)
  # Always include manually specified auxiliaries
  if (!is.null(additional_vars)) {
    selected_auxiliaries <- unique(c(selected_auxiliaries, additional_vars))
  }
  # Return results
  result$recommended_auxiliaries <- selected_auxiliaries
  result$correlation_details <- correlations
  result$missing_proportions <- missing_prop[selected_auxiliaries]
  return(result)
}
                                     
# Function 2: Select auxiliary variables with user review
# =============================================================================
select_auxiliary_variables <- function(data, target_vars, manual_auxiliary_vars = NULL, 
                                       threshold_cor = 0.5, threshold_miss = 0.05) {
  # Run function with additional auxiliary variables
  aux_results <- identify_auxiliary_vars(
    data, 
    target_vars, 
    threshold_cor = threshold_cor,
    threshold_miss = threshold_miss,
    additional_vars = manual_auxiliary_vars
  )
  # Print the final recommended auxiliary variables
  cat("Recommended auxiliary variables (", length(aux_results$recommended_auxiliaries), "):\n")
  print(aux_results$recommended_auxiliaries)
  # Print top correlations for reference
  if (nrow(aux_results$correlation_details) > 0) {
    cat("\nTop 20 correlations with target variables:\n")
    print(head(aux_results$correlation_details, 20))
  }
  return(aux_results)
}

# Function 3: Create datasets with one plausible value each
# =============================================================================
create_pv_datasets <- function(data, target_vars, survey_vars, auxiliary_vars, math_pv_vars) {
  # Check if I found any math PV variables
  if (length(math_pv_vars) == 0) {
    stop("No math PV variables found in the dataset")
  }
  # Determine how many sets of PVs we have (normally 10 for PISA 2022)
  pv_numbers <- unique(gsub(".*PV([0-9]+).*", "\\1", math_pv_vars))
  pv_numbers <- as.numeric(pv_numbers)
  # Create a unique student ID if not present
  if (!"STUDENT_ID" %in% names(data)) {
    data$STUDENT_ID <- 1:nrow(data)
  }
  # Variables to include in each dataset (excluding PVs)
  base_vars <- unique(c("STUDENT_ID", target_vars, survey_vars, auxiliary_vars))
  base_vars <- base_vars[base_vars %in% names(data)]
  # Create a list to store our datasets
  pv_datasets <- list()
  # For each PV number, create a dataset
  for (pv_num in pv_numbers) {
    # Find the variable for this PV number
    current_pv_var <- grep(paste0("PV", pv_num, ".*MATH"), math_pv_vars, value = TRUE)
    if (length(current_pv_var) > 0) {
      # If multiple matches, take the first one
      if (length(current_pv_var) > 1) {
        warning(paste("Multiple matches for PV", pv_num, ". Using:", current_pv_var[1]))
        current_pv_var <- current_pv_var[1]
      }
      # Create dataset with base variables + this PV
      vars_to_include <- c(base_vars, current_pv_var)
      temp_data <- data[, vars_to_include]
      # Rename the PV variable to a standard name
      names(temp_data)[names(temp_data) == current_pv_var] <- "MATH_PV"
      # Add to our list
      pv_datasets[[pv_num]] <- temp_data
      cat("Created dataset for PV", pv_num, "with dimensions:", 
          nrow(temp_data), "x", ncol(temp_data), "\n")
    }
  }
  return(pv_datasets)
}

# Function 4: Split datasets by country
# =============================================================================
split_by_country <- function(pv_datasets, countries) {
  country_pv_datasets <- list() # Initialize an empty list to store datasets split by country
  for (country in countries) { # Loop through each country of interest
    country_pv_datasets[[country]] <- list()
    for (pv_num in seq_along(pv_datasets)) { # Loop through each plausible value dataset
      # Extract data for this country
      country_data <- pv_datasets[[pv_num]] %>% 
        filter(CNT == country) # Filter rows where CNT matches the current country
      if (nrow(country_data) > 0) { # Check if the filtered dataset has any rows
        country_pv_datasets[[country]][[pv_num]] <- country_data # Store the filtered dataset in the nested list
        cat("Created dataset for country", country, "PV", pv_num, 
            "with dimensions:", nrow(country_data), "x", ncol(country_data), "\n")
      } else {
        warning(paste("No data found for country", country, "in PV", pv_num))
      }
    }
  }
  # Return the nested list of datasets split by country and PV
  return(country_pv_datasets) 
}

# Function 5: Create test subsets
# =============================================================================
create_test_subset <- function(country_pv_datasets, n_per_country = 100) {
  test_datasets <- list()
  for (country in names(country_pv_datasets)) {
    test_datasets[[country]] <- list()
    for (pv_num in seq_along(country_pv_datasets[[country]])) {
      # Get the current dataset
      current_data <- country_pv_datasets[[country]][[pv_num]]
      # Create test subset
      if (nrow(current_data) >= n_per_country) {
        # Set seed for reproducibility but vary by country and PV
        set.seed(as.numeric(factor(country)) * 100 + pv_num)
        test_subset <- current_data %>%
          slice_sample(n = n_per_country)
      } else {
        test_subset <- current_data
        warning(paste("Country", country, "PV", pv_num, 
                      "has fewer than", n_per_country, "observations"))
      }
      # Add to list
      test_datasets[[country]][[pv_num]] <- test_subset
      cat("Created test subset for country", country, "PV", pv_num, 
          "with dimensions:", nrow(test_subset), "x", ncol(test_subset), "\n")
    }
  }
  return(test_datasets)
}

# Function 6: Setup predictor matrix and imputation methods
# =============================================================================
setup_imputation_parameters <- function(data, survey_vars) {
  # Create predictor matrix 
  pred <- mice::make.predictorMatrix(data)
  # Exclude survey variables from being used as predictors or being imputed
  for (var in survey_vars) {
    if (var %in% colnames(pred)) {
      pred[var, ] <- 0  # Do not use this variable to predict others
      pred[, var] <- 0  # Do not predict this variable
    }
  }
  # Do not impute the math plausible value
  if ("MATH_PV" %in% colnames(pred)) {
    pred["MATH_PV", ] <- 0
  }
  # Setup imputation methods
  impMethod <- character(ncol(data))
  names(impMethod) <- colnames(data)
  # Default method for numeric variables is pmm
  for (var in names(impMethod)) {
    if (is.numeric(data[[var]])) {
      impMethod[var] <- "pmm"
    } else if (is.factor(data[[var]]) || is.character(data[[var]])) {
      impMethod[var] <- "polyreg"
    }
  }
  # Do not impute survey variables
  for (var in survey_vars) {
    if (var %in% names(impMethod)) {
      impMethod[var] <- ""
    }
  }
  # Do not impute MATH_PV
  if ("MATH_PV" %in% names(impMethod)) {
    impMethod["MATH_PV"] <- ""
  }
  return(list(pred = pred, impMethod = impMethod))
}

# Function 7: Perform imputation
# =============================================================================
perform_imputation <- function(data, pred, impMethod, m = 10, n.iter = 25, seed = NULL) {
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # Perform imputation
  imp <- mice(data, 
              m = m,                # Number of imputed datasets
              method = impMethod,   # Imputation methods
              predictorMatrix = pred, # Predictor matrix
              n.iter = n.iter,      # Iterations between imputations
              maxit = 25,
              chainMean = TRUE,     # Added for convergence checking
              printFlag = TRUE)     # Show progress
  return(imp)
}

# =============================================================================
# Main workflow
# =============================================================================

# Step 1: Identify auxiliary variables
# ------------------------------------

# Define manually selected auxiliary variables
manual_auxiliary_vars <- c("STRATUM", paste0("ST263Q", sprintf("%02d", seq(2, 8, by = 2)), "JA"),
                           paste0("ST267Q", sprintf("%02d", 1:8), "JA"), 
                           paste0("ST273Q", sprintf("%02d", 1:7), "JA"), 
                           paste0("ST290Q", sprintf("%02d", 1:9), "WA"),
                           paste0("ST291Q", sprintf("%02d", 1:10), "JA"),
                           paste0("ST293Q", sprintf("%02d", 1:8), "JA"),
                           "ANXMAT", "MATHPERS", "MATHEF21", "MATHEFF", 
                           "DISCLIM", "GROSAGR", "RELATST")

# Find auxiliary variables
result <- select_auxiliary_variables(Master_data, target_vars, manual_auxiliary_vars)
auxiliary_vars <- result$recommended_auxiliaries

# Step 2: Identify math PV variables
# ----------------------------------

math_pv_vars <- grep("PV.*MATH", names(Master_data), value = TRUE)
cat("Found", length(math_pv_vars), "math plausible value variables\n")

# Step 3: Create datasets with one PV each
# ----------------------------------------

pv_datasets <- create_pv_datasets(
  data = Master_data, 
  target_vars = target_vars,
  survey_vars = survey_vars,
  auxiliary_vars = auxiliary_vars,
  math_pv_vars = math_pv_vars
)

# Step 4: Split datasets by country
# ---------------------------------

country_pv_datasets <- split_by_country(pv_datasets, countries)

# =============================================================================
# Running test imputation
# =============================================================================

# Create test subsets
# ---------------------------

test_datasets <- create_test_subset(country_pv_datasets, n_per_country = 100)

# Perform test imputation
# ------------------------------

test_imputation_results <- list()

# For each country
for (country in countries) {
  test_imputation_results[[country]] <- list()
  # For each PV dataset
  for (pv_num in seq_along(test_datasets[[country]])) {
    # Get current test data
    test_data <- test_datasets[[country]][[pv_num]]
    # Setup imputation parameters
    imp_params <- setup_imputation_parameters(test_data, survey_vars)
    # Set seed based on country and PV for reproducibility
    seed_value <- as.numeric(factor(country)) * 100 + pv_num + 42
    cat("\nPerforming test imputation for country:", country, "- PV", pv_num, "\n")
    # Perform imputation
    imp <- perform_imputation(
      data = test_data,
      pred = imp_params$pred,
      impMethod = imp_params$impMethod,
      m = 5,           # 5 imputations for test run
      n.iter = 10,      # Iterations for test run
      maxit = 10,
      seed = seed_value
    )
    # Store results
    test_imputation_results[[country]][[pv_num]] <- imp
    cat("Completed test imputation for country:", country, "- PV", pv_num, "\n")
  }
}

# Save test imputation results
saveRDS(test_imputation_results, "pisa_2022_test_imputation_results.rds")
cat("Test imputation results saved to pisa_2022_test_imputation_results.rds\n")

# Load imputation results
test_imputation_results <- readRDS("pisa_2022_test_imputation_results.rds")

# =============================================================================
# Test imputation: Diagnostics
# =============================================================================

# Visual diagnostics for imputation quality
assess_imputation_quality <- function(imputation_results, key_variables = NULL, 
                                      plot_dir = "imputation_diagnostics") {
  # Create directory for plots if it doesn't exist
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  # Process each country
  for (country in names(imputation_results)) {
    cat("\nAssessing imputation quality for country:", country, "\n")
    # Create a PDF for each country
    pdf(file.path(plot_dir, paste0("imputation_quality_", country, ".pdf")))
    # Process each PV dataset
    for (pv_num in seq_along(imputation_results[[country]])) {
      # Get current imputation object
      imp <- imputation_results[[country]][[pv_num]]
      # Check if the object is a proper mice imputation
      if (!inherits(imp, "mids")) {
        warning(paste("Object for country", country, "PV", pv_num, "is not a valid mice imputation"))
        next
      }
      # If key variables are not specified, use all imputed variables
      if (is.null(key_variables)) {
        vars_to_check <- names(imp$imp)
      } else {
        # Only include variables that were actually imputed
        vars_to_check <- intersect(key_variables, names(imp$imp))
      }
      if (length(vars_to_check) == 0) {
        warning(paste("No variables to check for country", country, "PV", pv_num))
        next
      }
      # Get complete datasets
      complete_datasets <- mice::complete(imp, "all")
      # Add a title for this PV set
      plot.new()
      title(main = paste("Imputation Diagnostics -", country, "- PV", pv_num))
      # For each key variable, assess imputation quality
      for (var in vars_to_check) {
        # Skip variables with no imputations
        if (!var %in% names(imp$imp) || ncol(imp$imp[[var]]) == 0) next
        # 1. Compare original vs imputed distributions
        par(mfrow = c(1, 2))
        # Original data with missings
        original_data <- imp$data[[var]]
        hist(original_data, main = "Original Data (with NA)", 
             xlab = var, breaks = 20, col = "lightblue")
        # Density plots of all imputations
        plot(density(original_data, na.rm = TRUE), 
             main = "Density: Original vs Imputed", 
             xlab = var, col = "black", lwd = 2)
        # Add density curves for each imputation
        colors <- rainbow(length(complete_datasets))
        for (i in seq_along(complete_datasets)) {
          lines(density(complete_datasets[[i]][[var]], na.rm = TRUE), 
                col = colors[i], lwd = 1, lty = 2)
        }
        legend("topright", c("Original", paste0("Imp ", seq_along(complete_datasets))),
               col = c("black", colors), lwd = c(2, rep(1, length(colors))),
               lty = c(1, rep(2, length(colors))), cex = 0.7)
        # 2. Boxplots of all imputations
        par(mfrow = c(1, 1))
        # Prepare data for boxplot
        boxplot_data <- lapply(complete_datasets, function(ds) ds[[var]])
        names(boxplot_data) <- paste0("Imp.", seq_along(boxplot_data))
        boxplot(boxplot_data, main = paste("Boxplots of", var, "across imputations"),
                col = "lightblue", outpch = 20, outcex = 0.5)
        # 3. Compare means and standard deviations
        means <- sapply(boxplot_data, mean, na.rm = TRUE)
        sds <- sapply(boxplot_data, sd, na.rm = TRUE)
        par(mfrow = c(1, 2))
        # Plot means
        dotchart(means, main = paste("Means of", var), 
                 pch = 19, col = "blue", xlab = "Mean")
        abline(v = mean(original_data, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
        legend("topright", c("Original mean"), col = "red", lwd = 2, lty = 2, cex = 0.7)
        # Plot standard deviations
        dotchart(sds, main = paste("SDs of", var), 
                 pch = 19, col = "blue", xlab = "Standard Deviation")
        abline(v = sd(original_data, na.rm = TRUE), col = "red", lwd = 2, lty = 2)
        legend("topright", c("Original SD"), col = "red", lwd = 2, lty = 2, cex = 0.7)
      }
    }
    dev.off()
    cat("Created diagnostic plots for country:", country, "\n")
  }
}

# Run the assessment on test imputation results
assess_imputation_quality(
  imputation_results = test_imputation_results, 
  key_variables = key_vars
)

# =============================================================================
# Full imputation run
# =============================================================================

full_imputation_results <- list()

# For each country
for (country in countries) {
  full_imputation_results[[country]] <- list()
  # For each PV dataset 
  for (pv_num in seq_along(country_pv_datasets[[country]])) {
    # Get current data
    full_data <- country_pv_datasets[[country]][[pv_num]]
    # Setup imputation parameters
    imp_params <- setup_imputation_parameters(full_data, survey_vars)
    # Set seed based on country and PV for reproducibility
    seed_value <- as.numeric(factor(country)) * 100 + pv_num + 1000
    set.seed(seed_value)
    cat("\nPerforming full imputation for country:", country, "- PV", pv_num, "\n")
    # Perform imputation
    imp <- mice(full_data, 
                m = 10, # 10 imputations for full run
                method = imp_params$impMethod,
                predictorMatrix = imp_params$pred,
                n.iter = 25,        # Iterations for full run
                maxit = 25,
                chainMean = TRUE,
                printFlag = TRUE)
    # Store results
    full_imputation_results[[country]][[pv_num]] <- imp
    cat("Completed full imputation for country:", country, "- PV", pv_num, "\n")
    # Save intermediate results
    saveRDS(imp, paste0("pisa_2022_", country, "_PV", pv_num, "_imputation.rds"))
  }
}

# # Save full imputation results
saveRDS(full_imputation_results, "pisa_2022_full_imputation_results.rds")
cat("Full imputation results saved to pisa_2022_full_imputation_results.rds\n")

# # Load full imputation results
full_imputation_results <- readRDS("~/UiO/MAE4191/Coding/Imputation/pisa_2022_full_imputation_results.rds")

# =============================================================================
# Restructuring
# =============================================================================

# After full imputation is complete

imputation_results <- list()

# Retrieve and restructure imputation object
for (country in countries) {
  imputation_results[[country]] <- list()
  for (pv_num in 1:10) {  # 10 PVs
    # Get the mids object
    imp_object <- full_imputation_results[[country]][[pv_num]]
    # Extract all 10 complete datasets as a list
    complete_datasets <- lapply(1:10, function(i) {
      complete(imp_object, i)
    })
    # Store them
    imputation_results[[country]][[pv_num]] <- complete_datasets
  }
}

# Rename object to maintain overview
imputation_results_unreversed <- imputation_results

# Save full imputation results
saveRDS(imputation_results_unreversed, "pisa_2022_structured_imputation_results_unreversed.rds")
cat("Structured imputation results saved to pisa_2022_structured_imputation_results_unreversed.rds\n")

# # Load full imputation results
imputation_results_unreversed <- readRDS("~/UiO/MAE4191/Coding/Imputation/pisa_2022_structured_imputation_results_unreversed.rds")

# -----------------------------------------------------------------------------

# Function to reverse code Mathematics Anxiety in a dataset
reverse_code <- function(data) {
  # Reverse code ANXMAT (higher values = lower anxiety = higher achievement)
  data$ANXMAT <- -data$ANXMAT
  return(data)
}

# Apply reverse coding to all imputed datasets
for (country in names(imputation_results_unreversed)) {
  for (pv in 1:10) {
    for (imp in 1:10) {
      imputation_results_unreversed[[country]][[pv]][[imp]] <- reverse_code(imputation_results_unreversed[[country]][[pv]][[imp]])
    }
  }
}

# Rename object to maintain control
imputation_results <- imputation_results_unreversed

# -----------------------------------------------------------------------------

# Save full imputation results
saveRDS(imputation_results, "pisa_2022_ready_imputation_results.rds")
cat("Ready imputation results saved to pisa_2022_ready_imputation_results.rds\n")

# Load full imputation results
imputation_results <- readRDS("~/UiO/MAE4191/Coding/Imputation/pisa_2022_ready_imputation_results.rds")

# -----------------------------------------------------------------------------

# Function to count NAs per dataset
check_na <- lapply(imputation_results, function(country_data) {
  lapply(country_data, function(pv_data) {
    sapply(pv_data, function(df) sum(is.na(df)))
  })
})

# Print summary of NA counts
print(check_na)

# =============================================================================
# Full imputation: Diagnostics
# =============================================================================

# Define countries and variables
countries <- c("NOR", "SWE", "DNK", "FIN", "CHE")
variables <- c("ANXMAT", "MATHPERS", "MATHEFF", "MATHEF21", "DISCLIM", "GROSAGR", "RELATST", "ESCS", "IMMIG")

# Extract and combine all important variable means for thesis figures
create_combined_variable_means <- function(country) {
  # Create a dataframe to hold mean values for all variables across all PVs
  all_means <- data.frame()
  # For each variable
  for (var in variables) {
    var_means <- data.frame()
    # For each PV
    for (pv_num in 1:10) {
      imp_object <- full_imputation_results[[country]][[pv_num]]
      # Check if we can extract means from the imputation object
      if (inherits(imp_object, "mids")) {
        # Calculate the mean of the variable for each imputation in this PV
        for (m in 1:imp_object$m) {
          imputed_data <- mice::complete(imp_object, m)
          if (var %in% colnames(imputed_data)) {
            mean_val <- mean(imputed_data[[var]], na.rm = TRUE)
            var_means <- rbind(var_means, data.frame(
              variable = var,
              pv = pv_num,
              imputation = m,
              mean = mean_val
            ))
          }
        }
      }
    }
    # Add to overall means
    if (nrow(var_means) > 0) {
      all_means <- rbind(all_means, var_means)
    }
  }
  # If we have data, create summary plots
  if (nrow(all_means) > 0) {
    # Create boxplot of means across imputations by PV
    p1 <- ggplot(all_means, aes(x = as.factor(pv), y = mean, fill = as.factor(pv))) +
      geom_boxplot() +
      facet_wrap(~ variable, scales = "free_y") +
      labs(title = paste("Variable Means Across Imputations -", country),
           x = "Plausible Value",
           y = "Mean") +
      theme_minimal() +
      theme(legend.position = "none")
    # Save plot
    ggsave(paste0("variable_means_", country, ".pdf"), p1, width = 11, height = 8.5)
    ggsave(paste0("variable_means_", country, ".png"), p1, width = 11, height = 8.5, dpi = 300)
    # Create convergence line plot (mean across iterations)
    p2 <- ggplot(all_means, aes(x = imputation, y = mean, color = as.factor(pv), group = as.factor(pv))) +
      geom_line() +
      geom_point() +
      facet_wrap(~ variable, scales = "free_y") +
      labs(title = paste("Convergence of Means Across Imputations -", country),
           x = "Imputation Number",
           y = "Mean") +
      theme_minimal()
    # Save plot
    ggsave(paste0("convergence_means_", country, ".pdf"), p2, width = 11, height = 8.5)
    ggsave(paste0("convergence_means_", country, ".png"), p2, width = 11, height = 8.5, dpi = 300)
    cat("Created variable means plots for", country, "\n")
  }
}

# Create variable means plots for each country
for (country in countries) {
  create_combined_variable_means(country)
}

#######################################
#         MULTICOLLINEARITY           #
#######################################

calculate_vif_tolerance <- function(data, country_name) {
  # Create survey design using existing function
  design <- create_survey_design(data, country_name)
  # Create a formula with all predictors
  formula_str <- "MATH_PV ~ ANXMAT + MATHPERS + MATHEFF + DISCLIM + GROSAGR + RELATST" # + ESCS + IMMIG + MATHEF21"
  formula_obj <- as.formula(formula_str)
  # Fit survey-weighted linear model
  model <- svyglm(formula_obj, design = design)
  # Calculate VIF using car::vif
  vif_values <- car::vif(model)
  # Calculate tolerance values (1/VIF)
  tolerance_values <- 1 / vif_values
  # Return as list
  return(list(
    vif = vif_values,
    tolerance = tolerance_values
  ))
}

# -----------------------------------------------------------------------------

# Function to apply Rubin's Rules to VIF values
pool_vif_with_rubin <- function(vif_list) {
    # Extract variable names from the first result
      var_names <- names(vif_list[[1]]$vif)
        # Initialize result matrices
        n_vars <- length(var_names)
        n_imputations <- length(vif_list)
          # Matrix to store VIF values from each imputation
          vif_matrix <- matrix(NA, nrow = n_imputations, ncol = n_vars)
          colnames(vif_matrix) <- var_names
            # Matrix to store tolerance values from each imputation
            tolerance_matrix <- matrix(NA, nrow = n_imputations, ncol = n_vars)
            colnames(tolerance_matrix) <- var_names
              # Fill matrices with values from each imputation
              for (i in 1:n_imputations) {
                vif_matrix[i, ] <- vif_list[[i]]$vif
                tolerance_matrix[i, ] <- vif_list[[i]]$tolerance
                }
              # Calculate pooled results using Rubin's rules
              # For VIF and tolerance, we use the mean across imputations
              pooled_vif <- colMeans(vif_matrix)
                # Calculate within-imputation variance
                # For VIF we're not working with standard errors, so we calculate variances directly
                within_var_vif <- apply(vif_matrix, 2, var) * (n_imputations - 1) / n_imputations
                  # Calculate between-imputation variance
                  between_var_vif <- apply(vif_matrix, 2, function(x) var(x))
                    # Total variance using Rubin's formula
                    total_var_vif <- within_var_vif + between_var_vif + (between_var_vif / n_imputations)
                      # Standard errors
                      se_vif <- sqrt(total_var_vif)
                        # Repeat for tolerance
                        pooled_tolerance <- colMeans(tolerance_matrix)
                        within_var_tolerance <- apply(tolerance_matrix, 2, var) * (n_imputations - 1) / n_imputations
                        between_var_tolerance <- apply(tolerance_matrix, 2, function(x) var(x))
                        total_var_tolerance <- within_var_tolerance + between_var_tolerance + (between_var_tolerance / n_imputations)
                        se_tolerance <- sqrt(total_var_tolerance)
                          # Create results data frame
                          results <- data.frame(
                            Variable = var_names,
                            VIF = pooled_vif,
                            VIF_SE = se_vif,
                            Tolerance = pooled_tolerance,
                            Tolerance_SE = se_tolerance
                            )
                            return(results)
}

# -----------------------------------------------------------------------------

# Function to calculate VIF across all imputations for a country using Rubin's Rules
calculate_country_vif_with_rubins <- function(imputation_results, country) {
    all_vif_results <- list()
      # Loop through all combinations of PVs and imputations
      for (pv in 1:10) {
        for (imp in 1:10) {
          # Get the current imputed dataset
            data <- imputation_results[[country]][[pv]][[imp]]
              # Calculate VIF and tolerance
              vif_result <- calculate_vif_tolerance(data, country)
                # Store the result
                all_vif_results <- append(all_vif_results, list(vif_result))
                }
        }
      # Apply Rubin's Rules to pool results across imputations
      pooled_results <- pool_vif_with_rubin(all_vif_results)
        # Add country information
        pooled_results$Country <- country
          return(pooled_results)
}

# -----------------------------------------------------------------------------

# Function to run multicollinearity diagnostics for all countries
run_multicollinearity_diagnostics <- function(imputation_results, countries) {
    # Initialize list to store results
      vif_results_by_country <- list()
        # Run VIF analysis for each country
        for (country in countries) {
          vif_results_by_country[[country]] <- calculate_country_vif_with_rubins(imputation_results, country)
          # Print progress
            cat("Completed VIF analysis for", country, "\n")
          }
        # Combine all results into a single data frame
        all_results <- do.call(rbind, vif_results_by_country)
          # Create a summary table with interpretation
          summary_table <- all_results
          summary_table$Multicollinearity <- "No concern"
          summary_table$Multicollinearity[summary_table$VIF > 5] <- "Moderate concern"
          summary_table$Multicollinearity[summary_table$VIF > 10] <- "High concern"
            # Add confidence intervals (using 1.96 * SE for approximate 95% CI)
            summary_table$VIF_CI_Lower <- summary_table$VIF - 1.96 * summary_table$VIF_SE
            summary_table$VIF_CI_Upper <- summary_table$VIF + 1.96 * summary_table$VIF_SE
            summary_table$Tolerance_CI_Lower <- summary_table$Tolerance - 1.96 * summary_table$Tolerance_SE
            summary_table$Tolerance_CI_Upper <- summary_table$Tolerance + 1.96 * summary_table$Tolerance_SE
              # Return both detailed and summary results
              return(list(
                detailed_results = all_results,
                summary = summary_table
                ))
}

# -----------------------------------------------------------------------------

# Run the analysis
countries <- c("NOR", "SWE", "DNK", "FIN", "CHE")
multicollinearity_results <- run_multicollinearity_diagnostics(imputation_results, countries)

# Print summary results
print(multicollinearity_results$summary)

#######################################
#            MAIN ANALYSIS            #
#######################################

# =============================================================================
# Functions
# =============================================================================

# Function to create survey design
create_survey_design <- function(data, country) {
  country_data <- subset(data, CNT == country)
  if (nrow(country_data) == 0) {
    return(NULL)
  }
  # Identify strata and PSU distribution
  stratum_psu_counts <- table(country_data$STRATUM, country_data$CNTSCHID)
  # Determine which strata have insufficient PSUs
  strata_with_insufficient_psus <- rownames(stratum_psu_counts)[rowSums(stratum_psu_counts > 0) < 2]
  # Modify strata when necessary
  if (length(strata_with_insufficient_psus) > 0) {
    country_data$MODIFIED_STRATUM <- ifelse(
      country_data$STRATUM %in% strata_with_insufficient_psus, 
      "GENERIC_STRATUM", 
      as.character(country_data$STRATUM)
    )
    design <- svydesign(
      ids = ~CNTSCHID,
      weights = ~W_FSTUWT,
      strata = ~MODIFIED_STRATUM,
      data = country_data,
      nest = TRUE,
      repweights = "W_FSTURWT[1-80]+",
      type = "Fay",
      rho = 0.5,
      scale = 1,
      rscales = 1
    )
  } else {
    design <- svydesign(
      ids = ~CNTSCHID,
      weights = ~W_FSTUWT,
      strata = ~STRATUM,
      data = country_data,
      nest = TRUE,
      repweights = "W_FSTURWT[1-80]+",
      type = "Fay",
      rho = 0.5,
      scale = 1,
      rscales = 1
    )
  }
  return(design)
}

# -----------------------------------------------------------------------------

# Function to perform regression for a single dataset
run_regression <- function(design) {
  svyglm(
    MATH_PV ~ ANXMAT + MATHPERS + MATHEFF + 
      DISCLIM + GROSAGR + RELATST, # + MATHEF21 + IMMIG + ESCS) 
    design = design,
    family = gaussian()
  )
}

# -----------------------------------------------------------------------------

# Function to perform regression analysis
perform_regression_analysis <- function(imputation_results, countries) {
  regression_results <- list()
  for (country in countries) {
    all_regressions <- list()
    for (pv in 1:10) {
      for (imp in 1:10) {
        # Get imputed dataset
        data <- imputation_results[[country]][[pv]][[imp]]
        # Create survey design
        design <- create_survey_design(data, country)
        if (is.null(design)) next
        # Run regression
        reg <- run_regression(design)
        # Store regression results
        all_regressions <- append(all_regressions, list(reg))
      }
    }
    # Pool regression results using Rubin's rules
    pooled_reg <- MIcombine(all_regressions)
    # Save results
    regression_results[[country]] <- pooled_reg
  }
  return(regression_results)
}

# -----------------------------------------------------------------------------

# Function to format and summarize regression results
summarize_regressions <- function(regression_results) {
  summary_table <- data.frame()
  for (country in names(regression_results)) {
    reg_result <- regression_results[[country]]
    # Extract coefficients, SE, and p-values
    coefs <- coef(reg_result)
    se <- summary(reg_result)$se
    t_vals <- coefs / se
    p_vals <- 2 * pt(-abs(t_vals), reg_result$df[2])
    # Create country rows
    for (var in names(coefs)) {
      row <- data.frame(
        Country = country,
        Variable = var,
        Coefficient = coefs[var],
        SE = se[var],
        t_value = t_vals[var],
        p_value = p_vals[var],
        Significant = ifelse(p_vals[var] < 0.05, "*", "")
      )
      summary_table <- rbind(summary_table, row)
    }
  }
  return(summary_table)
}

# -----------------------------------------------------------------------------

# Function to perform model comparison
perform_model_comparison <- function(imputation_results, countries) {
  comparison_results <- list()
  for (country in countries) {
    # Store all null and full models for this country
    null_models <- list()
    full_models <- list()
    for (pv in 1:10) {
      for (imp in 1:10) {
        # Get imputed dataset
        data <- imputation_results[[country]][[pv]][[imp]]
        # Create survey design
        design <- create_survey_design(data, country)
        if (is.null(design)) next
        # Run null model (intercept only)
        null_model <- svyglm(MATH_PV ~ 1, design = design, family = gaussian())
        # Run full model (with all predictors)
        full_model <- svyglm(
          MATH_PV ~ ANXMAT + MATHPERS + MATHEFF + 
            DISCLIM + GROSAGR + RELATST + ESCS,
          design = design,
          family = gaussian()
        )
        # Store models
        null_models <- append(null_models, list(null_model))
        full_models <- append(full_models, list(full_model))
      }
    }
    # Pool null and full models using Rubin's rules
    pooled_null <- MIcombine(null_models)
    pooled_full <- MIcombine(full_models)
    # Calculate model fit statistics
    n_params_null <- length(coef(pooled_null))
    n_params_full <- length(coef(pooled_full))
    df_diff <- n_params_full - n_params_null
    # Calculate F-statistic using Wald test approach
    # Get the parameters that are in full but not in null model
    params_to_test <- setdiff(names(coef(pooled_full)), names(coef(pooled_null)))
    # Extract coefficients and variance-covariance matrix for those parameters
    coefs_to_test <- coef(pooled_full)[params_to_test]
    vcov_to_test <- vcov(pooled_full)[params_to_test, params_to_test]
    # Calculate Wald statistic (F-test)
    wald_stat <- t(coefs_to_test) %*% solve(vcov_to_test) %*% coefs_to_test
    f_stat <- wald_stat / df_diff
    # Calculate degrees of freedom for the denominator
    # Using Barnard and Rubin (1999) small sample adjustment
    df_complete <- pooled_full$df[2]  # Complete-data degrees of freedom
    # Calculate p-value
    p_value <- pf(f_stat, df_diff, df_complete, lower.tail = FALSE)
    # Store results
    comparison_results[[country]] <- list(
      F_statistic = f_stat,
      numerator_df = df_diff,
      denominator_df = df_complete,
      p_value = p_value,
      pooled_null = pooled_null,
      pooled_full = pooled_full
    )
  }
  return(comparison_results)
}

# -----------------------------------------------------------------------------

# Function to summarize model comparison results
summarize_model_comparison <- function(comparison_results) {
  summary_table <- data.frame(
    Country = character(),
    F_statistic = numeric(),
    numerator_df = numeric(),
    denominator_df = numeric(),
    p_value = numeric(),
    Significant = character(),
    stringsAsFactors = FALSE
  )
  for (country in names(comparison_results)) {
    result <- comparison_results[[country]]
    # Create row for this country
    row <- data.frame(
      Country = country,  
      F_statistic = as.numeric(result$F_statistic),
      numerator_df = result$numerator_df,
      denominator_df = result$denominator_df,
      p_value = as.numeric(result$p_value),
      Significant = ifelse(result$p_value < 0.05, "*", ""),
      stringsAsFactors = FALSE
    )
    # Add to summary table
    summary_table <- rbind(summary_table, row)
    rownames(summary_table) <- NULL
  }
  return(summary_table)
}

# =============================================================================
# Main Analysis
# =============================================================================

# Regression
# =============================================================================

# Define countries
countries <- c("NOR", "SWE", "DNK", "FIN", "CHE")

# Calculate regression results
regression_results <- perform_regression_analysis(imputation_results, countries)
regression_summary <- summarize_regressions(regression_results)

# Print regression results
print("Regression Results:")
print(regression_summary)

# Save results
saveRDS(regression_results, "pisa_regression_results.rds")
saveRDS(regression_summary, "pisa_regression_summary.rds")

# Model Comparison
# =============================================================================

# Perform model comparison
comparison_results <- perform_model_comparison(imputation_results, countries)
comparison_summary <- summarize_model_comparison(comparison_results)

# Print comparison results
print("Model Comparison Results:")
print(comparison_summary)

# Save results
saveRDS(comparison_results, "pisa_model_comparison_results.rds")
write.csv(comparison_summary, "pisa_model_comparison_summary.csv", row.names = FALSE)

# =============================================================================
# R2
# =============================================================================

# Function to calculate R2
calculate_r2_survey <- function(model) {
  # Check if it's a svyglm object
  if (!inherits(model, "svyglm")) {
    warning("Model is not a svyglm object, returning NA")
    return(NA)
  }
  # Extract the model summary to get deviances
  model_summary <- summary(model)
  # Calculate R2 as 1 - (residual deviance / null deviance)
  r2 <- 1 - (model$deviance / model$null.deviance)
  # Ensure R2 is within valid range [0,1]
  r2 <- max(0, min(1, r2))
  return(r2)
}

# -----------------------------------------------------------------------------

# Function to pool based on averages. 
# Rubin`s Rules do not apply for non-normally distributed estimates 
pool_r2 <- function(r2_values) {
  pooled_r2 <- mean(r2_values)
  variance <- var(r2_values)
  se <- sqrt(variance / length(r2_values))
  ci_lower <- pooled_r2 - 1.96 * se
  ci_upper <- pooled_r2 + 1.96 * se
  return(list(
    pooled_r2 = pooled_r2,
    ci_lower = max(0, ci_lower),
    ci_upper = min(1, ci_upper)
  ))
}

# -----------------------------------------------------------------------------

# Function to calculate and pool R2 
calculate_and_pool_r2 <- function(imputation_results, countries) {
  r2_results <- list()
  for (country in countries) {
    all_r2_imp <- numeric(0)  # Initialize as numeric vector
    for (imp in 1:10) {
      all_r2_pv <- numeric(0)  
      for (pv in 1:10) {
        # Get imputed dataset
        data <- imputation_results[[country]][[pv]][[imp]]
        if (is.null(data)) {
          message(paste("Missing data for", country, "PV", pv, "imp", imp))
          next
        }
        # Create survey design
        design <- create_survey_design(data, country)
        if (is.null(design)) {
          message(paste("Design failed for", country, "PV", pv, "imp", imp))
          next
        }
        # Run regression
        model <- tryCatch({
          run_regression(design)
        }, error = function(e) {
          message(paste("Model error:", country, "PV", pv, "imp", imp, "||", e$message))
          return(NULL)
        })
        # Calculate R2
        if (!is.null(model)) {
          r2 <- calculate_r2_survey(model)
          if (is.numeric(r2)) {
            all_r2_pv <- c(all_r2_pv, r2)  # Add to numeric vector
          }
        }
      }
      # Store average for this imputation only if we have values
      if (length(all_r2_pv) > 0) {
        avg_r2 <- mean(all_r2_pv)
        all_r2_imp <- c(all_r2_imp, avg_r2)  # Store in numeric vector
      } else {
        message(paste("No valid PVs for", country, "imp", imp))
      }
    }
    # Pool R2 values if we have any
    if (length(all_r2_imp) > 0) {
      r2_results[[country]] <- pool_r2(all_r2_imp)  # Call pool_r2 
    } else {
      message(paste("No valid imputations for", country))
      r2_results[[country]] <- list(pooled_r2 = NA, ci_lower = NA, ci_upper = NA)
    }
  }
  return(r2_results)
}

# -----------------------------------------------------------------------------

# Execute the analysis
countries <- c("NOR", "SWE", "DNK", "FIN", "CHE")
r2_results <- calculate_and_pool_r2(imputation_results, countries)

# Print results
for (country in countries) {
  cat("\nPooled R^2 for", country, ":\n")
  print(r2_results[[country]])
}

# Save results
saveRDS(r2_results, "pisa_r2_results.rds")

# =============================================================================
# Linear Model Assumptions
# =============================================================================

# Function to check normality of residuals for pooled MI results
check_normality_mi <- function(regression_results, imputation_results, countries) {
  normality_tests <- list()
  for (country in countries) {
    # Extract residuals from each imputation separately
    all_residuals <- numeric(0)
    # Loop through plausible values and imputations
    for (pv in 1:10) {
      for (imp in 1:10) {
        # Get the data
        data <- imputation_results[[country]][[pv]][[imp]]
        # Create design
        design <- create_survey_design(data, country)
        if (is.null(design)) next
        # Run the regression
        model <- svyglm(
          MATH_PV ~ ANXMAT + MATHPERS + MATHEFF +
            DISCLIM + GROSAGR + RELATST, # + MATHEF21 + ESCS + IMMIG,
          design = design,
          family = gaussian()
        )
        # Calculate residuals
        model_residuals <- residuals(model)
        all_residuals <- c(all_residuals, model_residuals)
      }
    }
    # Check if there are residuals to analyze
    if (length(all_residuals) > 0) {
      # Create a data frame to store residuals
      res_data <- data.frame(
        Residuals = all_residuals,
        Country = country
      )
      # Perform Shapiro-Wilk test if sample size permits
      if(length(all_residuals) <= 5000) {
        # Take a random sample if too many residuals
        sample_size <- min(5000, length(all_residuals))
        sample_idx <- sample(1:length(all_residuals), sample_size)
        sw_test <- shapiro.test(all_residuals[sample_idx])
        normality_tests[[country]] <- list(
          shapiro_p_value = sw_test$p.value,
          residuals_data = res_data
        )
      } else {
        # For large samples, just store the residuals for visual inspection
        normality_tests[[country]] <- list(
          shapiro_p_value = NA,
          residuals_data = res_data
        )
      }
    } else {
      # No residuals available
      normality_tests[[country]] <- NULL
    }
  }
  return(normality_tests)
}

# -----------------------------------------------------------------------------

# Function to check homoscedasticity for pooled MI results
check_homoscedasticity_mi <- function(regression_results, imputation_results, countries) {
  homoscedasticity_tests <- list()
  for (country in countries) {
    # extract fitted values and residuals from each imputation
    all_fitted <- numeric(0)
    all_residuals <- numeric(0)
    # Loop through plausible values and imputations
    for (pv in 1:10) {
      for (imp in 1:10) {
        # Get the data
        data <- imputation_results[[country]][[pv]][[imp]]
        # Create design
        design <- create_survey_design(data, country)
        if (is.null(design)) next
        # Run the regression
        model <- svyglm(
          MATH_PV ~ ANXMAT + MATHPERS + MATHEFF + 
            DISCLIM + GROSAGR + RELATST, # + MATHEF21 + IMMIG + ESCS,
          design = design,
          family = gaussian()
        )
        # Calculate fitted values and residuals
        model_fitted <- fitted(model)
        model_residuals <- residuals(model)
        all_fitted <- c(all_fitted, model_fitted)
        all_residuals <- c(all_residuals, model_residuals)
      }
    }
    # Check if we have data to analyze
    if (length(all_fitted) > 0 && length(all_residuals) > 0) {
      # Create data for residual vs fitted plot
      fitted_vs_residuals <- data.frame(
        Fitted = all_fitted,
        Residuals = all_residuals,
        Country = country
      )
      homoscedasticity_tests[[country]] <- list(
        fitted_vs_residuals = fitted_vs_residuals
      )
    } else {
      # No data available
      homoscedasticity_tests[[country]] <- NULL
    }
  }
  return(homoscedasticity_tests)
}

# -----------------------------------------------------------------------------

# Function to check linearity assumption for pooled MI results
check_linearity_mi <- function(imputation_results, countries) {
  linearity_checks <- list()
  for (country in countries) {
    # Use one imputation for visualization
    data <- NULL
    for (pv in 1:10) {
      for (imp in 1:10) {
        data <- imputation_results[[country]][[pv]][[imp]]
        if (!is.null(data) && nrow(data) > 0) break
      }
      if (!is.null(data) && nrow(data) > 0) break
    }
    if (!is.null(data) && nrow(data) > 0) {
      # Store predictor-outcome relationships
      predictor_outcome_data <- list()
      # Extract data for component+residual plots or partial residual plots
      predictors <- c("ANXMAT", "MATHPERS", "MATHEFF", "DISCLIM", 
                      "GROSAGR", "RELATST") # "MATHEF21", "IMMIG", "ESCS")
      for (predictor in predictors) {
        # Create data for scatter plots
        if (predictor %in% colnames(data)) {
          predictor_data <- data.frame(
            Predictor_Value = data[[predictor]],
            Outcome = data[["MATH_PV"]],
            Variable = predictor
          )
          # Remove NA values
          predictor_data <- predictor_data[!is.na(predictor_data$Predictor_Value) & 
                                             !is.na(predictor_data$Outcome), ]
          if (nrow(predictor_data) > 0) {
            predictor_outcome_data[[predictor]] <- predictor_data
          }
        }
      }
      linearity_checks[[country]] <- predictor_outcome_data
    } else {
      linearity_checks[[country]] <- NULL
    }
  }
  return(linearity_checks)
}

# -----------------------------------------------------------------------------

# Function to create diagnostic plots
create_diagnostic_plots_mi <- function(normality_tests, homoscedasticity_tests, linearity_checks) {
  for (country in names(normality_tests)) {
    if (is.null(normality_tests[[country]])) next
    # QQ plot for normality
    res_data <- normality_tests[[country]]$residuals_data
    if (!is.null(res_data) && nrow(res_data) > 0) {
      png(paste0("diagnostics_", country, "_qq_plot.png"), width = 800, height = 600)
      qqnorm(res_data$Residuals, main = paste("Normal Q-Q Plot for", country))
      qqline(res_data$Residuals)
      dev.off()
      # Histogram of residuals
      png(paste0("diagnostics_", country, "_residuals_hist.png"), width = 800, height = 600)
      hist(res_data$Residuals, main = paste("Histogram of Residuals for", country),
           xlab = "Residuals", breaks = 30)
      dev.off()
    }
    # Residual vs Fitted plot for homoscedasticity
    if (!is.null(homoscedasticity_tests[[country]])) {
      fitted_vs_res <- homoscedasticity_tests[[country]]$fitted_vs_residuals
      if (!is.null(fitted_vs_res) && nrow(fitted_vs_res) > 0) {
        png(paste0("diagnostics_", country, "_residuals_fitted.png"), width = 800, height = 600)
        plot(fitted_vs_res$Fitted, fitted_vs_res$Residuals, 
             main = paste("Residuals vs Fitted for", country),
             xlab = "Fitted values", ylab = "Residuals")
        abline(h = 0, col = "red", lty = 2)
        # Add loess smoothing line to check for patterns
        lines(lowess(fitted_vs_res$Fitted, fitted_vs_res$Residuals), col = "blue", lwd = 2)
        dev.off()
      }
    }
    # Scatter plots for linearity checks
    if (!is.null(linearity_checks[[country]])) {
      predictor_data <- linearity_checks[[country]]
      for (predictor in names(predictor_data)) {
        data <- predictor_data[[predictor]]
        if (!is.null(data) && nrow(data) > 0) {
          png(paste0("diagnostics_", country, "_", predictor, "_scatter.png"), width = 800, height = 600)
          plot(data$Predictor_Value, data$Outcome, 
               main = paste("Scatter plot of", predictor, "vs MATH_PV for", country),
               xlab = predictor, ylab = "MATH_PV")
          # Add loess smoothing line to check for non-linearity
          lines(lowess(data$Predictor_Value, data$Outcome), col = "blue", lwd = 2)
          dev.off()
        }
      }
    }
  }
}

# -----------------------------------------------------------------------------

# Check assumptions
normality_tests <- check_normality_mi(regression_results, imputation_results, countries)
homoscedasticity_tests <- check_homoscedasticity_mi(regression_results, imputation_results, countries)
linearity_checks <- check_linearity_mi(imputation_results, countries)

# Create diagnostic plots
create_diagnostic_plots_mi(normality_tests, homoscedasticity_tests, linearity_checks)

# Save diagnostic results
saveRDS(normality_tests, "pisa_normality_tests.rds")
saveRDS(homoscedasticity_tests, "pisa_homoscedasticity_tests.rds")
saveRDS(linearity_checks, "pisa_linearity_checks.rds")

#######################################
#              OUTLIERS               #
#######################################

# =============================================================================
# Before imputations - visual inspection
# =============================================================================

# List of variables to plot
variables <- c("ANXMAT", "MATHPERS", "MATHEFF", "DISCLIM",
               "GROSAGR", "RELATST") # "MATHEF21", "IMMIG", "ESCS")

# List of countries
countries <- c("NOR", "SWE", "DNK", "FIN", "CHE")

# Loop through each country
for (country in countries) {
  # Subset data for the current country
  country_data <- subset(Master_data, CNT == country)
  # Create survey design for the current country
  design <- svydesign(ids = ~CNTSCHID, weights = ~W_FSTUWT, strata = ~STRATUM, data = Master_data)
  # Loop through each variable and create a boxplot
  for (var in variables) {
    # Create a two-sided formula for svyboxplot
    formula <- as.formula(paste(var, "~ 1"))
    # Create the boxplot
    svyboxplot(formula, design, main = paste("Boxplot of", var, "in", country))
  }
}

# =============================================================================
# Before imputations - statistical inspection
# =============================================================================

# List of variables to analyze
variables <- c("ANXMAT", "MATHPERS", "MATHEFF", "DISCLIM",
               "GROSAGR", "RELATST") # "MATHEF21", "IMMIG", "ESCS")

# List of countries
countries <- c("NOR", "SWE", "DNK", "FIN", "CHE")

# Initialize the data frame to store pre-imputation outlier results
pre_imputation_outliers <- data.frame(
  Country = character(),
  Variable = character(),
  Count = integer(),
  Percentage = numeric(),
  Min_Z_Score = numeric(),
  Max_Z_Score = numeric(),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------

# Function to calculate z-scores and identify outliers
analyze_outliers <- function(data, variable) {
  # Skip if the variable does not exist in the data
  if (!variable %in% names(data)) {
    warning(paste("Variable", variable, "not found in dataset"))
    return(NULL)
  }
  # Skip if there are no non-NA values
  if (sum(!is.na(data[[variable]])) == 0) {
    warning(paste("No valid data for variable", variable))
    return(NULL)
  }
  # Calculate z-scores
  mean_val <- mean(data[[variable]], na.rm = TRUE)
  sd_val <- sd(data[[variable]], na.rm = TRUE)
  # Check for zero standard deviation
  if (sd_val == 0 || is.na(sd_val)) {
    warning(paste("Zero or NA standard deviation for variable", variable))
    return(NULL)
  }
  z_scores <- (data[[variable]] - mean_val) / sd_val
  # Identify outliers (|z| > 3)
  outliers <- which(abs(z_scores) > 3)
  # Return outlier indices and their z-scores
  if(length(outliers) > 0) {
    return(data.frame(
      index = outliers,
      value = data[[variable]][outliers],
      z_score = z_scores[outliers]
    ))
  } else {
    return(NULL)
  }
}

# -----------------------------------------------------------------------------

# Analysis for pre-imputation data
cat("\n=== Outlier Analysis Before Imputation ===\n")

for (country in countries) {
  cat("\nAnalyzing pre-imputation outliers for country:", country, "\n")
  # Subset data for the current country
  country_data <- subset(Master_data, CNT == country)
  # Loop through each variable
  for (var in variables) {
    # Perform outlier analysis
    outliers <- analyze_outliers(country_data, var)
    # Calculate percentage of outliers
    total_obs <- sum(!is.na(country_data[[var]]))
    outlier_count <- ifelse(is.null(outliers), 0, nrow(outliers))
    outlier_percentage <- ifelse(total_obs > 0, outlier_count / total_obs * 100, 0)
    # Add summary to the pre_imputation_outliers data frame
    pre_imputation_outliers <- rbind(pre_imputation_outliers, data.frame(
      Country = country,
      Variable = var,
      Count = outlier_count,
      Percentage = outlier_percentage,
      Min_Z_Score = ifelse(is.null(outliers), NA, min(outliers$z_score)),
      Max_Z_Score = ifelse(is.null(outliers), NA, max(outliers$z_score))
    ))
  }
}

# Print summary of pre-imputation outliers
print(pre_imputation_outliers)

# Save pre-imputation outlier results
write.csv(pre_imputation_outliers, "pre_imputation_outliers.csv", row.names = FALSE)

# =============================================================================
# After imputations - statistical inspection
# =============================================================================

# Function to compute outlier statistics for a single dataset
compute_outlier_stats <- function(data, variable) {
  # Skip if a variable does not exist or has no valid values
  if (!variable %in% names(data) || sum(!is.na(data[[variable]])) == 0) {
    return(list(
      count = 0,
      percentage = 0,
      min_z = NA,
      max_z = NA,
      variance = NA
    ))
  }
  # Calculate z-scores
  mean_val <- mean(data[[variable]], na.rm = TRUE)
  sd_val <- sd(data[[variable]], na.rm = TRUE)
  # Handle zero standard deviation
  if (sd_val == 0 || is.na(sd_val)) {
    return(list(
      count = 0,
      percentage = 0,
      min_z = NA,
      max_z = NA,
      variance = NA
    ))
  }
  z_scores <- (data[[variable]] - mean_val) / sd_val
  # Identify outliers (|z| > 3)
  outliers <- which(abs(z_scores) > 3)
  outlier_count <- length(outliers)
  # Calculate percentage
  total_obs <- sum(!is.na(data[[variable]]))
  outlier_percentage <- outlier_count / total_obs * 100
  # Calculate variance of the outlier count for Rubin's rules (binomial variance estimate)
  variance <- (outlier_percentage * (100 - outlier_percentage)) / total_obs
  # Return statistics
  return(list(
    count = outlier_count,
    percentage = outlier_percentage,
    min_z = ifelse(outlier_count > 0, min(abs(z_scores[outliers])), NA),
    max_z = ifelse(outlier_count > 0, max(abs(z_scores[outliers])), NA),
    variance = variance
  ))
}

# -----------------------------------------------------------------------------

# Initialize results storage
pooled_results <- data.frame(
  Country = character(),
  Variable = character(),
  Pooled_Count = numeric(),
  Pooled_Percentage = numeric(),
  Pooled_SE = numeric(),
  Min_Z = numeric(),
  Max_Z = numeric(),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------

# Outlier Analysis After Imputation with Rubin's Rules
all_imputation_results <- list()
for (country in countries) {
  cat("\nAnalyzing post-imputation outliers for country:", country, "\n")
  for (var in variables) {
    # Storage for this country-variable combination
    all_stats <- list()
    # For each PV 
    for (pv in 1:10) {
      # For each imputation
      for (imp in 1:10) {
        # Skip if the imputation dataset doesn't exist
        if (is.null(imputation_results[[country]][[pv]][[imp]])) {
          cat("Warning: Missing dataset for", country, "PV", pv, "Imputation", imp, "\n")
          next
        }
        # Extract the imputed dataset
        imputed_data <- imputation_results[[country]][[pv]][[imp]]
        # Compute outlier statistics
        stats <- compute_outlier_stats(imputed_data, var)
        # Store results
        all_stats[[length(all_stats) + 1]] <- stats
      }
    }
    # Apply Rubin's rules to pool results
    m <- length(all_stats)  # Number of imputations
    if (m == 0) {
      cat("Warning: No valid imputations for", country, "variable", var, "\n")
      next
    }
    # Extract statistics from all imputations
    counts <- sapply(all_stats, function(x) x$count)
    percentages <- sapply(all_stats, function(x) x$percentage)
    variances <- sapply(all_stats, function(x) x$variance)
    min_zs <- sapply(all_stats, function(x) x$min_z)
    max_zs <- sapply(all_stats, function(x) x$max_z)
    # Calculate pooled estimates using Rubin's rules
    pooled_percentage <- mean(percentages, na.rm = TRUE)
    # Within-imputation variance
    W <- mean(variances, na.rm = TRUE)
    # Between-imputation variance
    B <- var(percentages, na.rm = TRUE)
    # Total variance according to Rubin's rules
    T_var <- W + (1 + 1/m) * B
    # Standard error
    pooled_se <- sqrt(T_var)
    # Calculate pooled count
    pooled_count <- mean(counts, na.rm = TRUE)
    # Min and max z-scores (extremes across all imputations)
    min_z <- min(min_zs, na.rm = TRUE)
    max_z <- max(max_zs, na.rm = TRUE)
    # Add results to the pooled results dataframe
    pooled_results <- rbind(pooled_results, data.frame(
      Country = country,
      Variable = var,
      Pooled_Count = pooled_count,
      Pooled_Percentage = pooled_percentage,
      Pooled_SE = pooled_se,
      Min_Z = min_z,
      Max_Z = max_z
    ))
    # Store detailed results for debugging
    all_imputation_results[[paste(country, var, sep = "_")]] <- all_stats
  }
}

# -----------------------------------------------------------------------------

# Print pooled results
print(pooled_results)

# Save pooled results
write.csv(pooled_results, "pooled_outlier_analysis.csv", row.names = FALSE)

# -----------------------------------------------------------------------------

# Create comparison between pre and post imputation
comparison <- merge(pre_imputation_outliers, 
                    pooled_results, 
                    by = c("Country", "Variable"),
                    suffixes = c("_Pre", "_Post"))

# Calculate difference in percentage
comparison$Percentage_Diff <- comparison$Pooled_Percentage - comparison$Percentage

# Save comparison
write.csv(comparison, "outlier_comparison_pre_post_imputation.csv", row.names = FALSE)

# -----------------------------------------------------------------------------

# Visualizations
if (require(ggplot2)) {
  # Heat map of outlier percentages after imputation
  heatmap_post <- ggplot(pooled_results, 
                         aes(x = Variable, y = Country, fill = Pooled_Percentage)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    labs(title = "Percentage of Outliers After Imputation",
         x = "Variable", y = "Country", fill = "% Outliers")
  # Save plot
  ggsave("outlier_heatmap_post_imputation.png", heatmap_post, width = 10, height = 6)
  # Comparison plot
  if (nrow(comparison) > 0) {
    comparison_plot <- ggplot(comparison, aes(x = Variable)) +
      geom_point(aes(y = Percentage, color = "Pre-Imputation"), size = 3) +
      geom_point(aes(y = Pooled_Percentage, color = "Post-Imputation"), size = 3) +
      facet_wrap(~Country) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Comparison of Outlier Percentages Before and After Imputation",
           x = "Variable", y = "Percentage of Outliers", color = "Stage") +
      scale_color_manual(values = c("Pre-Imputation" = "blue", "Post-Imputation" = "red"))
    # Save plot
    ggsave("outlier_comparison_plot.png", comparison_plot, width = 12, height = 8)
  }
}

#######################################
#             STATISTICS              #
#######################################

# Means and SE
# =============================================================================

# Function to retrieve means and SE for all variables
get_predictor_means <- function(design) {
  # List of all predictor variables
  predictors <- c("ANXMAT", "MATHPERS", "MATHEFF", "DISCLIM",
                  "GROSAGR", "RELATST") # "MATHEF21", "IMMIG", "ESCS")
  # Create formula for svymean
  formula_str <- paste0("~", paste(predictors, collapse = "+"))
  pred_formula <- as.formula(formula_str)
  # Get means and variances
  means <- svymean(pred_formula, design)
  return(means)
}

# -----------------------------------------------------------------------------

# Function to retrieve dependent variable mean
get_outcome_mean <- function(design) {
  outcome_mean <- svymean(~MATH_PV, design)
  return(outcome_mean)
}

# -----------------------------------------------------------------------------

# Function to retrieve predictor means
pool_means_with_rubin <- function(means_list) {
  # Number of imputations
  m <- length(means_list)
  # Extract variable names from first imputation
  var_names <- names(means_list[[1]])
  # Initialize results
  pooled_means <- numeric(length(var_names))
  pooled_variances <- numeric(length(var_names))
  names(pooled_means) <- names(pooled_variances) <- var_names
  for (var in var_names) {
    # Extract estimates and variances for this variable
    estimates <- sapply(means_list, function(x) x[var])
    # Mean of the estimates (Q-bar)
    q_bar <- mean(estimates)
    # Between-imputation variance
    b <- var(estimates)
    # Within-imputation variance
    w <- mean(sapply(means_list, function(x) attr(x, "var")[var, var]))
    # Total variance according to Rubin's rules
    total_var <- w + b + (b/m)
    # Store results
    pooled_means[var] <- q_bar
    pooled_variances[var] <- total_var
  }
  return(list(means = pooled_means, variances = pooled_variances))
}

# -----------------------------------------------------------------------------

# Function to retrieve means, SE, and variances for all variables
get_descriptive_statistics <- function(imputation_results, countries) {
  descriptives <- list()
  for (country in countries) {
    all_predictor_means <- list()
    all_outcome_means <- list()
    all_predictor_se <- list()
    all_outcome_se <- list()
    for (pv in 1:10) {
      for (imp in 1:10) {
        # Get imputed dataset
        data <- imputation_results[[country]][[pv]][[imp]]
        # Create survey design
        design <- create_survey_design(data, country)
        if (is.null(design)) next
        # Calculate means and SE for all predictors
        predictor_means <- get_predictor_means(design)
        outcome_means <- get_outcome_mean(design)
        # Store results
        all_predictor_means <- append(all_predictor_means, list(coef(predictor_means)))
        all_outcome_means <- append(all_outcome_means, list(coef(outcome_means)))
        all_predictor_se <- append(all_predictor_se, list(sqrt(diag(attr(predictor_means, "var")))))
        all_outcome_se <- append(all_outcome_se, list(sqrt(diag(attr(outcome_means, "var")))))
      }
    }
    # Pool means and SE using Rubin's rules
    pooled_predictor_means <- apply(do.call(rbind, all_predictor_means), 2, mean)
    pooled_outcome_mean <- mean(unlist(all_outcome_means))
    pooled_predictor_se <- sqrt(apply(do.call(rbind, all_predictor_se)^2, 2, mean))
    pooled_outcome_se <- sqrt(mean(unlist(all_outcome_se)^2))
    # Store results
    descriptives[[country]] <- list(
      predictor_means = pooled_predictor_means,
      predictor_se = pooled_predictor_se,
      outcome_mean = pooled_outcome_mean,
      outcome_se = pooled_outcome_se
    )
  }
  return(descriptives)
}

# -----------------------------------------------------------------------------

# Function to format and summarize descriptive statistics
summarize_descriptives <- function(descriptives) {
  summary_table <- data.frame()
  for (country in names(descriptives)) {
    country_data <- descriptives[[country]]
    # Create country row for each predictor
    for (var in names(country_data$predictor_means)) {
      mean_val <- country_data$predictor_means[var]
      se_val <- country_data$predictor_se[var]
      row <- data.frame(
        Country = country,
        Variable = var,
        Mean = mean_val,
        SE = se_val,
        CI_Lower = mean_val - 1.96 * se_val,
        CI_Upper = mean_val + 1.96 * se_val
      )
      summary_table <- rbind(summary_table, row)
    }
    # Add row for outcome variable
    row <- data.frame(
      Country = country,
      Variable = "MATH_PV",
      Mean = country_data$outcome_mean,
      SE = country_data$outcome_se,
      CI_Lower = country_data$outcome_mean - 1.96 * country_data$outcome_se,
      CI_Upper = country_data$outcome_mean + 1.96 * country_data$outcome_se
    )
    summary_table <- rbind(summary_table, row)
  }
  return(summary_table)
}

# -----------------------------------------------------------------------------

# Define countries
countries <- c("NOR", "SWE", "DNK", "FIN", "CHE")

# Calculate descriptive statistics
descriptive_results <- get_descriptive_statistics(imputation_results, countries)
descriptive_summary <- summarize_descriptives(descriptive_results)

# Print descriptive statistics
print("Descriptive Statistics:")
print(descriptive_summary)

# Save results
saveRDS(descriptive_results, "pisa_descriptive_statistics.rds")
saveRDS(descriptive_summary, "pisa_descriptive_summary.rds")

# Skewness, kurtosis and correlations
# =============================================================================

# Function to calculate survey-weighted skewness
calculate_weighted_skewness <- function(design, variable) {
  # Calculate weighted moments
  m1 <- svymean(as.formula(paste0("~", variable)), design)
  m1 <- coef(m1)
  # Create formula for second moment
  formula_m2 <- as.formula(paste0("~I((", variable, " - ", m1, ")^2)"))
  m2 <- svymean(formula_m2, design)
  m2 <- coef(m2)
  # Create formula for third moment
  formula_m3 <- as.formula(paste0("~I((", variable, " - ", m1, ")^3)"))
  m3 <- svymean(formula_m3, design)
  m3 <- coef(m3)
  # Calculate skewness
  skewness <- m3 / (m2^(3/2))
  return(skewness)
}

# =============================================================================

# Function to calculate survey-weighted kurtosis
calculate_weighted_kurtosis <- function(design, variable) {
  # Calculate weighted moments
  m1 <- svymean(as.formula(paste0("~", variable)), design)
  m1 <- coef(m1)
  # Create formula for second moment
  formula_m2 <- as.formula(paste0("~I((", variable, " - ", m1, ")^2)"))
  m2 <- svymean(formula_m2, design)
  m2 <- coef(m2)
  # Create formula for fourth moment
  formula_m4 <- as.formula(paste0("~I((", variable, " - ", m1, ")^4)"))
  m4 <- svymean(formula_m4, design)
  m4 <- coef(m4)
  # Calculate kurtosis
  kurtosis <- m4 / (m2^2)
  return(kurtosis)
}

# =============================================================================

# Function to calculate survey-weighted correlations
calculate_weighted_correlation <- function(design, var1, var2) {
  # Create formula for correlation
  formula_cor <- as.formula(paste0("~", var1, "+", var2))
  # Calculate correlation
  cor_result <- svyvar(formula_cor, design)
  cor_matrix <- cov2cor(as.matrix(cor_result))
  return(cor_matrix[1, 2])
}

# =============================================================================

# Function to get skewness and kurtosis for continuous variables
get_skewness_kurtosis <- function(design, variables) {
  skewness_results <- numeric(length(variables))
  kurtosis_results <- numeric(length(variables))
  names(skewness_results) <- variables
  names(kurtosis_results) <- variables
  for (i in seq_along(variables)) {
    var <- variables[i]
    skewness_results[i] <- calculate_weighted_skewness(design, var)
    kurtosis_results[i] <- calculate_weighted_kurtosis(design, var)
  }
  return(list(
    skewness = skewness_results,
    kurtosis = kurtosis_results
  ))
}

# -----------------------------------------------------------------------------

# Function to get correlations for all variables
get_correlations <- function(design, variables) {
  n_vars <- length(variables)
  cor_matrix <- matrix(NA, nrow = n_vars, ncol = n_vars)
  rownames(cor_matrix) <- variables
  colnames(cor_matrix) <- variables
  for (i in 1:n_vars) {
    cor_matrix[i, i] <- 1
    for (j in 1:i) {
      if (i != j) {
        corr <- calculate_weighted_correlation(design, variables[i], variables[j])
        cor_matrix[i, j] <- corr
        cor_matrix[j, i] <- corr
      }
    }
  }
  return(cor_matrix)
}

# =============================================================================

# Main function to compute descriptive statistics with pooling
compute_descriptive_statistics <- function(imputation_results, countries) {
  results <- list()
  # Define continuous predictors
  continuous_vars <- c("ANXMAT", "MATHPERS", "MATHEFF", "DISCLIM",
                       "GROSAGR", "RELATST") # + MATHEF21 + IMMIG + ESCS) 
  # All variables for correlation
  all_vars <- c(continuous_vars)
  for (country in countries) {
    cat("Processing country:", country, "\n")
    # Initialize storage for pooling
    all_skewness <- list()
    all_kurtosis <- list()
    all_correlations <- list()
    all_proportions <- list()
    for (pv in 1:10) {
      cat("  Processing PV:", pv, "\n")
      for (imp in 1:10) {
        # Get imputed dataset
        data <- imputation_results[[country]][[pv]][[imp]]
        # Create survey design
        design <- create_survey_design(data, country)
        if (!is.null(design)) {
          # Calculate skewness and kurtosis for continuous variables
          sk_results <- get_skewness_kurtosis(design, continuous_vars)
          all_skewness[[length(all_skewness) + 1]] <- sk_results$skewness
          all_kurtosis[[length(all_kurtosis) + 1]] <- sk_results$kurtosis
          # Calculate correlations for all variables
          cor_results <- get_correlations(design, all_vars)
          all_correlations[[length(all_correlations) + 1]] <- cor_results
        }
      }
    }
    # For skewness
    skewness_array <- do.call(rbind, all_skewness)
    pooled_skewness <- colMeans(skewness_array, na.rm = TRUE)
    # For kurtosis
    kurtosis_array <- do.call(rbind, all_kurtosis)
    pooled_kurtosis <- colMeans(kurtosis_array, na.rm = TRUE)
    # For correlations (pooling each correlation separately)
    n_vars <- length(all_vars)
    pooled_correlations <- matrix(0, nrow = n_vars, ncol = n_vars)
    rownames(pooled_correlations) <- all_vars
    colnames(pooled_correlations) <- all_vars
    for (i in 1:n_vars) {
      pooled_correlations[i, i] <- 1
      for (j in 1:n_vars) {
        if (i != j) {
          # Extract all values for this specific correlation pair
          corr_values <- sapply(all_correlations, function(x) x[i, j])
          pooled_correlations[i, j] <- mean(corr_values, na.rm = TRUE)
        }
      }
    }
    # Store results for this country
    results[[country]] <- list(
      continuous_variables = list(
        skewness = pooled_skewness,
        kurtosis = pooled_kurtosis
      ),
      correlations = pooled_correlations
    )
  }
  return(results)
}

# -----------------------------------------------------------------------------

countries <- c("NOR", "SWE", "DNK", "FIN", "CHE")
desc_results <- compute_descriptive_statistics(imputation_results, countries)

# -----------------------------------------------------------------------------

# View results for a Norway

# View continuous variable statistics for Norway
desc_results$NOR$continuous_variables$skewness
desc_results$NOR$continuous_variables$kurtosis
# View correlation matrix for Norway
desc_results$NOR$correlations

# View results for a Sweden

# View continuous variable statistics for Sweden
desc_results$SWE$continuous_variables$skewness
desc_results$SWE$continuous_variables$kurtosis
# View correlation matrix for Sweden
desc_results$SWE$correlations

# View results for a Denmark

# View continuous variable statistics for Denmark
desc_results$DNK$continuous_variables$skewness
desc_results$DNK$continuous_variables$kurtosis
# View correlation matrix for Denmark
desc_results$DNK$correlations

# View results for a Finland

# View continuous variable statistics for Finland
desc_results$FIN$continuous_variables$skewness
desc_results$FIN$continuous_variables$kurtosis
# View correlation matrix for Finland
desc_results$FIN$correlations

# View results for a Switzerland

# View continuous variable statistics for Switzerland
desc_results$CHE$continuous_variables$skewness
desc_results$CHE$continuous_variables$kurtosis
# View correlation matrix for Switzerland
desc_results$CHE$correlations

# Sample sizes
# =============================================================================

# Load necessary libraries
library(dplyr)

# Define the dataset
data <- Master_data

# Define the variables
country_var <- "CNT"          # Country identifier
school_var <- "CNTSCHID"      # School identifier (PSU)
student_var <- "CNTSTUID"     # Student identifier

# Calculate the required statistics
summary_table <- data %>%
  group_by(!!sym(country_var), !!sym(school_var)) %>%  # Group by country and school
  summarise(
    n_students = n(),  # Count the number of students per school
    .groups = 'drop'   # Drop grouping after summarising
  ) %>%
  group_by(!!sym(country_var)) %>%  # Regroup by country
  summarise(
    n_schools = n(),  # Count the number of schools per country
    n_students_total = sum(n_students),  # Total students per country
    n_1_student = sum(n_students == 1),  # Schools with 1 student
    n_2_student = sum(n_students == 2),  # Schools with 2 student
    n_3_student = sum(n_students == 3)  # Schools with 3 student
  )

# Print the summary table
print(summary_table, width = Inf)

#######################################
#             THE END                 #
#######################################