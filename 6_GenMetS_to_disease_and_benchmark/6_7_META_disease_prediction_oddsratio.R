#################################################
# Pan Hong
# Sep 30, 2024
################################################

rm(list=ls())
library(readxl)
library(dplyr)
library(data.table)
library(meta)

OR_by_sex <- read_xlsx("../../0_data/OddsRatio_result/OR_by_sex.xlsx")
OR_by_sex <- OR_by_sex %>% rename(N=Case, Lower_CI=CI_Low, Upper_CI=CI_upper)
OR_by_sex$N <- OR_by_sex$N * 2 

# Define a function to conduct meta-analysis for each disease, incorporating sample size
run_meta_analysis <- function(disease_data) {
  # Calculate weights using inverse variance and sample size
  # Calculate the log odds ratios and their standard errors
  log_or <- log(disease_data$OR)
  se_log_or <- (log(disease_data$Upper_CI) - log(disease_data$Lower_CI)) / 3.92
  
  # Weighting factor: inverse variance (1 / variance) scaled by sample size (N)
  # sample size = case + control, and case=control
  weights <- disease_data$N / se_log_or^2
  
  # Meta-analysis using the metagen function, adding weights for sample sizes
  meta_analysis <- metagen(
    TE = log_or,                               # log odds ratio
    seTE = se_log_or,                          # standard error of log OR
    studlab = disease_data$Gender,     # study labels (men/women, cohort)
    sm = "OR",                                 # odds ratio scale
    weights = weights                          # adding weights based on sample size
  )
  
  # Extract OR, 95% CI, and p-value for the combined result
  combined_OR <- exp(meta_analysis$TE.random)                # Exponentiate to get OR
  lower_CI <- exp(meta_analysis$lower.random)                # Lower bound of 95% CI
  upper_CI <- exp(meta_analysis$upper.random)                # Upper bound of 95% CI
  p_value <- meta_analysis$pval.random                       # P-value
  
  return(data.frame(Disease = disease_data$Disease[1],
                    Combined_OR = combined_OR,
                    Lower_CI = lower_CI,
                    Upper_CI = upper_CI,
                    P_value = p_value,
                    Combined_N = sum(disease_data$N )))        # Add combined sample size
}

# Conduct the meta-analysis on SEX (men+women) for each disease for each cohort.
all_diseases <- unique(OR_by_sex$Disease)
all_cohorts <- unique(OR_by_sex$Cohort)
results_table <- data.frame()
for (cohort in all_cohorts){
  for (disease in all_diseases ) {
  cat("working disease=", disease, " cohort=", cohort, "\n")
  
  disease_data <- OR_by_sex %>% dplyr::filter(Disease==disease) %>% 
    dplyr::filter(Cohort==cohort)
  
  if(dim(disease_data)[1] > 0 ){
    result <- run_meta_analysis(disease_data)
    result$Cohort <- cohort
    # Append to the results table
    results_table <- rbind(results_table, result)
  }
}
}
print(results_table)


################META cross cohort
# Define a function to conduct meta-analysis for each disease, incorporating sample size
run_meta_analysis_cross_cohort <- function(disease_data) {
  # Calculate weights using inverse variance and sample size
  # Calculate the log odds ratios and their standard errors
  log_or <- log(disease_data$Combined_OR)
  se_log_or <- (log(disease_data$Upper_CI) - log(disease_data$Lower_CI)) / 3.92
  
  # Weighting factor: inverse variance (1 / variance) scaled by sample size (N)
  # sample size = case + control, and case=control
  weights <- disease_data$Combined_N / se_log_or^2
  
  # Meta-analysis using the metagen function, adding weights for sample sizes
  meta_analysis <- metagen(
    TE = log_or,                               # log odds ratio
    seTE = se_log_or,                          # standard error of log OR
    studlab = disease_data$Cohort,             # study labels (men/women, cohort)
    sm = "OR",                                 # odds ratio scale
    weights = weights                          # adding weights based on sample size
  )
  
  # Extract OR, 95% CI, and p-value for the combined result
  combined_OR <- exp(meta_analysis$TE.random)                # Exponentiate to get OR
  lower_CI <- exp(meta_analysis$lower.random)                # Lower bound of 95% CI
  upper_CI <- exp(meta_analysis$upper.random)                # Upper bound of 95% CI
  p_value <- meta_analysis$pval.random                       # P-value
  
  random_weights <- meta_analysis$w.random
  random_weights_percent <- (random_weights / sum(random_weights)) * 100
  
  
  output1 <- data.frame(Disease = disease_data$Disease[1],
                        Combined_OR = combined_OR,
                        Lower_CI = lower_CI,
                        Upper_CI = upper_CI,
                        P_value = p_value,
                        Combined_N = sum(disease_data$Combined_N),
                        Cohort= "META")
  
  output2 <- random_weights_percent       
  
  return(list(result_summary = output1, weights_percent = output2))
}

# Conduct the meta-analysis for cohorts each disease and store the results in a table
all_diseases <- unique(results_table$Disease)
results_table2 <- data.frame()

for (disease in all_diseases ) {
  cat("working disease=", disease, "\n")
  
  disease_data <- results_table %>% dplyr::filter(Disease==disease)
  results <- run_meta_analysis_cross_cohort(disease_data)
  
  meta_result <- results$result_summary; print(meta_result)
  weight_result <- round(results$weights_percent,2); print(weight_result)
  
  if(disease == "HTN" ){ #missing data in BBJ
    bbj_htn <- data.table(
      Disease = "HTN", 
      Combined_OR=NA,
      Lower_CI=NA,
      Upper_CI=NA,
      P_value = NA,
      Combined_N = NA,
      Cohort = "BBJ")
    
    result <- rbind(disease_data, bbj_htn, meta_result)
    
    # Access the second output (random-effects weights as percentages)
    weights_result <- c(round(results$weights_percent,2), NA, 100)
  }else{
    result <- rbind(disease_data, meta_result)
    weights_result <- c(round(results$weights_percent,2), 100)
  }
  # Append to the results table
  result$Weight <- weights_result
  results_table2 <- rbind(results_table2, result)
}

print(results_table2)

write.table(results_table2, "result_META_disease_prediction.csv", row.names = FALSE)
