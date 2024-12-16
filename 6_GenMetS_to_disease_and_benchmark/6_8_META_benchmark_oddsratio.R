#################################################
# Pan Hong
# Sep 30, 2024
################################################

rm(list=ls())
library(readxl)
library(dplyr)
library(data.table)
library(meta)

OR_by_sex <- read_xlsx("../../0_data/OddsRatio_result/OddsRatio_benchmark.xlsx")
OR_by_sex <- OR_by_sex %>% rename(N=case, OR = "OddsRatio", Lower_CI=`95%CI lower`, Upper_CI=`95%CI upper`, 
                                  Pvalue=`p-value`, Cohort = "cohort", Gender=Sex)
OR_by_sex$N <- OR_by_sex$N * 2 
OR_by_sex$Disease_and_predictor <- paste0(OR_by_sex$Disease, ";", OR_by_sex$predictor)
idx <- which(OR_by_sex$Disease == "Heart Failure")
OR_by_sex$Disease[idx]<- "Heart failure"

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
  
  return(data.frame(Combined_OR = combined_OR,
                    Lower_CI = lower_CI,
                    Upper_CI = upper_CI,
                    P_value = p_value,
                    Combined_N = sum(disease_data$N )))        # Add combined sample size
}

# Conduct the meta-analysis on SEX (men+women) for each disease, cohort and predictor
all_predictors <- unique(OR_by_sex$predictor)
all_diseases <- unique(OR_by_sex$Disease)
all_cohorts <- unique(OR_by_sex$Cohort)
results_table <- data.frame()

for(pred in all_predictors){
for (coh in all_cohorts){
  for (dis in all_diseases ) {
   #pred=all_predictor[1]; coh=all_cohorts[1]; dis=all_diseases[1];
   #cat("working disease=", dis, " cohort=", coh, "predictor=", pred, "\n")
  
  disease_data <- OR_by_sex %>% dplyr::filter(Disease==dis) %>% 
    dplyr::filter(Cohort==coh) %>%
    dplyr::filter(predictor==pred)
  
  if(dim(disease_data)[1] > 0 ){
    result <- run_meta_analysis(disease_data)
    result$Disease <- dis  
    result$Cohort <- coh    
    result$predictor <- pred   
    results_table <- rbind(results_table, result)
  }
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
  
  
  output1 <- data.frame(Combined_OR = combined_OR,
                        Lower_CI = lower_CI,
                        Upper_CI = upper_CI,
                        P_value = p_value,
                        Combined_N = sum(disease_data$Combined_N))
  
  output2 <- random_weights_percent       
  
  return(list(result_summary = output1, weights_percent = output2))
}

# Conduct the meta-analysis for cohorts each disease and store the results in a table
 
results_table2 <- data.frame()
for(pred in all_predictors){
for (dis in all_diseases ) {
   
  disease_data <- results_table %>% dplyr::filter(Disease==dis) %>%
    dplyr::filter(predictor == pred)
  
  if(dim(disease_data)[1] > 0 ){
  results <- run_meta_analysis_cross_cohort(disease_data)
  
  meta_result <- results$result_summary;
  meta_result$Disease <- dis;
  meta_result$Cohort <- "META";
  meta_result$predictor <- pred; 
  
  weight_result <- round(results$weights_percent,2); print(weight_result)
  
  result <- rbind(disease_data, meta_result)
  weights_result <- c(round(results$weights_percent,2), 100)
  
  result$Weight <- weights_result
  results_table2 <- rbind(results_table2, result)
  }
}
}

print(results_table2)

write.table(results_table2, "result_META_benchmark.csv", row.names = FALSE)
