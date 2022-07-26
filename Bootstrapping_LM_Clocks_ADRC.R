library(boot)
library(rsample)

library(readxl)
library(magrittr)
library(stringr)
library(tidyverse)


find_closest_age <- function(dataframe, grouped_df) {
  arrange(dataframe, abs(Age_at_Draw_Difference), ConnectivityZscore) %>% slice(1)
}

##############################################################################################
#################################-----Define Data Paths-----##################################
##############################################################################################


CHIP_Phen_Clocks_filepath <- "/Users/maurertm/Desktop/Chip_Phen_Clocks.csv"
proteome_filepath <- "/Users/maurertm/Desktop/dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021.csv"

##############################################################################################
###########################-----Import, Clean, and Filter Data-----###########################
##############################################################################################

CHIP_phen_clocks <- read_csv(CHIP_Phen_Clocks_filepath)
CHIP_phen_clocks_filtered <- CHIP_phen_clocks %>% filter(abs(Age_at_Draw_Difference) < 4) %>% filter(has_WGS_AND_clock_data)

CHIP_phen_clocks_closest <- CHIP_phen_clocks_filtered %>% group_by(Individual_ID) %>% group_modify(find_closest_age) 

##############################################################################################
###################################-----Resample Data-----####################################
##############################################################################################
boot <- CHIP_phen_clocks_closest %>% ungroup %>% bootstraps(1000) 

##############################################################################################
##################----Declare Function to Apply to Bootstrapped Data-----#####################
##############################################################################################
fit_lm_chip_class <- function(metadata_cases, metadata_df, metadata_controls) {
  metadata_controls_filter <- filter(metadata_controls, clock == metadata_df$clock)
  combined_df <- bind_rows(metadata_cases, metadata_controls_filter)
  combined_df$age_accel_resid <- lm(Estimated_Ages ~ Age_at_Prot_Draw, combined_df) %>% residuals
  #age_accel_lm <- lm(age_accel ~ Chronological_Age + Gender + has_chip, combined_df)
  age_accel_lm <- lm(age_accel_resid ~ Age_at_Prot_Draw + Gender + has_chip + WGS_Study + Study_Prot_Clocks + Storage + Storage_days, combined_df)
  age_accel_summary <- tidy(age_accel_lm, conf.int = T)
  age_accel_summary$n_chip <- nrow(metadata_cases)
  age_accel_summary$n_control <- nrow(metadata_controls_filter)
  age_accel_summary
}

run_lm <- function(dataframe){
  
  test_long <- as_tibble(dataframe) %>% pivot_longer(Adipose:Stomach, names_to = "clock", values_to = "Estimated_Ages")
  
  test_long_control <- filter(test_long, chip_class2 == "Control")
  test_long_cases <- filter(test_long, chip_class2 != "Control")
  
  chip_class_lm <- group_by(test_long_cases, clock, chip_class2) %>% 
    group_modify(fit_lm_chip_class, test_long_control) %>%
    filter(term != "(Intercept)")
  
  chip_class_lm
  
}

##############################################################################################
#########################-----Apply Function to Bootstrapped Data-----########################
##############################################################################################
limodels <- map(boot$splits, run_lm)

##############################################################################################
###################################-----Analyze Data-----#####################################
##############################################################################################
models <- bind_rows(limodels, .id = "Bootstrap") 
grouped_df_mean <- models %>%
  group_by(clock, chip_class2, term) %>% 
  summarise_at(vars("estimate"), mean) 

grouped_df_sd <- models %>%
  group_by(clock, chip_class2, term) %>% 
  summarise_at(vars("estimate"), sd) 

colnames(grouped_df_sd)[4] <- "sd_estimate"

final <- left_join(grouped_df_mean, grouped_df_sd)
final$ci_low <- final$estimate - (1.96*final$sd_estimate)
final$ci_high <- final$estimate + (1.96*final$sd_estimate)

final %>% filter(term == "has_chipTRUE") %>% filter(sign(ci_low) == sign(estimate) & sign(ci_high) == sign(estimate))

