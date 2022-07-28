library(readxl)
library(magrittr)
library(stringr)
library(tidyverse)

library(broom)


find_closest_age <- function(dataframe, grouped_df) {
  arrange(dataframe, abs(Age_at_Draw_Difference), ConnectivityZscore) %>% slice(1)
}

##############################################################################################
#################################-----Define Data Paths-----##################################
##############################################################################################

CHIP_Phen_Clocks_filepath <- "/Users/maurertm/Desktop/Chip_Phen_Clocks.csv"

##############################################################################################
###########################-----Import, Clean, and Filter Data-----###########################
##############################################################################################

CHIP_phen_clocks <- read_csv(CHIP_Phen_Clocks_filepath)
CHIP_phen_clocks_filtered <- CHIP_phen_clocks %>% filter(abs(Age_at_Draw_Difference) < 4) %>% filter(has_WGS_AND_clock_data)

CHIP_phen_clocks_closest <- CHIP_phen_clocks_filtered %>% group_by(Individual_ID) %>% group_modify(find_closest_age) 

##############################################################################################
################################-----Prep Metadata-----#######################################
##############################################################################################
metadata_CHIP_long <- pivot_longer(CHIP_phen_clocks_closest, Adipose:Stomach, names_to = "clock", values_to = "Estimated_Ages")
metadata_CHIP_long$age_accel <- metadata_CHIP_long$Age_at_Prot_Draw - metadata_CHIP_long$Estimated_Ages

##############################################################################################
##################################-----CHIP vs Age-----#######################################
##############################################################################################
chip_age_lm <- lm(Age_at_Prot_Draw ~ has_chip, CHIP_phen_clocks_closest) %>% tidy

##############################################################################################
###########################-----CHIP as a Whole vs Clocks-----################################
##############################################################################################

fit_lm_chip <- function(metadata_rows, metadata_df) {
  age_accel_lm <- lm(age_accel ~ Age_at_Prot_Draw + Gender + has_chip + WGS_Study + Study_Prot_Clocks + Storage + Storage_days, metadata_rows)
  age_accel_summary <- tidy(age_accel_lm, conf.int = T)
  age_accel_summary
}

has_chip_lm <- group_by(metadata_CHIP_long, clock) %>% group_modify(fit_lm_chip) %>% filter(term == "has_chipTRUE")
filter(has_chip_lm, term == "has_chipTRUE") %>% arrange(p.value)

##############################################################################################
############################-----CHIP Classes vs Clocks-----##################################
##############################################################################################

#############################---Compute Controls, Cases---####################################

metadata_CHIP_controls <- filter(metadata_CHIP_long, chip_class2 == "Control")
metadata_CHIP_cases <- filter(metadata_CHIP_long, chip_class2 != "Control")

#################################---Using Age Residuals---####################################
fit_lm_chip_class_residuals <- function(metadata_cases, metadata_df, metadata_controls) {
  metadata_controls_filter <- filter(metadata_controls, clock == metadata_df$clock)
  combined_df <- bind_rows(metadata_cases, metadata_controls_filter)
  combined_df$age_accel_resid <- lm(Estimated_Ages ~ Age_at_Prot_Draw, combined_df) %>% residuals %>% RankNorm()
  age_accel_lm <- lm(age_accel_resid ~ Age_at_Prot_Draw + Gender + has_chip + WGS_Study + Study_Prot_Clocks + Storage + Storage_days, combined_df)
  age_accel_summary <- tidy(age_accel_lm, conf.int = T)
  age_accel_summary$n_chip <- nrow(metadata_cases)
  age_accel_summary$n_control <- nrow(metadata_controls_filter)
  age_accel_summary
}

chip_class_lm <- group_by(metadata_CHIP_cases, clock, chip_class2) %>% group_modify(fit_lm_chip_class_residuals, metadata_CHIP_controls) %>%
  filter(term != "(Intercept)")
filter(chip_class_lm, term == "has_chipTRUE") %>% arrange(p.value)
chip_class_lm_filter <- filter(chip_class_lm, term == "has_chipTRUE")
chip_class_lm_filter$p_adjust <- p.adjust(chip_class_lm_filter$p.value, method = "fdr")
chip_class_lm_filter %>% arrange(p_adjust)
chip_class_lm_filter %>% arrange(p.value)

#################################---Using Age Difference---##################################
fit_lm_chip_class_difference <- function(metadata_cases, metadata_df, metadata_controls) {
  metadata_controls_filter <- filter(metadata_controls, clock == metadata_df$clock)
  combined_df <- bind_rows(metadata_cases, metadata_controls_filter)
  age_accel_lm <- lm(age_accel ~ Age_at_Prot_Draw + Gender + has_chip + WGS_Study + Study_Prot_Clocks + Storage + Storage_days, combined_df)
  age_accel_summary <- tidy(age_accel_lm, conf.int = T)
  age_accel_summary$n_chip <- nrow(metadata_cases)
  age_accel_summary$n_control <- nrow(metadata_controls_filter)
  age_accel_summary
}

chip_class_lm <- group_by(metadata_CHIP_cases, clock, chip_class2) %>% group_modify(fit_lm_chip_class_difference, metadata_CHIP_controls) %>%
  filter(term != "(Intercept)")
filter(chip_class_lm, term == "has_chipTRUE") %>% arrange(p.value)
chip_class_lm_filter <- filter(chip_class_lm, term == "has_chipTRUE")
chip_class_lm_filter$p_adjust <- p.adjust(chip_class_lm_filter$p.value, method = "fdr")
chip_class_lm_filter %>% arrange(p_adjust)
chip_class_lm_filter %>% arrange(p.value)
