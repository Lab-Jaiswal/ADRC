library(lme4)
library(lmerTest)

library(writexl)
library(readxl)
library(openxlsx)
library(magrittr)
library(tidyverse)

##############################################################################################
#################################-----Define Data Paths-----##################################
##############################################################################################
CHIP_Phen_Clocks_filepath <- "/Users/maurertm/Desktop/Chip_Phen_Clocks.csv"

##############################################################################################
###########################-----Import, Clean, and Subset Data-----###########################
##############################################################################################
CHIP_phen_clocks <- read_csv(CHIP_Phen_Clocks_filepath)
CHIP_phen_clocks_filtered <- CHIP_phen_clocks %>% filter(abs(Age_at_Draw_Difference) < 4) %>% filter(has_WGS_AND_clock_data)

#Convert Data to Long Form
metadata_CHIP_long <- pivot_longer(CHIP_phen_clocks_filtered, Adipose:Stomach, names_to = "clock", values_to = "Estimated_Ages")

#Compute age_accel column
metadata_CHIP_long$age_accel <- metadata_CHIP_long$Age_at_Prot_Draw - metadata_CHIP_long$Estimated_Ages

#define control and cases
metadata_CHIP_long_control <- filter(metadata_CHIP_long, chip_class2 == "Control")
metadata_CHIP_long_cases <- filter(metadata_CHIP_long, chip_class2 != "Control")

##############################################################################################
##################-----Discover and Remove Significant Hidden Factors-----####################
##############################################################################################
fit_lmer_chip_class <- function(metadata_cases, metadata_df, metadata_controls) {
  metadata_controls_filter <- filter(metadata_controls, clock == metadata_df$clock) %>%
    select(-chip_class2, -clock)
  combined_df <- bind_rows(metadata_cases, metadata_controls_filter)
  combined_df$age_accel_resid <- lm(Estimated_Ages ~ Age_at_Prot_Draw, combined_df) %>% residuals %>% RankNorm
  combined_df_select <- select(combined_df, -Estimated_Ages)
  age_accel_lmer <- lmer(age_accel_resid ~ (1 | Individual_ID) + Age_at_Prot_Draw + Gender + has_chip, combined_df_select)
  summary(age_accel_lmer) %>% pluck("coefficients") %>% as_tibble(rownames = "term")
}

chip_class_lmer <- group_by(metadata_CHIP_long_cases, clock, chip_class2) %>% 
  group_modify(fit_lmer_chip_class, metadata_CHIP_long_control) %>%
  filter(term != "(Intercept)")

chip_lmer_only <- filter(chip_class_lmer, term == "has_chipTRUE") %>% arrange(`Pr(>|t|)`)
chip_lmer_only$p_adjust <- p.adjust(chip_lmer_only$`Pr(>|t|)`, method="fdr")
chip_lmer_only %>% arrange(p_adjust)
