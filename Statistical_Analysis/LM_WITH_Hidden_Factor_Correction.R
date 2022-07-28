#library(limma)
#library(uwot)
#library(WGCNA)
#library(Cairo)
#library(vsn)

library(RNOmni)
library(PCAtools)
library(broom)

library(readxl)
library(magrittr)
library(stringr)
library(tidyverse)

find_closest_age <- function(dataframe, grouped_df) {
  arrange(dataframe, abs(Age_at_Draw_Difference), ConnectivityZscore) %>% slice(1)
}

chip_class_association <- function(data_rows, data_df) {
  eigenvector_lm <- lm(Eigenvector ~ chip_class2, data_rows)
  lm_summary <- anova(eigenvector_lm) %>% tidy
  lm_summary
}

chip_only_association <- function(data_rows, data_df) {
  eigenvector_lm <- lm(Eigenvector ~ has_chip, data_rows)
  lm_summary <- anova(eigenvector_lm) %>% tidy
  lm_summary
}

clock_association <- function(data_rows, data_df) {
  data_rows$age_accel_resid <- lm(Estimated_Ages ~ Age_at_Prot_Draw, data_rows) %>% residuals
  eigenvector_lm <- lm(Eigenvector ~ age_accel_resid, data_rows)
  lm_summary <- anova(eigenvector_lm) %>% tidy
  lm_summary
}

get_PCs <- function(grouped_rows, grouped_df){
  grouped_rows$PC
}

##############################################################################################
#################################-----Define Data Paths-----##################################
##############################################################################################

CHIP_Phen_Clocks_filepath <- "/Users/maurertm/Desktop/Chip_Phen_Clocks.csv"
proteome_filepath <- "/Users/maurertm/Desktop/dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021.csv"

##############################################################################################
###########################-----Import, Clean, and Subset Data-----###########################
##############################################################################################
CHIP_phen_clocks <- read_csv(CHIP_Phen_Clocks_filepath)
CHIP_phen_clocks_filtered <- CHIP_phen_clocks %>% filter(abs(Age_at_Draw_Difference) < 4) %>% filter(has_WGS_AND_clock_data)

CHIP_phen_clocks_closest <- CHIP_phen_clocks_filtered %>% group_by(Individual_ID) %>% group_modify(find_closest_age) 

#Convert Data to Long Form
metadata_CHIP_long <- pivot_longer(CHIP_phen_clocks_closest, Adipose:Stomach, names_to = "clock", values_to = "Estimated_Ages")

#Compute age_accel column
metadata_CHIP_long$age_accel <- metadata_CHIP_long$Age_at_Prot_Draw - metadata_CHIP_long$Estimated_Ages

#import and transpose the proteomics dataframe 
proteome <- read_csv(proteome_filepath)
proteome_mat <- select(proteome, -`...1`, -SampleType, -SampleId) %>% as.matrix %>% log2 %>% t 
colnames(proteome_mat) <- proteome$SampleId

#subset the proteomics dataframe to only contain values sharing a SampleId with metadata_youngest
proteome_subset <- proteome_mat[,CHIP_phen_clocks_closest$SampleId]

##############################################################################################
##################-----Discover and Remove Significant Hidden Factors-----####################
###############################---Associated with CHIP---#####################################
##############################################################################################

######################---Apply PCA to find PCs associated with CHIP---########################
vst_scaled <- apply(proteome_subset, 1, RankNorm) %>% t
pca_scaled <- pca(vst_scaled)
var_scaled <- as.vector(vst_scaled) %>% var
vst_npc <- chooseGavishDonoho(vst_scaled, var.explained = pca_scaled$sdev^2, noise = var_scaled)
vst_pcs <- pca_scaled$rotated[,1:vst_npc]
vst_pcs_scaled <- apply(vst_pcs, 2, RankNorm) %>% as_tibble(rownames = "SampleId")
CHIP_phen_clocks_closest_pcs <- left_join(CHIP_phen_clocks_closest, vst_pcs_scaled)
vst_pcs_scaled_long <- pivot_longer(CHIP_phen_clocks_closest_pcs, contains("PC"), names_to = "PC", values_to = "Eigenvector")

vst_pcs_lm_chip_class <- group_by(vst_pcs_scaled_long, PC) %>% group_modify(chip_class_association) %>% filter(term != "Residuals")
vst_pcs_lm_chip_class$p_adj <- p.adjust(vst_pcs_lm_chip_class$p.value, method = "bonferroni")
#remove p10

vst_pcs_lm_chip <- group_by(vst_pcs_scaled_long, PC) %>% group_modify(chip_only_association) %>% filter(term != "Residuals")
vst_pcs_lm_chip$p_adj <- p.adjust(vst_pcs_lm_chip$p.value, method = "bonferroni")
#keep p 1-20

######################---Add PC weights to CHIP_phen_clocks_closest---########################
CHIP_phen_clocks_closest_pcs <- left_join(CHIP_phen_clocks_closest, vst_pcs_scaled)
CHIP_phen_clocks_closest_pcs <- CHIP_phen_clocks_closest_pcs %>% ungroup() %>%
  select(Age_at_Prot_Draw, Gender, WGS_Study, Study_Prot_Clocks, Storage, Storage_days, PC1:PC21, Adipose:Stomach, chip_class2, has_chip)

##############################################################################################
##################-----Discover and Remove Significant Hidden Factors-----####################
###########################---Associated with Proteomic Clocks---#############################
##############################################################################################

#######################---Test for collinearity of clocks with PCs---#########################
vst_pcs_scaled_clocks_long <- pivot_longer(CHIP_phen_clocks_closest_pcs, contains("PC"), names_to = "PC", values_to = "Eigenvector") %>%
  pivot_longer(Adipose:Stomach, names_to = "clock", values_to = "Estimated_Ages")

vst_pcs_lm_clocks <- group_by(vst_pcs_scaled_clocks_long, PC, clock) %>% group_modify(clock_association) %>% filter(term != "Residuals")
vst_pcs_lm_clocks$p_adj <- p.adjust(vst_pcs_lm_clocks$p.value, method = "bonferroni")

####################---Find PCs associated with the proteomic clocks---#######################
##Get a list containing clock names and, within the named list, a list of PCs significantly associated with the proteomic clocks
list_of_clock_adj_pcs <- vst_pcs_lm_clocks %>% filter(p_adj < 0.05) %>% select(PC, clock) %>% group_by(clock) %>% group_map(get_PCs) 
names(list_of_clock_adj_pcs) <- unique(vst_pcs_lm_clocks$clock)

list_of_clock_unadj_pcs <- vst_pcs_lm_clocks %>% filter(p.value < 0.05) %>% select(PC, clock) %>% group_by(clock) %>% group_map(get_PCs) 
names(list_of_clock_unadj_pcs) <- unique(vst_pcs_lm_clocks$clock)

#########################-----Convert to Long, Controls, Cases-----###########################

CHIP_phen_clocks_closest_pcs_long <- pivot_longer(CHIP_phen_clocks_closest_pcs, Adipose:Stomach, names_to = "clock", values_to = "Estimated_Ages")

CHIP_phen_clocks_closest_pcs_controls <- filter(CHIP_phen_clocks_closest_pcs_long, chip_class2 == "Control")
CHIP_phen_clocks_closest_pcs_cases <- filter(CHIP_phen_clocks_closest_pcs_long, chip_class2 != "Control")

##############################################################################################
################################-----Create LM Models-----####################################
##############################################################################################

###########################-----CHIP as a Whole vs Clocks-----################################
#using age residuals
fit_whole_lm_chip_resid <- function(metadata_rows, metadata_df, clock_pc_list) {
  pcs_to_remove <- clock_pc_list[[metadata_df$clock]]
  metadata_rows$age_accel_resid <- lm(Estimated_Ages ~ Age_at_Prot_Draw, metadata_rows) %>% residuals %>% RankNorm
  metadata_rows_select <- select(metadata_rows, -Estimated_Ages, -chip_class2) %>% select(-one_of(pcs_to_remove))
  age_accel_lm <- lm(age_accel_resid ~ ., metadata_rows_select)
  age_accel_summary <- tidy(age_accel_lm, conf.int = T)
  age_accel_summary
}

#using age difference
fit_whole_lm_chip_diff <- function(metadata_rows, metadata_df, clock_pc_list) {
  pcs_to_remove <- clock_pc_list[[metadata_df$clock]]
  metadata_rows$age_accel <- metadata_rows$Estimated_Ages - metadata_rows$Age_at_Prot_Draw
  metadata_rows_select <- select(metadata_rows, -Estimated_Ages, -chip_class2) %>% select(-one_of(pcs_to_remove))
  age_accel_lm <- lm(age_accel ~ ., metadata_rows_select)
  age_accel_summary <- tidy(age_accel_lm, conf.int = T)
  age_accel_summary
}

############################-----CHIP Classes vs Clocks-----##################################
#using age residuals
fit_lm_chip_class_resid <- function(metadata_cases, metadata_df, metadata_controls, clock_pc_list) {
  pcs_to_remove <- clock_pc_list[[metadata_df$clock]]
  metadata_controls_filter <- filter(metadata_controls, clock == metadata_df$clock) %>%
    select(-chip_class2, -clock)
  combined_df <- bind_rows(metadata_cases, metadata_controls_filter)
  combined_df$age_accel_resid <- lm(Estimated_Ages ~ Age_at_Prot_Draw, combined_df) %>% residuals %>% RankNorm
  combined_df_select <- select(combined_df, -Estimated_Ages) %>% select(-one_of(pcs_to_remove))
  age_accel_lm <- lm(age_accel_resid ~ ., combined_df_select) # this removes one assc with
  age_accel_summary <- tidy(age_accel_lm, conf.int = T)
  age_accel_summary$n_chip <- nrow(metadata_cases)
  age_accel_summary$n_control <- nrow(metadata_controls_filter)
  age_accel_summary
}

#using age difference
fit_lm_chip_class_difference <- function(metadata_cases, metadata_df, metadata_controls, clock_pc_list) {
  pcs_to_remove <- clock_pc_list[[metadata_df$clock]]
  metadata_controls_filter <- filter(metadata_controls, clock == metadata_df$clock) %>%
    select(-chip_class2, -clock)
  combined_df <- bind_rows(metadata_cases, metadata_controls_filter)
  combined_df$age_accel <- combined_df$Estimated_Ages - combined_df$Age_at_Prot_Draw
  combined_df_select <- select(combined_df, -Estimated_Ages) %>% select(-one_of(pcs_to_remove))
  age_accel_lm <- lm(age_accel ~ ., combined_df_select)
  age_accel_summary <- tidy(age_accel_lm, conf.int = T)
  age_accel_summary$n_chip <- nrow(metadata_cases)
  age_accel_summary$n_control <- nrow(metadata_controls_filter)
  age_accel_summary
}


##############################################################################################
##########################-----Analysis with Unadj Clocks-----################################
##############################################################################################

###########################-----CHIP as a Whole vs Clocks-----################################
#using age residuals
chip_clocks_lm <- group_by(CHIP_phen_clocks_closest_pcs_long, clock) %>% group_modify(fit_whole_lm_chip_resid, list_of_clock_unadj_pcs) 
has_chip_lm <- filter(chip_clocks_lm, term == "has_chipTRUE") %>% arrange(p.value)
has_chip_lm$p_adjust <- p.adjust(has_chip_lm$p.value, method = "fdr")
chip_whole_unajd_resid <- has_chip_lm %>% arrange(p_adjust)
chip_whole_unajd_resid

#using age difference
chip_clocks_lm <- group_by(CHIP_phen_clocks_closest_pcs_long, clock) %>% group_modify(fit_whole_lm_chip_diff, list_of_clock_unadj_pcs) 
has_chip_lm <- filter(chip_clocks_lm, term == "has_chipTRUE") %>% arrange(p.value)
has_chip_lm$p_adjust <- p.adjust(has_chip_lm$p.value, method = "fdr")
chip_whole_unajd_diff <- has_chip_lm %>% arrange(p_adjust)
chip_whole_unajd_diff

############################-----CHIP Classes vs Clocks-----##################################
#using age residuals

#list_of_clock_unadj_pcs can be swapped for list_of_clock_adj_pcs to only remove PCs
#that are significant when the p.value is adjusted
chip_class_lm <- group_by(CHIP_phen_clocks_closest_pcs_cases, clock, chip_class2) %>% 
  group_modify(fit_lm_chip_class_resid, CHIP_phen_clocks_closest_pcs_controls, list_of_clock_unadj_pcs) %>%
  filter(term != "(Intercept)")
filter(chip_class_lm, term == "has_chipTRUE") %>% arrange(p.value)
chip_class_lm_filter <- filter(chip_class_lm, term == "has_chipTRUE")
chip_class_lm_filter$p_adjust <- p.adjust(chip_class_lm_filter$p.value, method = "fdr")
chip_classes_unadj_age_resid <- chip_class_lm_filter %>% arrange(p_adjust)
chip_classes_unadj_age_resid


#using age difference

#list_of_clock_unadj_pcs can be swapped for list_of_clock_adj_pcs to only remove PCs
#that are significant when the p.value is adjusted
chip_class_lm <- group_by(CHIP_phen_clocks_closest_pcs_cases, clock, chip_class2) %>% 
  group_modify(fit_lm_chip_class_difference, CHIP_phen_clocks_closest_pcs_controls, list_of_clock_unadj_pcs) %>%
  filter(term != "(Intercept)")
filter(chip_class_lm, term == "has_chipTRUE") %>% arrange(p.value)
chip_class_lm_filter <- filter(chip_class_lm, term == "has_chipTRUE")
chip_class_lm_filter$p_adjust <- p.adjust(chip_class_lm_filter$p.value, method = "fdr")
chip_classes_unadj_age_diff <- chip_class_lm_filter %>% arrange(p_adjust)
chip_classes_unadj_age_diff

##############################################################################################
############################-----Analysis with Adj Clocks-----################################
##############################################################################################

###########################-----CHIP as a Whole vs Clocks-----################################
#using age residuals
chip_clocks_lm <- group_by(CHIP_phen_clocks_closest_pcs_long, clock) %>% group_modify(fit_whole_lm_chip_resid, list_of_clock_adj_pcs) 
has_chip_lm <- filter(chip_clocks_lm, term == "has_chipTRUE") %>% arrange(p.value)
has_chip_lm$p_adjust <- p.adjust(has_chip_lm$p.value, method = "fdr")
chip_whole_unajd_resid <- has_chip_lm %>% arrange(p_adjust)
chip_whole_unajd_resid

#using age difference
chip_clocks_lm <- group_by(CHIP_phen_clocks_closest_pcs_long, clock) %>% group_modify(fit_whole_lm_chip_diff, list_of_clock_adj_pcs) 
has_chip_lm <- filter(chip_clocks_lm, term == "has_chipTRUE") %>% arrange(p.value)
has_chip_lm$p_adjust <- p.adjust(has_chip_lm$p.value, method = "fdr")
chip_whole_ajd <- has_chip_lm %>% arrange(p_adjust)
chip_whole_ajd

############################-----CHIP Classes vs Clocks-----##################################
#using age residuals

#list_of_clock_unadj_pcs can be swapped for list_of_clock_adj_pcs to only remove PCs
#that are significant when the p.value is adjusted
chip_class_lm <- group_by(CHIP_phen_clocks_closest_pcs_cases, clock, chip_class2) %>% 
  group_modify(fit_lm_chip_class_resid, CHIP_phen_clocks_closest_pcs_controls, list_of_clock_adj_pcs) %>%
  filter(term != "(Intercept)")
filter(chip_class_lm, term == "has_chipTRUE") %>% arrange(p.value)
chip_class_lm_filter <- filter(chip_class_lm, term == "has_chipTRUE")
chip_class_lm_filter$p_adjust <- p.adjust(chip_class_lm_filter$p.value, method = "fdr")
chip_classes_adj_age_resid <- chip_class_lm_filter %>% arrange(p_adjust)
chip_classes_adj_age_resid


#using age difference

#list_of_clock_unadj_pcs can be swapped for list_of_clock_adj_pcs to only remove PCs
#that are significant when the p.value is adjusted
chip_class_lm <- group_by(CHIP_phen_clocks_closest_pcs_cases, clock, chip_class2) %>% 
  group_modify(fit_lm_chip_class_difference, CHIP_phen_clocks_closest_pcs_controls, list_of_clock_adj_pcs) %>%
  filter(term != "(Intercept)")
filter(chip_class_lm, term == "has_chipTRUE") %>% arrange(p.value)
chip_class_lm_filter <- filter(chip_class_lm, term == "has_chipTRUE")
chip_class_lm_filter$p_adjust <- p.adjust(chip_class_lm_filter$p.value, method = "fdr")
chip_classes_adj_age_diff <- chip_class_lm_filter %>% arrange(p_adjust)
chip_classes_adj_age_diff

