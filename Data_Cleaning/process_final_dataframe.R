library(readxl)
library(magrittr)
library(stringr)
library(tidyverse)

##############################################################################################
#########################-----Helpful Processing Functions-----###############################
##############################################################################################
find_mutation_by_type <- function(dataframe, mutation) {
  dataframe_mutation <- filter(CHIP_phen_clocks, chip_class == mutation)
  dataframe_multiple <- filter(dataframe, chip_class == "Multiple")
  
  dataframe_multiple_mutation <- filter(dataframe_multiple, Hugo_Symbol == mutation | 
                                          Hugo_Symbol2 == mutation | Hugo_Symbol3 == mutation | Hugo_Symbol4 == mutation)
  
  with_multiple <- dataframe_multiple_mutation$filename
  with_single <- dataframe_mutation$filename
  
  with_mutation <- append(with_multiple, with_single)
  filtered_dataframe <- dataframe %>% filter(is_in(filename, with_mutation))
}

find_closest_age <- function(dataframe, grouped_df) {
  arrange(dataframe, abs(Age_at_Draw_Difference), ConnectivityZscore) %>% slice(1)
}

##############################################################################################
##########################-----Load and Filter Data-----######################################
##############################################################################################
CHIP_phen_clocks_filename <- "ADRC_Project/Processed_DataFrames/Chip_Phen_Clocks_12_19_at_16_14.csv"
CHIP_phen_clocks <- read_csv(CHIP_phen_clocks_filename)

CHIP_phen_clocks_filtered <- CHIP_phen_clocks %>% filter(abs(Age_at_Draw_Difference) < 4) %>% filter(has_WGS_AND_clock_data)
CHIP_phen_clocks_closest <- CHIP_phen_clocks_filtered %>% group_by(filename) %>% group_modify(find_closest_age) 

CHIP_phen_clocks_closest_select <- CHIP_phen_clocks_closest %>% select(filename, Age_at_WGS_Draw, Age_at_Prot_Draw, Age_at_Draw_Difference, Visit_Consensus, Gender, Diagnosis_group,
                                                                       Diagnosis_consensus, Storage_days, WGS_Study, PlateId, PlateRunDate,
                                                                       ScannerID, PlatePosition, SlideId, Subarray, PercentDilution, SampleMatrix,
                                                                       ConnectivityZscore, Fed.Fasted, Freeze.Thaw, Storage, Time.To.Decant, Time.To.Spin,
                                                                       Time.To.Ship, Time.To.Freeze, Study_Prot_Clocks, Adipose, Artery, Brain, CognitionBrain,
                                                                       Esophagus, Heart, Immune, Intestine, Kidney, Liver, Lung, Muscle, Organismal, Pancreas,
                                                                       Pituitary, Salivary, Stomach, has_clock_data, has_WGS_AND_clock_data, has_WGS_data,
                                                                       Largest_VAF, Multiple_VAF, Multiple_t_ref, Multiple_t_alt, contains("Hugo"), contains("Variant"),
                                                                       contains("Protein_Change"), contains("tumor_f"), contains("t_ref"), contains("t_alt"), chip_class,
                                                                       chip_class2, has_chip, vcf_FILTER, contains("AS_FilterStatus"), Chromosome, Start_Position, End_Position, Reference_Allele,
                                                                       contains("Tumor_Seq"), Transcript_Exon, Transcript_Position, cDNA_Change, Codon_Change, gc_content, DP, AS_SB_TABLE,
                                                                       ECNT, GERMQ, MPOS, POPAF, TLOD, OREGANNO_ID, OREGANNO_Values, Other_Transcripts, ref_context)
                                                                       
filepath <- 'ADRC_Project/Processed_DataFrames/Chip_Phen_Clocks_filtered.csv'
time <- paste0(sub('\\..*', '', filepath), format(Sys.time(),'_%m_%d_at_%H_%M'), '.csv')
write_csv(CHIP_phen_clocks_closest_select, time)                                                                    
                                                                       
                                                               
                                                                    
                                                                    



