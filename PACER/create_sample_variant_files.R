library(readxl)
library(magrittr)
library(stringr)
library(tidyverse)

set_colnames <- function(df_given, column_names){
  colnames(df_given) <- column_names
  df_given
}

split_by_gene <- function(mutation_number){
  if (mutation_number == 1) {
    number=""
  } else {
    number=mutation_number
  }
  
  dataframe1 <- dataframe %>% select(filename)
  hugo_column <- paste("Hugo_Symbol", number, sep="")
  t_ref_column <- paste("t_ref_count", number, sep ="")
  t_alt_column <- paste("t_alt_count", number, sep="")
  tumor_f_column <- paste("tumor_f", number, sep="")
  
  dataframe1 <- dataframe %>% select(filename, t_ref_column, t_alt_column, hugo_column, tumor_f_column, WGS_Study, CHIP, Gender, vcf_sample_header) 
  dataframe1$AD <- do.call(paste, c(dataframe1[t_ref_column], ",", Multiple_carriers[t_alt_column], sep = ""))
  dataframe1$Mutation_number <- mutation_number
  dataframe1 <- dataframe1 %>% select(filename, hugo_column, AD, tumor_f_column, CHIP, WGS_Study, Gender, Mutation_number, vcf_sample_header) %>% 
    set_colnames(c("filename", "Gene", "AD", "VAF", "haschip", "STUDY", "Gender", "Mutation_number", "vcf_sample_header"))
  dataframe1
} 

combined_split_genes <- function(dataframe, number){
  dataframe <- dataframe
  range_number <- 1:4
  multiple_list <- map(range_number, split_by_gene)
  multiple_list
}
  

##############################################################################################
################################-----Declare Filepaths-----###################################
##############################################################################################
metadata_filepath <- "/Users/maurertm/Desktop/Chip_Phen_Clocks.csv"
PACER_vcf_list_filepath <- "/Users/maurertm/Desktop/PACER_vcfs.tsv"
ADRC_key_filepath <- "/Users/maurertm/Desktop/New_ADRC_Key.xlsx"

##############################################################################################
#############################-----Load and Clean Metadata-----################################
##############################################################################################
metadata <- read_csv(metadata_filepath) %>% filter(has_WGS_AND_clock_data == 1) %>% filter(CHIP == 1) %>% select(filename, vcf_sample_header, Gender, contains("Hugo_"), CHIP, WGS_Study, Study_Prot_Clocks, Largest_VAF, has_chip, chip_class, chip_class2, Largest_VAF, Multiple_VAF, contains("tumor_f"), contains("t_ref"), contains("t_alt"))

Multiple_carriers <- filter(metadata, chip_class == "Multiple")
count_max_genes <- metadata %>% select(contains("Hugo")) %>% ncol

dataframe <- Multiple_carriers

Multiple_carriers_long <- combined_split_genes(Multiple_carriers, 4) %>%
  bind_rows(.id = "Mutation_Number") %>% select(-Mutation_number) %>%
  filter(!is.na(VAF))

Single_carriers <- filter(metadata, chip_class != "Multiple")
Single_carriers$Gene <- Single_carriers$chip_class
Single_carriers$AD <- paste(Single_carriers$t_ref_count, ",", Single_carriers$t_alt_count, sep = "")
Single_carriers$VAF <- Single_carriers$tumor_f
Single_carriers$Mutation_Number <- NA
Single_carriers <- Single_carriers %>% select(filename, Hugo_Symbol, AD, VAF, CHIP, WGS_Study, Gender, Mutation_Number, vcf_sample_header)%>% 
  set_colnames(c("filename", "Gene", "AD", "VAF", "haschip", "STUDY", "Gender", "Mutation_Number", "vcf_sample_header"))

combined_metadata <- bind_rows(Multiple_carriers_long, Single_carriers) 
combined_metadata$INFERRED_SEX <- ifelse(combined_metadata$Gender == "M", 1, 2)
combined_metadata$haschip <- ifelse(combined_metadata$haschip, 1, 2)

combined_metadata_final <- combined_metadata %>% select(vcf_sample_header, INFERRED_SEX, Gene, haschip, STUDY, AD, VAF, Mutation_Number) %>%
  set_colnames(c("Sample", "INFERRED_SEX", "Gene", "haschip", "STUDY", "AD", "VAF", "Mutation_number"))

##############################################################################################
##########################-----Create Metadata DFs 1 and 2-----###############################
##############################################################################################
sample_df <- combined_metadata_final %>% filter(is.na(Mutation_number) | Mutation_number == 1) %>% select(Sample, INFERRED_SEX, Gene, haschip)
sample_df$Study <- "ADRC"
colnames(sample_df) <- c("Sample", "INFERRED_SEX", "Gene", "haschip", "STUDY")
sample_df <- sample_df %>% filter(!is_in(Sample, c("Stanford_0563", "Stanford_0583"))) %>% distinct

variant_df <- combined_metadata_final %>% select(Sample, Gene, AD, VAF)
colnames(variant_df) <- c("Sample", "Gene", "AD", "VAF")
variant_df <- variant_df %>% filter(!is_in(Sample, c("Stanford_0563", "Stanford_0583"))) %>% distinct

##############################################################################################
##############################-----Output Metadata DFs-----###################################
##############################################################################################
write_tsv(sample_df, "/Users/maurertm/Desktop/metadata_df_1.txt" )
write_tsv(variant_df, "/Users/maurertm/Desktop/metadata_df_2.txt" )
