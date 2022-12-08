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
PACER_list_filepath <- "/Users/maurertm/Desktop/Projects/ADRC/PACER/WGS_files.csv"
ADRC_key_filepath <- "/Users/maurertm/Desktop/Projects/ADRC/PACER/New_ADRC_Key.xlsx"

##############################################################################################
#############################-----Load and Clean Metadata-----################################
##############################################################################################
in_wgs_missing_phen <- c("584","465", "931", "415", "583", "565", "474","610", "563", "84478")
metadata <- read_csv(metadata_filepath) %>% filter(!is.na(filename))
  
PACER <- read_csv(PACER_list_filepath)
colnames(PACER) <- "filepath"
PACER$filename <- gsub("\\..*", "", PACER$filepath) 
PACER$filename_zero_removed<- str_replace(PACER$filename, "^0+" ,"")
PACER$filename_underscore_removed <- str_replace(PACER$filename_zero_removed, "_.*" ,"")
PACER <- PACER %>% filter(!is_in(filename_underscore_removed, in_wgs_missing_phen))

incorrect_ids <- setdiff(PACER$filename_underscore_removed, metadata$filename)
PACER_correct_ids <- PACER %>% filter(!is_in(filename_underscore_removed, incorrect_ids))
PACER_correct_ids$metadata_filename <- PACER_correct_ids$filename_underscore_removed
PACER_correct_ids$directory_filename <- PACER_correct_ids$filename
PACER_correct_ids <- PACER_correct_ids %>% select(metadata_filename, directory_filename)

##the mystery BAM/BAI ID is useless to use as it matches to nothing in the metadata
##so I created a "key" datarame containing the ADRC_IDs and their corresponding BAM/ BAI IDs
new_ADRC_table <- read_excel(ADRC_key_filepath)
colnames(new_ADRC_table) <- c("ADRC_ID", "BAM_BAI_ID", "first", "filename_underscore_removed", "last" )
ADRC_key <- select(new_ADRC_table, ADRC_ID, filename_underscore_removed)
ADRC_key <- transform(ADRC_key,  filename_underscore_removed= as.character(filename_underscore_removed))

PACER_incorrect_ids <- PACER %>% filter(is_in(filename_underscore_removed, incorrect_ids))
PACER_incorrect_ids <- left_join(PACER_incorrect_ids, ADRC_key, on = "filename_underscore_removed")
PACER_incorrect_ids$metadata_filename <- PACER_incorrect_ids$ADRC_ID
PACER_incorrect_ids$directory_filename <- PACER_incorrect_ids$filename
PACER_incorrect_ids <- PACER_incorrect_ids %>% select(metadata_filename, directory_filename)

PACER_data <- rbind(PACER_correct_ids, PACER_incorrect_ids)
PACER_data <- transform(PACER_data,  metadata_filename= as.numeric(metadata_filename))

intersect(PACER_data$metadata_filename, metadata$filename) %>% length
#525
setdiff(PACER_data$metadata_filename, metadata$filename)
#character(0)
colnames(PACER_data) <- c("filename", "directory_filename")

metadata_filenames<- PACER_data$filename
metadata_filtered <- metadata %>% filter(is_in(filename, metadata_filenames))

metadata_PACER_joined <- left_join(PACER_data, metadata_filtered, on="filename") %>% distinct

metadata_PACER_joined <-  metadata_PACER_joined %>% select(directory_filename, filename, vcf_sample_header, Gender, contains("Hugo_"), CHIP, WGS_Study, Study_Prot_Clocks, Largest_VAF, has_chip, chip_class, chip_class2, Largest_VAF, Multiple_VAF, contains("tumor_f"), contains("t_ref"), contains("t_alt"))

##############################################################################################
######################-----Split into Single and Multiple Carriers-----#######################
##############################################################################################
#first, split up into who has chip and who doesn't
no_chip <- metadata_PACER_joined %>% filter(has_chip == FALSE) 
no_chip$Gene <- NA
no_chip$AD <- NA
no_chip$VAF <- NA
no_chip$Mutation_Number <- NA
no_chip <- no_chip %>% select(directory_filename, filename, Hugo_Symbol, AD, VAF, CHIP, Gender, Mutation_Number)%>% 
  set_colnames(c("directory_filename", "filename", "Gene", "AD", "VAF", "haschip", "Gender", "Mutation_Number")) %>% distinct()

has_chip <- metadata_PACER_joined %>% filter(has_chip == TRUE)

intersect(no_chip$filename, has_chip$filename)
#expected: numeric(0)

#then, split up into multiple and single chip carriers
Multiple_carriers <- filter(has_chip, chip_class == "Multiple")
count_max_genes <- has_chip %>% select(contains("Hugo")) %>% ncol
#expected: 4

Multiple_carriers_long <- combined_split_genes(Multiple_carriers, 4) %>%
  bind_rows(.id = "Mutation_Number") %>% select(-Mutation_number) %>%
  filter(!is.na(VAF)) %>% unique

Multiple_carriers_long <- Multiple_carriers_long %>% 
  select(Mutation_Number, filename, Gene, AD, VAF, haschip, Gender, Mutation_Number, vcf_sample_header) %>% 
  left_join(PACER_data, on = filename) %>% 
  select(directory_filename, filename, Gene, AD, VAF, haschip, Gender, Mutation_Number)

Single_carriers <- filter(has_chip, chip_class != "Multiple")
Single_carriers$Gene <- Single_carriers$chip_class
Single_carriers$AD <- paste(Single_carriers$t_ref_count, ",", Single_carriers$t_alt_count, sep = "")
Single_carriers$VAF <- Single_carriers$tumor_f
Single_carriers$Mutation_Number <- NA
Single_carriers <- Single_carriers %>% select(directory_filename, filename, Hugo_Symbol, AD, VAF, CHIP, Gender, Mutation_Number)%>% 
  set_colnames(c("directory_filename", "filename", "Gene", "AD", "VAF", "haschip", "Gender", "Mutation_Number")) %>% distinct()

intersect(Single_carriers$filename, Multiple_carriers$filename)
#expected: numeric(0)

#then, combine all 3 dataframes
combined_metadata <- bind_rows(Multiple_carriers_long, Single_carriers) 
combined_metadata <- bind_rows(combined_metadata, no_chip) 
combined_metadata$filename %>% unique %>% length
#expected:25

combined_metadata$INFERRED_SEX <- ifelse(combined_metadata$Gender == "M", 1, 2)
combined_metadata$haschip <- ifelse(combined_metadata$haschip, 1, 2)
combined_metadata$Mutation_Number <- ifelse(is.na(combined_metadata$Mutation_Number), 0, combined_metadata$Mutation_Number)


combined_metadata_final <- combined_metadata %>% select(directory_filename, INFERRED_SEX, Gene, haschip, AD, VAF, Mutation_Number) %>%
  set_colnames(c("Sample", "INFERRED_SEX", "Gene", "haschip", "AD", "VAF", "Mutation_number"))

##############################################################################################
##########################-----Create Metadata DFs 1 and 2-----###############################
##############################################################################################
#metadata1: Sample, INFERRED_SEX, Gene, haschip, STUDY
sample_df <- combined_metadata_final %>% filter(Mutation_number < 2) %>% 
                                         select(Sample, INFERRED_SEX, Gene, haschip) %>% unique()

sample_df_multiple_sex <- sample_df %>% filter(is_in(Sample, c("0423", "39393")))
#39393" is male according to redcap (1)
#423 is female according to redcap (2)
sample_df_multiple_sex <- sample_df_multiple_sex[-c(1,4),]

sample_df_single_sex <- sample_df %>% filter(!is_in(Sample, c("0423", "39393")))

sample_df_final <- rbind(sample_df_multiple_sex, sample_df_single_sex)
sample_df_final %>% nrow
#expected: 525
sample_df_final$STUDY <- "ADRC"
colnames(sample_df_final) <- c("Sample", "INFERRED_SEX", "Gene", "haschip", "STUDY")


#metadata2: Sample, Gene, AD, VAF
variant_df <- combined_metadata_final %>% distinct
variant_df_multiple_sex <- variant_df %>% filter(is_in(Sample, c("0423", "39393")))
variant_df_multiple_sex <- variant_df_multiple_sex[-c(1,4),]
variant_df_single_sex <- variant_df %>% filter(!is_in(Sample, c("0423", "39393")))
variant_df_combined <- rbind(variant_df_single_sex, variant_df_multiple_sex)
variant_df_combined %>% unique %>% filter(Mutation_number < 2) %>% group_by(Sample) %>% summarise(count = n_distinct(AD)) %>% arrange(desc(count)) %>% print(n=10)
variant_df %>% filter(is_in(Sample, c("3900004", "3900032", "3900152", "3900157", "3900246"))) %>% arrange(desc(Sample)) %>% print(n=30)
metadata %>% filter(is_in(filename, c("3900004", "3900032", "3900152", "3900157", "3900246"))) %>% select(contains("Protein_Change"), filename, contains("tumor_f")) %>% arrange(desc(filename)) %>% print(n=30)

variant_df_multiple_mutations_1_gene <- variant_df %>% filter(is_in(Sample, c("3900004", "3900032", "3900152", "3900157", "3900246"))) %>% arrange(desc(Sample))
variant_df_multiple_mutations_1_gene$Mutation_number <- c(1,2,3,4, 1,2,3,4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 5, 6)

variant_df_one_mutations_1_gene <- variant_df %>% filter(!is_in(Sample, c("3900004", "3900032", "3900152", "3900157", "3900246")))

variant_df_combined_2 <- rbind(variant_df_multiple_mutations_1_gene, variant_df_one_mutations_1_gene)
variant_df_combined_2 %>% unique %>% filter(Mutation_number < 2) %>% group_by(Sample) %>% summarise(count = n_distinct(AD)) %>% arrange(desc(count)) %>% print(n=10)
#expected: 525 rows

variant_df_final <- variant_df_combined_2 %>% select(Sample, Gene, AD, VAF)
colnames(variant_df_final) <- c("Sample", "Gene", "AD", "VAF")
variant_df_final <- variant_df_final %>% distinct

##############################################################################################
##############################-----Output Metadata DFs-----###################################
##############################################################################################
write_tsv(sample_df_final, "/Users/maurertm/Desktop/metadata_df_1.txt" )
write_tsv(variant_df_final, "/Users/maurertm/Desktop/metadata_df_2.txt" )
