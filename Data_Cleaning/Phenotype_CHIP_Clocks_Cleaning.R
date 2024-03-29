library(readxl)
library(magrittr)
library(stringr)
library(arrow)
library(tidyverse)

choose_best <- function(metadata_rows, metadata_df) { 
  arrange(metadata_rows, desc(ConnectivityZscore)) %>% slice(1)
}
#get the unqiue value in each row that matches to an scg/ WGS id (if there is one)
get_unique_value <- function(row, row_df) {
  unique_ids <- c(row$PIDN_files, row$SampleID_files, row$ADRC_ID_files) %>% na.omit %>% unique
  row$filename <- str_c(unique_ids, collapse = ";")
  row
}

widen_df <- function(column_name, subject_rows) {
  subject_rows_wider <- arrange(subject_rows, desc(tumor_f)) %>% select(ID, !!column_name) %>% 
    pivot_wider(names_from = "ID", values_from = column_name) 
  colnames(subject_rows_wider) <- str_c(column_name, c("", 2:ncol(subject_rows_wider)))
  subject_rows_wider
}

choose_largest_vaf <- function(subject_rows, subject_df) {
  arrange(subject_rows, desc(tumor_f)) %>% slice(1)
}

squash_multiple_variants <- function(subject_rows, subject_df) {
  temp_ids <- str_c(subject_df$Sample[1], 1:nrow(subject_rows), sep = "_")
  subject_rows$ID <- temp_ids
  colnames_widen <- select(subject_rows, -ID) %>% colnames
  subject_rows_wider <- map(colnames_widen, widen_df, subject_rows) %>% bind_cols
  subject_rows_wider
}

filter_results <- function(dataframe){
  pass <- dataframe %>%
    mutate(CT = if_else(
      (REF=="C" & ALT =="T") | (REF=="T" & ALT =="C") | (REF=="G" & ALT =="A") | (REF=="A" & ALT =="G"), true="yes", false="no"))
  
  pass_filter <- pass %>%
    filter (CT == "yes") %>%
    filter (TLOD > 10) %>%
    filter (VAF < 0.3) %>%
    filter (VAF > 0.0) %>%
    filter (AD>1) %>%
    count(NWD_ID)
  
  pass_filter
  
}

##############################################################################################
#################################-----Define Data Paths-----##################################
##############################################################################################
#Import clock data from the 4 different ADRC studies
Wagner_filepath <- "ADRC_Project/Source_Data/Aging_Clock_Data/df_all_tissue_bootstrap_lasso_ENRICHED_COGNITIONBRAIN_ABOVE60_WAGNER_RAW_Age_prediction_mean_results.csv"
Poston_filepath <- "ADRC_Project/Source_Data/Aging_Clock_Data/df_all_tissue_bootstrap_lasso_ENRICHED_COGNITIONBRAIN_ABOVE60_POSTON_RAW_Age_prediction_mean_results.csv"
Kercher_filepath <- "ADRC_Project/Source_Data/Aging_Clock_Data/df_all_tissue_bootstrap_lasso_ENRICHED_COGNITIONBRAIN_ABOVE60_KERCHNER_RAW_Age_prediction_mean_results.csv"
ADRC_filepath <- "ADRC_Project/Source_Data/Aging_Clock_Data/df_all_tissue_bootstrap_lasso_ENRICHED_COGNITIONBRAIN_ABOVE60_ADRC_RAW_Age_prediction_mean_results.csv"

phen_metadata_filepath <- "ADRC_Project/Source_Data/Metadata/Clinical_Sample_Level/Plasma_metadata_FINAL_052021_ADRC_additionalQC_new.csv"
chip_carriers_filepath <- "ADRC_Project/Source_Data/mutect2_Results/Chip_non_Chip_Lists/chip_carriers.txt"
non_chip_carriers_filepath <- "ADRC_Project/Source_Data/mutect2_Results/Chip_non_Chip_Lists/non_chip_carriers.txt"
new_ADRC_table_filepath <- "ADRC_Project/Source_Data/Metadata/Filename_Correction_Keys/New_ADRC_Key.xlsx"

Mutect_results_with_CHIP_filepath <- "ADRC_Project/Source_Data/mutect2_Results/mutect_somatic_042822.csv"

Yanns_metadata_filepath <- "ADRC_Project/Source_Data/Metadata/Phenotype/Incomplete_Phenotype_MetaData2.csv"

vcf_header_incorrect_names <- "ADRC_Project/Source_Data/Metadata/Filename_Correction_Keys/whole_exome_normal_headers_list.txt"
vcf_header_normal_names <- "ADRC_Project/Source_Data/Metadata/Filename_Correction_Keys/whole_exome_incorrect_headers_list.txt"

parquet_fp <- "/ADRC_Project/Results/PACER/singletons_2023_02_03.parquet"
parquet_old_fp <- "/ADRC_Project/Results/PACER/pacer_output.tsv"
parquet_100_fp <- "/ADRC_Project/Results/PACER/singletons_2023_02_05_100.parquet"
parquet_200_fp <- "/ADRC_Project/Results/PACER/singletons_2023_02_05_200.parquet"
parquet_300_fp <- "/ADRC_Project/Results/PACER/singletons_2023_02_06_300.parquet"
parquet_400_fp<- "/ADRC_Project/Results/PACER/singletons_2023_02_06_400.parquet"

##############################################################################################
####################-########----Load and Combine Clock Data-----#############################
##############################################################################################
Wagner <- read_csv(Wagner_filepath)
Poston <- read_csv(Poston_filepath)
Kercher <- read_csv(Kercher_filepath)
ADRC <- read_csv(ADRC_filepath)

#Join the files
#expect and ignore warning, it's fine
Wagner <- as.tibble(Wagner) %>% mutate(Study = "Wagner")
Poston <- as.tibble(Poston) %>% mutate(Study = "Poston")
Kercher <- as.tibble(Kercher) %>% mutate(Study = "Kercher")
ADRC <- as.tibble(ADRC) %>% mutate(Study = "ADRC")

Wagner_Poston <- rbind(Wagner, Poston) 
Wagner_Poston_Kercher <- rbind(Wagner_Poston, Kercher)
clock_combined <- rbind(Wagner_Poston_Kercher, ADRC) %>% 
  select(-dage_resid, -dage_resid_zscored, -yhat)
pivot_clock <- pivot_wider(clock_combined, names_from = "tissue", values_from = "Pred_Age")

##############################################################################################
########################-----Import and Clean Phenotype Metadata-----#########################
##############################################################################################
phen_metadata <- read_csv(phen_metadata_filepath)

#three people have multiple genders listed. I went to redcap and pulled their genetic sex
#it is unlikely a study of this size would have three transgender participants, but regardless, we should confirm with Jarod that all 'genders' recorded are actaully genetic sexes
wrong_gender_adrc <- c("3900044", "3900050", "3900055")
phen_metadata[is_in(phen_metadata$ADRC_ID, wrong_gender_adrc),]$Gender <- "F"

##############################################################################################
##########################-----Combine Clock and Phenotype Data-----##########################
##############################################################################################
setdiff(phen_metadata$Barcode, pivot_clock$Barcode)
#expect:character(0)
missing_in_metadata <- setdiff(pivot_clock$Barcode, phen_metadata$Barcode)
#expect: "S1228665" "S1228593" "S1228625" "S1228535" "S1228669" "S1229627" "S1227123"
clock_pivot_filtered <- filter(pivot_clock, !is_in(Barcode, missing_in_metadata))
colnames(clock_pivot_filtered) <- c("Barcode", "Age_at_Prot_Draw", "Study_Prot_Clocks", "Adipose", "Artery", "Brain", "CognitionBrain", "Esophagus", "Heart",         
                                    "Immune", "Intestine", "Kidney", "Liver", "Lung", "Muscle", "Organismal", "Pancreas","Pituitary",     
                                    "Salivary", "Stomach")

#join pivot_filtered and phen_metadata on barcodes
join_on_barcodes <- left_join(phen_metadata, clock_pivot_filtered, on="Barcode") %>% filter(!is.na(Adipose))

#check that there are no rows with missing clock data (this returns an empty tibble)
left_join(phen_metadata, clock_pivot_filtered, on="Barcode") %>% 
  filter(is.na(Adipose))

#find which Individual IDs contain multiple barcodes (and thus multiple proteomic samples)
barcode_tally <- join_on_barcodes %>% 
  group_by(Individual_ID) %>% 
  tally() %>% 
  arrange(desc(n)) %>% 
  filter(n > 1)

#split the dataset into those with multiple barcodes/IID and those without
#we are only looking at "Barcode" and "Indivdual_ID" so we can make a "key" of sorts
dup_barcodes <- join_on_barcodes %>% filter(is_in(Individual_ID, barcode_tally$Individual_ID)) %>% select(Individual_ID, Barcode, ConnectivityZscore)
single_barcodes <- join_on_barcodes %>% filter(!is_in(Individual_ID, barcode_tally$Individual_ID)) %>% select(Individual_ID, Barcode, ConnectivityZscore)

duplicate_barcodes_best <- group_by(dup_barcodes, Individual_ID) %>% group_modify(choose_best)

best_barcodes <- rbind(single_barcodes, duplicate_barcodes_best)

#now remove any barcodes with a Connectivityscore less than -2
valid_barcodes <- filter(best_barcodes, ConnectivityZscore >= -2)
filter(best_barcodes, ConnectivityZscore <= -2)

#filter metadata for is in valid_barcodes
barcodes_better <- join_on_barcodes %>% filter(is_in(Barcode, valid_barcodes$Barcode))  %>% 
  select(Barcode, Individual_ID)

add_individual_id_df <- filter(clock_pivot_filtered,is_in(Barcode, barcodes_better$Barcode)) %>% 
  left_join(barcodes_better) %>% 
  select(-Barcode)

add_individual_id_df %>% group_by(Individual_ID) %>% tally() %>% arrange(desc(n)) %>% filter(n > 1)
#^ this gives 0 and is a test to show that each Individual_ID only contains 1 barcode

#now join to phen_metadata by the Individual_ID in add_individual_id_df
phenotype_clocks <- left_join(phen_metadata, add_individual_id_df, by = "Individual_ID") %>% filter(!is.na(Adipose))

#make a columns called has_clock_data, which contains TRUE and FALSE for if we have proteomic data on the individual
phenotype_clocks$has_clock_data <- !is.na(phenotype_clocks$Adipose)

#change some column names
colnames(phenotype_clocks)[5] <- "Age_at_WGS_Draw"
colnames(phenotype_clocks)[16] <- "WGS_Study"

##############################################################################################
########################-----Get list of Samples with WGS Data-----###########################
##############################################################################################
#the bam fiels are located in 2 different folders. Here are their contents (with .hg38.bam and .hg38.bai removed)
chip_carriers <- read_tsv(chip_carriers_filepath) %>% as.tibble
non_chip_carriers <- read_tsv(non_chip_carriers_filepath) %>% as.tibble

incorrect_ids <- read_tsv(vcf_header_incorrect_names, col_names = FALSE)
incorrect_ids[c('X1', 'X2')] <- str_split_fixed(incorrect_ids$X1, ' ', 2)
incorrect_ids <- as.data.frame(sapply(incorrect_ids, function(x) gsub("\"", "", x)))

normal_ids <- read_tsv(vcf_header_normal_names, col_names = FALSE)

filenames_vcf_ids <- bind_rows(incorrect_ids, normal_ids)
filenames_vcf_ids <- filenames_vcf_ids %>% mutate(Sample = gsub(".hg38_GRCh38_mutect2.vcf", "", X1))

colnames(chip_carriers) <- "Sample"
colnames(non_chip_carriers) <- "Sample"

combined_sequenced_list <- rbind(chip_carriers, non_chip_carriers) %>% distinct

setdiff(as.numeric(combined_sequenced_list$Sample), as.numeric(filenames_vcf_ids$Sample))
#expect: (and Warning Message)
#"0428"    "0449"    "0881"    "0894"    "0919"    "1065"    "1095"    "3900049" "3900131" "3900227" "3900294"
#"3900295" "3900318" "3900357" "3900358" "3900379" "3900385" "3900388" "3900393" "3900394" "3900400" "3900402"
setdiff(filenames_vcf_ids$Sample, combined_sequenced_list$Sample)
#expect: character(0)

complete_scg_id_df <- left_join(combined_sequenced_list, filenames_vcf_ids, by = "Sample")

#some of the files are not labeled by any recognized id in the phenotype or clock data (I call these the "BAM/BAI mystery ids")
#these files originally  came insets of folders, where the folders was named with the subjects ADRC id. Within each folder, the bam/ bai files were named with the "BAM/BAI mystery id"
##the mystery BAM/BAI ID is useless to use as it matches to nothing in the metadata
##so I created a "key" datarame containing the ADRC_IDs and their corresponding BAM/ BAI IDs
new_ADRC_table <- read_excel(new_ADRC_table_filepath)
colnames(new_ADRC_table) <- c("ADRC_ID", "BAM_BAI_ID", "first", "Sample", "last" )
ADRC_key <- select(new_ADRC_table, ADRC_ID, Sample)
ADRC_key$Sample <- as.numeric(ADRC_key$Sample)
mystery_bam_bais <- ADRC_key %>% pull(Sample)
new_bams_ADRC_IDs <- ADRC_key %>% select(ADRC_ID)
colnames(new_bams_ADRC_IDs) <- "Sample"

#split the combined sequence
new_ids <- filter(complete_scg_id_df, is_in(Sample, mystery_bam_bais)) %>% distinct 
new_ids$Sample <- as.numeric(new_ids$Sample)
setdiff(new_ids$Sample, ADRC_key$Sample)
#expect: numeric(0)
setdiff(ADRC_key$Sample, new_ids$Sample)
#expect: numeric(0)

new_ids_joined <- left_join(new_ids, ADRC_key, by="Sample") %>% select(-Sample)
colnames(new_ids_joined) <- c("name_of_file", "vcf_sample_header", "Sample")

old_ids <- filter(complete_scg_id_df, !is_in(Sample, mystery_bam_bais)) 
colnames(old_ids) <- c("Sample", "name_of_file", "vcf_sample_header")

#replace the old mystery ids with their corresponding ADRC_IDs
updated_ids <- rbind(old_ids, new_ids_joined) %>% distinct

#remove the underscore and leading zeros
updated_ids$no_leading_zero <- str_remove_all(updated_ids$Sample, "^0")

samples_unique <-  mutate(updated_ids, no_underscore = sub("\\_.*", "", no_leading_zero)) 

genotypes <- dplyr::select(samples_unique, no_underscore, name_of_file, vcf_sample_header)
genotypes %>% group_by(no_underscore) %>% tally() %>% arrange(desc(n)) %>% filter(n>1)
#expect:  868          2
problematic <- genotypes %>% filter(no_underscore == "868")
# A tibble: 2 × 3
#filename name_of_vcf                  vcf_sample_header 
#1 868      0868.hg38_GRCh38_mutect2.vcf "   Stanford_0868"
#2 868      0868.hg38_GRCh38_mutect2.vcf "Stanford_0868
problematic <- problematic[2,]
non_problematic <- genotypes %>% filter(no_underscore != "868")
genotypes <- rbind(problematic, non_problematic)

colnames(genotypes) <- c("filename", "name_of_vcf", "vcf_sample_header")
WGS_sample_list <- genotypes %>% select(filename) %>% as.list

##############################################################################################
#############-----Get Column with SCG_ID (aka filename) in phenotype_clocks-----##############
##############################################################################################
#clean PIDN, SampleID, and ADRC to match the patterns in genotype (no leading 0's, etc.)
phenotype_clocks$PIDN_cleaned <- str_remove_all(phenotype_clocks$PIDN, "^0") %>%
  str_remove_all(" .*$") %>% str_remove_all("_.*$") %>% str_remove_all("[A-Z]")

phenotype_clocks$SampleID_cleaned <- str_remove_all(phenotype_clocks$SampleID, "^0") %>%
  str_remove_all(" .*$") %>% str_remove_all("_.*$") %>% str_remove_all("[A-Z]")

phenotype_clocks$ADRC_ID_cleaned <- str_remove_all(phenotype_clocks$ADRC_ID, "^0") %>%
  str_remove_all(" .*$") %>% str_remove_all("_.*$") %>% str_remove_all("[A-Z]")
#create columns called PIDN_files, SampleID_files, ADRC_files that contain the id in scg (WGS id)
#if these columns are NA it means that ID pattern is not present in the scg/ WGS ids
phenotype_clocks$PIDN_files <- phenotype_clocks$PIDN_cleaned
phenotype_clocks$PIDN_files[!is_in(phenotype_clocks$PIDN_files, genotypes$filename)] <- NA
intersect(genotypes$filename, phenotype_clocks$PIDN_cleaned)
#expect: 163 results

phenotype_clocks$SampleID_files <- phenotype_clocks$SampleID_cleaned
phenotype_clocks$SampleID_files[!is_in(phenotype_clocks$SampleID_files, genotypes$filename)] <- NA
intersect(genotypes$filename, phenotype_clocks$SampleID_cleaned)
#expect: 147 results

phenotype_clocks$ADRC_ID_files <- phenotype_clocks$ADRC_ID_cleaned
phenotype_clocks$ADRC_ID_files[!is_in(phenotype_clocks$ADRC_ID_files, genotypes$filename)] <- NA
intersect(genotypes$filename, phenotype_clocks$ADRC_ID_cleaned)
#expect: 384 results

#these are the data in scg (the WGS data) that is not in the metadata file Jarod provided
#notice there are about 11- this is what we expect
genotypes$missing_temp <- !is_in(genotypes$filename, phenotype_clocks$PIDN_cleaned) & !is_in(genotypes$filename, phenotype_clocks$SampleID_cleaned) & !is_in(genotypes$filename, phenotype_clocks$ADRC_ID_cleaned)

missing_WGS_in_metadata <- filter(genotypes, missing_temp == TRUE) %>% pull(filename)

#create a column called filename that contains the scg/ WGS id
phenotype_clocks_with_filenames <- rowwise(phenotype_clocks) %>% group_map(get_unique_value) %>% bind_rows
phenotype_clocks_with_filenames$filename[nchar(phenotype_clocks_with_filenames$filename) == 0] <- NA

#chack that no Individual id is linked to multiple filenames, then link filename to Individual_ID so the CHIP data is linked to the person instead of the visit
phenotype_clocks_with_filenames %>% select(Individual_ID, filename) %>% filter(!is.na(filename)) %>% group_by(Individual_ID) %>% distinct %>% tally() %>% arrange(desc(n))
phenotype_clocks_with_IID_filenames <- phenotype_clocks_with_filenames %>% select(Individual_ID, filename) %>% filter(!is.na(filename)) %>% group_by(Individual_ID) %>% distinct 
phenotype_clocks_wo_filenames <- phenotype_clocks_with_filenames %>% select(-filename)
phenotype_clocks_with_filenames <- left_join(phenotype_clocks_wo_filenames, phenotype_clocks_with_IID_filenames, on = "Individual_ID")

##############################################################################################
##################-----Add has_WGS_data and has_clocks_data Columns-----######################
##############################################################################################
#make a column called has_WGS_data which contains TRUE and FALSE for if we have WGS data on the individual
phenotype_clocks_with_filenames$has_WGS_data <- !is.na(phenotype_clocks_with_filenames$filename)
phenotype_clocks_with_filenames$has_clock_data <- !is.na(phenotype_clocks_with_filenames$Adipose)

phenotype_clocks_with_filenames %>% filter(has_WGS_data == TRUE & has_clock_data == TRUE) %>% select(Individual_ID) %>%unique
#551 match

#make a columne called has_WGS_AND_proteomic_data, which contains TRUE and FALSE for if we have both data on the individual
phenotype_clocks_with_filenames$has_WGS_AND_clock_data <- phenotype_clocks_with_filenames$has_WGS_data == TRUE & phenotype_clocks_with_filenames$has_clock_data == TRUE

phenotype_clocks_with_filenames %>% filter(has_WGS_AND_clock_data) %>% pull(Individual_ID) %>% unique()
#551 match

phenotype_clocks_with_genotypes <- left_join(phenotype_clocks_with_filenames, genotypes, by = "filename")

##############################################################################################
############################-----Import and Clean Mutect Data-----############################
##############################################################################################
#the samples in the file below have CHIP. Other files the WGS sequencing was completed on were determined in the steps above
Mutect_results_with_CHIP <- read_csv(Mutect_results_with_CHIP_filepath)
Mutect_results_with_CHIP$Sample %<>% str_remove_all("\\.hg38") %<>% str_remove_all("^0") %<>% str_remove_all("_2$")

#we need to transform the "mystery bam bai" ids in this dataframe into their corresponding ADRC_ID
old_ids <- filter(Mutect_results_with_CHIP, !is_in(Sample, mystery_bam_bais)) %>%  transform(Sample = as.numeric(Sample))
new_ids <- filter(Mutect_results_with_CHIP, is_in(Sample, mystery_bam_bais)) %>%  transform(Sample = as.numeric(Sample))

ADRC_key_has_CHIP <- filter(ADRC_key, is_in(Sample, new_ids$Sample))
updated_new_ids_with_CHIP <- new_ids %>% left_join(ADRC_key_has_CHIP, by = "Sample") %>% select(-Sample)

colnames(updated_new_ids_with_CHIP)[33] <- "Sample"

has_CHIP_df <- rbind(updated_new_ids_with_CHIP, old_ids)
has_CHIP_df$CHIP <- 1

no_CHIP <- setdiff(genotypes$filename, has_CHIP_df$Sample)
no_CHIP_df <- tibble(Sample = no_CHIP, CHIP = 0) %>% transform(Sample = as.numeric(Sample))

mutect_results <- bind_rows(has_CHIP_df, no_CHIP_df)

##############################################################################################
#######################-----Fix Multiple Mutation Category-----###############################
##############################################################################################
#results labled somatic have been judged  by Sidd to have CHIP
results_wo_CHIP <- mutect_results %>% filter(CHIP == 0)

results_with_CHIP <- mutect_results %>% filter(CHIP == 1) %>% select(Sample, Hugo_Symbol, Variant_Classification, Protein_Change, tumor_f, t_ref_count, t_alt_count) %>% distinct

#need to deal with if people have multiple mutatiosn
results_n_mutations <- group_by(results_with_CHIP, Sample) %>% tally %>% arrange(desc(n))
results_n_multiple <- filter(results_n_mutations, n > 1)

results_multiple <- filter(mutect_results, is_in(Sample, results_n_multiple$Sample)) %>% arrange(Sample) 

results_multiple_largest_vaf <- results_multiple %>% group_by(Sample) %>% group_modify(choose_largest_vaf) %>%
  select("Sample", "Hugo_Symbol", "tumor_f", "t_ref_count", "t_alt_count")

colnames(results_multiple_largest_vaf) <- c("Sample", "Largest_VAF", "Multiple_VAF", "Multiple_t_ref", "Multiple_t_alt")

results_multiple_squashed <- group_by(results_multiple, Sample) %>% group_modify(squash_multiple_variants) %>%
  select(Sample, contains("Hugo_Symbol"), contains("Variant_Classification"), contains("Protein_Change"), contains("tumor_f"), contains("t_ref"), contains("t_alt"))
results_multiple_squashed$chip_class <- "Multiple"
results_multiple_squashed$CHIP <- 1
#results_multiple_squashed[3,26] <- "TET2"

results_multiple_squashed_joined <- left_join(results_multiple_largest_vaf, results_multiple_squashed, by = "Sample")

results_single <- filter(mutect_results, !is_in(Sample, results_n_multiple$Sample)) %>% arrange(Sample) 
results_single$chip_class <- results_single$Hugo_Symbol

mutect_results_all <- bind_rows(results_multiple_squashed_joined, results_single)
mutect_results_all$filename <- mutect_results_all$Sample %>% str_remove_all("^0") %>% str_remove_all("_.*$") 
mutect_results_all <- transform(mutect_results_all, filename = as.numeric(filename))

table(mutect_results_all$chip_class) %>% sort

##############################################################################################
####################-----BIND MUTECT DATA WITH PHEN/CLOCK DATA-----###########################
##############################################################################################
setdiff(phenotype_clocks_with_genotypes$filename, mutect_results_all$filename)
#NA
missing_CHIP_filenames <- setdiff(mutect_results_all$filename, phenotype_clocks_with_genotypes$filename)
#expect:  [1]     583     415     449     465     474     563     565     584     610     931 3900308

setdiff(missing_CHIP_filenames, phenotype_clocks_with_genotypes$ADRC_ID_cleaned)
#expect:  [1]     583     415     449     465     474     563     565     584     610     931 3900308
setdiff(missing_CHIP_filenames, phenotype_clocks_with_genotypes$SampleID_cleaned)
#expect:  [1]     583     415     449     465     474     563     565     584     610     931 3900308
setdiff(missing_CHIP_filenames, phenotype_clocks_with_genotypes$PIDN_cleaned)
#expect:  [1]     583     415     449     465     474     563     565     584     610     931 3900308

#we will filter out the missing ids and then add them in seperatly (as they are joining on a different column)
mutect_results_filtered <- filter(mutect_results_all, !is_in(filename, missing_CHIP_filenames))
mutect_results_filtered <- transform(mutect_results_filtered, filename = as.character(filename))

CHIP_phen_clocks <- left_join(phenotype_clocks_with_genotypes, mutect_results_filtered, by = "filename")

CHIP_phen_clocks_missing_filename <- CHIP_phen_clocks %>% filter(is.na(filename))
intersect(as.character(CHIP_phen_clocks_missing_filename$SampleId), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$SampleID), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$ADRC_ID), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$PIDN), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$Individual_ID), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$Barcode), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$ADRC_ID_cleaned), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$SampleID_cleaned), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$PIDN_cleaned), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$PIDN_files), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$SampleID_files), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$ADRC_ID_files), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$filename), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$name_of_vcf), as.character(mutect_results_all$Sample))
#expect: character(0)
intersect(as.character(CHIP_phen_clocks_missing_filename$Sample), as.character(mutect_results_all$Sample))
#expect: character(0)

#check that the number of rows in the new file is unchanged from proteomics_metadata_filename
nrow(CHIP_phen_clocks) - nrow(phenotype_clocks_with_filenames)


##############################################################################################
##################################-----Processing-----########################################
##############################################################################################

################################-----Mutation Type Split-----#################################
CHIP_phen_clocks$chip_class2 <- NA
CHIP_phen_clocks$chip_class2[CHIP_phen_clocks$chip_class == "DNMT3A"] <- "DNMT3A"
CHIP_phen_clocks$chip_class2[CHIP_phen_clocks$chip_class == "TET2"] <- "TET2"
CHIP_phen_clocks$chip_class2[CHIP_phen_clocks$chip_class == "ASXL1"] <- "ASXL1"
CHIP_phen_clocks$chip_class2[CHIP_phen_clocks$chip_class == "Multiple"] <- "Multiple"
CHIP_phen_clocks$chip_class2[is.na(CHIP_phen_clocks$chip_class)] <- "Control"
CHIP_phen_clocks$chip_class2[is.na(CHIP_phen_clocks$chip_class2)] <- "Other"
CHIP_phen_clocks$has_chip <- CHIP_phen_clocks$chip_class2 != "Control"
CHIP_phen_clocks$chip_class2 %<>% factor(levels = c("Control", "DNMT3A", "TET2", "Multiple", "ASXL1", "Other")) 

################################-----Age at Draw-----#################################
CHIP_phen_clocks$Age_at_Draw_Difference <- CHIP_phen_clocks$Age_at_WGS_Draw - CHIP_phen_clocks$Age_at_Prot_Draw

############################-----Singletons Resids-----################################
parquet <- read_parquet("/Users/maurertm/Desktop/ADRC_Project/Results/PACER/singletons_2023_02_03.parquet")
parquet_old <- read_tsv("/Users/maurertm/Desktop/ADRC_Project/Results/PACER/pacer_output.tsv")
parquet_100 <- read_parquet("/Users/maurertm/Desktop/ADRC_Project/Results/PACER/singletons_2023_02_05_100.parquet")
parquet_200 <- read_parquet("/Users/maurertm/Desktop/ADRC_Project/Results/PACER/singletons_2023_02_05_200.parquet")
parquet_300 <- read_parquet("/Users/maurertm/Desktop/ADRC_Project/Results/PACER/singletons_2023_02_06_300.parquet")
parquet_400 <- read_parquet("/Users/maurertm/Desktop/ADRC_Project/Results/PACER/singletons_2023_02_06_400.parquet")

tallied_sample_ids <- parquet %>% filter_results()
colnames(tallied_sample_ids) <- c("NWD_ID", "all_analysis_n")

tallied_sample_ids_old <- parquet_old  %>% filter_results() 
colnames(tallied_sample_ids_old) <- c("NWD_ID", "old_analysis_n")

tallied_sample_ids_100 <- parquet_100  %>% filter_results()
colnames(tallied_sample_ids_100) <- c("NWD_ID", "100_analysis_n")

tallied_sample_ids_200 <- parquet_200  %>% filter_results()
colnames(tallied_sample_ids_200) <- c("NWD_ID", "200_analysis_n")

tallied_sample_ids_300 <- parquet_300  %>% filter_results()
colnames(tallied_sample_ids_300) <- c("NWD_ID", "300_analysis_n")

tallied_sample_ids_400 <- parquet_400  %>% filter_results()
colnames(tallied_sample_ids_400) <- c("NWD_ID", "400_analysis_n")

joined_old_new <- left_join(tallied_sample_ids, tallied_sample_ids_400) %>% 
  left_join(tallied_sample_ids_old) %>%
  left_join(tallied_sample_ids_300) %>%
  left_join(tallied_sample_ids_200) %>%
  left_join(tallied_sample_ids_100) 

colnames(joined_old_new) <- c("vcf_sample_header", "all_analysis_n", "400_analysis_n", "old_analysis_n", "300_analysis_n", "200_analysis_n", "100_analysis_n")
CHIP_phen_clocks_pacer <- left_join(CHIP_phen_clocks, joined_old_new, by="vcf_sample_header")


############################-----Save Output-----################################
filepath <- 'ADRC_Project/Processed_DataFrames/Chip_Phen_Clock.csv'
time <- paste0(sub('\\..*', '', filepath), format(Sys.time(),'_%m_%d_at_%H_%M'), '.csv')
write_csv(CHIP_phen_clocks, time)
