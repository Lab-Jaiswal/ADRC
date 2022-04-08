library(writexl)
library(biomaRt)
library(rtracklayer)
library(GenomicRanges)
library(readxl)
library(magrittr)
library(tidyverse)
library(stringr)


##############################################################################################################
################################-------REARRANGE CLOCK METADATA FILE-------###################################
##############################################################################################################
#the clock metadata file has three columns: tissue, id (428 unique), and mean_dage_resid_zscored
#I want to have one row per sample id, which contians all of the mean residual zscores for each tissue (each in a seperate column)
clock_metadata_file <-"filepath"
clock_metadata <- read_csv(clock_metadata_file)
proteomics <- pivot_wider(clock_metadata, names_from = "tissue", values_from = "mean_dage_resid_zscored")
colnames(proteomics)[1] <- "Barcode"

##############################################################################################################
##################################-------ADD CLOCK INFO TO METADATA-------####################################
##############################################################################################################
metadata_file <-"/Users/maurertm/Desktop/Projects/ADRC/Metadata/Plasma_metadata_FINAL_052021_ADRC_additionalQC_new.csv"
metadata <- read_csv(metadata_file)

metadata_proteomics <- left_join(metadata, proteomics)


##############################################################################################################
################################-------GET LIST WITH ALL GENOTYPED IDS-------#################################
##############################################################################################################
sample_list <- "//Users/maurertm/Desktop/Projects/ADRC/file_lists/files.txt"
sample_list2 <- "/Users/maurertm/Desktop/Projects/ADRC/file_lists/files_R.tsv"


samples<- read_tsv(sample_list, col_names=FALSE)
samples_unique <- samples[!duplicated(samples[ , "X1"]), ] 

samples2<- read_tsv(sample_list2, col_names=FALSE)
samples2 <- samples2[2:502, 2]
samples_unique2 <- samples2 %>% unique()
colnames(samples_unique2) <- "X1"

setdiff(samples_unique, samples_unique2) %>% unique()

#remove the underscore
#samples_unique$X2 <- str_remove_all(samples_unique$X1, "^0")
samples_no_underscore <- samples_unique %>% mutate(ID = sub("\\_.*", "", X1)) %>% dplyr::select(ID)

#remove the 0 before numbers
samples_no_underscore$ID_nospace <- str_remove_all(samples_no_underscore$ID, "^0")

genotypes <- dplyr::select(samples_no_underscore, ID_nospace)
colnames(genotypes) <- "ID"

metadata_proteomics$PIDN_cleaned <- str_remove_all(metadata_proteomics$PIDN, "^0") %>%
  str_remove_all(" .*$") %>% str_remove_all("_.*$") %>% str_remove_all("[A-Z]")

metadata_proteomics$SampleID_cleaned <- str_remove_all(metadata_proteomics$SampleID, "^0") %>%
  str_remove_all(" .*$") %>% str_remove_all("_.*$") %>% str_remove_all("[A-Z]")

metadata_proteomics$ADRC_ID_cleaned <- str_remove_all(metadata_proteomics$ADRC_ID, "^0") %>%
  str_remove_all(" .*$") %>% str_remove_all("_.*$") %>% str_remove_all("[A-Z]")

metadata_ids$match <- as.numeric(metadata_ids$ID %in% genotypes_filtered_list_unnammed)

metadata_proteomics_no_plus <- filter(metadata_proteomics, !str_detect(PIDN_cleaned, "\\+"))
metadata_proteomics_plus <- filter(metadata_proteomics, str_detect(PIDN_cleaned, "\\+"))

fix_id <- function(row, row_df) {
  duplicated <- bind_rows(row, row)
  fixed_pidn <- str_split_fixed(duplicated$PIDN_cleaned[1], "\\+", 2)
  fixed_sampleid <- str_split_fixed(duplicated$SampleID_cleaned[1], "\\+", 2)
  duplicated$PIDN_cleaned <- as.vector(fixed_pidn)
  duplicated$SampleID_cleaned <- as.vector(fixed_pidn)
  duplicated
}

metadata_proteomics_fixed <- rowwise(metadata_proteomics_plus) %>% group_map(fix_id) %>% bind_rows
metadata_proteomics_recombined <- bind_rows(metadata_proteomics_no_plus, metadata_proteomics_fixed)

metadata_proteomics_recombined$PIDN_files <- metadata_proteomics_recombined$PIDN_cleaned
metadata_proteomics_recombined$PIDN_files[!is_in(metadata_proteomics_recombined$PIDN_files, genotypes$ID)] <- NA

metadata_proteomics_recombined$SampleID_files <- metadata_proteomics_recombined$SampleID_cleaned
metadata_proteomics_recombined$SampleID_files[!is_in(metadata_proteomics_recombined$SampleID_files, genotypes$ID)] <- NA

metadata_proteomics_recombined$ADRC_ID_files <- metadata_proteomics_recombined$ADRC_ID_cleaned
metadata_proteomics_recombined$ADRC_ID_files[!is_in(metadata_proteomics_recombined$ADRC_ID_files, genotypes$ID)] <- NA

get_unique_value <- function(row, row_df) {
  unique_ids <- c(row$PIDN_files, row$SampleID_files, row$ADRC_ID_files) %>% na.omit %>% unique
  row$filename <- str_c(unique_ids, collapse = ";")
  row
}

metadata_proteomics_filenames <- rowwise(metadata_proteomics_recombined) %>% group_map(get_unique_value) %>% bind_rows
metadata_proteomics_filenames$filename[nchar(metadata_proteomics_filenames$filename) == 0] <- NA
filter(metadata_proteomics_filenames, str_detect(unique_id, ";"))

tissues <- filter(metadata_proteomics_filenames, !is.na(Esophagus)) 
tissues_tally <- group_by(tissues, Individual_ID) %>% tally
#nine have 2- which ones do we use?

#metadata_proteomics_filenames this is the version that has all of the information
#we are now going to filter out the rows that do not have tissue data

#what columns do we care about- how do we combine things such as dx that change over multiple visits?
group_and_concat <- annotated_metadata %>%
  dplyr::select(match, SCG_ID, SampleId,	Individual_ID,	Age,	ADRC_ID,PIDN,	Visit,	Gender,	Diagnosis_group, Diagnosis_consensus,Date.of.draw,	Study) %>% 
  group_by(Individual_ID) %>%
  mutate(SampleIds = paste(SampleId, collapse = ", ")) %>%
  mutate(Ages = paste(Age, collapse = ", ")) %>%
  mutate(ADRC_IDs = paste(ADRC_ID, collapse = ", ")) %>% 
  mutate(PIDNs = paste(PIDN, collapse = ", ")) %>% 
  mutate(Visits = paste(Visit, collapse = ", ")) %>% 
  mutate(Diagnoses_consensus = paste(Diagnosis_consensus, collapse = ", ")) %>% 
  mutate(Diagnosis_groups = paste(Diagnosis_group, collapse = ", ")) %>% 
  mutate(Date.of.draws = paste(Date.of.draw, collapse = " | ")) %>% 
  mutate(Studies = paste(Study, collapse = ", ")) %>% 
  mutate(Ages = paste(SCG_ID, collapse = ", ")) %>%
  mutate(matches = paste(match, collapse = "")) 


group_and_concats <- group_and_concat %>% dplyr::select(matches, Individual_ID, SampleIds, ADRC_IDs, PIDNs, Ages, Gender, Visits, Diagnosis_group, Diagnoses_consensus, Date.of.draws, Studies)


##############################################################################################################
####################################-------FIND DIFF BTW THE LISTS-------#####################################
##############################################################################################################
individual_id <- setdiff(genotypes$ID, metadata_proteomics_filenames$filename)
#"415" "424" "465" "474" "563" "565" "583" "584" "610" -> these all belong to the "IndividualId?
missing_genotypes <- c(415, 424, 465, 474, 563, 565, 583, 584, 610)

genotypes_individual <- subset(genotypes, ID %in% individual_id)
genotypes_filtered <- genotypes %>% filter(!ID %in% missing_genotypes) 
genotypes_filtered_list <- dplyr::select(genotypes_filtered, ID) %>% as.list()
genotypes_filtered_list_unnammed <- unname(genotypes_filtered_list) %>% unlist()

setdiff(genotypes_filtered, metadata_id_df)
#none, which is what we expect


##############################################################################################################
####################################-------ADD COLUMN TO METADATA------#######################################
####################################-------SAYS IF IS IN GENOTYPES------######################################
####################################-------DOES NOT INCLUDE INDIV ID------#####################################
##############################################################################################################
write_tsv(data.frame(annotated_metadata), "/Users/maurertm/Desktop/Projects/ADRC/annotated_metadata_April_7.tsv")





CHIP_calls_file <- "/Users/maurertm/Downloads/mutect_aggregated_noWL_Apr_5.tsv"
CHIP_calls <-  read_tsv(CHIP_calls_file)


CHIP_calls$filename <- str_remove_all(CHIP_calls$Sample, "^0") %>% str_remove_all("\\.hg38") %>%
  str_remove_all(" .*$") %>% str_remove_all("_.*$") %>% str_remove_all("[A-Z]") 

CHIP_calls_filtered <- filter(CHIP_calls, t_ref_count > 0 & 
                              tumor_f <= 0.35 & longest_repeat < 5 &
                              is_in(vcf_FILTER, c("PASS", "weak_evidence", "clustered_events", "clustered_events;weak_evidence")))

CHIP_proteomics <- left_join(metadata_proteomics_filenames, CHIP_calls)

