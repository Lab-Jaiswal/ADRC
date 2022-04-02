library(writexl)

library(biomaRt)
library(rtracklayer)
library(GenomicRanges)
library(readxl)
library(magrittr)
library(tidyverse)
library(stringi)

##############################################################################################################
################################-------GET COLUMN WITH ALL METADATA IDS-------################################
##############################################################################################################
metadata_file <-"/Users/maurertm/Desktop/Projects/ADRC/Metadata/Plasma_metadata_FINAL_052021_ADRC_additionalQC(2).csv"
metadata_file <-"/Users/maurertm/Desktop/Plasma_metadata_FINAL_052021_ADRC_additionalQCs.csv"
metadata_file <- "/Users/maurertm/Downloads/updated_metadata_file.csv"

metadata <- read_csv(metadata_file)
metadata$index <- rownames(metadata)

metadata_sample <- metadata %>% 
  dplyr::select(SampleId, ADRC_ID, PIDN) %>% 
  mutate(SampleId_nospace = sub(" .*", "", SampleId)) 

#1. three dataframes for the three different types of IDs
metadata_sample$index <- rownames(metadata)
ADRC_ID <- metadata_sample %>% dplyr::select(ADRC_ID, index)
PIDN_ID <- metadata_sample %>% dplyr::select(PIDN, index)
SI <- metadata_sample %>% dplyr::select(SampleId, index)

#2. Clean the data frames
#SI- remove the space after the numbers and put into a seperate column
SI_cleaned1 <- separate(SI, SampleId, into = c("SampleId", "visit_details"), sep = " (?=[^ ]+$)")
SI_cleaned2 <- separate(SI_cleaned1, SampleId, into = c("SampleId", "visit_details2"), sep = " (?=[^ ]+$)")
SI_cleaned2$details_visit <- str_c(SI_cleaned2$visit_details, " ", SI_cleaned2$visit_details2)
SI_cleaned2$index_val <- str_c(SI_cleaned2$index, "_", "SI")
#remove the 0 before the numbers
SI_cleaned2$ID_nozero <- str_remove_all(SI_cleaned2$SampleId, "^0")
SI_cleaned2 <- SI_cleaned2[order(SI_cleaned2$SampleId), ]
SI_cleaned_final <- SI_cleaned2 %>% dplyr::select(ID_nozero, index, index_val, details_visit)

colnames(SI_cleaned_final) <- c("ID", "index", "index_val", "details_visit")


#ADRC_ID- no cleaning needed, except add key to index and details visit column
ADRC_ID_cleaned_final <- ADRC_ID
ADRC_ID_cleaned_final$index_val <- str_c(ADRC_ID_cleaned_final$index, "_", "ADRC")
ADRC_ID_cleaned_final$details_visit <- NA
colnames(ADRC_ID_cleaned_final) <- c("ID", "index", "index_val", "details_visit")

#PIDN_ID
#I need to convert all of the values that can be numbers to numbers
PIDN_numbers <- filter(PIDN_ID, !is.na(as.numeric(PIDN)))
PIDN_letters <- filter(PIDN_ID, is.na(as.numeric(PIDN)))

PIDN_numbers[, 1] <- sapply(PIDN_numbers[, 1], as.numeric)

PIDN_cleaned_final <- rbind(PIDN_numbers, PIDN_letters)
PIDN_cleaned_final$index_val <- str_c(PIDN_cleaned_final$index, "_", "PIDN")
PIDN_cleaned_final$details_visit <- NA
colnames(PIDN_cleaned_final) <- c("ID", "index", "index_val", "details_visit")


#3. Combine all three dataframes 
metadata_ids <- rbind(PIDN_cleaned_final, ADRC_ID_cleaned_final, SI_cleaned_final)
metadata_id_df <- dplyr::select(metadata_ids, ID)

##############################################################################################################
################################-------GET LIST WITH ALL METADATA IDS-------################################
##############################################################################################################
metadata_id_list <- dplyr::select(metadata_ids, ID) %>% as.list()
metadata_id_unnammed <- unname(metadata_id_list) %>% unlist()


##############################################################################################################
################################-------GET LIST WITH ALL GENOTYPED IDS-------#################################
##############################################################################################################
#sample_list <- "/Users/maurertm/Downloads/files.txt"
sample_list2 <- "/Users/maurertm/Downloads/files_R.tsv"

samples<- read_tsv(sample_list, col_names=FALSE)
samples_unique <- samples[!duplicated(samples[ , "X1"]), ] 

samples2<- read_tsv(sample_list2, col_names=FALSE)
samples2 <- samples2[2:502, 2]
samples_unique2 <- samples2 %>% unique()
colnames(samples_unique2) <- "X1"

setdiff(samples_unique, samples_unique2)

#remove the underscore
#samples_unique$X2 <- str_remove_all(samples_unique$X1, "^0")
samples_no_underscore <- samples_unique %>% mutate(ID = sub("\\_.*", "", X1)) %>% dplyr::select(ID)

#remove the 0 before numbers
samples_no_underscore$ID_nospace <- str_remove_all(samples_no_underscore$ID, "^0")

genotypes <- dplyr::select(samples_no_underscore, ID_nospace)
colnames(genotypes) <- "ID"

##############################################################################################################
####################################-------FIND DIFF BTW THE LISTS-------#####################################
##############################################################################################################
individual_id <- setdiff(genotypes, metadata_id_df)
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
metadata_ids$match <- as.numeric(metadata_ids$ID %in% genotypes_filtered_list_unnammed)

PIDN_matched <- metadata_ids %>% filter(match == 1 & grepl("PIDN", index_val))
PIDN_matched <- PIDN_matched %>% 
  mutate(PIDN_ID = ifelse(!grepl("3900", ID), paste0("0", ID), ID))
PIDN_no_matched <- metadata_ids %>% filter(match == 0 & grepl("PIDN", index_val))
PIDN_no_matched$PIDN_ID <- NA
PIDN_combined <- rbind(PIDN_matched, PIDN_no_match) %>% select(index, PIDN_ID)


ADRC_matched <- metadata_ids %>% filter(match == 1 & grepl("ADRC", index_val))
ADRC_matched <- ADRC_matched %>% 
  mutate(ADRC_ID = ifelse(!grepl("3900", ID), paste0("0", ID), ID))
ADRC_no_matched <- metadata_ids %>% filter(match == 0 & grepl("ADRC", index_val))
ADRC_no_matched$ADRC_ID <- NA
ADRC_combined <- rbind(ADRC_matched, ADRC_no_matched) %>% select(index, ADRC_ID)


SI_matched <- metadata_ids %>% filter(match == 1 & grepl("SI", index_val))
SI_matched <- SI_matched %>% 
  mutate(SI_ID = ifelse(!grepl("3900", ID), paste0("0", ID), ID))
SI_no_matched <- metadata_ids %>% filter(match == 0 & grepl("SI", index_val))
SI_no_matched$SI_ID <- NA
SI_combined <- rbind(SI_matched, SI_no_matched) %>% select(index, SI_ID)

PIDN_ADRC <- merge(PIDN_combined, ADRC_combined, by="index")

complete_ID <- merge(PIDN_ADRC, SI_combined, by="index")
#complete_ID[is.na(complete_ID)] <- 0

grouped_metadata_index <- aggregate(match ~ index, metadata_index, sum)

cols <- c('PIDN_ID' , 'ADRC_ID' , 'SI_ID')
complete_ID$SCG_IDs<- apply( complete_ID[ , cols ] , 1 , paste , collapse = " " )
complete_merged<- complete_ID %>% 
  mutate(new_id = ifelse(grepl("[1-9]", SCG_IDs), gsub("NA","",as.character(complete_ID$SCG_IDs)), NA)) %>% 
  mutate(sequenced = ifelse(grepl("[1-9]", new_id), "Yes", "No")) %>% 
  mutate(SCG_ID = new_id) %>%
  select(index, SCG_ID, sequenced) %>%
  merge(grouped_metadata_index, by = "index")


annotated_metadata <- merge(complete_merged, metadata, by = "index")

write_tsv(data.frame(annotated_metadata), "/Users/maurertm/Desktop/Projects/ADRC/annotated_metadata.tsv")

group_and_concat <- annotated_metadata %>%
  dplyr::select(match, SampleId,	Individual_ID,	Age,	ADRC_ID,PIDN,	Visit,	Gender,	Diagnosis_group, Diagnosis_consensus,Date.of.draw,	Study) %>% 
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

group_and_concats$genotyped <- ifelse(grepl("1|2", group_and_concats$matches), "yes", "no")
group_and_concats <- mutate(group_and_concats, Num_Visits = str_count(Visits, ',') + 1)

nrow(unique(group_and_concats))
unique_grouped <- unique(group_and_concats)

group_and_concats <- group_and_concat %>% dplyr::select(Individual_ID, SampleIds, ADRC_IDs, PIDNs, genotyped, Ages, Gender, Num_Visits, Diagnoses_consensus, Date.of.draw, Studies)
write_tsv(data.frame(unique_grouped), "/Users/maurertm/Desktop/Projects/ADRC/annotated_grouped_metadata.tsv")

C_numeric <- lapply(samples_ADRC, as.numeric)
samples_ADRC_numeric_unnammed <- unname(samples_ADRC_numeric) %>% unlist()


