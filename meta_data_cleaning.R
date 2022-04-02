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

library("xlsx")
write_xlsx(unique_grouped, "/Users/maurertm/Desktop/Projects/ADRC/annotated_grouped_metadata.xlsx")
matched <- annotated_metadata %>% filter(match > 0)

matched$Individual_ID %>% unique() %>% length()
#489
samples_unique$X1 %>% length()
#497
##############################################################################################################
####################################-------ADD COLUMN TO METADATA------#######################################
####################################-------IF IN INDIV ID------############################################3##
##############################################################################################################
#metadata_indiv <- metadata %>% 
#  dplyr::select(Individual_ID) 

#metadata_indiv$index <- rownames(metadata_indiv)
#metadata_indiv_no_ind <- metadata_indiv %>% mutate_all(~gsub("Ind", "", .))
#metadata_indiv_no_ind$details_visit <- NA
#metadata_indiv_no_ind$index_val <- str_c(metadata_indiv_no_ind$index, "_", "Indiv")
#colnames(metadata_indiv_no_ind) <- c("ID", "index", "details_visit", "index_val")
#metadata_indiv_no_ind <- dplyr::select(metadata_indiv_no_ind, ID, index, index_val, details_visit)


#metadata_ids_new <- rbind(ADRC_ID_cleaned_final, metadata_indiv_no_ind)
#metadata_ids_test <- dplyr::select(metadata_ids_new, ID)
#setdiff(genotypes_filtered, metadata_ids_test)

#remove ind from the id


#metadata_indiv_no_index <- dplyr::select(metadata_indiv_no_ind, -index)
#colnames(metadata_indiv_no_index) <- "ID"

#setdiff(genotypes, metadata_indiv_no_index)
#intersect(genotypes_filtered, metadata_indiv_no_index)
#does the 0 actually mean something

#metadata_sample <- metadata %>% 
 # dplyr::select(Individual_ID) %>% 
  #mutate(SampleId_nospace = sub(" .*", "", SampleId)) 

##############################################################################################################
####################################-------TESTING STUFF------#######################################
####################################-------IF IN INDIV ID------############################################3##
##############################################################################################################
matched <- annotated_metadata %>% filter(match==1)
matched_ADRC <- filter(matched, nchar(ADRC_ID) > 4, grepl("^3900",ADRC_ID))
ADRC_PIDN <- matched_ADRC %>% dplyr::select(PIDN)
#405
PIDN_list <- c(416, 423, 428, 449, 447, 460,463, 475, 456, 518, 530, 482, 527, 470, 560, 448, 550, 503, 468)
PIND_not_ADRC_matched <- subset(annotated_metadata, PIDN %in% PIDN_list)
#all of these have ADRCID=NA
individual_id <- c("Ind415", "Ind424", "Ind465", "Ind474", "Ind563", "Ind565", "Ind583", "Ind584", "Ind610")
IndividualIDDF <- subset(annotated_metadata, Individual_ID %in% individual_id)
IndividualIDDF$PIDN
#"1179" "1180" "1159" "990"  "1056" "999"  "1043" "1208" "1157" "990" 
IndividualIDDF$ADRC_ID
#3900341 3900342 3900323 3900217 3900255      NA      NA 3900359 3900321 3900217

filter(matched_ADRC, ADRC_ID==3900341) %>% dplyr::select(ADRC_ID, PIDN, SampleID, index)
#  ADRC_ID PIDN SampleID index
# 3900341 1179     1179   132
filter(matched_ADRC, ADRC_ID==3900342) %>% dplyr::select(ADRC_ID, PIDN, SampleID, index)
#ADRC_ID PIDN SampleID index
#1 3900342 1180     1180   133
filter(matched_ADRC, ADRC_ID==3900323) %>% dplyr::select(ADRC_ID, PIDN, SampleID, index)
#not present
filter(matched_ADRC, ADRC_ID==3900217) %>% dplyr::select(ADRC_ID, PIDN, SampleID, index)
#  ADRC_ID PIDN SampleID index
#1 3900217  990      990   292
#2 3900217  990  0990 Y2    97
filter(matched_ADRC, ADRC_ID==3900255) %>% dplyr::select(ADRC_ID, PIDN, SampleID, index)
ADRC_ID PIDN SampleID index
#1 3900255 1056     1056   455
filter(matched_ADRC, ADRC_ID==3900359) %>% dplyr::select(ADRC_ID, PIDN, SampleID, index)
#not present
filter(matched_ADRC, ADRC_ID==3900321) %>% dplyr::select(ADRC_ID, PIDN, SampleID, index)
ADRC_ID PIDN SampleID index
#1 3900321 1157     1157    95



# "1159"  "999"  "1043" "1208" 
 filter(metadata_annotated, PIDN==1159) %>% dplyr::select(ADRC_ID, PIDN, SampleID, index)
 #ADRC_ID PIDN SampleID index
 #1 3900323 1159     1159   168
 #2 3900323 1159     1159   167
 #3 3900323 1159     1159  1029
 
 filter(metadata_annotated, PIDN==999) %>% dplyr::select(ADRC_ID, PIDN, SampleID, index)
# ADRC_ID PIDN SampleID index
# 1      NA  999      999   714
# 2      NA  999      999   685
# 3      NA  999      999   983
 
 filter(metadata_annotated, PIDN==1043) %>% dplyr::select(ADRC_ID, PIDN, SampleID, index)
 #ADRC_ID PIDN SampleID index
 #1      NA 1043     1043   722
 #2      NA 1043     1043   693
 #3      NA 1043     1043  1007
 
 filter(metadata_annotated, PIDN==1208) %>% dplyr::select(ADRC_ID, PIDN, SampleID, index)
# ADRC_ID PIDN SampleID index
# 1 3900359 1208     1208   958
# 2 3900359 1208     1208   900
# 3 3900359 1208     1208   251
 
matched_PIDN <- matched %>% dplyr::select(PIDN)
setdiff(ADRC_PIDN, matched_PIDN)
setdiff(matched_PIDN, ADRC_PIDN)

grep("415", metadata_id_unnammed, value=TRUE)
#3900415
grep("424", metadata_id_unnammed, value=TRUE)
#"C424"    "3900424"
grep("465", metadata_id_unnammed, value=TRUE)
character(0)
grep("474", metadata_id_unnammed, value=TRUE)
character(0)
grep("563", metadata_id_unnammed, value=TRUE)
character(0)
grep("565", metadata_id_unnammed, value=TRUE)
character(0)
grep("583", metadata_id_unnammed, value=TRUE)
character(0)
grep("584", metadata_id_unnammed, value=TRUE)
character(0)
grep("610", metadata_id_unnammed, value=TRUE)
character(0)



samples_ADRC <- filter(samples_no_underscore, nchar(ID) > 4, grepl("^3900",ID)) %>% dplyr::select(ID) %>% as.list()
samples_ADRC_numeric <- lapply(samples_ADRC, as.numeric)
samples_ADRC_numeric_unnammed <- unname(samples_ADRC_numeric) %>% unlist()









#2. create a new column with the extra info removed from SampleId (everything after the space and leading zeros)
#the ones with pluses are removed because those are not removed in the PI
metadata_plus <- metadata_nospace %>% filter(str_detect(SampleId_nospace, "\\+"))
metadata_noplus <- metadata_nospace[-c(109, 495, 514, 996), ]

metadata_nospace <- metadata_sample %>% mutate(SampleId_nospace = sub("\\ .*", "", SampleId_nospace)) 
metadata_nospace$SampleId_nozero <- str_remove_all(metadata_nospace$SampleId_nospace, "^0")

#find cases where the PIDN and SampleId are different

metadata_PIDN_SI <- filter(metadata_nospace, PIDN == SampleId_nozero)

metadata_no_PIDN_SI <- filter(metadata_nospace, PIDN != SampleId_nozero)

metadata_no_PIDN_SI$id <- filter(metadata_nospace, PIDN != SampleId_nozero)

metadata_combined <- rbind(metadata_noplus, metadata_plus_before, metadata_plus_after)




metadata_sample$index <- rownames(metadata_sample)


metadata_combined_NA_ADRC_ID <- filter(metadata_combined, is.na(ADRC_ID))
metadata_combined_NA_ADRC_ID$combined_id <- str_c("NA", "_", metadata_combined_NA_ADRC_ID$SampleId_nospace)


metadata_combined_no_NA <- filter(metadata_combined, !is.na(ADRC_ID))
metadata_combined_Sample_ADRC_Equal <- filter(metadata_combined_no_NA, ADRC_ID == SampleId_nospace)
metadata_combined_Sample_ADRC_Equal$combined_id <- str_c(metadata_combined_Sample_ADRC_Equal$ADRC_ID, "_", "NA")

metadata_combined_no_NA_equals <- filter(metadata_combined_no_NA, ADRC_ID != SampleId_nospace)
metadata_combined_no_NA_equals$combined_id <- str_c(metadata_combined_no_NA_equals$ADRC_ID, "_", metadata_combined_no_NA_equals$SampleId_nospace)

metadata <- rbind(metadata_combined_NA_ADRC_ID, metadata_combined_Sample_ADRC_Equal, metadata_combined_no_NA_equals)










sample_list <- "/Users/maurertm/Downloads/sample_list2"
samples<- read_tsv(sample_list, col_names=FALSE)
samples_unique <- samples[!duplicated(samples[ , "X1"]), ] 
#samples_unique$X2 <- str_remove_all(samples_unique$X1, "^0")
samples_no_underscore <- samples_unique %>% mutate(ID = sub("\\_.*", "", X1)) %>% dplyr::select(ID)

samples_ADRC <- filter(samples_no_underscore, nchar(ID) > 4, grepl("^3900",ID)) %>% dplyr::select(ID) %>% as.list()
samples_ADRC_numeric <- lapply(samples_ADRC, as.numeric)
samples_ADRC_numeric_unnammed <- unname(samples_ADRC_numeric) %>% unlist()


ADRC_list <- dplyr::select(metadata, ADRC_ID) %>% as.list()
ADRC_list_unnammed <- unname(ADRC_list) %>% unlist()


ADRC_shared <- intersect(samples_ADRC_numeric_unnammed, ADRC_list_unnammed)
setdiff(ADRC_list_unnammed, samples_ADRC_numeric_unnammed) 


samples_non_ADRC <- filter(samples_no_underscore, !grepl("^3900",ID)) %>% as.list()
samples_nonADRC_numeric <- lapply(samples_non_ADRC, as.numeric)
samples_nonADRC_numeric_unnammed <- unname(samples_nonADRC_numeric) %>% unlist()

non_ADRC_list <- dplyr::select(metadata, SampleId_nospace) %>% as.list()
non_ADRC_list_unnammed <- unname(non_ADRC_list) %>% unlist() %>% sort()
nonADRC_numeric <- lapply(non_ADRC_list_unnammed, as.numeric) %>% unlist()
#samples_unique$X2 <- str_remove_all(samples_unique$X1, "^0")

non_ADRC_shared <- intersect(samples_nonADRC_numeric_unnammed, non_ADRC_list_unnammed)
setdiff(samples_nonADRC_numeric_unnammed, nonADRC_numeric)

metadata$nonADRC_match <- as.numeric(metadata$SampleId_nospace %in% non_ADRC_shared)

metadata$ADRC_match <- as.numeric(metadata$ADRC_ID %in% ADRC_shared)



metadata$ADRC_ID_match <- samples_ADRC_numeric  %in% ADRC_list




metadata_noADRCmatch <- filter(metadata, ADRC_ID_match==0)
metadata_ADRCmatch <- filter(metadata, ADRC_ID_match==1)


metadata_noADRCmatch$non_ADRC_ID_match <- as.integer(metadata_noADRCmatch$SampleId_nospace %in% samples_non_ADRC$ID)
metadata_noADRCmatch_sampleID_match <- filter(metadata_noADRCmatch, non_ADRC_ID_match==1)


metadata$match <- metadata$non_ADRC_ID_match + metadata$ADRC_ID_match

matching <- filter(metadata, match == 1)


metadata_filtered_ADRC_ID <- filter(metadata, is_in(ADRC_ID, samples$X2))
metadata_filtered_SampleID <- filter(metadata, is_in(SampleID, samples$X2))

metadata_filtered <- bind_rows(metadata_filtered_ADRC_ID, metadata_filtered_SampleID)

metadata_baseline_filtered <- filter(metadata_filtered, Visit_Consensus == "Baseline")

metadata_baseline_filtered %>% nrow()

setdiff(unique(metadata_filtered$ADRC_ID), unique(metadata_baseline_filtered$ADRC_ID))
#3900330 3900028 3900350-- all three are missing baseline visits

metadata_unsequenced <- filter(metadata, !is_in(ADRC_ID, samples$X2) & !is_in(SampleID, samples$X2))
metadata_sequenced <- filter(metadata, is_in(ADRC_ID, samples$X2) | is_in(SampleID, samples$X2))
metadata_unsequenced_ADRC_ID <- filter(metadata_unsequenced, !is.na(ADRC_ID) & !duplicated(ADRC_ID))
metadata_unsequenced_SampleID <- filter(metadata_unsequenced, is.na(ADRC_ID) & !duplicated(SampleID) & !str_detect(SampleID, "Wave"))



setdiff(unique(metadata_filtered$SampleID), unique(metadata_baseline_filtered$SampleID))
setdiff(samples$X1, metadata$SampleID)
col_names = FALSE
name <- "/Users/maurertm/Downloads/sample_list.tsv"
write_tsv(data.frame(sampels), name)

write_tsv()
