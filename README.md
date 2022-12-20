# ADRC Analysis
## Data Sources
* Once you have been approved to recieve the data from the ADRC, @maurermaggie will share the approved files with you <br />
* All of the following files have been de-identified
#### Aging Clock Data
* All have filepaths "ADRC_Project/Source_Data/Aging_Clock_Data/df_all_tissue_bootstrap_lasso_ENRICHED_COGNITIONBRAIN_ABOVE60" + : <br />
     * "_WAGNER_RAW_Age_prediction_mean_results.csv"
     * "_POSTON_RAW_Age_prediction_mean_results.csv"
     * "_KERCHNER_RAW_Age_prediction_mean_results.csv"
     * "_ADRC_RAW_Age_prediction_mean_results.csv"
#### Phenotype Metadata
* Source: Greicius Lab, Wyss-Coray Lab, RedCap
      * File: "ADRC_Project/Source_Data/Metadata/Phenotype/Incomplete_Phenotype_MetaData2.csv"
#### Clinical and Sample Level Metadata
* Source: Montgomery Lab
      * File: "ADRC_Project/Source_Data/Metadata/Clinical_Sample_Level/Plasma_metadata_FINAL_052021_ADRC_additionalQC_new.csv"
#### Filename correction keys
* We recieved WGS files from multiple sources, resulting in two problems: <br />
     1. In some cases, the vcf header name does not match the file name
          * To fix this, the following two files are necesssary: 
               * "ADRC_Project/Source_Data/Metadata/Filename_Correction_Keys/whole_exome_incorrect_headers_list.txt"
               * "ADRC_Project/Source_Data/Metadata/Filename_Correction_Keys/whole_exome_normal_headers_list.txt"
     2. In others, the file name does not match any identifiers in the provided phenotype metadata with using a "key" given from the ADRC
          * "ADRC_Project/Source_Data/Metadata/Filename_Correction_Keys/New_ADRC_Key.xlsx"
#### Mutect2 Results
* "ADRC_Project/Source_Data/mutect2_Results/mutect_somatic_042822.csv"
* "ADRC_Project/Source_Data/mutect2_Results/Chip_non_Chip_Lists/chip_carriers.txt"
* "ADRC_Project/Source_Data/mutect2_Results/Chip_non_Chip_Lists/non_chip_carriers.txt"
#### Combined dataframe containing Metadata (Phenotype, Clinical, Sample Level), Aging Clock Scores, and Mutect2 results
* created using `Phenotype_CHIP_Clocks_Cleaning.R`
* "ADRC_Project/Processed_DataFrames/Chip_Phen_Clocks_12_19_at_12_54.csv"
#### Blood Plasma Proteome Data
* "ADRC_Project/Source_Data/Metadata/Blood_Plasma_Data/dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021.csv"
#### CSF Plasma Proteome Data
* "ADRC_Project/Source_Data/Metadata/CSV_Plasma_Data/datProt_SS-205083_v4_CSF.hybNorm.medNormInt.plateScale.calibrate.medNormRefSMP.csv"


## Data Processing
### ID Mishaps 
* Because of the multi-study nature of this project, people are connected to multiple IDs. 
    * Some IDs are *individual specific*
    * Some IDs are *study specific*
    * Some IDs, due to incorrect clinical reporting (and the fact that people were enrolled in multiple studies at once), are not sample or study specific

| Sample_ID1        | Sample_ID2           | Sample_ID3  | Study | 
|  :---:            | :---:                | :---:       | :---: | 
| 100               | 0A50                 | 100         | Alpha |
| 100               | 0B75                 | 7776        | Beta  |
| 200               | 0A51                 | NA          | Alpha |
| 200               | 0B76                 | 8998        | Beta  |

Here, Sample_ID1 is *individual specific*, Sample_ID2 is *study specific*, and Sample_ID3 is neither. Unfortunately, most of the IDs given resembled Sample_ID3 and none adequately captured all of the study participants. Because of this, *I reccomend using the column labled `filename`*, an *individual level ID* I created in `Phenotype_CHIP_Clocks_Cleaning.R` through an extensive ID matching and filtering process. 
* I included the other IDs incase you want to used them; however, if you do please examine my logic in `Phenotype_CHIP_Clocks_Cleaning.R`
* You will notice there are 232 rows in Chip_Phen_Clocks_12_19_at_12_54.csv with NA as filename. Unfortunately that is becuase there are 232 rows, likley  representing ~100 people whose aging clock data cannot be linked to the WGS data provided by the ADRC. I have spoken to the ADRC and our collaborators in the Montgomery, Greicius, and Wyss-Coray labs, combed through the ADRC's RedCap stores, and we have come to the conclusion that those 100 people never recieved whole genome sequencing. If you have any ideas concerning their identity, please let us know!

### Subsetting Data
#### If you are looking at patterns between the WGS and plasma proteome data, we reccomend filtering the data in the following ways:
1. Filter out all samples were the `Age_at_Draw_Difference`, a column with the difference between the age at which the proteome sample and WGS sample were drawn, is greater than 4 or less than -4. This is in line with other aging studies. 
2. Filter out all samples where `has_WGS_AND_clock_data` is FALSE (this column has TRUE or FALSE for all rows and represents if the sample has both WGS and clock data
3. Many analyses will require you to subset the data to one blood draw per person. 
     * Using the find_closest_age function in `process_final_dataframe.R`, will allow you to select the sample for each individual that is closest by Age_at_Draw_Difference and, if tied, then by connectivity z score. 
     
### Data Oddities
* The Diagnosis_Group and Diagnosis_Consensus columns might not match each other within the same row
* The Diagnosis_Group is not static and changes over time for patients
* Some patients have multiple sexes recorded. We are currently working on calculating genetic sex from the WGS data.
