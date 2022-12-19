# ADRC Analysis
## Data Sources
Once you have been approved to recieve the data from the ADRC, @maurermaggie will share the approved files with you. We are currently performing analyses on the following files <br />
#### Aging Clock Data
* All have filepaths "df_all_tissue_bootstrap_lasso_ENRICHED_COGNITIONBRAIN_ABOVE60" + : <br />
     * "_WAGNER_RAW_Age_prediction_mean_results.csv"
     * "_POSTON_RAW_Age_prediction_mean_results.csv"
     * "_KERCHNER_RAW_Age_prediction_mean_results.csv"
     * "_ADRC_RAW_Age_prediction_mean_results.csv"
#### Phenotype Metadata
##### WGS Metadata
* Source: Montgomery Lab
      * File: "Plasma_metadata_FINAL_052021_ADRC_additionalQC_new.csv"
* Source: Greicius Lab
      * File: "Greicius_MetaData2.csv"
* Lists of CHIP and non-CHIP carriers:
      * Source: mutect2 analysis
      * Files: "chip_carriers.txt", "non_chip_carriers.txt"
##### Filename corrections
* We recieved WGS files from multiple sources, resulting in two problems:
      1. In some cases, the vcf header name does not match the file name
           * To fix this, the following two files are necesssary: whole_exome_incorrect_headers_list.txt and whole_exome_normal_headers_list.txt
      2. In others, the file name does not match any identifiers in the provided phenotype metadata with using a "key" given from the ADRC
           * New_ADRC_Key.xlsx
#### Mutect2 Results
* mutect_somatic_042822.csv

#### Blood Plasma Proteome Data
* dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021.csv

#### CSF Plasma Proteome Data
* datProt_SS-205083_v4_CSF.hybNorm.medNormInt.plateScale.calibrate.medNormRefSMP.csv
