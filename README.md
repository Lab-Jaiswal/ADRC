# ADRC Analysis
## Data Sources
* Once you have been approved to recieve the data from the ADRC, @maurermaggie will share the approved files with you <br />
* All of the following files have been de-identified
#### Aging Clock Data
* All have filepaths "df_all_tissue_bootstrap_lasso_ENRICHED_COGNITIONBRAIN_ABOVE60" + : <br />
     * "_WAGNER_RAW_Age_prediction_mean_results.csv"
     * "_POSTON_RAW_Age_prediction_mean_results.csv"
     * "_KERCHNER_RAW_Age_prediction_mean_results.csv"
     * "_ADRC_RAW_Age_prediction_mean_results.csv"
#### Phenotype Metadata
* Source: Greicius Lab, RedCap
      * File: "Greicius_MetaData2.csv"
#### Clinical and Sample Level Metadata
* Source: Montgomery Lab
      * File: "Plasma_metadata_FINAL_052021_ADRC_additionalQC_new.csv"
#### Mutect2 Results
* mutect_somatic_042822.csv
#### Blood Plasma Proteome Data
* dataProt_SS-205063.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP_ADRC_Feb2021.csv
#### CSF Plasma Proteome Data
* datProt_SS-205083_v4_CSF.hybNorm.medNormInt.plateScale.calibrate.medNormRefSMP.csv
#### Filename corrections
* We recieved WGS files from multiple sources, resulting in two problems: <br />
     1. In some cases, the vcf header name does not match the file name
          * To fix this, the following two files are necesssary: whole_exome_incorrect_headers_list.txt and whole_exome_normal_headers_list.txt
     2. In others, the file name does not match any identifiers in the provided phenotype metadata with using a "key" given from the ADRC
          * New_ADRC_Key.xlsx

## Data Processing
* Because of the multi-study nature of this project, many people have multiple IDs. 
    * Some IDs are *sample specific*
    * Some IDs are *study specific*
    * Some IDs, due to incorrect clinical reporting (and the fact that people were enrolled in multiple studies at once), are not sample or study specific

| Sample_ID1        | Sample_ID2           | Sample_ID3  | Study | 
|  :---:            | :---:                | :---:       | :---: | 
| 100               | 0A50                 | 100         | Alpha |
| 100               | 0B75                 | 7776        | Beta  |
| 200               | 0A51                 | 200         | Alpha |
| 200               | 0B76                 | 8998        | Beta  |

Here, Sample_ID1 is *sample specific*, Sample_ID2 is *study specific*, and Sample_ID3 is neither. Unfortunately, most of the IDs given resembled Sample_ID3 and none adequately captured all of the study participants. Because of this, I reccomend using the column labled `filename`. If you insist on using the other IDs in the CHIP
