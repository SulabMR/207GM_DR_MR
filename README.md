## Identification of causal gut microbes for diabetic retinopathy using systematic Mendelian Randomization

**This README contains:**

- project background and aims
- project work/code stored in this repo


### Project Background 

Diabetic retinopathy (DR) is a leading causes of vision loss in adults aged 20–74 years with increasing prevalence in multiple countries worldwide. While increasing evidence has highlighted the importance of gut microbiome in etiology and development of DR, the causal link of gut microbes with DR and DR risk factors remains unknown. In this study, we applied two-sample Mendelian Randomization (MR) to estimate causal effect of 207 gut microbes on DR and assessed mediation effects by DR risk factors, using summary statistics data from genome-wide association studies (GWAS). We found causal association of 8 bacterial taxa with DR (FDR < 0.05), among which genus Collinsella and species Collinsella_aerofaciens were associated with increased risk of DR while species Bacteroides_faecis, Burkholderiales bacterium_1_1_47, Ruminococcus_torques, Streptococcus_salivarius, genus Burkholderiales_noname and family Burkholderiales_noname were protective for DR.


### This repository

Main analysis scripts and metadata (see details below):

```
├── 01_data_clump.R
├── 02_MR_main.R
├── 03_MR_Risk-DR.R
├── 04_MR_GM-Risk.R
├── 05_MVMR.sh
├── 06_my_presso.R
├── functions.R
├── combined_results.R
│   └── Supplementary table 1.xlsx

```



### How to reproduce this analysis, i.e. Main Workflow


1. R script `01_data_clump.R` is required for processing data that comes as text files (i.e. GWAS summary stats) from the IEU GWAS pipeline or other sources. This script has to be run to convert raw data into `outcome` data frames and to extract instruments from each GWAS (in `exposure` format) and save them to be used directly in MR analysis in subsequent scripts. 

2. R script `02_MR_main.R` is performed two-sample MR analysis from `TwoSampleMR` package to evaluate causal relationship between gut microbes and DR by defining gut microbe as exposures and DR as outcome.

3. R script `03_MR_Risk-DR.R` is performed two-sample MR analysis from `TwoSampleMR` package to evaluate causal relationship between risk factors and DR by defining risk factors as exposures and DR as outcome.

4. R script `04_MR_GM-Risk.R` is performed two-sample MR analysis from `TwoSampleMR` package to evaluate causal relationship between gut microbes and risk factors by defining gut microbe as exposures and risk factors as outcome.

5. shell script `05_MVMR.sh` is performed multivariable Mendelian randomization to estimate the effect of DR risk factor on DR (Gut microbes → DR risk factors), adjusting for the effect of corresponding gut microbe.

6. R script `06_my_presso.R` is performed MR-PRESSO test using R package `MRPRESSO`.

7. R script `functions.R` contains supplementary functions to display the results.

7. R script `combined_results.R` displays the results as a table to make it easy to copy values.


_Supplementary data_

`Supplementary table 1.xlsx` is a table that contains information about GWAS datasets.
