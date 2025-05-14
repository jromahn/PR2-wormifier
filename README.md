# PR2-wormifier

## Introduction

## Input

As Species List the pipeline needs a comma-separated table with the columns: "genus","species","taxon".
The taxon column is important to control if the downloaded and found species really belongs to the same taxonomic group, which is essential to check due to non-unique genus names. The taxonomic group added to the column should fit WoRMS taxonomy. If you prefer to not define the taxonomy, write NA.

As example file check: Formated_species_list.csv 
This is a species list for the Baltic Sea and was used to create the reference database.

IMPORTANT:

You need to have access to Algaebase database to run the script, you need an API key which you can request with the help of the following homepage: [Algaebase API link](https://www.algaebase.org/api/)
The key, the email address for NCBI requests and the filename and path of the species list can be saved in ```00_login_data.ini```. Therefore, adaption of the signle R scripts are not neceassary.

## Dependencies

To run the pipeline, you need to have access to Algaebase. For this write a line to:

The pipeline is tested for the following versions: 


How to install all packages
```
options(repos = c(CRAN = "https://cloud.r-project.org/"))
# dependencies for pr2database:
 install.packages("devtools")
 install.packages("BiocManager")
 BiocManager::install("Biostrings")
devtools::install_github("pr2database/pr2database")

#the rest
instal.packages(c("tidyverse","readxl","stringi","ini","worrms","taxonomizr","phylotools","jsonlite","curl","utils","rentrez","gapminder","treemapify","ggpubr"))

```
Essential R packages: 
```
require(tidyverse)
require(stringi) # try next time without
require(readxl)
require(ini)

#databases and taxonomy
require(worrms)
require(pr2database)
require(taxonomizr)


#fasta file handling
require(phylotools)

# extract from the internet
library(jsonlite)
library(curl)
require(utils)
require(rentrez)

# for plots
require(gapminder)
require(treemapify)
require(ggpubr)


```


## Scripts

- 00_Function_Library.R
- 01_Clean_PR2.R
- 02_Clean_Input.R
- 03_Download_from_PR2__modJR.R
- 04_Identify_missing_species__modJR.R
- 05_Check_NCBI_with_biopython__advanced.py
- 06_Clean_Sort_NCBI_Downloads__modJR.R
- 06_PrepareDatabase.R
- 07_Download_Rest_PR2_modJR.R
- 08_Combine_Sort_Files__modJR.R
- 09_Taxonomy__modJR.R
- 10_Sort_Fasta__modJR.R
- 11_Create_Database_Versions.R
- 12_Convert_header_Crabs.R
-062024_bash_crabs.sh

## Benchmark Scripts

- 13_Assignment_Stats_Crabs.R
- 14_Treemaps_db_results.R
- 15_DB_Stats.R

## Metadata:

| Taxonomy | |
| ------------- | ------------- |
| 1 | WoRMS, AphiaID in PR2 |
| 2 | WoRMS, without AphiaID in PR2, but a sequence of the same genus has|
| 3 | WoRMS, without AphiaID of any sequence of that genus in PR |
| 4 | WoRMS, without AphiaID and non unique genus name|
| 5 | WoRMS, not 1-4|
| 6 | Algaebase |

## References
