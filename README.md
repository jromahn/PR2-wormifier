# PR2-wormifier

The manuscript of this pipeline will be soon submitted to Molecular Ecology Resources. For a detailed information please check the manuscript, will be linked as soon as possible.

## Introduction
THe PR2-wormifier builds a customized, curated 18S reference database by combining sequences from PR² and NCBI with validated taxonomy from WoRMS and AlgaeBase. It filters out low-quality or non-nuclear sequences and standardizes species names and metadata. The result is a more comprehensive database that improves taxonomic resolution in metabarcoding studies.

it consists out of 10 scripts, 9 written in R and 1 written in Python

![ObiMAGIC pipeline overview](PR2-wormifier.png)

## Overview

[Input](#Input)  
[Dependencies](#Dependencies)   
[Pipeline Scripts](#Pipeline-Scripts)  
[Scripts Associated with Benchmarking in folder: 13_Benchmark_Script](#Scripts-Associated-with-Benchmarking-in-folder:-13-Benchmark_Script)
[Metadata](#Metadata)  
[References](#References)  

## Input

As Species List the pipeline needs a comma-separated table with the columns: "genus","species","taxon".
The taxon column is important to control if the downloaded and found species really belongs to the same taxonomic group, which is essential to check due to non-unique genus names. The taxonomic group added to the column should fit WoRMS taxonomy. If you prefer to not define the taxonomy, write NA.

As example file check: `Formated_species_list.csv`
This is a species list for the Baltic Sea and was used to create the reference database.

For benchmarking a sedaDNA metabarcoding dataset was used. It is available under the following DOI in FigShare: [10.6084/m9.figshare.28457489](https://doi.org/10.6084/m9.figshare.28457489)

IMPORTANT:

You need to have access to Algaebase database to run the script, you need an API key which you can request with the help of the following homepage: [Algaebase API link](https://www.algaebase.org/api/)
The key, the email address for NCBI requests and the filename and path of the species list can be saved in ```00_login_data.ini```. Therefore, adaption of the signle R scripts are not neceassary.

## Dependencies

To run the pipeline, you need to have access to Algaebase. For this write a line to:

The pipeline is tested for the following versions: 
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

###  Versions used

All R scripts were conducted in R version 4.5.0. The following R packages were used: `pr2database` (v5.1.0; Vaulot, 2025), `dplyr` (v1.1.4; Wickham et al., 2023), `tidyr` (v1.3.1;  Wickham et al., 2024), `worrms` (v0.4.3; Chamberlain & Vanhoorne, 2023), `stringr` (v1.5.1; Wickham, 2023), `stringi` (v1.8.7; Gagolewski, 2022), `jsonlite` (v2.0.0; Cooley, 2022), `curl` (v6.2.2; Ooms, 2025), `Biostrings` (v2.76.0; Pagès et al., 2025), `rentrez` (v1.2.3; Winter, 2017), `taxonomizr` (v0.11.1; Sherrill-Mix, 2025), and `phylotools` (v0.2.2;Zhang, 2017) .
The Python script was executed using Python version 3.8.3. The following external libraries were used: `numpy` (v.1.23.4; C. R. Harris et al., 2020) and `Biopython` (v1.76, including `Bio.Seq`, `Bio.SeqIO`, and `Bio.Entrez`; Cock et al., 2009).


The comparison was conducted in R with the following additional R packages: `ggplot2` (v3.5.1; Kassambara, 2023) `tibble` (v3.2.1; Müller & Wickham, 2023), `gapminder` (v1.0.0; Bryan, 2023), `treemapify` (v2.5.6; Wilkins, 2023), `ggpubr` (v0.6.0; Kassambara, 2023), `ghibli` (v0.3.4; Henderson, 2024).

Note: `dplyr`, `tidyr`, `stringr`, `tibble` and `ggplot2` are all part of the R package `tidyverse`

### How to install all packages
```
options(repos = c(CRAN = "https://cloud.r-project.org/"))
# dependencies for pr2database:
 install.packages("devtools")
 install.packages("BiocManager")
 BiocManager::install("Biostrings")
devtools::install_github("pr2database/pr2database")

#the rest
install.packages(c("tidyverse","stringi","ini","worrms","taxonomizr","phylotools","jsonlite","curl","utils","rentrez","gapminder","treemapify","ggpubr", "ghibli"))
```

## Pipeline Scripts

Part of all scripts is that the environment is saved, which allow easy repetition of the scripts. The numbers of the output refer always to the script which created them. Next to the folders in which the sequences will be saved, three folders will be created: 1.) workspace 2.) intermediate results (Script 1-8) 3.) final results (9-10).

E.g. `2.9_F_species_FINAL_withAlgaebase.csv` is an output file of the script `02_Clean_Input.R`.

Pipeline Scripts:
- 00_Function_Library.R
- 01_Clean_PR2.R
- 02_Clean_Input.R
- 03_Download_from_PR2.R
- 04_Identify_missing_species.R
- 05_Check_NCBI.py
- 06_Clean_Sort_NCBI_Downloads.R
- 07_Download_Rest_PR2.R
- 08_Combine_Sort_Files.R
- 09_Taxonomy.R
- 10_Sort_Fasta.R

### 00_Function_Library.R

Function: R script containing all functions and also filter options

### 01_Clean_PR2.R

Function: Filters the PR² reference database to retain only nuclear 18S sequences identified to genus or species level, removing organelle-derived and low-resolution entries.

Input:
- 00_login_data.ini

Output:
- 1.12_F_Cleaned_pr2_database_wAlgbase.tsv
- 1.12_F_Cleaned_pr2_database_wAlgbase

### 02_Clean_Input.R

Function: Prepares and formats the user-provided species list to ensure compatibility with the cleaned PR² reference database.

Input:
- 00_login_data.ini
- Formated_species_list.csv
  
Output:
- 2.4_F_species_FINAL.csv # includes only WoRMS results 
- 2.8_F_Algaebase_specieslist.csv  # includes only Algaebase results 
- 2.9_F_species_FINAL_withAlgaebase.csv

### 03_Download_from_PR2.R

Function: Searches the cleaned PR² database for user-specified species and related taxa, downloads matching sequences, and generates a detailed metadata file including taxonomic and sequence information.

Input:
- 1.12_F_Cleaned_pr2_database_wAlgbase.tsv
- 2.9_F_species_FINAL_withAlgaebase.csv

Output:
- 3.1_F_Overview_Species.cs
- PR2_Sequences/Search/SEQUENCES.fasta


### 04_Identify_missing_species.R

Function: Compares the user's species list with downloaded PR² data to identify missing species and prepare input for NCBI searches.

Input:
- 2.9_F_species_FINAL_withAlgaebase.csv
- 3.1_F_Overview_Species.csv

Output:
- 4.1_Missing_Species.csv
- 4.2_Present_Species.csv


### 05_Check_NCBI.py

Function:  Searches NCBI for nuclear 18S sequences of missing or related species, filters results, and downloads relevant sequences. The created folder and metadata will contain the data when the downloading was started.

Input:
- 4.1_Missing_Species.csv
- 4.2_Present_Species.csv
  
Output:
- 4.2_Present_Species_completed_withNCBI.csv
- 02_NCBI_${DATE}/SEQUENCES.fasta
- 02_NCBI_${DATE}_results/


### 06_Clean_Sort_NCBI_Downloads.R

Function: Cleans and filters downloaded NCBI sequences, removing duplicates and non-nuclear 18S entries, and standardizes metadata.

Input:
- 05_Missing_PR2_Downloaded_species_NCBI_info__${DATE}.tsv
- 1.12_F_Cleaned_pr2_database_wAlgbase
- 2.9_F_species_FINAL_withAlgaebase.csv
- 4.1_Missing_Species.csv

Output:
- 6.3_Species_NCBI.csv

### 07_Download_Rest_PR2.R

Function:  Downloads any remaining PR² sequences not yet retrieved, ensuring database completeness.

Input:
- 1.12_F_Cleaned_pr2_database_wAlgbase.tsv
- 3.1_F_Overview_Species.csv
  
Output:
- 7.1_Overview_PR2_Rest.csv
- PR2_Sequences/Rest/SEQUENCES.fasta

### 08_Combine_Sort_Files.R

Function: Merges cleaned NCBI and PR² sequences, assigns metadata, and organizes files to avoid duplication.

Input:
- 3.1_F_Overview_Species.csv
- 6.3_Species_NCBI.csv
- 7.1_Overview_PR2_Rest.csv
- 
Output:
- PR2_Sequences/NCBI/SEQUENCES.fasta
- 8.1_Overview_Sequences_ALL.csv


### 09_Taxonomy.R

Function: Assigns a standardized nine-level WoRMS taxonomy to all sequences, resolving ambiguities and excluding entries without valid taxonomy.

Input:
- 00_login_data.ini
- 1.7_F_Cleanded_pr2_taxonomy.tsv
- 8.1_Overview_Sequences_ALL.csv

Output:
- PR2_metadata_Algaebase.csv
- 9.5_Taxonomy_FINAL.tax
- 9.5_Taxonomy_FINAL2.tax
- 9.6_Overview_Sequences_FINAL.csv


### 10_Sort_Fasta.R

Function: Produces the final fasta file with only taxonomically validated nuclear 18S reference sequences.

Input:
- 9.5_Taxonomy_FINAL2.tax
- 9.6_Overview_Sequences_FINAL.csv

Output:
- 10.1_Sequences_FINAL.fasta
- 10.2_Taxonomy_FINAL_detail.tax
- 10.2_Taxonomy_FINAL.tax
- 10.3__Overview_species_lacking_taxonomy.csv
- 10.4__Overview_species_lacking_taxonomy_detailed.csv

##  Scripts Associated with Benchmarking

These are located in 13_Benchmark_Script and should be executed from the parent directory of `13_Benchmark_Script`, meaning `PATH/TO/FOLDER`, not from within `PATH/TO/FOLDER/13_Benchmark_Script`.

To execute in the following order:

- 11_Create_Database_Versions.R
- 12_Convert_header_Crabs_part1.R
- 12_Run_CRABS_DB.sh
- 12_Convert_header_Crabs_part2.R
- 13_combine_mothur.sh
- 13_README_mothur.txt
- 13_Assignment_Stats_Crabs.R
- 14_Treemaps_db_results.R
- 15_DB_Stats.R

  Community files were splitted and concatinated after mothur assignment due to RAM issues.


## Metadata

Explanation for the different taxanomy assignment strategies:
| Taxonomy | |
| ------------- | ------------- |
| 1 | WoRMS, AphiaID in PR2 |
| 2 | WoRMS, without AphiaID in PR2, but a sequence of the same genus has|
| 3 | WoRMS, without AphiaID of any sequence of that genus in PR |
| 4 | WoRMS, without AphiaID and non unique genus name|
| 5 | WoRMS, not 1-4|
| 6 | Algaebase |

## References
