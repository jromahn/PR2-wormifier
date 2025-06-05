# PR2-wormifier


The manuscript of this pipeline will be soon submitted to Molecular Ecology Resources.

## Introduction


![Graphical overview of the pipeline](PR2-wormifier.png)
## Input

As Species List the pipeline needs a comma-separated table with the columns: "genus","species","taxon".
The taxon column is important to control if the downloaded and found species really belongs to the same taxonomic group, which is essential to check due to non-unique genus names. The taxonomic group added to the column should fit WoRMS taxonomy. If you prefer to not define the taxonomy, write NA.

As example file check: Formated_species_list.csv 
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


##  Scripts Associated with Benchmarking in folder: 00_Analysis

To execute in the following order:
The scripts are written, that they are executed from the path of their location.

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

## Detailed overview of the pipeline scripts

For a detailed overview of the function of the script is included in the manuscript, here a overview of the general fnction, the input and the output is provided.

### 00_Function_Library.R

Function: 
Input:
Output:

### 01_Clean_PR2.R

Function: Filters the PR² reference database to retain only nuclear 18S sequences identified to genus or species level, removing organelle-derived and low-resolution entries.
Input:
Output:

### 02_Clean_Input.R

Function: Prepares and formats the user-provided species list to ensure compatibility with the cleaned PR² reference database.
Input:
Output:


### 03_Download_from_PR2.R

Function: Searches the cleaned PR² database for user-specified species and related taxa, downloads matching sequences, and generates a detailed metadata file including taxonomic and sequence information.
Input:
Output:


### 04_Identify_missing_species.R

Function: Compares the user's species list with downloaded PR² data to identify missing species and prepare input for NCBI searches.
Input:
Output:


### 05_Check_NCBI.py

Function:  Searches NCBI for nuclear 18S sequences of missing or related species, filters results, and downloads relevant sequences.
Input:
Output:


### 06_Clean_Sort_NCBI_Downloads.R

Function: Cleans and filters downloaded NCBI sequences, removing duplicates and non-nuclear 18S entries, and standardizes metadata.
Input:
Output:



### 07_Download_Rest_PR2.R

Function:  Downloads any remaining PR² sequences not yet retrieved, ensuring database completeness.
Input:
Output:


### 08_Combine_Sort_Files.R

Function: Merges cleaned NCBI and PR² sequences, assigns metadata, and organizes files to avoid duplication.
Input:
Output:


### 09_Taxonomy.R

Function: Assigns a standardized nine-level WoRMS taxonomy to all sequences, resolving ambiguities and excluding entries without valid taxonomy.
Input:
Output:


### 10_Sort_Fasta.R

Function: Produces the final fasta file with only taxonomically validated nuclear 18S reference sequences.
Input:
Output:


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
