##########################################################################################
#This script is supposed to take the pr2 version of the database and the final database
#version of the pipeline and create the "in-between" versions (both .fasta and .tax file)
##########################################################################################

#setwd("/PATH/TO/02_Scripts_folder") # to path above the Script folder
source("00_Function_Library.R") # read file with functions
read_simple_ini("00_login_data.ini")

require(phylotools)
require(dplyr)
require(stringr)


#define paths
path_results <- file.path(path_to_output,"10_FINAL_results")
path_pr2 <- "/PATH/TO/PR2_DATABASE"

output_path <- "11_Benchmark_Database_Versions_Euka02"
primer="Euka02"


pr2_fasta<-"pr2_version_5.1.1_SSU_mothur.fasta"
pr2_tax <- "pr2_version_5.1.1_SSU_mothur.tax"

########################################################################################
path_creation( output_path)


print("Read pr2 and new database files")
#Read in files (taxonomy and sequences of both the first and last version of the database + overview file)
#pr2
taxonomy_1 <- read.table(file.path(path_pr2,pr2_tax))
sequences_1 <- read.fasta(file.path(path_pr2,pr2_fasta))
#own
taxonomy_4 <- read.table(file.path(path_results, "10.2_Taxonomy_FINAL.tax"))
sequences_4 <- read.fasta(file.path(path_results,"10.1_Sequences_FINAL.fasta")) %>% mutate(accession = str_split_i(seq.name, " ", 1))
overview <- read.table(file.path(path_results,"9.6_Overview_Sequences_FINAL.csv"), header =  TRUE, sep=",") %>% select( -taxon, -Rep_new)


#colnames(overview)
#quit()

print("Create intermediate databases")
#Get information about what sequences are in the versions 2 and 3 of the database
overview_4 <- overview %>% filter(Taxonomy != 5) 
overview_3_2 <- overview_4 %>% filter(Database == "pr2")

#Get all the sequences that are in version 2 and 3
sequences_3_2 <- sequences_4 %>% filter(accession %in% overview_3_2$pr2_accession)

#Write the reduced taxonomy for versions 2 and 3
taxonomy_3 <- taxonomy_4 %>% filter(V1 %in% sequences_3_2$accession)
taxonomy_2 <- taxonomy_1 %>% filter(V1 %in% sequences_3_2$accession)

#remove dummy column
sequences_3_2 <- sequences_3_2 %>% select(-accession)
sequences_4 <- sequences_4 %>% select(-accession)

#Qrite all database versions into seperate file
print("Write taxonomy files")
write.table(taxonomy_1 %>% arrange(V1), file.path(output_path,paste("DB_V1",primer, "Taxonomy.tax", sep="_")), quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(taxonomy_2 %>% arrange(V1), file.path(output_path,paste("DB_V2",primer, "Taxonomy.tax", sep="_")), quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(taxonomy_3 %>% arrange(V1), file.path(output_path,paste("DB_V3",primer, "Taxonomy.tax", sep="_")), quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(taxonomy_4 %>% arrange(V1), file.path(output_path,paste("DB_V4",primer, "Taxonomy.tax", sep="_")), quote = FALSE, col.names = FALSE, row.names = FALSE)

print("Write fasta files")
dat2fasta(sequences_1 %>% arrange(seq.name), file.path(output_path,paste("DB_V1",primer, "Sequences.fasta", sep="_")))
dat2fasta(sequences_3_2 %>% arrange(seq.name), file.path(output_path,paste("DB_V2",primer, "Sequences.fasta", sep="_")))
dat2fasta(sequences_3_2 %>% arrange(seq.name), file.path(output_path,paste("DB_V3",primer, "Sequences.fasta", sep="_")))
dat2fasta(sequences_4 %>% arrange(seq.name), file.path(output_path,paste("DB_V4",primer, "Sequences.fasta", sep="_")))
