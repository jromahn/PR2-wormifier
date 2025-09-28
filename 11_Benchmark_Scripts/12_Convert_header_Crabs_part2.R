################################################################################
#This script is supposed to convert the sequence headers of the reference 
#databases with "fake" one to be compatible with CRABS. It the can be used to 
#convert back the fake header of the amplicons to their original ones
################################################################################


#setwd("/PATH/TO/02_Scripts_folder") # to path above the Script folder

require(dplyr)
require(stringr)
require(phylotools)

path <- "11_Benchmark_Database_Versions"
primer <- "Euka02"

CRAB_file_db1 <- "DB_V1_Amplicons_fake_header_L500_P0.85.fasta"
CRAB_file_db2 <- "DB_V2_Amplicons_fake_header_L500_P0.85.fasta"
CRAB_file_db3 <- "DB_V3_Amplicons_fake_header_L500_P0.85.fasta"
CRAB_file_db4 <- "DB_V4_Amplicons_fake_header_L500_P0.85.fasta"


########################################################################################

### read in original sequences
sequences_1 <- read.fasta(file.path(path,"DB_V1_Sequences.fasta")) %>%
  mutate(fake_header = paste("AB", row_number(), ".1", sep = "")) %>%
  mutate(dummy = 10 - nchar(fake_header)) %>%
  mutate(row_num = row_number()) %>%
  rowwise() %>%
  mutate(fake_header = paste("AB", paste(rep("0", dummy), collapse = ""), row_num, ".1", sep = "")) %>%
  select(-c(dummy, row_num))

sequences_2 <- read.fasta(file.path(path,"DB_V2_Sequences.fasta")) %>%
  mutate(fake_header = paste("AB", row_number(), ".1", sep = "")) %>%
  mutate(dummy = 10 - nchar(fake_header)) %>%
  mutate(row_num = row_number()) %>%
  rowwise() %>%
  mutate(fake_header = paste("AB", paste(rep("0", dummy), collapse = ""), row_num, ".1", sep = "")) %>%
  select(-c(dummy, row_num))

sequences_3 <- read.fasta(file.path(path,"DB_V3_Sequences.fasta")) %>%
  mutate(fake_header = paste("AB", row_number(), ".1", sep = "")) %>%
  mutate(dummy = 10 - nchar(fake_header)) %>%
  mutate(row_num = row_number()) %>%
  rowwise() %>%
  mutate(fake_header = paste("AB", paste(rep("0", dummy), collapse = ""), row_num, ".1", sep = "")) %>%
  select(-c(dummy, row_num))

sequences_4 <- read.fasta(file.path(path,"DB_V4_Sequences.fasta")) %>%
  mutate(fake_header = paste("AB", row_number(), ".1", sep = "")) %>%
  mutate(dummy = 10 - nchar(fake_header)) %>%
  mutate(row_num = row_number()) %>%
  rowwise() %>%
  mutate(fake_header = paste("AB", paste(rep("0", dummy), collapse = ""), row_num, ".1", sep = "")) %>%
  select(-c(dummy, row_num))

#after run in-silico PCR with CRABS
################################################################################

### Part 2: Convert header back
## Read header and convert back

amplicons_1 <- read.fasta(file.path(path,CRAB_file_db1)) %>%
  rename(fake_header = "seq.name", amplicon = "seq.text") %>%
  left_join(sequences_1, by = "fake_header") %>%
  select(-fake_header, -seq.text) %>%
  relocate(seq.name, .before = amplicon) %>%
  arrange(seq.name) %>%
  rename(seq.text = amplicon)

amplicons_2 <- read.fasta(file.path(path,CRAB_file_db2)) %>%
  rename(fake_header = "seq.name", amplicon = "seq.text") %>%
  left_join(sequences_2, by = "fake_header") %>%
  select(-fake_header, -seq.text) %>%
  relocate(seq.name, .before = amplicon) %>%
  arrange(seq.name) %>%
  rename(seq.text = amplicon)

amplicons_3 <- read.fasta(file.path(path,CRAB_file_db3)) %>%
  rename(fake_header = "seq.name", amplicon = "seq.text") %>%
  left_join(sequences_3, by = "fake_header") %>%
  select(-fake_header, -seq.text) %>%
  relocate(seq.name, .before = amplicon) %>%
  arrange(seq.name) %>%
  rename(seq.text = amplicon)

amplicons_4 <- read.fasta(file.path(path,CRAB_file_db4)) %>%
  rename(fake_header = "seq.name", amplicon = "seq.text") %>%
  left_join(sequences_4, by = "fake_header") %>%
  select(-fake_header, -seq.text) %>%
  relocate(seq.name, .before = amplicon) %>%
  arrange(seq.name) %>%
  rename(seq.text = amplicon)

## Write tables with amplified sequences and original header

dat2fasta(amplicons_1, file.path(path,paste("11",primer, "database_V1.fasta", sep="_")))
dat2fasta(amplicons_2, file.path(path,paste("11",primer, "database_V2.fasta", sep="_")))
dat2fasta(amplicons_3, file.path(path,paste("11",primer, "database_V3.fasta", sep="_")))
dat2fasta(amplicons_4, file.path(path,paste("11",primer, "database_V4.fasta", sep="_")))



## Filter taxonomy to keep only sequences that have been amplified by in-silico PCR
#Version 1
amplicons_1_accessions <- amplicons_1 %>%
  mutate(accession = str_split_i(seq.name, " ", 1)) %>%
  select(-seq.name, -seq.text)

taxonomy_1 <- read.table(file.path(path,"DB_V1_Taxonomy.tax")) %>%
  filter(V1 %in% amplicons_1_accessions$accession) %>%
  arrange(V1)

#Version 2
amplicons_2_accessions <- amplicons_2 %>%
  mutate(accession = str_split_i(seq.name, " ", 1)) %>%
  select(-seq.name, -seq.text)

taxonomy_2 <- read.table(file.path(path,"DB_V2_Taxonomy.tax")) %>%
  filter(V1 %in% amplicons_2_accessions$accession) %>%
  arrange(V1)

#Version 3
amplicons_3_accessions <- amplicons_3 %>%
  mutate(accession = str_split_i(seq.name, " ", 1)) %>%
  select(-seq.name, -seq.text)

taxonomy_3 <- read.table(file.path(path,"DB_V3_Taxonomy.tax")) %>%
  filter(V1 %in% amplicons_3_accessions$accession) %>%
  arrange(V1)

#Version 4
amplicons_4_accessions <- amplicons_4 %>%
  mutate(accession = str_split_i(seq.name, " ", 1)) %>%
  select(-seq.name, -seq.text)

taxonomy_4 <- read.table(file.path(path,"DB_V4_Taxonomy.tax")) %>%
  filter(V1 %in% amplicons_4_accessions$accession) %>%
  arrange(V1)


## Write reduced taxonomy
write.table(taxonomy_1, file.path(path,paste("11", primer, "Taxonomy_DB1.tax",sep="_")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(taxonomy_2, file.path(path,paste("11", primer, "Taxonomy_DB2.tax",sep="_")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(taxonomy_3, file.path(path,paste("11", primer, "Taxonomy_DB3.tax",sep="_")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(taxonomy_4, file.path(path,paste("11", primer, "Taxonomy_DB4.tax",sep="_")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


