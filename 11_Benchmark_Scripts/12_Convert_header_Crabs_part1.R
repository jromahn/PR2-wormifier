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


########################################################################################
### Part 1: Create "fake" header

## Read in sequences and replace header with "fake" accession numbers to be compliant with crabs
#Version 1
sequences_1 <- read.fasta(file.path(path,"DB_V1_Sequences.fasta")) %>%
  mutate(fake_header = paste("AB", row_number(), ".1", sep = "")) %>%
  mutate(dummy = 10 - nchar(fake_header)) %>%
  mutate(row_num = row_number()) %>%
  rowwise() %>%
  mutate(fake_header = paste("AB", paste(rep("0", dummy), collapse = ""), row_num, ".1", sep = "")) %>%
  select(-c(dummy, row_num))

new_sequences_1 <- sequences_1 %>%
  select(-seq.name) %>%
  relocate(fake_header, .before = "seq.text") %>%
  rename(seq.name = fake_header)


#Version 2
sequences_2 <- read.fasta(file.path(path,"DB_V2_Sequences.fasta")) %>%
  mutate(fake_header = paste("AB", row_number(), ".1", sep = "")) %>%
  mutate(dummy = 10 - nchar(fake_header)) %>%
  mutate(row_num = row_number()) %>%
  rowwise() %>%
  mutate(fake_header = paste("AB", paste(rep("0", dummy), collapse = ""), row_num, ".1", sep = "")) %>%
  select(-c(dummy, row_num))

new_sequences_2 <- sequences_2 %>%
  select(-seq.name) %>%
  relocate(fake_header, .before = "seq.text") %>%
  rename(seq.name = fake_header)


#Version 3
sequences_3 <- read.fasta(file.path(path,"DB_V3_Sequences.fasta")) %>%
  mutate(fake_header = paste("AB", row_number(), ".1", sep = "")) %>%
  mutate(dummy = 10 - nchar(fake_header)) %>%
  mutate(row_num = row_number()) %>%
  rowwise() %>%
  mutate(fake_header = paste("AB", paste(rep("0", dummy), collapse = ""), row_num, ".1", sep = "")) %>%
  select(-c(dummy, row_num))

new_sequences_3 <- sequences_3 %>%
  select(-seq.name) %>%
  relocate(fake_header, .before = "seq.text") %>%
  rename(seq.name = fake_header)


#Version 4
sequences_4 <- read.fasta(file.path(path,"DB_V4_Sequences.fasta")) %>%
  mutate(fake_header = paste("AB", row_number(), ".1", sep = "")) %>%
  mutate(dummy = 10 - nchar(fake_header)) %>%
  mutate(row_num = row_number()) %>%
  rowwise() %>%
  mutate(fake_header = paste("AB", paste(rep("0", dummy), collapse = ""), row_num, ".1", sep = "")) %>%
  select(-c(dummy, row_num))

new_sequences_4 <- sequences_4 %>%
  select(-seq.name) %>%
  relocate(fake_header, .before = "seq.text") %>%
  rename(seq.name = fake_header)

## Write files with sequences with "fake" header
dat2fasta(new_sequences_1, file.path(path,"DB_V1_Sequences_fake_header.fasta"))
dat2fasta(new_sequences_2, file.path(path,"DB_V2_Sequences_fake_header.fasta"))
dat2fasta(new_sequences_3, file.path(path,"DB_V3_Sequences_fake_header.fasta"))
dat2fasta(new_sequences_4, file.path(path,"DB_V4_Sequences_fake_header.fasta"))
