
rm(list=ls())
#setwd("/PATH/TO/11_Benchmark_Scripts") # to path above the Script folder NOT including it
source("00_Function_Library.R") # read file with functions

require(dplyr)
require(tidyr)
require(Biostrings)
require(stringr)
require(rentrez)

#Import files
#Import files
output_path <- "01_intermediate_results" # results of a script
workspace_path <- "00_workspace"

spec_pr2 <- read.table(file.path(output_path, "3.1_F_Overview_Species.csv"), header = TRUE, sep=",") %>% mutate(file = NA)
spec_ncbi <- read.table(file.path(output_path, "6.3_Species_NCBI.csv"), header = TRUE)
rest_pr2 <- read.table(file.path(output_path, "7.1_Overview_PR2_Rest.csv"), header = TRUE) %>% mutate(file = NA)

################################################################################
##output
out_path_seq <- "PR2_Sequences/NCBI"
path_creation(c(out_path_seq, output_path))

#Combine all downloaded sequences
overview_A <- bind_rows(spec_pr2, spec_ncbi, rest_pr2) %>%
  group_by(Clean_Name) %>%
  nest() %>%
  rowwise() %>%
  mutate(num = nrow(data)) %>%
  mutate(no_Aphia = list(subset(data, is.na(data$AphiaID)))) %>%
  mutate(n_NA = nrow(no_Aphia)) 

overview_B <- overview_A %>%
  filter(n_NA != 0) %>%
  filter(n_NA != num) %>%
  unnest(cols = data) %>%
  group_by(Clean_Name) %>%
  mutate(AphiaID = AphiaID[!is.na(AphiaID)][1]) %>%
  mutate(Acc_Name = Acc_Name[!is.na(Acc_Name)][1]) %>%
  bind_rows(overview_A %>% unnest(cols = data) %>% filter(n_NA == 0 | n_NA == num)) %>%
  select(-num, -n_NA, -no_Aphia) %>%
  mutate(dummy = Acc_Name, dummy = ifelse(! is.na(dummy),dummy, gsub( " sp\\.", "", Clean_Name, perl=T)))
  #mutate(dummy = Acc_Name, dummy = replace(dummy, is.na(dummy), gsub( " sp\\.", "", Clean_Name, perl=T))) #str_remove(Clean_Name, " sp.")


#Renumber all sequences and remove duplicate NCBI downloads
overview_cleaned <- overview_B %>%
  filter(Database == "NCBI") %>%
  distinct(genbank_accession, .keep_all = TRUE) %>% 
  bind_rows(overview_B %>% filter(Database == "pr2")) %>%
  group_by(dummy) %>%
  arrange(desc(Database)) %>%
  arrange(dummy) %>%
  mutate(Rep_new = row_number()) %>%
  mutate(Replicate = Rep_new) %>%
  arrange(Replicate) %>%
  arrange(dummy)


#Write NCBI downloads that should be kept into new files
print(paste("Start processing NBCI: ",nrow(overview_cleaned)))

pb <- txtProgressBar(min = 1, max = nrow(overview_cleaned), style = 3)

for (i in c(1:nrow(overview_cleaned))) {
  if (overview_cleaned[i,]$Database == "NCBI") {
    
    group <- str_split_i(overview_cleaned[i,]$file, "_", 3)
    file <- overview_cleaned[i,]$file
    
    seqs <- readDNAStringSet(file)
    seqs_names <- names(seqs)
    names(seqs) <- str_split_i(names(seqs), " ", 1)
    
    seq <- seqs[names(seqs) == overview_cleaned[i,]$genbank_accession]
    num_entry <- which(names(seqs) == names(seq))
    names(seq) <- seqs_names[num_entry]
    
    if(length(seq)>0){ # check if seq is not deleted
      
      overview_cleaned$sequence_length[i] <- width(seq)
      writeXStringSet(seq, sprintf(paste(out_path_seq,"/%s%s", sep=""), str_replace_all(overview_cleaned[i,]$dummy, " ", "_"), paste("_", overview_cleaned$Replicate[i], sep = "")))
    }
  } 
  setTxtProgressBar(pb, i)
};close(pb); rm(pb, seq, seqs, group, file, num_entry, i, seqs_names)


overview_FINAL <- overview_cleaned %>% ungroup %>% select(-file, -dummy)
write.table(overview_FINAL, file.path(output_path, "8.1_Overview_Sequences_ALL.csv"), col.names = TRUE, row.names = FALSE)


save.image(file.path(workspace_path,"8.2_Workspace"))
