
rm(list=ls())
#setwd("/PATH/TO/11_Benchmark_Scripts") # to path above the Script folder NOT including it
source("00_Function_Library.R") # read file with functions

require(dplyr)
require(stringr)
require(tidyr)
require(Biostrings)

#Import files
output_path <- "01_intermediate_results" # results of a script
workspace_path <- "00_workspace"
pr2 <- read.table(file.path(output_path,"1.12_F_Cleaned_pr2_database_wAlgbase.tsv"), header = TRUE, sep="\t")
spec_pr2 <- read.table(file.path( output_path,"3.1_F_Overview_Species.csv"), header = TRUE, sep=",")

output_path_seq <- "PR2_Sequences/Rest"

path_creation(c(output_path_seq, output_path))

################################################################################
#Sequences not yet downloaded
print("Create list for downloading")


pr2_Rest <- pr2 %>%
  filter(!pr2_accession %in% spec_pr2$pr2_accession) %>%
  rowwise() %>%
  mutate(dummy = Acc_Name, dummy = replace(dummy, is.na(dummy), str_remove(Clean_Name, " sp.")))

pr2_rest_overview <- c()

#Download sequences with downloading bar
print(paste("Start downloading: ",length(unique(pr2_Rest$dummy))), "species")

pb <- txtProgressBar(min = 1, max = length(unique(pr2_Rest$dummy)), style = 3)
counter=0
for (i in unique(pr2_Rest$dummy)) {
  sequences <- pr2_Rest %>% filter(dummy == i)
  
  for(i in c(1:nrow(sequences))) {
    
    seq <- DNAStringSet(sequences[i,]$sequence)
    names(seq) <- paste(sequences[i,]$pr2_accession, sequences[i,]$dummy, sep = " ")
    writeXStringSet(seq, sprintf(paste(output_path_seq,"/%s%s", sep=""), str_replace_all(sequences[i,]$dummy, " ", "_"), paste("_", i, sep = "")))
    
    
    info_1 <- sequences[i,] %>% select(pr2_accession, genbank_accession, AphiaID, Acc_Name, Clean_Name, species, sequence_length)
    info_2 <- tibble(i, paste("pr2")); colnames(info_2) <- c("Replicate", "Database")
    info <- bind_cols(info_1, info_2)
    
    pr2_rest_overview <- bind_rows(pr2_rest_overview, info)
  }
  counter = counter + 1
  setTxtProgressBar(pb, counter)
};close(pb); rm(pb, counter, sequences, seq, info, info_1, info_2)



#Write overview file
write.table(pr2_rest_overview, file.path(output_path, "7.1_Overview_PR2_Rest.csv"))

save.image(file.path(workspace_path, "7.2_Workspace"))
