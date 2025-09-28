
rm(list=ls())
#setwd("/PATH/TO/11_Benchmark_Scripts") # to path above the Script folder NOT including it
source("00_Function_Library.R") # read file with functions



require(dplyr)
require(Biostrings)
require(stringr)

#input
output_path <- "01_intermediate_results" # results of a script
workspace_path <- "00_workspace"
######

##output
output_seq <- "PR2_Sequences/Search"


path_creation(c(output_path, "PR2_Sequences", output_seq))


pr2 <- read.table(file.path(output_path,"1.12_F_Cleaned_pr2_database_wAlgbase.tsv"), header = TRUE, sep="\t",)
species_list <- read.table(file.path(output_path,"2.9_F_species_FINAL_withAlgaebase.csv"), header = TRUE)


# filter for accepted name 
spec <- species_list %>%
  filter(status == "aktuell") %>%
  mutate(acc_Name = replace(acc_Name, is.na(acc_Name), "NONE"))

# filter for overlapping genera from pr2 sequences
spec_relevant_sequences <- pr2 %>%
  rowwise() %>%
  mutate(dummy = Acc_Name, dummy = replace(dummy, is.na(dummy), str_remove(Clean_Name, " sp."))) %>%
  mutate(Genus = str_split_i(dummy, " ", 1)) %>%
  filter(Genus %in% spec$genus)


spec_overview <- c()

for (i in unique(spec_relevant_sequences$dummy)) {
  sequences <- spec_relevant_sequences %>% filter(dummy == i)
  #print(unique(sequences$Acc_Name))
  for(i in c(1:nrow(sequences))) {
    
    seq <- DNAStringSet(sequences[i,]$sequence)
    names(seq) <- paste(sequences[i,]$pr2_accession, sequences[i,]$dummy, sep=" ")
    writeXStringSet(seq, sprintf(paste(output_seq,"/%s%s",sep=""), str_replace_all(sequences[i,]$dummy, " ", "_"), paste("_", i, sep = "")))
    
    
    info_1 <- sequences[i,] %>% select(pr2_accession, genbank_accession, AphiaID, Acc_Name, Clean_Name, species, sequence_length)
    info_2 <- tibble(i, paste("pr2")); colnames(info_2) <- c("Replicate", "Database")
    info <- bind_cols(info_1, info_2)
    
    spec_overview <- bind_rows(spec_overview, info)
  }
}; rm(sequences, seq, info, info_1, info_2)



write.table(spec_overview, file.path(output_path,"3.1_F_Overview_Species.csv"), row.names = FALSE, sep=",")
save.image(file.path(workspace_path, "3.2_Workspace"))
load(file.path(workspace_path, "3.2_Workspace"))
