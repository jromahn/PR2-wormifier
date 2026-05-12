
rm(list=ls())
#setwd("/PATH/TO/11_Benchmark_Scripts") # to path above the Script folder NOT including it
source("00_Function_Library.R") # read file with functions
read_simple_ini("00_login_data.ini")


require(dplyr)
require(Biostrings)
require(stringr)

##################################################################################
## Aim: Downloading the species list entries from PR2
## Function: Searches the cleaned PR² database for user-specified species and related taxa, 
######  downloads matching sequences, and generates a detailed metadata file including taxonomic and sequence information.
######
# 1.) Filtering
# 2.) Sequence export
# 3.) Metadata collection
##################################################################################

#input ────────────────────────────────────────────────────────────────
output_path <- file.path(path_to_output,"01_intermediate_results" ) # results of a script
workspace_path <- file.path(path_to_output, "00_workspace")
######

##output ────────────────────────────────────────────────────────────────
output_seq <- file.path(path_to_output,"PR2_Sequences/Search" )
path_creation(c(output_path, file.path(path_to_output,"PR2_Sequences"), output_seq))

pr2 <-readRDS(file.path(output_path,"1.12_F_Cleaned_pr2_database_wAlgbase.RDS"))
species_list <- readRDS(file.path(output_path,"2.9_F_species_FINAL_withAlgaebase.RDS"))


# filter for accepted name  ────────────────────────────────────────────────────────────────
spec <- species_list %>%
  filter(status == "aktuell") %>%
  mutate(acc_Name = replace(acc_Name, is.na(acc_Name), "NONE"))

# filter for overlapping genera from pr2 sequences ────────────────────────────────────────────────────────────────
spec_relevant_sequences <- pr2 %>%
  rowwise() %>%
  mutate(dummy = Acc_Name, dummy = replace(dummy, is.na(dummy), str_remove(Clean_Name, " sp."))) %>%
  mutate(Genus = str_split_i(dummy, " ", 1))%>%
  ungroup() %>%
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
save.image(file.path(workspace_path, "3.2_Workspace.RData"))
load(file.path(workspace_path, "3.2_Workspace.RData"))
