
rm(list=ls())
#setwd("/home/jromahn/2024april_refdb_stefanie")
source("00_Function_Library.R") # read file with functions
read_simple_ini("00_login_data.ini")


require(dplyr)
require(stringr)

##################################################################################
## Aim: Identify the missing species and genera in the PR2 database
## Function: Compares the user's species list with downloaded PR² data to identify missing species and prepare input for NCBI searches.
######
# 1.) Missing species & genera identification
# 2.) NCBI query preparation for missing taxa (species & genus)
# 3.) NCBI context list for existing taxa (species & genus)
##################################################################################

#input
output_path <- file.path(path_to_output, "01_intermediate_results" )# results of a script
workspace_path <- file.path(path_to_output, "00_workspace")

#load(file.path(workspace_path, "4.3_Workspace"))

#Import files
species_list <- read.table(file.path(output_path,"2.9_F_species_FINAL_withAlgaebase.csv"), header = TRUE)
spec_overview <- read.table(file.path(output_path,"3.1_F_Overview_Species.csv"), header = TRUE, sep=",")



species_list_akt <- species_list %>%
  filter(status == "aktuell") 

spec_missing <- species_list_akt %>%
  mutate(acc_Name = str_replace_na(acc_Name, "NONE")) %>%
  mutate(species = str_replace_all(species, "sp$", "sp.")) %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  filter(!acc_Name %in% spec_overview$Acc_Name) %>%
  filter(!species %in% spec_overview$species)

#################################################  ────────────────────────────────────────────────────────────────
spec_present <- species_list_akt %>%
  mutate(acc_Name = str_replace_na(acc_Name, "NONE")) %>%
  mutate(species = str_replace_all(species, "sp$", "sp.")) %>%
  mutate(species = str_replace_all(species, " ", "_")) %>%
  filter(acc_Name %in% spec_overview$Acc_Name | species %in% spec_overview$species)

spec_for_NCBI <- spec_missing %>%
  select(genus, species, acc_Name, taxon) %>%
  rowwise() %>%
  mutate(acc_Name = replace(acc_Name, acc_Name == "NONE", species)) %>%
  mutate(acc_Name = str_remove_all(acc_Name, "_sp\\.")) %>%
  mutate(acc_Name = str_replace_all(acc_Name, "_", " "))%>%
  ungroup() %>%
  select(-species) %>%
  dplyr::rename(species_name = acc_Name) %>%
  distinct()

spec_for_NCBI_FINAL <- spec_for_NCBI %>%
  filter(str_detect(species_name, "^\\w+$")) %>%
  mutate(species_name = paste(species_name, "ABCXYZ", sep = " ")) %>%
  bind_rows(spec_for_NCBI %>% filter(str_detect(species_name, "^\\w+$", negate = TRUE))) %>%
  arrange(species_name)

write.table(spec_for_NCBI_FINAL, file.path(output_path, "4.1_Missing_Species.csv"), row.names = FALSE, sep = ",")

#################################################  ────────────────────────────────────────────────────────────────
# create list of species of a missing genus which are already represented so they don't have to be downloaded again
spec_for_NCBI_gen <- spec_present %>%
  select(genus, species, acc_Name) %>%
  rowwise() %>%
  mutate(acc_Name = replace(acc_Name, acc_Name == "NONE", species)) %>%
  mutate(acc_Name = str_remove_all(acc_Name, "_sp\\.")) %>%
  mutate(acc_Name = str_replace_all(acc_Name, "_", " "))%>%
  ungroup() %>%
  select(-species) %>%
  dplyr::rename(species_name = acc_Name) %>%
  distinct()

spec_for_NCBI_gen_FINAL <- spec_for_NCBI_gen %>%
  filter(str_detect(species_name, "^\\w+$")) %>%
  mutate(species_name = paste(species_name, "ABCXYZ", sep = " ")) %>%
  bind_rows(spec_for_NCBI_gen %>% filter(str_detect(species_name, "^\\w+$", negate = TRUE))) %>%
  arrange(species_name) 

write.table(spec_for_NCBI_gen_FINAL, file.path(output_path, "4.2_Present_Species.csv"), row.names = FALSE, sep = ",")
save.image(file.path(workspace_path, "4.3_Workspace.RData"))
