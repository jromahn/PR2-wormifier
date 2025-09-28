
rm(list=ls())
#setwd("/PATH/TO/11_Benchmark_Scripts") # to path above the Script folder NOT including it
source("00_Function_Library.R") # read file with functions
read_simple_ini("00_login_data.ini")
####################
## Aim: Pre-cleaning of pr2 database
##
######################

require(utils)
require(pr2database)
require(dplyr)
require(tidyr)
require(worrms)
require(stringr)# try next time without
require(stringi) # try next time without
library(jsonlite)
library(curl)
#require(purrr)

####################################
# data storage
output_path <- "01_intermediate_results" # results of a script
workspace_path <- "00_workspace"


path_creation(c(workspace_path, output_path))
####################################


#Import Pr2-database and Pr2-taxonomy
pr2 <- pr2_database()
pr2_tax <- pr2_taxonomy()


######
sink("1.9_F_Overview_and_statistics.txt")
print(paste("Downloaded pr2_database & taxonomy:", format(Sys.time(), "%a %b %d %X %Y")))
sink()
unlink("1.9_F_Overview_and_statistics.txt")
######

## Taxonomy ##

## update 13 May 2024 : update of pr2, no longer worms_id is missing using pr2_taxonomy() adding via the pr2 dataframe
if (! "worms_id" %in% names(pr2_tax)) {
  pr2_tax <- pr2_tax %>% left_join(pr2 %>% select(species_url,worms_id )%>% unique(), by = "species_url")
}


#Clean taxonomy
pr2_tax_A <- clean_pr2_taxonomy(pr2_tax) 
#colnames(pr2_tax)
#colnames(pr2)
#quit()
#############################################################################

#Create columns for accepted names and corresponding AphiaIDs
pr2_tax_B <- pr2_tax_A %>%
  mutate(Acc_Name = NA) %>%
  mutate(Acc_Aphia = NA)

print("Start pr2_tax_C")
#Get accepted name and AphiaID where AphiaID is already known
pr2_tax_C <- pr2_tax_B %>%
  filter(!is.na(worms_id)) %>%
  rowwise() %>%
  mutate(Acc_Name = wm_record(worms_id)$valid_name ) %>%
  mutate(Acc_Aphia = wm_record(worms_id)$valid_AphiaID) %>%
  rbind(pr2_tax_B %>% filter(is.na(worms_id)))

write.table(pr2_tax_C, file.path(workspace_path,"1.1_pr2_tax_C"), row.names = FALSE)
pr2_tax_C <- read.table(file.path(workspace_path,"1.1_pr2_tax_C"), header = TRUE)



#Get Worms records for entries where AphiaID is not already known
print("Start pr2_tax_D")
pr2_tax_D <- pr2_tax_C %>%
  mutate(Clean_Name = str_replace_all(Clean_Name, "_", " ")) %>%
  rowwise() %>%
  mutate(Worms_record = list(tibble()))


pb <- txtProgressBar(min = 1, max = nrow(pr2_tax_D), style = 3)

############## needs long
for (i in c(1:nrow(pr2_tax_D))) {
  if (is.na(pr2_tax_D$Acc_Name[i])) {
    pr2_tax_D$Worms_record[i] <- try_worms(pr2_tax_D$Clean_Name[i], marine_only = FALSE)
  }
  if (i %% 100 == 0) {
    saveRDS(pr2_tax_D, file.path(workspace_path, "1.2_pr2_tax_D"))
    write.table(i, file.path(workspace_path, "1.3_last_save"), row.names = FALSE, col.names = FALSE)
  }
  setTxtProgressBar(pb, i)
}

close(pb); rm(pb)


saveRDS(pr2_tax_D, file.path(output_path,"1.2_pr2_tax_D"))
pr2_tax_D <-readRDS(file.path(output_path,"1.2_pr2_tax_D"))


#Rerun last step to make sure no entries are missed due to connection problems
print("Start pr2_tax_E")
pr2_tax_E <- pr2_tax_D

#pr2_tax_E <- readRDS(file.path(workspace_path,"1.4_pr2_tax_E"))

pb <- txtProgressBar(min = 1, max = nrow(pr2_tax_E), style = 3)

for (i in c(1:nrow(pr2_tax_E))) {
  if (is.na(pr2_tax_E$Acc_Name[i]) && sum(!is.na(pr2_tax_E$Worms_record[[i]])) == 0) {
    pr2_tax_E$Worms_record[i] <- try_worms(pr2_tax_E$Clean_Name[i], marine_only = FALSE)
  }
  if (i %% 100 == 0) {
    saveRDS(pr2_tax_E, file.path(workspace_path,"1.4_pr2_tax_E"))
    write.table(i, file.path(workspace_path,"1.5_last_save"), row.names = FALSE, col.names = FALSE)
  }
  setTxtProgressBar(pb, i)
}

close(pb); rm(pb)


saveRDS(pr2_tax_E, file.path(workspace_path,"1.4_pr2_tax_E"))
pr2_tax_E <- readRDS(file.path(workspace_path,"1.4_pr2_tax_E"))

#Split up results: entries that already have AphiaIDs (F_1), new entries with AphiaID (F_2), rest (F_3)
print("Start pr2_tax_F")
pr2_tax_F_1 <- pr2_tax_E %>%
  filter(!is.na(Acc_Name)) %>%
  select(-Worms_record)

pr2_tax_F_2 <- pr2_tax_E %>%
  filter(is.na(Acc_Name)) %>%
  filter(sum(!is.na(Worms_record[[1]])) > 0) %>%
  mutate(n_records = nrow(Worms_record)) %>%
  unnest(cols = Worms_record, names_sep = "_") %>%
  #remove records that are not exact matches for entries on species level
  filter(!(Worms_record_rank == "Genus" & Worms_record_match_type != "exact")) 


pr2_tax_F_2.1 <- filter_for_right_taxa(pr2_tax_F_2 )%>%  group_by(across(domain:Acc_Aphia)) %>% # mod by JR

  nest() %>%
  rowwise() %>%
  
  #choose newest AphiaID for duplicate Worms entries
  mutate(data = list(subset(data, data$Worms_record_AphiaID == max(data$Worms_record_AphiaID)))) %>%
  unnest(cols = data) %>%
  
  mutate(Acc_Name = Worms_record_valid_name) %>%
  mutate(Acc_Aphia = Worms_record_valid_AphiaID) %>%
  select(-starts_with("Worms_record_"), -n_records)

## filtered out information
pr2_tax_F_2.2 <- pr2_tax_F_2 %>% filter(!species %in%pr2_tax_F_2.1$species ) %>%
  relocate(Worms_record_kingdom, .after=species) %>%
  relocate(Worms_record_phylum, .after=Worms_record_kingdom) %>%
  relocate(Worms_record_order, .after=Worms_record_phylum)%>%
  relocate(Worms_record_class, .after=Worms_record_order)
  
  
pr2_tax_F_3 <- pr2_tax_E %>%
  filter(!species %in% pr2_tax_F_1$species) %>%
  filter(!species %in% pr2_tax_F_2.1$species) %>%
  select(-Worms_record)

#Combine all entries
pr2_tax_FINAL <- bind_rows(pr2_tax_F_1, pr2_tax_F_2.1, pr2_tax_F_3)


## Sequence database ##

#Remove sequences without taxonomy
print("Start pr2_A")
pr2_A <- pr2 %>%
  filter(species %in% pr2_tax_FINAL$species)

#Add information: Cleaned names, AphiaID, accepted name
pr2_B <- pr2_A %>%
  dplyr::rename(spec = species) %>%
  rowwise() %>%
  mutate(Acc = list(subset(pr2_tax_FINAL, pr2_tax_FINAL$species == spec))) %>%
  relocate(Acc, .before = pr2_accession) %>%
  unnest(cols = Acc, names_sep = "_") %>%
  dplyr::rename(species = spec) %>%
  dplyr::rename(Clean_Name = Acc_Clean_Name) %>%
  dplyr::rename(Name_accepted = Acc_Acc_Name) %>%
  dplyr::rename(AphiaID = Acc_Acc_Aphia) %>%
  select(-starts_with("Acc_"))

saveRDS(pr2_B, file.path(workspace_path,"1.6_pr2_B"))
pr2_B <- readRDS(file.path(workspace_path,"1.6_pr2_B"))

pr2_FINAL <- pr2_B %>%
  dplyr::rename(Acc_Name = Name_accepted)

write.table(pr2_tax_FINAL, file=file.path(output_path,"1.7_F_Cleanded_pr2_taxonomy"), row.names = FALSE)
write.table(pr2_FINAL, file=file.path( output_path,"1.8_F_Cleanded_pr2_database"), row.names = FALSE)


## Overview ##

perc_taxa <- round((nrow(pr2_tax)-nrow(pr2_tax_FINAL))/nrow(pr2_tax) * 100, 2)
perc_seq <- round((nrow(pr2)-nrow(pr2_FINAL))/nrow(pr2) * 100, 2)
num_clean <- nrow(pr2_tax_A %>% filter(species != Clean_Name))

aphia_start <- nrow(pr2_tax %>% filter(!is.na(worms_id))) 
aphia_end <- nrow(pr2_tax_FINAL %>% filter(!is.na(Acc_Aphia)))

print(getwd())
sink("1.9_F_Overview_and_statistics.txt")

print(paste("Finalised:", format(Sys.time(), "%a %b %d %X %Y")))
print(paste("Before cleaning: \t", nrow(pr2), "sequences of", nrow(pr2_tax), "different taxa \n"))
print(paste("After cleaning: \t", nrow(pr2_FINAL), "sequences of", nrow(pr2_tax_FINAL), "different taxa \n"))
print(paste("Removed: \t\t ", perc_taxa, "% of taxa and ", perc_seq, "% of sequences\n", sep = ""))
print("\n")
print(paste("Names cleaned for", num_clean, "species \n")) 
print(paste("AphiaID added for", aphia_end - aphia_start, "species \n"))
print(paste("AphiaID available for", aphia_end, "species \n"))

## sequences of. cleaning steps




## remove plastid etc
pr2_taxt_non18S <- pr2_tax %>%
  filter(str_detect(species, ":mito")  | str_detect(species, ":nucl") |
           str_detect(species, ":plas")| str_detect(species, ":chro")|
           str_detect(species, ":apic")) 
nrow(pr2_taxt_non18S) # records
nrow(pr2 %>% filter(species %in% pr2_taxt_non18S$species )) # sequences

## bad taxonomy
pr2_badTax <- pr2_tax %>% filter( ! species %in% pr2_taxt_non18S$species & ! species %in% pr2_tax_A$species )
nrow(pr2_badTax) # records
nrow(pr2 %>% filter(species %in% pr2_badTax$species )) # sequences

#terminates the connection 
sink()
unlink("1.9_F_Overview_and_statistics.txt")
save.image(file.path(workspace_path, "1.10_Workspace"))
load(file.path(workspace_path, "1.10_Workspace"))


### Check phytoplankton  with algaebaese
print("Start algaebaese")

pr2_phyto <- filter_pr2_for_algae(pr2_FINAL)
pr2_non_phyto <- pr2_FINAL %>% filter(!pr2_accession %in%pr2_phyto$pr2_accession )


######### query
pr2_phyto <-pr2_phyto %>%mutate(Full_search= ifelse(is.na(Acc_Name), Clean_Name,Acc_Name )) %>%
  mutate(Alg.genus=gsub("^([\\w\\-]+) .*", "\\1", Full_search, perl = T), 
         Alg.epithel= gsub("^[\\w\\-]+ ?([\\w\\-\\.]*).*", "\\1", Full_search, perl = T))



search <- unique(pr2_phyto[,c("Alg.genus", "Alg.epithel", "class", "Full_search")]) %>% filter(Alg.epithel !="")
algaebase_species <- check_species_against_Algaebase(search, key, password)
save(algaebase_species, file= file.path(workspace_path, "PR2__algaebase_info.RData"))


load( file= file.path(workspace_path, "PR2__algaebase_info.RData"))

algaebase.species.df <- algaebase_species[[1]]
algaebase.species.df <- unique(algaebase.species.df)
non.algaebase.df <- algaebase_species[[2]]
rm(non.algaebase.df)

## add orginal search
search <- search%>% mutate(request = paste(Alg.genus, Alg.epithel) )%>%
                    mutate( Alg.genus= NULL, Alg.epithel = NULL)%>%
                    dplyr::rename("taxon"=class)

algaebase.species.df <- merge( algaebase.species.df, search, by= c("request", "taxon"))
rm(search)

algaebase.species.df_cleaned <- algaebase_species_clean_up(algaebase.species.df)
write.table(algaebase.species.df_cleaned, file.path( output_path,"1.11_F_Algaebase_pr2_database.csv"), row.names = FALSE ,sep=",")



algaebase.species.df_short <- algaebase.species.df_cleaned[,c("Algaebase_search", "acc_Name", "Algaebase_entry", "taxon")] %>% 
  unique()%>%
  dplyr::rename("class" = taxon)

#combien with non phytoplankton data
colnames(pr2_phyto)
pr2_phyto_F <- pr2_phyto %>%
    dplyr::rename("Algaebase_search"=Full_search)

pr2_phyto_F <-algaebase.species.df_short %>% 
              right_join( pr2_phyto_F,  by = c("Algaebase_search", "class")) %>% 
              mutate(Acc_Name=ifelse(is.na(acc_Name),Acc_Name,  acc_Name))%>% 
              mutate(  acc_Name= NULL)



pr2_FINAL2 <- plyr::rbind.fill(pr2_non_phyto,pr2_phyto_F)
write.table(pr2_FINAL2, file.path( output_path,"1.12_F_Cleaned_pr2_database_wAlgbase.tsv"), row.names = FALSE, sep="\t", quote=T)
#write.table(pr2_FINAL2, file.path( output_path,"1.12_F_Cleaned_pr2_database_wAlgbase"), row.names = FALSE)
save.image(file.path(workspace_path, "1.13_Workspace"))
load(file.path(workspace_path, "1.13_Workspace"))



