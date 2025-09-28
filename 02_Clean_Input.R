rm(list=ls())
#setwd("/PATH/TO/11_Benchmark_Scripts") # to path above the Script folder NOT including it
source("00_Function_Library.R") # read file with functions
read_simple_ini("00_login_data.ini")

####################
## Aim: Pre-cleaning of pr2 database
##
## written by Stefanie Knell
## adapted by Juliane Romahn - 29.jan 2024
####
# We can list the libraries that are actually loaded doing
(.packages())
unloadNamespace( "pr2database" )
require(dplyr)
require(stringr)
require(worrms)
require(tidyr)
library(jsonlite) #needed for algaebase
library(curl) # needed for algaebase


##################
# data storage
output_path <- "01_intermediate_results" # results of a script
workspace_path <- "00_workspace"
##

path_creation(c(workspace_path, output_path))

species_list <- read.table(file.path(input_path, input_table), header = TRUE, sep = ",")

#########
spec_A <- species_list %>%
  unique()%>%
  mutate(species = str_replace_all(species, "_", " ")) %>%
  mutate(worms_record = NA) %>%
  mutate(n_records = NA) %>%
  mutate(NR= c(1:length(species)), status= "aktuell")

#spec_A<- spec_A[c(1:100, 4400:4500),]

#check if species has worm record
pb <- txtProgressBar(min = 1, max = nrow(spec_A), style = 3)

for (i in c(1:nrow(spec_A))) {
  spec_A$worms_record[i] <- try_worms(spec_A$species[i])
  setTxtProgressBar(pb, i)
  
  if (sum(!is.na(spec_A$worms_record[[i]])) > 0) {
    spec_A$n_records[i] <- nrow(spec_A$worms_record[[i]])
  }else{
    spec_A$n_records[i] <- 0
  }
  
}
close(pb); rm(pb)

saveRDS(spec_A,file = file.path(workspace_path, "2.1_spec_A"))
spec_A <- readRDS(file.path(workspace_path,"2.1_spec_A"))

# remove non-target taxon or keep all if taxon is not set
print("Start spec_B")
spec_B <- spec_A %>%
  filter(n_records > 0) %>%
  rowwise() %>%
  mutate(worms_record = list(subset(worms_record, mapply(grepl, taxon,  paste(kingdom,phylum,class,order,family , sep=";"))|is.na(taxon)))) %>%
  mutate(n_records = nrow(worms_record)) %>%
  bind_rows(spec_A %>% filter(n_records == 0))

#keep most uptodate aphiaID 
spec_B2 <- spec_B %>%
  filter(n_records > 0) %>%
  rowwise() %>%
  ## accepted or most up to date
  mutate(worms_record = list(subset(worms_record, worms_record$status == "accepted"|
                                      as.Date(worms_record$modified) == as.Date(worms_record$modified)[which.max(as.Date(worms_record$modified))]))) %>%
  #if several accepted keep most upto date
  mutate(worms_record = list(subset(worms_record, 
                                    as.Date(worms_record$modified) == as.Date(worms_record$modified)[which.max(as.Date(worms_record$modified))])))%>%
  # if still several entries keep the highest aphia id
  mutate(worms_record = list(subset(worms_record, worms_record$AphiaID == max(worms_record$AphiaID)))) %>%
  mutate(n_records = nrow(worms_record)) %>%
  bind_rows(spec_B %>% filter(n_records == 0))

## add two rows
spec_C <- spec_B2 %>%
  mutate(acc_Name = NA) %>%
  mutate(AphiaID = NA) 

## write valid names from df to db
spec_D <- spec_C %>%
  filter(n_records > 0) %>%
  rowwise() %>%
  mutate(AphiaID = worms_record$valid_AphiaID) %>%
  mutate(acc_Name = worms_record$valid_name) %>%
  mutate(acc_Name = worms_record$valid_name) %>%
  mutate(taxonomy = paste(worms_record$kingdom,worms_record$phylum,worms_record$class, sep=";")) %>%
  bind_rows(spec_C %>% filter(n_records == 0)) %>%
  select(-worms_record, - n_records) %>%
  arrange(NR)


saveRDS(spec_D,file = file.path(workspace_path, "2.1_spec_D"))
spec_D <- readRDS(file.path(workspace_path,"2.1_spec_D"))
taxonomy <- spec_D[,c("NR", "genus", "species", "acc_Name","taxon", "taxonomy")]
spec_D$taxonomy <- NULL

write.table(spec_D, file.path(workspace_path,"2.2_spec_D.csv"), row.names = FALSE)


## create column to keep synoynms
print("Start spec_E")
spec_E <- spec_D %>%
  group_by(NR) %>%
  nest() %>%
  rowwise() %>%
  mutate(akt = list(subset(data, data$status == "aktuell"))) %>%
  mutate(syn = list(subset(data, data$status == "synonym"))) %>%
  mutate(n_syn = nrow(syn)) %>%
  mutate(syn = list(subset(syn, syn$acc_Name == akt$acc_Name)))

pb <- txtProgressBar(min = 1, max = nrow(spec_E), style = 3)

##search for synoymns
spec_E <- spec_E %>% mutate(new_syn =NA)
for (i in c(1:nrow(spec_E))) {
  spec_E$new_syn[i] <- try_synonyms(spec_E[i,]$akt[[1]]$AphiaID) 
  setTxtProgressBar(pb, i)
}

close(pb); rm(pb)

saveRDS(spec_E, file.path(workspace_path,"2.3_spec_E"))
spec_E <- readRDS(file.path(workspace_path,"2.3_spec_E"))

print("Start spec_FINAL")
spec_FINAL <- spec_E %>%
  unnest(cols = new_syn) %>%
  mutate(Rec_syn = try_record(new_syn)) %>%
  unnest(cols = Rec_syn) %>%
  select(NR, data, akt, syn,  scientificname) %>%
  rowwise() %>%
  mutate(genus = NA) %>%
  dplyr::rename(species = scientificname) %>%
  mutate(status = "synonym") %>%
  mutate(acc_Name = akt$acc_Name) %>%
  mutate(AphiaID = akt$AphiaID) %>%
  relocate(genus, .before = species) %>%
  dplyr::rename(old = data) %>%
  group_by(NR, old, akt, syn) %>%
  nest() %>%
  rowwise() %>%
  mutate(new = list(bind_rows(akt, syn, data))) %>%
  mutate(n_all = nrow(new)) %>%
  mutate(new = list(distinct(new, species, .keep_all = TRUE))) %>%
  mutate(n_short = nrow(new)) %>%
  ungroup() %>%
  select(NR, new) %>%
  unnest(cols = new) %>%
  filter(!is.na(species))

write.table( spec_FINAL, file.path(output_path,"2.4_F_species_FINAL.csv"))

########### double check phytoplankton species with algaebase
#Chlorophyta #Rhodophyta Charophyta Ochrophyta Myzozoa Cryptophyta Haptophyta Cyanobacteria Tracheophyta Bacillariophyta Euglenozoa
print("Start Algaebase")
  
species_phyta <- taxonomy %>% filter(
  mapply(grepl, "Chlorophyta|Rhodophyta|Charophyta|Ochrophyta|Myzozoa|Cryptophyta|Haptophyta|Cyanobacteria|Tracheophyta|Bacillariophyta|Euglenozoa",  taxonomy)|
  mapply(grepl, "Chlorophyta|Rhodophyta|Charophyta|Ochrophyta|Myzozoa|Cryptophyta|Haptophyta|Cyanobacteria|Tracheophyta|Bacillariophyta|Euglenozoa",  taxon) )


######### prepare query
species_phyta <-species_phyta %>%mutate(Full_search= ifelse(is.na(acc_Name), species,acc_Name )) %>%
  mutate(Alg.genus=gsub("^([\\w-]+) .*", "\\1", Full_search, perl = T), 
         Alg.epithel= gsub("^[\\w-]+ ?([\\w\\-\\.]*).*", "\\1", Full_search, perl = T))%>%
  mutate(Search.taxon= ifelse(is.na(taxon), gsub("^\\w+;(\\w+);.*", "\\1", taxonomy, perl = T), taxon))

#query
search <- unique(species_phyta[,c("Alg.genus", "Alg.epithel", "Search.taxon", "Full_search", "NR")]) %>% filter(Alg.epithel !="")
algaebase_species <- check_species_against_Algaebase(search, key, password)
save(algaebase_species, file= file.path(output_path, "SpeciesList__algaebase_info.RData"))


algaebase.species.df <- algaebase_species[[1]]
algaebase.species.df <- unique(algaebase.species.df)

search$request <- paste(search$Alg.genus, search$Alg.epithel)
search$taxon<- search$Search.taxon
search <- search[,c("request","taxon", "Full_search", "NR")]
algaebase.species.df <- merge(algaebase.species.df, search, by=c("request", "taxon"), all =T)

#cleaning results
algaebase.species.df <-algaebase_species_clean_up(algaebase.species.df)
write.table(algaebase.species.df, file.path( output_path,"2.8_F_Algaebase_specieslist.csv"), row.names = FALSE)

## combine with specFINAL
print("Start spec_FINAL_new")
spec_FINAL_algaebase <- algaebase.species.df[,c("NR","acc_Name")]
spec_FINAL_new <- merge(spec_FINAL,spec_FINAL_algaebase, by = "NR" , all =T)


spec_FINAL_new<- spec_FINAL_new %>% mutate(acc_Name= ifelse(is.na(acc_Name.y), acc_Name.x, acc_Name.y))

new_synoyms <- spec_FINAL_new %>%filter(acc_Name.y!=acc_Name.x & status=="aktuell")
new_synoyms <-new_synoyms %>% mutate(genus= NA, taxon =NA,AphiaID= NA, status = "synonym", species= acc_Name.x)

spec_FINAL_new <- rbind(spec_FINAL_new, new_synoyms)%>% arrange(NR) %>% mutate(acc_Name.x=NULL, acc_Name.y = NULL)
write.table( spec_FINAL_new, file.path(output_path,"2.9_F_species_FINAL_withAlgaebase.csv"))

save.image(file.path(workspace_path, "2.10_Workspace"))
load(file.path(workspace_path, "2.10_Workspace"))
