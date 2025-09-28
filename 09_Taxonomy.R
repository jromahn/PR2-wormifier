
rm(list=ls())
#setwd("/PATH/TO/11_Benchmark_Scripts") # to path above the Script folder NOT including it
source("00_Function_Library.R") # read file with functions
read_simple_ini("00_login_data.ini")

require(dplyr)
require(stringr)
require(tidyr)
require(worrms)
require(taxonomizr)
library(jsonlite)
library(curl)
library(pr2database)


####################
## Aim: Create sequence overview including the taxonomy
####################

##################
# data storage
output_path <- "01_intermediate_results" # results of a script
final_path <- "10_FINAL_results"
workspace_path <- "00_workspace"


#Import sequences and Pr2-Taxonomy
pr2_tax <- read.table(file.path(output_path, "1.7_F_Cleanded_pr2_taxonomy.tsv"), header = TRUE, sep="\t")
sequences <- read.table(file.path(output_path,"8.1_Overview_Sequences_ALL.csv"), header = TRUE)



################################################################################
path_creation(c(workspace_path, output_path,final_path))

sequences2 <- sequences %>% 
  mutate(  Clean_Name =gsub("^([\\w\\-ë]+) ([\\w\\-ë\\.\\-]+).*", "\\1 \\2", Clean_Name, perl=T))


#Prepare Taxonomy File
tax_A <- sequences %>%
  select(-sequence_length, -Database, -Replicate) %>%
  rowwise() %>%
  mutate(accession_number = pr2_accession) %>%
  mutate(accession_number = replace(accession_number, is.na(accession_number), genbank_accession)) %>%
  select(-pr2_accession, -genbank_accession) %>%
  relocate(accession_number, .before = Clean_Name) %>%
  mutate(Clean_Genus = str_split_i(Acc_Name, " ", 1)) %>%
  mutate(Clean_Genus = replace(Clean_Genus, is.na(Clean_Genus), str_split_i(Clean_Name, " ", 1)))

#prepare table for taxonomy
tax_B <- tax_A %>%
  group_by(Clean_Name, Clean_Genus, AphiaID) %>%
  nest() %>%
  rowwise() %>%
  mutate(species = data$species[1]) %>%
  mutate(Worms_record = NA)

# create new column
tax_C <- tax_B %>%
  select(-Worms_record) %>%
  mutate(Kingdom = NA) %>%
  mutate(Phylum = NA) %>%
  mutate(Subphylum = NA) %>%
  mutate(Class = NA) %>%
  mutate(Subclass = NA) %>%
  mutate(Order = NA) %>%
  mutate(Family = NA) %>%
  mutate(Genus = NA) %>%
  mutate(Species = NA)


### create dataframe with cnbi stored taxonomy, so it has not be downloadd later
taxo_ncbi  <- sequences %>%
  select(-sequence_length, -Database, -Replicate) %>%
  rowwise() %>%
  mutate(accession_number = pr2_accession) %>%
  mutate(accession_number = replace(accession_number, is.na(accession_number), genbank_accession)) %>%
  select(-pr2_accession, -genbank_accession) %>%
  relocate(accession_number, .before = Clean_Name) %>%
  mutate(Clean_Genus = str_split_i(Acc_Name, " ", 1)) %>%
  mutate(Clean_Genus = replace(Clean_Genus, is.na(Clean_Genus), str_split_i(Clean_Name, " ", 1)))%>%
  group_by(Clean_Name, Clean_Genus, AphiaID, ncbi.taxonomy) %>%
  nest() %>%
  rowwise() %>%
  mutate(species = data$species[1]) %>%
  mutate(Worms_record = NA)%>% rowwise() %>%
  mutate(accession = data$accession_number[1]) %>%
  filter(!is.na(ncbi.taxonomy))


################ Step 1
print("Step 1 - Get taxonomy for sequences that already have an AphiaID")
#Get taxonomy for sequences that already have an AphiaID
tax_1C <- tax_C 
ranks <- colnames(tax_1C)[6:14]
print(ranks) # should be "Kingdom"   "Phylum"    "Subphylum" "Class"     "Subclass"  "Order"     "Family"    "Genus"     "Species"  

pb <- txtProgressBar(min = 1, max = nrow(tax_1C), style = 3)
 
#i = grep("Miamiensis_avidus", tax_C$species)

for (i in c(1:nrow(tax_1C))) {
  if(! is.na(tax_1C$AphiaID[i])){
    tax <- try_taxonomy(tax_1C$AphiaID[i])
    
    tax$rank <- str_replace_all(tax$rank, "Phylum \\(Division\\)", "Phylum")
    tax$rank <- str_replace_all(tax$rank, "Subphylum \\(Subdivision\\)", "Subphylum")
    
    for (k in ranks){
      if (k %in% tax$rank) {
        index <- match(k, tax$rank)
        tax_1C[i, k] <- tax$scientificname[index]
      }
    }
  }
  setTxtProgressBar(pb, i)
}; close(pb); rm(pb, i, k, index, tax)

saveRDS(tax_1C, file.path(workspace_path, "9.1_tax_1C"))
tax_1C <- readRDS(file.path(workspace_path,"9.1_tax_1C"))

## fill out ranks with downloaded taxonomy
tax_1D <- tax_1C %>%
  rowwise() %>%
  mutate(Phylum = replace(Phylum, is.na(Phylum), paste(Kingdom, "_X", sep = ""))) %>%
  mutate(Subphylum = replace(Subphylum, is.na(Subphylum), paste(Phylum, "_X", sep = ""))) %>%
  mutate(Class = replace(Class, is.na(Class), paste(Subphylum, "_X", sep = ""))) %>%
  mutate(Subclass = replace(Subclass, is.na(Subclass), paste(Class, "_X", sep = ""))) %>%
  mutate(Order = replace(Order, is.na(Order), paste(Subclass, "_X", sep = ""))) %>%
  mutate(Family = replace(Family, is.na(Family), paste(Order, "_X", sep = ""))) %>%
  mutate(Genus = replace(Genus, is.na(Genus), paste(Family, "_X", sep = ""))) %>%
  mutate(Species = replace(Species, is.na(Species), paste(Genus, "sp.", sep = " ")))


tax_1E_FINAL <- tax_1D %>%
  unnest(cols = data, names_sep = "_") %>%
  ungroup() %>%
  select(-data_Acc_Name, -Clean_Name, -Clean_Genus, -AphiaID, -data_species, -species) %>%  
  dplyr::rename("accession_number" = data_accession_number)

saveRDS(tax_1E_FINAL, file.path(workspace_path, "9.1_tax_1E_FINAL"))

#Find genera with double names
pr2_multi <- pr2_tax %>%
  rowwise() %>%
  mutate(Clean_Genus = str_split_i(Acc_Name, " ", 1)) %>%
  mutate(Clean_Genus = replace(Clean_Genus, is.na(Clean_Genus), str_split_i(Clean_Name, " ", 1))) %>%
  group_by(Clean_Genus) %>%
  mutate(multi = n_distinct(class)) %>%
  filter(multi > 1)

pr2_exp_tax <- pr2_tax %>%
  rowwise() %>%
  mutate(Clean_Genus = str_split_i(Acc_Name, " ", 1)) %>%
  mutate(Clean_Genus = replace(Clean_Genus, is.na(Clean_Genus), str_split_i(Clean_Name, " ", 1))) %>%
  group_by(Clean_Genus) %>%
  mutate(multi = n_distinct(class)) %>%
  select(Clean_Genus, division, subdivision, class) %>%unique()%>% 
  dplyr::rename(pr2_division = division,
                pr2_subdivision = subdivision,
                pr2_class = class)

saveRDS(pr2_multi, file.path(workspace_path, "9.1_pr2_multi_genera"))
pr2_multi <- readRDS(file.path(workspace_path,"9.1_pr2_multi_genera"))

################ Step 2
print("Step 2 - Sequences that do not have an AphiaID yet, but a sequence of the same genus has")
#Sequences that do not have an AphiaID yet, but a sequence of the same genus has
tax_2C <- tax_C %>%
  filter(is.na(AphiaID)) %>%
  filter(Clean_Genus %in% tax_1E_FINAL$Genus) %>%
  filter(!Clean_Genus %in% pr2_multi$Clean_Genus)

pb <- txtProgressBar(min = 1, max = nrow(tax_2C), style = 3)

for (i in c(1:nrow(tax_2C))) {
  tax <- subset(tax_1E_FINAL, tax_1E_FINAL$Genus == tax_2C$Clean_Genus[i], select = -accession_number)[1,]

  for (k in ranks) {
    tax_2C[i, k] <- replace(tax[[k]], is.null(tax[[k]]), NA)
  }

  setTxtProgressBar(pb, i)
}; close(pb); rm(pb, i, k, tax)

saveRDS(tax_2C, file.path(workspace_path, "9.1_tax_2C"))
tax_2C <- readRDS(file.path(workspace_path,"9.1_tax_2C"))

tax_2D_FINAL <- tax_2C %>%
  unnest(cols = data, names_sep = "_") %>%
  ungroup() %>%
  mutate(Species = Clean_Name) %>%
  select(-data_Acc_Name, -Clean_Name, -Clean_Genus, -AphiaID, -data_species, -species) %>%
  dplyr::rename(accession_number = data_accession_number)

saveRDS(tax_2D_FINAL, file.path(workspace_path, "9.1_tax_2D_FINAL"))

########### Step 3
#Sequences that do not have an AphiaID; Search on genus level
print("Step 3 - Sequences that do not have an AphiaID; Search on genus level")
tax_3C <- tax_C %>%
  filter(is.na(AphiaID)) %>%
  filter(!Clean_Genus %in% tax_1E_FINAL$Genus) %>%
  filter(!Clean_Genus %in% pr2_multi$Clean_Genus)

tax_3D <- tax_3C %>%
  mutate(Worms_record = NA) %>%
  rowwise() %>%
  mutate(accession = data$accession_number[1]) %>%
  mutate(database = "NCBI") %>%
  mutate(database = replace(database, str_detect(accession, "\\..+\\."), "pr2")) %>%
  group_by(AphiaID, Clean_Genus, Kingdom, Phylum, Subphylum, Class, Subclass, Order, Family, Genus, Species, Worms_record) %>%
  nest() %>%
  rowwise() %>%
  mutate(accession = data$accession[1]) %>%
  mutate(database = data$database[1]) %>%
  mutate(species = data$species[1])
  

#check for genus information in worms
pb <- txtProgressBar(min = 1, max = nrow(tax_3D), style = 3)
for (i in c(1:nrow(tax_3D))) {
  tax_3D$Worms_record[i] <- list(try_names(tax_3D$Clean_Genus[i]))
  setTxtProgressBar(pb, i)
}; close(pb); rm(pb, i)

saveRDS(tax_3D, file.path(workspace_path,"9.2_tax_3D"))
tax_3D <- readRDS(file.path(workspace_path,"9.2_tax_3D"))


#repeat check for genus information in worms
tax_3E <- tax_3D

pb <- txtProgressBar(min = 1, max = nrow(tax_3E), style = 3)
for (i in c(1:nrow(tax_3E))) {
  if(is.na(tax_3E$Worms_record[i])) {
    tax_3E$Worms_record[i] <- list(try_names(tax_3E$Clean_Genus[i]))
  }
  setTxtProgressBar(pb, i)
}; close(pb); rm(pb, i)


#repeat check for genus information in worms again
tax_3F <- tax_3E 

pb <- txtProgressBar(min = 1, max = nrow(tax_3F), style = 3)
for (i in c(1:nrow(tax_3F))) {
  if(is.na(tax_3F$Worms_record[i])) {
    tax_3F$Worms_record[i] <- list(try_names(tax_3F$Clean_Genus[i]))
  }
  setTxtProgressBar(pb, i)
}; close(pb); rm(pb, i)

saveRDS(tax_3F, file.path(workspace_path,"9.3_tax_3F"))
tax_3F <- readRDS(file.path(workspace_path,"9.3_tax_3F"))

#check success
test1 <- tax_3D %>% 
  filter(!is.na(Worms_record[[1]][1]))
test2 <- tax_3E %>%
  filter(!is.na(Worms_record[[1]][1]))
test3 <- tax_3F %>%
  filter(!is.na(Worms_record[[1]][1]))

length(test1$AphiaID)==length(test2$AphiaID)
length(test1$AphiaID)==length(test3$AphiaID)

rm(test1, test2, test3)

#####
#remove entries without worms entry
tax_3G <- tax_3F %>%
  filter(!is.na(Worms_record[[1]][1])) %>%
  rowwise() %>%
  mutate(n_record = nrow(Worms_record))

## get ncbi taxonomy and just keep if worms entry match expected taxonomy

# Define parser for NCBI taxonomy based on fixed positions
parse_taxonomy <- function(taxonomy_str, species_name) {
  tax <- unlist(strsplit(taxonomy_str, ";"))
  out <- list(
    superkingdom = ifelse(length(tax) >= 1, tax[1], NA),
    phylum = ifelse(length(tax) >= 4, tax[4], NA),
    class = ifelse(length(tax) >= 6, tax[6], NA),
    order = ifelse(length(tax) >= 10, tax[10], NA),
    family = ifelse(length(tax) >= 13, tax[13], NA),
    genus = ifelse(length(tax) >= 14, tax[14], NA),
    species = species_name,
    stringsAsFactors = FALSE
  )
  return(out)
}

print("Control NCBI sequences")
colnames(tax_3G)
tax_3H <- tax_3G %>%
  filter(database == "NCBI") %>%
  left_join(taxo_ncbi %>% select(accession,ncbi.taxonomy), by = "accession")%>%
  rowwise()%>%
  #not neceassary anymore these data are downloaded with ncbi
  #mutate(taxID = accessionToTaxa(accession, accession_ncbi)) %>% #accessionToTaxa(accession[1], accession_ncbi)
  #mutate(tax_NCBI = getTaxonomy(taxID, sqlFile = accession_ncbi)) %>%
  #new version
  mutate(tax_NCBI = list(parse_taxonomy(ncbi.taxonomy, species))) %>%
  #continued
  mutate(superkingdom = tax_NCBI$superkingdom,
        phylum       = tax_NCBI$phylum,
        class        = tax_NCBI$class,
        order        = tax_NCBI$order,
        family       = tax_NCBI$family,
        genus        = tax_NCBI$genus,
        species      = tax_NCBI$species)%>%
  select(-tax_NCBI, -ncbi.taxonomy) %>%
  dplyr::rename(phy = phylum) %>%
  dplyr::rename(cla = class) %>%
  mutate(Worms_record = list(subset(Worms_record, Worms_record$phylum == phy | Worms_record$class == cla))) %>%
  mutate(n_record = nrow(Worms_record)) %>%
  mutate(Worms_record = list(Worms_record[1,])) %>%
  mutate(n_record = nrow(Worms_record)) %>%
  #select(-taxID)%>%
  select( -superkingdom, -phy, -cla, -order, -family, -genus, -species) %>%
  bind_rows(tax_3G %>% filter(database == "pr2"))

#extract taxonomy of ncbi downloaded species from pr2
pr2_tax_short <- pr2_tax %>%
  filter(species %in% tax_3H$species) %>%
  select(domain, supergroup, division, subdivision, class, order, family, genus, species)

# filter  and just keeps entries which match with worms detected taxonomy
tax_3I_int  <- tax_3H %>%
  filter(database == "pr2") %>%
  mutate(PR2 = list(pr2_tax_short)) %>%
  dplyr::rename(spec = species) %>%
  mutate(PR2 = list(subset(PR2, PR2$species == spec))) %>%
  select(-spec) %>%
  unnest(cols = PR2) %>%
  unnest(cols = Worms_record, names_sep = "_") 


tax_3I <- filter_for_right_taxa(tax_3I_int)%>%
  dplyr::rename(da = data) %>%
  group_by(AphiaID, Clean_Genus, Kingdom, Phylum, Subphylum, Class, Subclass, Order, Family, Genus, Species, da, accession, database, n_record, domain, supergroup, 
           division, subdivision, class, order, family, genus, species) %>%
  
  nest() %>%
  dplyr::rename(Worms_record = data) %>%
  dplyr::rename(data = da) %>%
  ungroup() %>%
  select(-domain, -supergroup, -division, -subdivision, -class, -order, -family, -genus) %>%
  relocate(Worms_record, .before = data) 

rm(pr2_tax_short)

# add species which are not in 3I & are original pr2 and don't have a worm record & ncbi species
print("add species which are not in 3I & are original pr2 and don't have a worm record & ncbi species")
tax_3J <- tax_3I %>%
  bind_rows(tax_3H %>% filter(!Clean_Genus %in% tax_3I$Clean_Genus) %>% filter(database == "pr2") %>% mutate(Worms_record = NA)) %>%
  rowwise() %>%
  mutate(dummy = list(Worms_record$Worms_record_AphiaID[[1]]))

# add intermediate step because of error message
tax_3J2 <- tax_3H %>% filter(database == "NCBI") %>%
  mutate(dummy=NA)%>%
  rowwise() %>%
  mutate(dummy = list(replace(dummy, database == "NCBI", Worms_record$AphiaID))) 

#combine
tax_3J <- tax_3J %>%
  bind_rows(tax_3J2)%>%
  mutate(dummy = replace(dummy, is.null(dummy), NA)) 


tax_3K <- tax_3J %>%
  filter(!is.na(dummy)) %>%
  select(data, Kingdom, Phylum, Subphylum, Class, Subclass, Order, Family, Genus, Species, dummy)


#extract taxonomy
pb <- txtProgressBar(min = 1, max = nrow(tax_3K), style = 3)
for (i in c(1:nrow(tax_3K))) {
  tax <- try_taxonomy(tax_3K$dummy[i])
  
  tax$rank <- str_replace_all(tax$rank, "Phylum \\(Division\\)", "Phylum")
  tax$rank <- str_replace_all(tax$rank, "Subphylum \\(Subdivision\\)", "Subphylum")
  
  for (k in ranks){
    if (k %in% tax$rank) {
      index <- match(k, tax$rank)
      tax_3K[i, k] <- tax$scientificname[index]
    }
  }
  setTxtProgressBar(pb, i)
}; close(pb); rm(pb, i, k, index, tax)

saveRDS(tax_3K, file.path(workspace_path,"9.3_tax_3K"))
tax_3K <- readRDS(file.path(workspace_path,"9.3_tax_3K"))

## replace unknown taxonomy levels with _X
tax_3L <- tax_3K %>%
  rowwise() %>%
  mutate(Phylum = replace(Phylum, is.na(Phylum), paste(Kingdom, "_X", sep = ""))) %>%
  mutate(Subphylum = replace(Subphylum, is.na(Subphylum), paste(Phylum, "_X", sep = ""))) %>%
  mutate(Class = replace(Class, is.na(Class), paste(Subphylum, "_X", sep = ""))) %>%
  mutate(Subclass = replace(Subclass, is.na(Subclass), paste(Class, "_X", sep = ""))) %>%
  mutate(Order = replace(Order, is.na(Order), paste(Subclass, "_X", sep = ""))) %>%
  mutate(Family = replace(Family, is.na(Family), paste(Order, "_X", sep = ""))) %>%
  mutate(Genus = replace(Genus, is.na(Genus), paste(Family, "_X", sep = "")))

# clean up
tax_3M_FINAL <- tax_3L %>%
  select(-dummy) %>%
  unnest(cols = data) %>%
  mutate(Species = Clean_Name) %>%
  select(-Clean_Name, -accession, -species, -database) %>%
  unnest(cols = data) %>%
  select(-Acc_Name, -species)

saveRDS(tax_3M_FINAL, file.path(workspace_path,"9.4_tax_3FINAL"))
tax_3FINAL <- readRDS(file.path(workspace_path,"9.4_tax_3FINAL"))
  

############# Step 4 - Sequences with double names
print("Step 4 - Check Sequences with double names")

pr2_multi <- readRDS(file.path(workspace_path,"9.1_pr2_multi_genera"))

tax_4C <- tax_C %>%
  filter(is.na(AphiaID)) %>%
  filter(Clean_Genus %in% pr2_multi$Clean_Genus) %>%
  mutate(Worms_record = NA) 


pb <- txtProgressBar(min = 1, max = nrow(tax_4C), style = 3)
for (i in c(1:nrow(tax_4C))) {
  tax_4C$Worms_record[i] <- list(try_names(tax_4C$Clean_Genus[i]))
  
  setTxtProgressBar(pb, i)
}; close(pb); rm(pb, i)

# get pr2 taxonomy
pr2_tax_short <- pr2_tax %>%
  filter(species %in% tax_4C$species) %>%
  select(domain, supergroup, division, subdivision, class, order, family, genus, species)

tax_4D_int <- tax_4C %>%
  mutate(PR2 = list(pr2_tax_short)) %>%
  dplyr::rename(spec = species) %>%
  mutate(PR2 = list(subset(PR2, PR2$species == spec))) %>%
  unnest(cols = PR2) %>%
  unnest(cols = Worms_record, names_sep = "_") 



## check if found worms taxonomy if overlapping with expected pr2 tax
tax_4D <- filter_for_right_taxa(tax_4D_int)%>%  
  dplyr::rename(da = data) %>%
  group_by(Clean_Name, AphiaID, Clean_Genus, da, spec, Kingdom, Phylum, Subphylum, Class, Subclass, Order, Family, Genus, Species, domain, supergroup, 
           division, subdivision, class, order, family, genus, species) %>%
  
  nest() %>%
  dplyr::rename(Worms_record = data) %>%
  dplyr::rename(data = da) %>%
  ungroup() %>%
  select(-domain, -supergroup, -division, -subdivision, -class, -order, -family, -genus, -spec) %>%
  relocate(Worms_record, .before = data) 

rm(pr2_tax_short)

tax_4E <- tax_4D %>% 
  rowwise() %>%
  mutate(num = nrow(Worms_record)) %>%
  rowwise() %>%
  mutate(Worms_record = list(subset(Worms_record, Worms_record$Worms_record_AphiaID == max(Worms_record$Worms_record_AphiaID)))) %>%
  mutate(dummy = Worms_record$Worms_record_AphiaID) %>%
  select(Clean_Name, data, Kingdom, Phylum, Subphylum, Class, Subclass, Order, Family, Genus, Species, dummy)


tax_4F <- tax_4E
print(" check if found worms taxonomy if overlapping with expected pr2 tax")

pb <- txtProgressBar(min = 1, max = nrow(tax_4F), style = 3)
for (i in c(1:nrow(tax_4F))) {
  tax <- try_taxonomy(tax_4F$dummy[i])
  
  tax$rank <- str_replace_all(tax$rank, "Phylum \\(Division\\)", "Phylum")
  tax$rank <- str_replace_all(tax$rank, "Subphylum \\(Subdivision\\)", "Subphylum")
  
  for (k in ranks){
    if (k %in% tax$rank) {
      index <- match(k, tax$rank)
      tax_4F[i, k] <- tax$scientificname[index]
    }
  }
  setTxtProgressBar(pb, i)
}; close(pb); rm(pb, i, k, index, tax)

tax_4G_FINAL <- tax_4F %>%
  mutate(Species = Clean_Name) %>%
  select(-Clean_Name, -dummy) %>%
  unnest(cols = data) %>%
  select(-Acc_Name, -species) %>%
  rowwise() %>%
  mutate(Phylum = replace(Phylum, is.na(Phylum), paste(Kingdom, "_X", sep = ""))) %>%
  mutate(Subphylum = replace(Subphylum, is.na(Subphylum), paste(Phylum, "_X", sep = ""))) %>%
  mutate(Class = replace(Class, is.na(Class), paste(Subphylum, "_X", sep = ""))) %>%
  mutate(Subclass = replace(Subclass, is.na(Subclass), paste(Class, "_X", sep = ""))) %>%
  mutate(Order = replace(Order, is.na(Order), paste(Subclass, "_X", sep = ""))) %>%
  mutate(Family = replace(Family, is.na(Family), paste(Order, "_X", sep = ""))) %>%
  mutate(Genus = replace(Genus, is.na(Genus), paste(Family, "_X", sep = "")))
  
saveRDS(tax_4G_FINAL, file.path(workspace_path,"9.3_tax_4G_FINAL"))
save.image(file.path(workspace_path,"9.7_Workspace_before_algaebase"))


########################## Check Species with Algaebase ########################## 
print("Step 6 - Check Species with Algaebase ")

# important note - Clean_Genus doesn't have to be identical with genus of Clean_Name!!!
tax_6C <- tax_C %>%
  rowwise() %>%
  mutate(accession = data$accession_number[1]) %>%
  mutate(database = "NCBI") %>%
  mutate(database = replace(database, str_detect(accession, "\\..+\\."), "pr2"))%>%
  mutate(Clean_Name2 = gsub("\\(.*\\) ","", Clean_Name , perl=T)) %>%
  mutate(Clean_Epithel= gsub("^([\\w\\-ë]+) ([\\w\\-ë\\.]+).*", "\\2", Clean_Name2, perl=T)) %>%
  mutate(taxon = data$taxon[1])

########################
#check all species names again if they are not sp.
tax_6D <- tax_6C %>% select(Clean_Genus, Clean_Epithel, taxon, accession, Clean_Name, Clean_Name2)%>% 
          filter(Clean_Epithel !="sp.")%>% filter(Clean_Epithel !="") 

# 1.) check species name
tax_6D_search <- tax_6D%>% select(Clean_Genus, Clean_Epithel, taxon)%>% unique()
tax_6D.list <- check_species_against_Algaebase(tax_6D_search, key, password)
saveRDS(tax_6D.list, file= file.path(workspace_path, "tax_6D_algaebase"))

# clean up the results
algaebase.species.df <- tax_6D.list[[1]] %>% unique()
tax_6E <- tax_6D %>% mutate(request = paste(Clean_Genus, Clean_Epithel), taxon= taxon)  %>%
      right_join(algaebase.species.df, by=c("request", "taxon"),relationship = "many-to-many")
tax_6F <- algaebase_species_clean_up(tax_6E)


# 2.) search taxonomy of genus for known species names
tax_6F_search <- tax_6F %>% group_by(taxon, acc_Genus) %>% select(taxon, acc_Genus) %>% unique()
tax_6F.list <- check_genus_n_taxonomy_algaebase(tax_6F_search, key, password)
saveRDS(tax_6F.list, file= file.path(workspace_path, "tax_6F_algaebase_known"))

# just keep accepted data because we are using valid species names
tax_6G <- tax_6F.list[[1]] %>% unique()
tax_6G.accepted <- tax_6G %>% 
  filter(request_genus == `dwc:genus`)  %>% 
  filter(`dwc:taxonomicStatus` =="currently accepted taxonomically") %>%
  dplyr::rename("acc_Genus"=request_genus,
         "taxon" = request_taxon,
         "algaebase.scientificName.genus" =  `dwc:scientificName`)%>% 
  unique()

#if Neotype & Holotype exists keep accepted
tax_6G.accepted <- tax_6G.accepted %>%
  group_by(taxon, acc_Genus)%>%
  mutate(n_records =n())

tax_6G.accepted <- tax_6G.accepted%>%
  filter(n_records >1)%>%
  group_by(taxon, acc_Genus)%>%
  filter(`dwc:nomenclaturalStatus`=="nom. et typ. cons.")%>%
  bind_rows(tax_6G.accepted %>% filter(n_records == 1))%>%
  mutate(n_records=NULL)

tax_6H <- tax_6G.accepted %>% 
  mutate(`dwc:genus` =NULL,`dwc:taxonRank` =NULL,`dwc:taxonomicStatus`=NULL, `dwc:acceptedNameUsage` = NULL, `dwc:acceptedNameUsageID`=NULL)
saveRDS(tax_6H, file= file.path(workspace_path, "tax_6H_algaebase_known"))

## transform dataframe & just keep algae entries & combine with species df
tax_6H <- tax_6H %>% unique()
tax_6H.algae <- filter_algaebase_for_algae(tax_6H)


tax_6I <- tax_6H.algae %>% left_join( tax_6F, by =c("taxon", "acc_Genus"), suffix=c(".genus", ".species") )

# filter if taxon hits expected pr2 taxonomy
pr2_short_tax <- pr2_tax %>%  select(division, subdivision, class, Clean_Name)%>%
  dplyr::rename(pr2_division =division, 
                pr2_subdivision =subdivision, 
                pr2_class=class)%>%
  filter(Clean_Name %in% tax_6I$Clean_Name)%>%
  unique()%>% group_by(Clean_Name)%>%
  filter(n()==1)
tax_6I_pr2 <-tax_6I %>%   left_join(pr2_short_tax, by= "Clean_Name")
tax_6I.FINAL <- filter_for_right_taxa_algaebase(tax_6I_pr2)

 
## create genus df for further analysis
tax_6H_avd <- tax_6F %>% 
  group_by(taxon, acc_Genus, Clean_Genus) %>%
  select(taxon, acc_Genus, Clean_Genus) %>% 
  unique()%>% 
  right_join( tax_6H, by =c("taxon", "acc_Genus") )

## identify multi renamed
tax_6H_avd.multi <- tax_6H_avd  %>%group_by(Clean_Genus) %>% mutate( taxon = NULL)
tax_6H_avd.multi <- tax_6H_avd.multi %>%group_by(Clean_Genus) %>%
  mutate(renamed = n())  %>%  filter(Clean_Genus != acc_Genus) %>% unique() %>% filter(renamed > 1)


######################## 
## species with known genus but not multi genus and not multi renames
tax_6Inon <-  tax_6C %>% 
  select(Clean_Genus, Clean_Epithel, taxon, accession, Clean_Name, Clean_Name2) %>%
  filter(Clean_Genus %in% tax_6F$Clean_Genus)%>%
  filter(! Clean_Name %in% tax_6I$Clean_Name)%>%
  filter(! Clean_Genus %in% pr2_multi$Clean_Genus)%>%
  filter(! Clean_Genus %in% tax_6H_avd.multi$Clean_Genus)

tax_6J_non <- merge(tax_6H_avd, tax_6Inon, by = c("taxon", "Clean_Genus")) %>%
  mutate( acc_Name = paste(acc_Genus, Clean_Epithel))


# filter if taxon hits expected pr2 taxonomy
pr2_short_tax <- pr2_tax %>%  select(division, subdivision, class, Clean_Name)%>%
  dplyr::rename(pr2_division =division, 
                pr2_subdivision =subdivision, 
                pr2_class=class)%>%
  filter(Clean_Name %in% tax_6J_non$Clean_Name)%>%
  unique()%>% group_by(Clean_Name)%>%
  filter(n()==1)

tax_6J_non <-tax_6J_non %>%   left_join(pr2_short_tax, by= "Clean_Name")
tax_6J_non.FINAL <- filter_for_right_taxa_algaebase(tax_6J_non)

#filter for algae
tax_6J_non.FINAL.algae <- filter_algaebase_for_algae(tax_6J_non.FINAL) 

########################
## check for other species with so far unknown genera if genera is known in algaebase
tax_6J <-  tax_6C %>% 
  select(Clean_Genus, Clean_Epithel, taxon, accession, Clean_Name, Clean_Name2) %>%
  filter(! Clean_Genus %in% tax_6F$Clean_Genus)%>%
  filter(! Clean_Genus %in% pr2_multi$Clean_Genus)

# search on algaebase
tax_6J.genus <- tax_6J %>% select( taxon, Clean_Genus ) %>% unique()
tax_6J.list <- check_genus_n_taxonomy_algaebase(tax_6J.genus, key, password)
saveRDS(tax_6J.list, file= file.path(workspace_path, "tax_6J_algaebase_known"))

# filter out non algae and non phytoplankton
tax_6K <- tax_6J.list[[1]] %>% unique()%>% 
  filter(request_genus == `dwc:genus`) %>%
  dplyr::rename(Clean_Genus=request_genus, taxon = request_taxon)

#check if its matching own expected taxonomy
tax_6K <- tax_6K %>% filter( is.na(taxon)|
                            mapply( grepl, taxon , paste(`dwc:phylum`, `dwc:class`, `dwc:order`, sep=";"))|
                            (taxon=="Ochrophyta"& `dwc:phylum` =="Heterokontophyta"))

# check if the found taxonomy match with the expected pr2 taxonomy
tax_6K_pr2 <-merge(tax_6K, pr2_exp_tax, by = "Clean_Genus") 
tax_6K_pr2 <- filter_for_right_taxa_algaebase(tax_6K_pr2) %>% 
                            unique()%>%
                            mutate(pr2_division =NULL, pr2_subdivision =NULL, pr2_class =NULL)



# filter only for algae & phytoplankton
tax_6K.algae <- filter_algaebase_for_algae(tax_6K_pr2) 


##change to accapted names
tax_6K.accepted <- tax_6K.algae %>% 
  filter(Clean_Genus == `dwc:genus`)  %>% 
  filter(`dwc:taxonomicStatus` =="currently accepted taxonomically") %>%
  mutate(acc_Genus=Clean_Genus) %>%
  dplyr::rename("algaebase.scientificName.genus" =  `dwc:scientificName`)%>% 
  unique()

#if Neotype & Holotype exists keep accepted
tax_6K.accepted <- tax_6K.accepted %>%
  group_by(taxon, acc_Genus)%>%
  mutate(n_records =n())

tax_6K.accepted <- tax_6K.accepted%>%
  filter(n_records >1)%>%
  group_by(taxon, acc_Genus)%>%
  filter(`dwc:nomenclaturalStatus`=="nom. et typ. cons.")%>%
  bind_rows(tax_6K.accepted %>% filter(n_records == 1))%>%
  mutate(n_records=NULL)

## nonaccepted
tax_6K.nonacc<- tax_6K.algae %>% 
  filter(`dwc:taxonomicStatus` !="currently accepted taxonomically") %>%
  filter(! Clean_Genus %in% tax_6K.accepted$acc_Genus)%>% 
  mutate(acc_Genus = ifelse(is.na(`dwc:acceptedNameUsage`), `dwc:genus`, `dwc:acceptedNameUsage`),
         algaebase.scientificName.genus = ifelse(is.na(`dwc:scientificName`), "Uncertain", `dwc:scientificName`))%>%
  mutate(`dwc:scientificName` =NULL)

tax_6K.cleaned <- rbind(tax_6K.accepted, tax_6K.nonacc)


## merge
tax_6L <- merge(tax_6J, tax_6K.cleaned, by =c("taxon","Clean_Genus"))%>%
  mutate( acc_Name = paste(acc_Genus, Clean_Epithel))

########################
# what to do with multi renamed species/genera??

tax_6M <- tax_6C %>% 
  select(Clean_Genus, Clean_Epithel, taxon, accession, Clean_Name, Clean_Name2) %>%
  filter(! Clean_Name %in% tax_6F$Clean_Name)%>%
  filter( Clean_Genus %in% pr2_multi$Clean_Genus)

pr2_short_tax <- pr2_tax %>%  select(division, subdivision, class, Clean_Name)%>%
  dplyr::rename(pr2_division =division, 
                pr2_subdivision =subdivision, 
                pr2_class=class)%>%
  filter(Clean_Name %in% tax_6M$Clean_Name)%>%
  unique()%>% group_by(Clean_Name)%>%
  filter(n()==1)

#filter for algae
tax_6M.pr2 <- tax_6M %>%select(Clean_Genus, Clean_Epithel,taxon,Clean_Name, accession)%>%unique() %>% 
  left_join(pr2_short_tax, by ="Clean_Name") %>%
  filter(pr2_division =="Rhodophyta" | pr2_division=="Streptophyta"| pr2_subdivision=="Gyrista")

#search
tax_6M.pr2.search <- tax_6M.pr2 %>% select(taxon, Clean_Genus) %>% unique()
tax_6M.list <- check_genus_n_taxonomy_algaebase(tax_6M.pr2.search, key, password)

#check if fitting expected
tax_6N <- tax_6M.list[[1]] %>% unique()%>%
  filter(`dwc:genus`==request_genus)%>%
  dplyr::rename(Clean_Genus=request_genus,
                taxon = request_taxon)

tax_6O <- merge(tax_6N, tax_6M.pr2, by=c("taxon", "Clean_Genus"))
tax_6O.pr2 <- filter_for_right_taxa_algaebase(tax_6O )%>%
                unique()

#change to accepted names

tax_6O.accepted <- tax_6O.pr2 %>% 
  filter(Clean_Genus == `dwc:genus`)  %>% 
  filter(`dwc:taxonomicStatus` =="currently accepted taxonomically") %>%
  mutate(acc_Genus=Clean_Genus) %>%
  dplyr::rename("algaebase.scientificName.genus" =  `dwc:scientificName`)%>% 
  unique()

#if Neotype & Holotype exists keep accepted
tax_6O.accepted <- tax_6O.accepted %>%
  group_by(taxon, acc_Genus)%>%
  mutate(n_records =n())

tax_6O.accepted <- tax_6O.accepted%>%
  filter(n_records >1)%>%
  group_by(taxon, acc_Genus)%>%
  filter(`dwc:nomenclaturalStatus`=="nom. et typ. cons.")%>%
  bind_rows(tax_6O.accepted %>% filter(n_records == 1))%>%
  mutate(n_records=NULL)

## unaccepted
tax_6O.nonacc<- tax_6O.pr2 %>% 
  filter(`dwc:taxonomicStatus` !="currently accepted taxonomically") %>%
  filter(! Clean_Genus %in% tax_6O.accepted$acc_Genus)%>% 
  mutate(acc_Genus = as.character(ifelse(is.na(`dwc:acceptedNameUsage`), `dwc:genus`, `dwc:acceptedNameUsage`)),
         algaebase.scientificName.genus = as.character(ifelse(is.na(`dwc:scientificName`), "Uncertain", `dwc:scientificName`)))%>%
  mutate(`dwc:scientificName` =NULL)

#final
tax_6O.cleaned <- rbind(tax_6O.accepted, tax_6O.nonacc)%>% 
  mutate(pr2_division =NULL, pr2_subdivision=NULL, pr2_class =NULL)%>%
  mutate(acc_Name= paste(acc_Genus, Clean_Epithel))


### final datasets
colnames_of_interest <- c("taxon" ,"Clean_Genus", "Clean_Epithel", "Clean_Name" ,"Clean_Name2","accession",
                          "dwc:kingdom","dwc:phylum", "dwc:class","dwc:order" ,"dwc:family" ,
                           "acc_Genus" , "acc_Name","dwc:taxonRank", "typeSpeciesId","date.genus" ,
                          "dwc:isFossil","dwc:isFreshwater","dwc:isMarine","dwc:isTerrestrial"  )

tax_6P <- plyr::rbind.fill(tax_6O.cleaned,tax_6L )
tax_6P <- plyr::rbind.fill(tax_6P,tax_6J_non.FINAL.algae )
tax_6P <- plyr::rbind.fill(tax_6P,tax_6I.FINAL )

tax_6P  <- tax_6P %>% 
  unique()%>%
  select(all_of(colnames_of_interest))%>% 
  left_join(tax_6C %>% select(Clean_Genus, Clean_Epithel, taxon, Clean_Name, accession, data), 
            by= c("Clean_Genus", "Clean_Epithel", "taxon", "Clean_Name","accession"))%>%
  mutate(Clean_Genus=NULL,Clean_Epithel=NULL, taxon=NULL,Clean_Name2=NULL) %>% unique()


## have to unpack all accession numbers & reformat to style of Stefanie
tax_6Q <-filter_algaebase_for_algae(tax_6P) %>%
  mutate(Species = acc_Name,
         Genus= acc_Genus) %>%
  select(-Clean_Name, -acc_Name, -acc_Genus ) %>%
  unnest(cols = data, names_sep = "_") %>%
  select(-data_Acc_Name, -data_species) %>%
  dplyr::rename("Kingdom" = `dwc:kingdom` ,
                "Phylum"= `dwc:phylum`,
                "Class" = `dwc:class`,
                "Order"= `dwc:order`, 
                "Family"= `dwc:family`)%>%
  filter(is.na(data_ncbi.taxonomy)|mapply( grepl, Class , data_ncbi.taxonomy))
  
tax_6R <- tax_6Q %>%
  mutate(Subphylum= paste(Phylum, "X", sep="_"),
         Subclass= paste(Class, "X", sep="_"))%>%
  dplyr::rename(accession_number = data_accession_number)%>%
  mutate(accession=NULL)%>%
  relocate(accession_number, , .before = everything())%>%
  relocate( data_genus, .after = accession_number)%>%
  relocate( data_ncbi.taxonomy, .after=data_genus)%>%
  relocate( data_header , .after=data_ncbi.taxonomy)%>%
  relocate(data_taxon , .after=data_header)%>%
  relocate( Subphylum, .after = Phylum) %>% 
  relocate(Subclass, .after= Class)

write.table(tax_6R, file = file.path(final_path, "PR2_metadata_Algaebase.csv"), row.names = F)

#####

tax_6R_FINAL <-tax_6R %>% select(all_of(colnames(tax_1E_FINAL))) 
saveRDS(tax_6R_FINAL, file.path(workspace_path, "9.1_tax_6R_FINAL"))


##########################################################################################   
print("Final steps")
############# Sequences without taxonomy
save.image(file.path(workspace_path,"9.6_Workspace"))
#load(file.path(workspace_path,"9.6_Workspace"))


# so far all species without AphiaID are still included
tax_1E_FINAL <- tax_1E_FINAL %>% filter(!is.na(Kingdom))

tax_5B_FINAL <- tax_A %>%
  filter(!accession_number %in% tax_1E_FINAL$accession_number) %>% # already have an AphiaID
  filter(!accession_number %in% tax_2D_FINAL$accession_number) %>% # not have an AphiaID yet, but a sequence of the same genus has
  filter(!accession_number %in% tax_3M_FINAL$accession_number) %>% # do not have an AphiaID; Search on genus level
  filter(!accession_number %in% tax_4G_FINAL$accession_number) %>%    # Sequences with double names
  filter(!accession_number %in% tax_6R_FINAL$accession_number) 

############# Subtract Algaebase hits
tax_1E_FINAL.2 <-tax_1E_FINAL %>% filter(!accession_number %in% tax_6R_FINAL$accession_number) 
tax_2D_FINAL.2 <-tax_2D_FINAL %>% filter(!accession_number %in% tax_6R_FINAL$accession_number) 
tax_3M_FINAL.2 <-tax_3M_FINAL %>% filter(!accession_number %in% tax_6R_FINAL$accession_number) 
tax_4G_FINAL.2 <-tax_4G_FINAL %>% filter(!accession_number %in% tax_6R_FINAL$accession_number) 

############# Combine
tax_FINAL <- bind_rows(tax_1E_FINAL.2, tax_2D_FINAL.2, tax_3M_FINAL.2, tax_4G_FINAL.2,tax_6R_FINAL, tax_5B_FINAL)
tax_EXIST <- bind_rows(tax_1E_FINAL.2, tax_2D_FINAL.2, tax_3M_FINAL.2, tax_4G_FINAL.2, tax_6R_FINAL)

############# Write Taxonomy File

taxonomy <- tax_EXIST %>%
  select(- ncbi.taxonomy, -genus, -header)%>%
  unite(taxonomy, Kingdom, Phylum, Subphylum, Class, Subclass, Order, Family, Genus, Species, sep = ";") %>%
  mutate(taxonomy = paste(taxonomy, ";", sep = "")) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "X_X", "XX")) %>%
  mutate(taxonomy = str_replace_all(taxonomy, " ", "_")) %>%
  arrange(accession_number)

write.table(taxonomy, file.path(output_path,"9.5_Taxonomy_FINAL.tax"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(taxonomy, file.path(output_path,"9.5_Taxonomy_FINAL2.tax"), col.names = T, row.names = FALSE, sep = "\t", quote = T)

############# Add Taxonomy information to overview file
sequences_FINAL <- sequences %>%
  select(- ncbi.taxonomy, -genus, -header)%>%
  mutate(Taxonomy = NA) %>%
  mutate(Taxonomy = replace(Taxonomy, (pr2_accession %in% tax_1E_FINAL.2$accession_number | genbank_accession %in% tax_1E_FINAL.2$accession_number), 1)) %>%
  mutate(Taxonomy = replace(Taxonomy, (pr2_accession %in% tax_2D_FINAL.2$accession_number | genbank_accession %in% tax_2D_FINAL.2$accession_number), 2)) %>%
  mutate(Taxonomy = replace(Taxonomy, (pr2_accession %in% tax_3M_FINAL.2$accession_number | genbank_accession %in% tax_3M_FINAL.2$accession_number), 3)) %>%
  mutate(Taxonomy = replace(Taxonomy, (pr2_accession %in% tax_4G_FINAL.2$accession_number | genbank_accession %in% tax_4G_FINAL.2$accession_number), 4)) %>%
  mutate(Taxonomy = replace(Taxonomy, (pr2_accession %in% tax_5B_FINAL$accession_number | genbank_accession %in% tax_5B_FINAL$accession_number), 5)) %>%
  mutate(Taxonomy = replace(Taxonomy, (pr2_accession %in% tax_6R_FINAL$accession_number | genbank_accession %in% tax_6R_FINAL$accession_number), 6)) %>%
  rowwise() %>%
  mutate(dummy = Acc_Name) %>%
  mutate(dummy = replace(dummy, is.na(dummy), Clean_Name)) %>%
  arrange(Replicate) %>%
  arrange(dummy) %>%
  select(-dummy)

write.table(sequences_FINAL,file.path(final_path, "9.6_Overview_Sequences_FINAL.csv"), col.names = TRUE, row.names = FALSE)
  
save.image(file.path(workspace_path,"9.7_Workspace"))


diff <- sequences_FINAL %>%filter( Clean_Name != Acc_Name) %>%select(Clean_Name, Acc_Name, species) %>% unique()
