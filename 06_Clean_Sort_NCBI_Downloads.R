
rm(list=ls())
#setwd("/PATH/TO/11_Benchmark_Scripts") # to path above the Script folder NOT including it
source("00_Function_Library.R") # read file with functions


require(dplyr)
require(stringr)
require(tidyr)
require(worrms)
require(taxonomizr)

####################
## Aim: Control is the downloaded NCBI sequences  represent the expected taxonomy
####################

#input
input_path_ncbi <- "02_NCBI_15_May_2025_results" # results of a script

# output from python script 05_Check_NCBI_with_biopython__advanced.py
ncbi_downloaded_table <- "05_Missing_PR2_Downloaded_species_NCBI_info__15_May_2025.tsv" 

##output
output_path <- "01_intermediate_results" # results of a script
workspace_path <- "00_workspace"

################################################################################
# load(file.path(workspace_path,"6.4_Workspace"))

#Import Files
pr2 <- read.table(file.path(output_path,"1.12_F_Cleaned_pr2_database_wAlgbase.tsv"), header = TRUE, sep="\t")
spec <- read.table(file.path(output_path,"2.9_F_species_FINAL_withAlgaebase.csv"), header = TRUE)
missing_spec<- read.table(file=file.path(output_path, "4.1_Missing_Species.csv"), header=T, sep=",")
                  
#ncbi
spec_new <- read.table(file.path(input_path_ncbi,ncbi_downloaded_table), sep = "\t", header=F)
colnames(spec_new) <- c("search", "species", "accession_number", "file", "ncbi.taxonomy")



#remove entries already existing in p2
spec_A <- spec_new %>%
  filter(!str_remove_all(accession_number, "\\..") %in% pr2$genbank_accession)

# clean up species names
spec_B <- spec_A %>%
  mutate(Clean_Name = str_replace_all(species, "sp\\..+", "sp.")) %>%
  mutate(Clean_Name = str_replace_all(Clean_Name, "cf$", "sp.")) %>%
  mutate(Clean_Name = str_replace_all(Clean_Name, "cf\\..+", "sp.")) %>%
  mutate(Clean_Name = str_replace_all(Clean_Name, "aff\\..+", "sp.")) %>% 
  mutate(Clean_Name = str_remove_all(Clean_Name, "uncultured ")) %>%
  rowwise() %>%
  mutate(Clean_Name = str_remove_all(Clean_Name, " [[:upper:]].+")) %>%
  mutate(Clean_Name = str_remove_all(Clean_Name, " \\d.?")) %>%
  mutate(Clean_Name = replace(Clean_Name, !str_detect(Clean_Name, " "), paste(Clean_Name, "sp.")))

#add header of fasta to dataframe
spec_X <- spec_B %>%
  group_by(search, file) %>%
  nest() %>%
  #head(15) %>%
  mutate(header = NA)

i <- 1
spec_X[i,]
for (i in c(1:nrow(spec_X))) {
  fasta <- readLines(file.path(spec_X$file[i]))
  ids = grepl(">", fasta)
  spec_X$header[i] <- list(data.frame(fasta[ids], fasta[!ids]))
}
spec_Y <- spec_X %>%
  unnest(cols = data) %>%
  unnest(cols = header) %>%
  rowwise() %>%
  mutate(length = str_length(fasta..ids.)) %>%
  dplyr::rename(header = fasta.ids.) %>%
  filter(str_detect(header, accession_number)) %>%
  select(-fasta..ids.)

#filter ncbi for non 18S seq
spec_Z <- spec_Y %>%
  filter(str_detect(header, "18S") | str_detect(header, "small")) %>%
  filter(!str_detect(header, "protein")) %>%
  filter(!str_detect(header, "\\w+ase\\W")) %>%
  filter(!str_detect(header, "peptide")) %>%
  filter(!str_detect(header, "transporter")) %>%
  filter(!((str_detect(header, "ITS") | str_detect(header, "internal transcribed spacer")) & length < 1000)) %>%
  filter(!((str_detect(header, "28S") | str_detect(header, "large")) & length < 2000)) %>%
  select(-length)

removed_spec <- spec_Y %>% filter(! accession_number %in%spec_Z$accession_number )
write.table(removed_spec, file.path(output_path,"6.3_Accessions_NCBI_removed.tsv"), sep="\t", row.names = F)

spec_B2 <- spec_Z

# unqie species and save all ncbi related info into data
spec_C <- spec_B2 %>%
  group_by(Clean_Name) %>%
  nest() %>%
  mutate(acc_Name = NA) %>%
  mutate(AphiaID = NA) %>%
  rowwise() %>%
  mutate(acc_Name = replace(acc_Name, Clean_Name %in% spec$acc_Name, Clean_Name)) %>%
  rowwise() %>%
  mutate(AphiaID = replace(AphiaID, Clean_Name %in% spec$acc_Name, subset(spec, spec$acc_Name == Clean_Name)$AphiaID[1]))

#extract species with missing worms information
print("Start spec_D")
spec_D <- spec_C %>%
  mutate(worms_record = NA)


# try to extract worms species names information
pb <- txtProgressBar(min = 1, max = nrow(spec_D), style = 3)

for (i in c(1:nrow(spec_D))) {
  if (is.na(spec_D$AphiaID[i])) {
    spec_D$worms_record[i] <- try_worms(spec_D$Clean_Name[i]) 
  }
  setTxtProgressBar(pb, i)
};close(pb); rm(pb)

saveRDS(spec_D, file.path(workspace_path,"6.2_spec_D"))
spec_D <- readRDS(file.path(workspace_path,"6.2_spec_D"))

## modJR get original taxonomy data

# Identify species with duplicated entries & taxon is partially known
print("Start duplicated_species")
duplicated_species <- missing_spec %>%
  group_by(genus,species_name) %>%
  filter(n() > 1) %>%
  filter(is.na(taxon))

# Remove rows where species is duplicated and one occurrence has not a not assigned
filtered_dataframe <- missing_spec %>%
  anti_join(duplicated_species, by =c("genus", "species_name", "taxon")) 


missing_gen <- unique(filtered_dataframe[,c("genus", "taxon")])

spec_D2 <- spec_D %>% mutate(genus = gsub("^([\\w+\\-]+).*", "\\1", Clean_Name, perl=T))%>% unique()
spec_D3 <- merge(spec_D2, missing_gen, by ="genus")


##need to change code, remove entries with non matching taxonomy
spec_E <- spec_D3 %>% 
  rowwise() %>%
  mutate(n_records = sum(!is.na(worms_record[[1]]))) %>%
  filter(n_records > 0) %>%
  # check if taxonomy fits
  mutate(worms_record = list(subset(worms_record, worms_record$kingdom == taxon |   worms_record$phylum == taxon |worms_record$class == taxon |
                                      worms_record$order == taxon |worms_record$family == taxon | is.na(worms_record$phylum)| is.na(taxon)))) %>%
  mutate(worms_record = list(subset(worms_record,  !is.na(worms_record$genus)))) %>%
  mutate(n_records = nrow(worms_record)) 


spec_E2 <- spec_E %>%
  filter(n_records > 0) %>%
  rowwise() %>%
  #keep most uptodate entry
  mutate(worms_record = list(subset(worms_record, worms_record$status == "accepted"|
                                      as.Date(worms_record$modified) == as.Date(worms_record$modified)[which.max(as.Date(worms_record$modified))]))) %>%
  #if several accepted keep most upto date
  mutate(worms_record = list(subset(worms_record, 
                                    as.Date(worms_record$modified) == as.Date(worms_record$modified)[which.max(as.Date(worms_record$modified))])))%>%
  # if still several entries keep the highest aphia id
  mutate(worms_record = list(subset(worms_record, worms_record$AphiaID == max(worms_record$AphiaID)))) %>%
  mutate(n_records = nrow(worms_record)) %>%
  bind_rows(spec_D3 %>% rowwise() %>% mutate(n_records = sum(!is.na(worms_record[[1]]))) %>% filter(n_records == 0))

spec_F <- spec_E2 %>%
  filter(n_records > 0) %>%
  rowwise() %>%
  # change to accepted name
  mutate(acc_Name = worms_record$valid_name) %>%
  mutate(AphiaID = worms_record$valid_AphiaID) %>%
  # add all data
  bind_rows(spec_E2 %>% filter(n_records == 0)) %>%
  arrange(acc_Name) %>%
  #remove worms records
  select(-worms_record, -n_records) %>%
  unnest(cols = data) %>%
  mutate(pr2_accession = NA) %>%
  relocate(pr2_accession, .before = species) %>%
  rowwise() %>%
  # calulcate replicate number for species
  mutate(dummy = acc_Name) %>%
  mutate(dummy = replace(dummy, is.na(dummy), species)) %>%
  group_by(dummy) %>%
  mutate(Replicate = row_number()) %>%
  ungroup() %>%
  select(-dummy, -search) %>%
  mutate(Database = "NCBI") %>%
  #rename and reorder columns
  dplyr::rename(genbank_accession = accession_number) %>%
  relocate(genbank_accession, .after = pr2_accession) %>%
  relocate(AphiaID, .after = genbank_accession) %>%
  dplyr::rename(Acc_Name = acc_Name) %>%
  relocate(Acc_Name, .after = AphiaID) %>%
  relocate(file, .after = Database) %>%
  mutate(sequence_length = NA) %>%
  relocate(sequence_length, .after = species) %>%
  relocate(Clean_Name, .before = species)


## remove sequences which are not having the excepted taxonomy in ncbi 
print("Start spec_F_FINAL")
spec_F_FINAL <- spec_F %>% rowwise() %>% filter(is.na(taxon)|grepl(taxon, ncbi.taxonomy))
spec_F_removed <- spec_F %>% rowwise() %>% filter(!grepl(taxon, ncbi.taxonomy))


## double check removed if taxon is just not known by ncbi
spec_F_FINAL2 <- data.frame()
if( nrow(spec_F_removed)>0){
  spec_F_removed <-spec_F_removed%>% mutate(taxon2 = taxon )%>%
    mutate(taxon2 = gsub("Bacillariophyceae", "Bacillariophyta",taxon2))%>%
    mutate(taxon2 = gsub("Charophyta", "Viridiplantae",taxon2))%>%
    mutate(taxon2 = gsub("Cyanobacteria", "Cyanobacteriota",taxon2))%>%
    mutate(taxon2 = gsub("Chlorophyta", "Chlorophyceae",taxon2))%>%
    mutate(taxon2 = gsub("Cryptophyta", "Cryptophyceae",taxon2))
  spec_F_NOTremoved <- spec_F_removed %>% rowwise() %>% filter(grepl(taxon2, ncbi.taxonomy)) %>% mutate(taxon2 =NULL)
  spec_F_removed2 <- spec_F_removed %>% rowwise() %>% filter(!grepl(taxon2, ncbi.taxonomy))
  
  spec_F_FINAL2 <- rbind(spec_F_FINAL, spec_F_NOTremoved)
}else{
  spec_F_FINAL2 <- spec_F_FINAL
}


write.table(spec_F_FINAL2, file.path(output_path,"6.3_Species_NCBI.csv"))


save.image(file.path(workspace_path, "6.4_Workspace"))

