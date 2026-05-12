
rm(list=ls())
#setwd("/PATH/TO/11_Benchmark_Scripts") # to path above the Script folder NOT including it

source("00_Function_Library.R") # read file with functions
read_simple_ini("00_login_data.ini")


require(dplyr)
require(stringr)
require(phylotools)


################################################################################
## Aim: Create sequence overview including the taxonomy
## Function: Produces the final fasta file with only taxonomically validated nuclear 18S reference sequences.
######
# 1.) FASTA loading 
# 2.) Taxonomy filtering
# 3.) NCBI subset export
# 4.) Final FASTA & taxonomy output
################################################################################

##################
# data storage
input_path <- file.path(path_to_output,"01_intermediate_results" ) # results of a script
workspace_path <- file.path(path_to_output, "00_workspace")
output_path <- file.path(path_to_output,"10_FINAL_results")

#load(file.path(workspace_path,"10.3_Workspace.RData"))

path_creation(c(workspace_path, output_path))

#Read all fasta into dataframe format   ────────────────────────────────────────────────────────────────
NCBI_files <- list.files(path=file.path(path_to_output,"PR2_Sequences/NCBI/"), full.names = T)
REST_files <- list.files(path=file.path(path_to_output,"PR2_Sequences/Rest"), full.names = T)
INTEREST_files <- list.files(path=file.path(path_to_output,"PR2_Sequences/Search/"), full.names = T)


sequences <- data.frame()  
for (files in c(NCBI_files,REST_files ,INTEREST_files )){
  seq_all <- read.fasta(files)
  colnames(seq_all) <-  c("header", "sequence")
  sequences <- sequences %>%bind_rows( seq_all)
}


save.image(file.path(workspace_path,"10.1_Workspace.RData"))
#load(file.path(workspace_path,"10.1_Workspace.RData"))

#Import all sequences and taxonomy file   ────────────────────────────────────────────────────────────────
taxonomy <- read.table(file.path(input_path, "9.5_Taxonomy_FINAL2.tax"), header =  T, sep = "\t")%>%
  mutate(accession_number= iconv(accession_number, from = "", to = "ASCII//TRANSLIT"))

#Filter out sequences not in taxonomy file   ────────────────────────────────────────────────────────────────
seq_A <- sequences %>%
  rowwise() %>%
  mutate(accession = str_split_i(header, " ", 1)) %>%
  mutate(accession= iconv(accession, from = "", to = "ASCII//TRANSLIT"))%>%
  ungroup() %>%
  filter(accession %in% taxonomy$accession_number |
           sub("\\.\\d+$", "", accession) %in% taxonomy$accession_number)%>% # V1
  arrange(accession)%>%
  relocate(accession, .before = header)%>%
  unique()

######################   ────────────────────────────────────────────────────────────────
###### extract the taxonomy of ncbi entries and of those genus which include a ncbi genus
ncbi_taxonomy <- taxonomy %>%
  mutate(genus = sub(";.*", "", sub("^([^;]+;){7}", "", taxonomy)))%>%
  filter(!is.na(data_ncbi.taxonomy))

# Extract genera present in ncbi_taxonomy   ────────────────────────────────────────────────────────────────
ncbi_genera <- unique(ncbi_taxonomy$genus)

# Subset: sequences in ncbi_taxonomy OR sharing the same genus   ────────────────────────────────────────────────────────────────
taxonomy_subset <- taxonomy %>%
  mutate(genus = sub(";.*", "", sub("^([^;]+;){7}", "", taxonomy))) %>%
  filter(!is.na(data_ncbi.taxonomy) | genus %in% ncbi_genera)%>%
  left_join(seq_A %>% mutate(
    accession_number= case_when(grepl("\\.\\d+$", accession)~ gsub("\\.\\d+$", "", accession),
                                TRUE ~ accession)), by = "accession_number")

taxonomy_subset_unique <- taxonomy_subset %>%
  group_by(genus, sequence) %>%
  filter(n() == 1 | row_number() == 1) %>%
  ungroup()

# save those as fasta 
seq_NCBI <- taxonomy_subset_unique %>% 
    mutate(header_new = paste(accession, taxonomy, sep =" " ))%>% 
    select(header_new, sequence)
           
colnames(seq_NCBI) <-c("seq.name","seq.text" ) # for dat2fasta
dat2fasta(seq_NCBI, outfile = file.path(output_path,"10.05_Sequences_closeNCBI_sequences.fasta"))

taxon_list <- c("Plantae", "Chromista", "Animalia", "Bacteria","Protozoa", "Fungi", "Archaea")
counter=0
for (tax in taxon_list){
  seq_subset <- seq_NCBI %>% filter(mapply(grepl,tax , seq_NCBI$seq.name))
  counter= counter + nrow(seq_subset)
  if(nrow(seq_subset) > 0){
    dat2fasta(seq_subset, outfile = file.path(output_path,paste0("10.05_Sequences_closeNCBI_sequences__",tax,".fasta")))
  }
}; rm(seq_subset)
print(counter)

######################
### filter out taxonomy for species without sequence ###   ────────────────────────────────────────────────────────────────
tax_FINAL <- taxonomy %>% 
  inner_join(seq_A %>% select(accession) %>% 
               mutate(accession_number= gsub("\\.\\d+$", "", accession)), by="accession_number" )%>%
  relocate(accession_number, .before = accession)%>%
  select(-accession_number)%>%
  dplyr::rename(accession_number= accession)%>%
  arrange(accession_number)%>%
  unique()



print("Identical rownumber between seq and taxonomy?")
print(nrow(tax_FINAL)==nrow(seq_A)) 
### Should be fixed in Skript 3 ###


seq_B_FINAL <- seq_A %>%
  select(-header)


#Write final versions of sequence and taxonomy file    ────────────────────────────────────────────────────────────────
colnames(seq_B_FINAL) <-c("seq.name","seq.text" ) # for dat2fasta
dat2fasta(seq_B_FINAL, outfile = file.path(output_path,"10.1_Sequences_FINAL.fasta"))
write.table(tax_FINAL, file.path(output_path,"10.2_Taxonomy_FINAL_detail.tax"), col.names = F, row.names = F, sep = "\t", quote = F)

tax_FINAL_essential <- tax_FINAL%>%
  select(accession_number, taxonomy)%>%
  mutate( accession_number = iconv(accession_number, from = "", to = "ASCII//TRANSLIT"))
write.table(tax_FINAL_essential, file.path(output_path,"10.2_Taxonomy_FINAL.tax"), col.names = F, row.names = F, sep = "\t", quote = F)

save.image(file.path(workspace_path,"10.3_Workspace.RData"))

sequences_OVERVIEW <- read.table(file.path(output_path, "9.6_Overview_Sequences_FINAL.csv"), header = TRUE, sep=",")

ncbi <- sequences_OVERVIEW %>% filter(Database =="NCBI")

seq_A %>%
  filter(sapply(accession, function(x) any(grepl(x, ncbi$genbank_accession))))

seq_unknown <- sequences %>%
  rowwise() %>%
  mutate(accession = str_split_i(header, " ", 1))%>%
  mutate( accession = iconv(accession, from = "", to = "ASCII//TRANSLIT"))%>%
  ungroup() %>%
  filter(! accession %in% taxonomy$accession_number) %>% # V1
  arrange(accession)%>%
  unique() %>%
  left_join(sequences_OVERVIEW, by = c("accession"="pr2_accession"))

seq_unknown__species <- seq_unknown%>% select(Clean_Name, species, Database, Taxonomy) %>%
  unique()

write.table(seq_unknown, file.path(output_path,"10.3__Overview_species_lacking_taxonomy.tsv"), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(seq_unknown, file.path(output_path,"10.4__Overview_species_lacking_taxonomy_detailed.tsv"), col.names = T, row.names = F, sep = "\t", quote = F)
