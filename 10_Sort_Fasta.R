
rm(list=ls())
#setwd("/PATH/TO/11_Benchmark_Scripts") # to path above the Script folder NOT including it
source("00_Function_Library.R") # read file with functions

require(dplyr)
require(stringr)
require(phylotools)

##################
# data storage
output_path <- "10_FINAL_results"
workspace_path <- "00_workspace"
input_path <- "01_intermediate_results" # results of a script

path_creation(c(workspace_path, output_path))

#Read all fasta into dataframe format
NCBI_files <- list.files(path="PR2_Sequences/NCBI/", full.names = T)
REST_files <- list.files(path="PR2_Sequences/Rest", full.names = T)
INTEREST_files <- list.files(path="PR2_Sequences/Search/", full.names = T)


sequences <- data.frame()  
for (files in c(NCBI_files,REST_files ,INTEREST_files )){
  seq_all <- read.fasta(files)
  colnames(seq_all) <-  c("header", "sequence")
  sequences <- rbind(sequences, seq_all)
}

save.image(file.path(workspace_path,"10.1_Workspace"))
load(file.path(workspace_path,"10.1_Workspace"))

#Import all sequences and taxonomy file
taxonomy <- read.table(file.path(input_path, "9.5_Taxonomy_FINAL2.tax"), header =  T, sep = "\t")%>%
  mutate(accession_number= iconv(accession_number, from = "", to = "ASCII//TRANSLIT"))

#Filter out sequences not in taxonomy file
seq_A <- sequences %>%
  rowwise() %>%
  mutate(accession = str_split_i(header, " ", 1)) %>%
  mutate(accession= iconv(accession, from = "", to = "ASCII//TRANSLIT"))%>%
  filter(accession %in% taxonomy$accession_number) %>% # V1
  arrange(accession)%>%
  relocate(accession, .before = header)%>%
  unique()

#rint(colnames(seq_A))
#print(head(seq_A))
#quit()

### filter out taxonomy for species without sequence ###
tax_FINAL <- taxonomy %>% 
  filter(accession_number %in% seq_A$accession)%>%
  arrange(accession_number)%>%
  unique()

print("Identical rownumber between seq and taxonomy?")
print(nrow(tax_FINAL)==nrow(seq_A)) 
### Should be fixed in Skript 3 ###


seq_B_FINAL <- seq_A %>%
  select(-header)



#Write final versions of sequence and taxonomy file
colnames(seq_B_FINAL) <-c("seq.name","seq.text" ) # for dat2fasta
dat2fasta(seq_B_FINAL, outfile = file.path(output_path,"10.1_Sequences_FINAL.fasta"))
write.table(tax_FINAL, file.path(output_path,"10.2_Taxonomy_FINAL_detail.tax"), col.names = F, row.names = F, sep = "\t", quote = F)

tax_FINAL_essential <- tax_FINAL%>%
  select(accession_number, taxonomy)%>%
  mutate( accession_number = iconv(accession_number, from = "", to = "ASCII//TRANSLIT"))
write.table(tax_FINAL_essential, file.path(output_path,"10.2_Taxonomy_FINAL.tax"), col.names = F, row.names = F, sep = "\t", quote = F)

save.image(file.path(workspace_path,"10.3_Workspace"))


sequences_OVERVIEW <- read.table(file.path(output_path, "9.6_Overview_Sequences_FINAL.csv"), header = TRUE)

seq_unknown <- sequences %>%
  rowwise() %>%
  mutate(accession = str_split_i(header, " ", 1))%>%
  mutate( accession = iconv(accession, from = "", to = "ASCII//TRANSLIT")) %>%
  filter(! accession %in% taxonomy$accession_number) %>% # V1
  arrange(accession)%>%
  unique() %>%
  left_join(sequences_OVERVIEW, by = c("accession"="pr2_accession"))

seq_unknown__species <- seq_unknown%>% select(Clean_Name, species, Database, Taxonomy) %>%
  unique()

write.table(seq_unknown, file.path(output_path,"10.3__Overview_species_lacking_taxonomy.csv"), col.names = T, row.names = F, sep = "\t", quote = F)
write.table(seq_unknown, file.path(output_path,"10.4__Overview_species_lacking_taxonomy_detailed.csv"), col.names = T, row.names = F, sep = "\t", quote = F)
