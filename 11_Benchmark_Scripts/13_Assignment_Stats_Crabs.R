################################################################################
#This script can be used to asses how the different database versions performed
#regarding taxonomic identification
################################################################################
#Run after mothur
################################################################################
rm(list = ls())

require(tidyverse)
require(stringr)
require(ggplot2)
require(tibble)
library(gapminder) 

options(scipen = 999) # stop scientific notation
#setwd("/PATH/TO/02_Scripts_folder") # to path above the Script folder
setwd("/Users/juliane/Documents/00_Work_SGN/00_PhytoArk/XX_PAPERS/2025_stefanie_wormifier/00_code_new")
source("00_Function_Library.R") # read file with functions

## create output folder
input_path <- "12_MOTHUR_ASSIGNMENTS"
output_path <- "13_Analyses_2025"
if (!dir.exists(output_path)){ dir.create(output_path) }

cleaned_community_file <- "00_input/Supplementary_Table_5_CommunityMatrix_Cleaned.csv"
mothur_results_db1 <- file.path(input_path, "9_PhytoArk_Euka02_MUC__Taxonomy_amplicons__DB1.wang.taxonomy")
mothur_results_db2 <- file.path(input_path,"9_PhytoArk_Euka02_MUC__Taxonomy_amplicons__DB2.wang.taxonomy")
mothur_results_db3 <- file.path(input_path,"9_PhytoArk_Euka02_MUC__Taxonomy_amplicons__DB3.wang.taxonomy")
mothur_results_db4 <- file.path(input_path,"9_PhytoArk_Euka02_MUC__Taxonomy_amplicons__DB4.wang.taxonomy")


read_simple_ini("00_login_data.ini")

## list of Baltic Species
baltic_list <- file.path(input_path, input_table)
baltic_ciliates <- read.table(baltic_list, header = T, sep=",") %>% filter(taxon=="Ciliophora")

#######################################################################################
### Import counttable (cleaned version)

counttable <- read.table(cleaned_community_file, header = TRUE, sep = ",") %>%
  select(-c(sample_id, tag, station, depth, replicate, total_reads)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  rowwise() %>%
  mutate(total_reads = sum(across(starts_with("V")))) %>%
  select(ASV, total_reads) %>%
  dplyr::rename(V1 = ASV)


###
assignment_stats <- data.frame()

### (1) Assignment done with original PR2-database and PR2 taxonomy
#ASVs and Taxonomic Assignment
results_1_mothur <- read.table(mothur_results_db1)

results_1_clean <- results_1_mothur %>%
  mutate(Kingdom = str_split_i(V2, ";", 1)) %>%
  mutate(King_perc = str_split_i(Kingdom, "\\(", 2)) %>%
  mutate(Kingdom = str_remove_all(Kingdom, "\\(.+\\)")) %>%
  mutate(King_perc = str_remove_all(King_perc, "\\)")) %>%
  
  mutate(Phylum = str_split_i(V2, ";", 2)) %>%
  mutate(Phyl_perc = str_split_i(Phylum, "\\(", 2)) %>%
  mutate(Phylum = str_remove_all(Phylum, "\\(.+\\)")) %>%
  mutate(Phyl_perc = str_remove_all(Phyl_perc, "\\)")) %>%
  
  mutate(Subphylum = str_split_i(V2, ";", 3)) %>%
  mutate(Subp_perc = str_split_i(Subphylum, "\\(", 2)) %>%
  mutate(Subphylum = str_remove_all(Subphylum, "\\(.+\\)")) %>%
  mutate(Subp_perc = str_remove_all(Subp_perc, "\\)")) %>%
  
  mutate(Class = str_split_i(V2, ";", 4)) %>%
  mutate(Clas_perc = str_split_i(Class, "\\(", 2)) %>%
  mutate(Class = str_remove_all(Class, "\\(.+\\)")) %>%
  mutate(Clas_perc = str_remove_all(Clas_perc, "\\)")) %>%
  
  mutate(Subclass = str_split_i(V2, ";", 5)) %>%
  mutate(Subc_perc = str_split_i(Subclass, "\\(", 2)) %>%
  mutate(Subclass = str_remove_all(Subclass, "\\(.+\\)")) %>%
  mutate(Subc_perc = str_remove_all(Subc_perc, "\\)")) %>%
  
  mutate(Order = str_split_i(V2, ";", 6)) %>%
  mutate(Orde_perc = str_split_i(Order, "\\(", 2)) %>%
  mutate(Order = str_remove_all(Order, "\\(.+\\)")) %>%
  mutate(Orde_perc = str_remove_all(Orde_perc, "\\)")) %>%
  
  mutate(Family = str_split_i(V2, ";", 7)) %>%
  mutate(Fami_perc = str_split_i(Family, "\\(", 2)) %>%
  mutate(Family = str_remove_all(Family, "\\(.+\\)")) %>%
  mutate(Fami_perc = str_remove_all(Fami_perc, "\\)")) %>%
  
  mutate(Genus = str_split_i(V2, ";", 8)) %>%
  mutate(Genu_perc = str_split_i(Genus, "\\(", 2)) %>%
  mutate(Genus = str_remove_all(Genus, "\\(.+\\)")) %>%
  mutate(Genu_perc = str_remove_all(Genu_perc, "\\)")) %>%
  
  mutate(Species = str_split_i(V2, ";", 9)) %>%
  mutate(Spec_perc = str_split_i(Species, "[^_]\\(", 2)) %>%
  mutate(Species = str_remove_all(Species, "\\(\\d+\\)$")) %>%
  mutate(Spec_perc = str_remove_all(Spec_perc, "\\)")) %>%
  
  select(-V2)

#Classified Ciliates (all, species level, genus level)
classified_1_Cil <- results_1_clean %>%
  filter(Class == "Ciliophora") %>%
  select(V1, Species, Spec_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_1_Cil_spec <- results_1_clean %>%
  filter(Class == "Ciliophora") %>%
  filter(str_detect(Species, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Species, "sp.", negate = TRUE)) %>%
  filter(str_detect(Species, "sp", negate = TRUE)) %>%
  filter(str_detect(Species, "_X", negate = TRUE)) %>%
  select(V1, Species, Spec_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_1_Cil_gen <- results_1_clean %>%
  filter(Class == "Ciliophora") %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE)) %>%
  select(V1, Genus, Genu_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_1_Cil_fam <- results_1_clean %>%
  filter(Class == "Ciliophora") %>%
  filter(str_detect(Family, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Family, "_X", negate = TRUE)) %>%
  select(V1, Family, Fami_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

#Classified (all ASVs classified to species/genus level)
classified_1_spec <- results_1_clean %>%
  filter(str_detect(Species, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Species, "sp.", negate = TRUE)) %>%
  filter(str_detect(Species, "sp", negate = TRUE)) %>%
  filter(str_detect(Species, "_X", negate = TRUE)) %>%
  select(V1, Species, Spec_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)
  
classified_1_gen <- results_1_clean %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE)) %>%
  select(V1, Genus, Genu_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

stats <- data.frame(database="v1") %>% 
    mutate(genus= nrow(classified_1_gen),
           species= nrow(classified_1_spec),
           ciliate= nrow(classified_1_Cil),
           cil_fam = nrow(classified_1_Cil_fam),
           cil_gen_no = length(unique(classified_1_Cil_gen$Genus)),
           cil_gen = nrow(classified_1_Cil_gen),
           cil_spec = length(unique(classified_1_Cil_spec$Species)),
           cil_ASV = nrow(classified_1_Cil_spec),
           measurement= "ASV")%>%
    bind_rows( data.frame(database="v1", measurement="Reads") %>%
                 mutate(genus= sum(classified_1_gen$total_reads),
                        species= sum(classified_1_spec$total_reads),
                        ciliate= sum(classified_1_Cil$total_reads),
                        cil_fam = sum(classified_1_Cil_fam$total_reads),
                        cil_gen = sum(classified_1_Cil_gen$total_reads),
                        cil_spec = sum(classified_1_Cil_spec$total_reads)))
assignment_stats <- rbind(assignment_stats,stats)

### (2) Assignment with original PR2 taxonomy but only sequences where a WoRMS entry exist
#ASVs and Taxonomic Assignment
results_2_mothur <- read.table(mothur_results_db2)

results_2_clean <- results_2_mothur %>%
  mutate(Kingdom = str_split_i(V2, ";", 1)) %>%
  mutate(King_perc = str_split_i(Kingdom, "\\(", 2)) %>%
  mutate(Kingdom = str_remove_all(Kingdom, "\\(.+\\)")) %>%
  mutate(King_perc = str_remove_all(King_perc, "\\)")) %>%
  
  mutate(Phylum = str_split_i(V2, ";", 2)) %>%
  mutate(Phyl_perc = str_split_i(Phylum, "\\(", 2)) %>%
  mutate(Phylum = str_remove_all(Phylum, "\\(.+\\)")) %>%
  mutate(Phyl_perc = str_remove_all(Phyl_perc, "\\)")) %>%
  
  mutate(Subphylum = str_split_i(V2, ";", 3)) %>%
  mutate(Subp_perc = str_split_i(Subphylum, "\\(", 2)) %>%
  mutate(Subphylum = str_remove_all(Subphylum, "\\(.+\\)")) %>%
  mutate(Subp_perc = str_remove_all(Subp_perc, "\\)")) %>%
  
  mutate(Class = str_split_i(V2, ";", 4)) %>%
  mutate(Clas_perc = str_split_i(Class, "\\(", 2)) %>%
  mutate(Class = str_remove_all(Class, "\\(.+\\)")) %>%
  mutate(Clas_perc = str_remove_all(Clas_perc, "\\)")) %>%
  
  mutate(Subclass = str_split_i(V2, ";", 5)) %>%
  mutate(Subc_perc = str_split_i(Subclass, "\\(", 2)) %>%
  mutate(Subclass = str_remove_all(Subclass, "\\(.+\\)")) %>%
  mutate(Subc_perc = str_remove_all(Subc_perc, "\\)")) %>%
  
  mutate(Order = str_split_i(V2, ";", 6)) %>%
  mutate(Orde_perc = str_split_i(Order, "\\(", 2)) %>%
  mutate(Order = str_remove_all(Order, "\\(.+\\)")) %>%
  mutate(Orde_perc = str_remove_all(Orde_perc, "\\)")) %>%
  
  mutate(Family = str_split_i(V2, ";", 7)) %>%
  mutate(Fami_perc = str_split_i(Family, "\\(", 2)) %>%
  mutate(Family = str_remove_all(Family, "\\(.+\\)")) %>%
  mutate(Fami_perc = str_remove_all(Fami_perc, "\\)")) %>%
  
  mutate(Genus = str_split_i(V2, ";", 8)) %>%
  mutate(Genu_perc = str_split_i(Genus, "\\(", 2)) %>%
  mutate(Genus = str_remove_all(Genus, "\\(.+\\)")) %>%
  mutate(Genu_perc = str_remove_all(Genu_perc, "\\)")) %>%
  
  mutate(Species = str_split_i(V2, ";", 9)) %>%
  mutate(Spec_perc = str_split_i(Species, "[^_]\\(", 2)) %>%
  mutate(Species = str_remove_all(Species, "\\(\\d+\\)$")) %>%
  mutate(Spec_perc = str_remove_all(Spec_perc, "\\)")) %>%
  
  select(-V2)

#Classified Ciliates (all, species level, genus level)
classified_2_Cil <- results_2_clean %>%
  filter(Class == "Ciliophora") %>%
  select(V1, Species, Spec_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_2_Cil_spec <- results_2_clean %>%
  filter(Class == "Ciliophora") %>%
  filter(str_detect(Species, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Species, "sp.", negate = TRUE)) %>%
  filter(str_detect(Species, "sp", negate = TRUE)) %>%
  filter(str_detect(Species, "_X", negate = TRUE)) %>%
  select(V1, Species, Spec_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_2_Cil_gen <- results_2_clean %>%
  filter(Class == "Ciliophora") %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE)) %>%
  select(V1, Genus, Genu_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_2_Cil_fam <- results_2_clean %>%
  filter(Class == "Ciliophora") %>%
  filter(str_detect(Family, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Family, "_X", negate = TRUE)) %>%
  select(V1, Family, Fami_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)



#Classified (all ASVs classified to species/genus level)
classified_2_spec <- results_2_clean %>%
  filter(str_detect(Species, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Species, "sp.", negate = TRUE)) %>%
  filter(str_detect(Species, "sp", negate = TRUE)) %>%
  filter(str_detect(Species, "_X", negate = TRUE)) %>%
  select(V1, Species, Spec_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_2_gen <- results_2_clean %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE)) %>%
  select(V1, Genus, Genu_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)


stats <- data.frame(database="v2") %>% 
  mutate(genus= nrow(classified_2_gen),
         species= nrow(classified_2_spec),
         ciliate= nrow(classified_2_Cil),
         cil_fam = nrow(classified_2_Cil_fam),
         cil_gen_no = length(unique(classified_2_Cil_gen$Genus)),
         cil_gen = nrow(classified_2_Cil_gen),
         cil_spec = length(unique(classified_2_Cil_spec$Species)),
         cil_ASV = nrow(classified_2_Cil_spec),
         measurement= "ASV")%>%
  bind_rows( data.frame(database="v2", measurement="Reads") %>%
               mutate(genus= sum(classified_2_gen$total_reads),
                      species= sum(classified_2_spec$total_reads),
                      ciliate= sum(classified_2_Cil$total_reads),
                      cil_fam = sum(classified_2_Cil_fam$total_reads),
                      cil_gen = sum(classified_2_Cil_gen$total_reads),
                      cil_spec = sum(classified_2_Cil_spec$total_reads)))
assignment_stats <- rbind(assignment_stats,stats)


### (3) Assignment with cleaned PR2 database and WoRMS Taxonomy
#ASVs and Taxonomic Assignment
results_3_mothur <- read.table(mothur_results_db3)

results_3_clean <- results_3_mothur %>%
  mutate(Kingdom = str_split_i(V2, ";", 1)) %>%
  mutate(King_perc = str_split_i(Kingdom, "\\(", 2)) %>%
  mutate(Kingdom = str_remove_all(Kingdom, "\\(.+\\)")) %>%
  mutate(King_perc = str_remove_all(King_perc, "\\)")) %>%
  
  mutate(Phylum = str_split_i(V2, ";", 2)) %>%
  mutate(Phyl_perc = str_split_i(Phylum, "\\(", 2)) %>%
  mutate(Phylum = str_remove_all(Phylum, "\\(.+\\)")) %>%
  mutate(Phyl_perc = str_remove_all(Phyl_perc, "\\)")) %>%
  
  mutate(Subphylum = str_split_i(V2, ";", 3)) %>%
  mutate(Subp_perc = str_split_i(Subphylum, "\\(", 2)) %>%
  mutate(Subphylum = str_remove_all(Subphylum, "\\(.+\\)")) %>%
  mutate(Subp_perc = str_remove_all(Subp_perc, "\\)")) %>%
  
  mutate(Class = str_split_i(V2, ";", 4)) %>%
  mutate(Clas_perc = str_split_i(Class, "\\(", 2)) %>%
  mutate(Class = str_remove_all(Class, "\\(.+\\)")) %>%
  mutate(Clas_perc = str_remove_all(Clas_perc, "\\)")) %>%
  
  mutate(Subclass = str_split_i(V2, ";", 5)) %>%
  mutate(Subc_perc = str_split_i(Subclass, "\\(", 2)) %>%
  mutate(Subclass = str_remove_all(Subclass, "\\(.+\\)")) %>%
  mutate(Subc_perc = str_remove_all(Subc_perc, "\\)")) %>%
  
  mutate(Order = str_split_i(V2, ";", 6)) %>%
  mutate(Orde_perc = str_split_i(Order, "\\(", 2)) %>%
  mutate(Order = str_remove_all(Order, "\\(.+\\)")) %>%
  mutate(Orde_perc = str_remove_all(Orde_perc, "\\)")) %>%
  
  mutate(Family = str_split_i(V2, ";", 7)) %>%
  mutate(Fami_perc = str_split_i(Family, "\\(", 2)) %>%
  mutate(Family = str_remove_all(Family, "\\(.+\\)")) %>%
  mutate(Fami_perc = str_remove_all(Fami_perc, "\\)")) %>%
  
  mutate(Genus = str_split_i(V2, ";", 8)) %>%
  mutate(Genu_perc = str_split_i(Genus, "\\(", 2)) %>%
  mutate(Genus = str_remove_all(Genus, "\\(.+\\)")) %>%
  mutate(Genu_perc = str_remove_all(Genu_perc, "\\)")) %>%
  
  mutate(Species = str_split_i(V2, ";", 9)) %>%
  mutate(Spec_perc = str_split_i(Species, "[^_]\\(", 2)) %>%
  mutate(Species = str_remove_all(Species, "\\(\\d+\\)$")) %>%
  mutate(Spec_perc = str_remove_all(Spec_perc, "\\)")) %>%
  
  select(-V2)

#Classified Ciliates (all, species level, genus level)
classified_3_Cil <- results_3_clean %>%
  filter(Phylum == "Ciliophora") %>%
  select(V1, Species, Spec_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_3_Cil_spec <- results_3_clean %>%
  filter(Phylum == "Ciliophora") %>%
  filter(str_detect(Species, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Species, "sp.", negate = TRUE)) %>%
  filter(str_detect(Species, "sp", negate = TRUE)) %>%
  filter(str_detect(Species, "_X", negate = TRUE)) %>%
  select(V1, Species, Spec_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_3_Cil_gen <- results_3_clean %>%
  filter(Phylum == "Ciliophora") %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE)) %>%
  select(V1, Genus, Genu_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_3_Cil_fam <- results_3_clean %>%
  filter(Phylum == "Ciliophora") %>%
  filter(str_detect(Family, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Family, "_X", negate = TRUE)) %>%
  select(V1, Family, Fami_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)


#Classified (all ASVs classified to species/genus level)
classified_3_spec <- results_3_clean %>%
  filter(str_detect(Species, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Species, "sp.", negate = TRUE)) %>%
  filter(str_detect(Species, "sp", negate = TRUE)) %>%
  filter(str_detect(Species, "_X", negate = TRUE)) %>%
  select(V1, Species, Spec_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_3_gen <- results_3_clean %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE)) %>%
  select(V1, Genus, Genu_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

stats <- data.frame(database="v3") %>% 
  mutate(genus= nrow(classified_3_gen),
         species= nrow(classified_3_spec),
         ciliate= nrow(classified_3_Cil),
         cil_fam = nrow(classified_3_Cil_fam),
         cil_gen_no = length(unique(classified_3_Cil_gen$Genus)),
         cil_gen = nrow(classified_3_Cil_gen),
         cil_spec = length(unique(classified_3_Cil_spec$Species)),
         cil_ASV = nrow(classified_3_Cil_spec),
         measurement= "ASV")%>%
  bind_rows( data.frame(database="v3", measurement="Reads") %>%
               mutate(genus= sum(classified_3_gen$total_reads),
                      species= sum(classified_3_spec$total_reads),
                      ciliate= sum(classified_3_Cil$total_reads),
                      cil_fam = sum(classified_3_Cil_fam$total_reads),
                      cil_gen = sum(classified_3_Cil_gen$total_reads),
                      cil_spec = sum(classified_3_Cil_spec$total_reads)))
assignment_stats <- rbind(assignment_stats,stats)

### (4) Assignment done with updated PR2 database and WoRMS taxonomy
#ASVs and Taxonomic Assignment
results_4_mothur <- read.table(mothur_results_db4)

results_4_clean <- results_4_mothur %>%
  mutate(Kingdom = str_split_i(V2, ";", 1)) %>%
  mutate(King_perc = str_split_i(Kingdom, "\\(", 2)) %>%
  mutate(Kingdom = str_remove_all(Kingdom, "\\(.+\\)")) %>%
  mutate(King_perc = str_remove_all(King_perc, "\\)")) %>%
  
  mutate(Phylum = str_split_i(V2, ";", 2)) %>%
  mutate(Phyl_perc = str_split_i(Phylum, "\\(", 2)) %>%
  mutate(Phylum = str_remove_all(Phylum, "\\(.+\\)")) %>%
  mutate(Phyl_perc = str_remove_all(Phyl_perc, "\\)")) %>%
  
  mutate(Subphylum = str_split_i(V2, ";", 3)) %>%
  mutate(Subp_perc = str_split_i(Subphylum, "\\(", 2)) %>%
  mutate(Subphylum = str_remove_all(Subphylum, "\\(.+\\)")) %>%
  mutate(Subp_perc = str_remove_all(Subp_perc, "\\)")) %>%
  
  mutate(Class = str_split_i(V2, ";", 4)) %>%
  mutate(Clas_perc = str_split_i(Class, "\\(", 2)) %>%
  mutate(Class = str_remove_all(Class, "\\(.+\\)")) %>%
  mutate(Clas_perc = str_remove_all(Clas_perc, "\\)")) %>%
  
  mutate(Subclass = str_split_i(V2, ";", 5)) %>%
  mutate(Subc_perc = str_split_i(Subclass, "\\(", 2)) %>%
  mutate(Subclass = str_remove_all(Subclass, "\\(.+\\)")) %>%
  mutate(Subc_perc = str_remove_all(Subc_perc, "\\)")) %>%
  
  mutate(Order = str_split_i(V2, ";", 6)) %>%
  mutate(Orde_perc = str_split_i(Order, "\\(", 2)) %>%
  mutate(Order = str_remove_all(Order, "\\(.+\\)")) %>%
  mutate(Orde_perc = str_remove_all(Orde_perc, "\\)")) %>%
  
  mutate(Family = str_split_i(V2, ";", 7)) %>%
  mutate(Fami_perc = str_split_i(Family, "\\(", 2)) %>%
  mutate(Family = str_remove_all(Family, "\\(.+\\)")) %>%
  mutate(Fami_perc = str_remove_all(Fami_perc, "\\)")) %>%
  
  mutate(Genus = str_split_i(V2, ";", 8)) %>%
  mutate(Genu_perc = str_split_i(Genus, "\\(", 2)) %>%
  mutate(Genus = str_remove_all(Genus, "\\(.+\\)")) %>%
  mutate(Genu_perc = str_remove_all(Genu_perc, "\\)")) %>%
  
  mutate(Species = str_split_i(V2, ";", 9)) %>%
  mutate(Spec_perc = str_split_i(Species, "[^_]\\(", 2)) %>%
  mutate(Species = str_remove_all(Species, "\\(\\d+\\)$")) %>%
  mutate(Spec_perc = str_remove_all(Spec_perc, "\\)")) %>%
  
  select(-V2)

#Classified Ciliates (all, species level, genus level)
classified_4_Cil <- results_4_clean %>% 
  filter(Phylum == "Ciliophora") %>%
  select(V1, Species, Spec_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_4_Cil_spec <- results_4_clean %>%
  filter(Phylum == "Ciliophora") %>%
  filter(str_detect(Species, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Species, "sp.", negate = TRUE)) %>%
  filter(str_detect(Species, "sp", negate = TRUE)) %>%
  filter(str_detect(Species, "_X", negate = TRUE)) %>%
  select(V1, Species, Spec_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_4_Cil_gen <- results_4_clean %>%
  filter(Phylum == "Ciliophora") %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE)) %>%
  select(V1, Genus, Genu_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_4_Cil_fam <- results_4_clean %>%
  filter(Phylum == "Ciliophora") %>%
  filter(str_detect(Family, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Family, "_X", negate = TRUE)) %>%
  select(V1, Family, Fami_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

#Classified (all ASVs classified to species/genus level)
classified_4_spec <- results_4_clean %>%
  filter(str_detect(Species, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Species, "sp.", negate = TRUE)) %>%
  filter(str_detect(Species, "sp", negate = TRUE)) %>%
  filter(str_detect(Species, "_X", negate = TRUE)) %>%
  select(V1, Species, Spec_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)

classified_4_gen <- results_4_clean %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE)) %>%
  select(V1, Genus, Genu_perc) %>%
  inner_join(counttable, by = "V1") %>%
  filter(total_reads > 0)


stats <- data.frame(database="v4") %>% 
  mutate(genus= nrow(classified_4_gen),
         species= nrow(classified_4_spec),
         ciliate= nrow(classified_4_Cil),
         cil_fam = nrow(classified_4_Cil_fam),
         cil_gen_no = length(unique(classified_4_Cil_gen$Genus)),
         cil_gen = nrow(classified_4_Cil_gen),
         cil_spec = length(unique(classified_4_Cil_spec$Species)),
         cil_ASV = nrow(classified_4_Cil_spec),
         measurement= "ASV")%>%
  bind_rows( data.frame(database="v4", measurement="Reads") %>%
               mutate(genus= sum(classified_4_gen$total_reads),
                      species= sum(classified_4_spec$total_reads),
                      ciliate= sum(classified_4_Cil$total_reads),
                      cil_fam = sum(classified_4_Cil_fam$total_reads),
                      cil_gen = sum(classified_4_Cil_gen$total_reads),
                      cil_spec = sum(classified_4_Cil_spec$total_reads)))
assignment_stats <- rbind(assignment_stats,stats)


##compare families
stats_family1 <- classified_1_Cil_fam %>% full_join(classified_2_Cil_fam, by = "V1", suffix = c(".db1",".db2"))
stats_family2 <- classified_3_Cil_fam %>% full_join(classified_4_Cil_fam, by = "V1", suffix = c(".db3",".db4"))
stats_family <- stats_family1 %>% full_join(stats_family2, by = "V1")
write.table(stats_family, file=file.path(output_path,"13_Ciliate_family_list_RefDBcomparison.tsv"), sep="\t", row.names = F)
rm(stats_family1, stats_family2)

#investigating
new_family_ASV_db2 <-stats_family %>% filter( is.na(Family.db1)& !is.na(Family.db2))%>% pull(V1)
results_1_clean %>% filter(V1 %in%new_family_ASV_db2 )%>% select(Class, Subclass, Order, Family)%>%unique

new_family_ASV_db3 <-stats_family %>% filter( is.na(Family.db2)& !is.na(Family.db3))%>% pull(V1)
results_2_clean %>% filter(V1 %in%new_family_ASV_db3 )%>% select(Class, Subclass, Order, Family)%>%unique()


### baltic genera
baltic_4_Cil_gen <- classified_4_Cil_gen %>% filter(Genus %in% baltic_ciliates$genus)
balticNon_4_Cil_gen <- classified_4_Cil_gen %>% filter(!Genus %in% baltic_ciliates$genus)
nrow(classified_4_Cil_gen)
nrow(baltic_4_Cil_gen)


baltic_4_Cil_spec <- classified_4_Cil_spec %>% mutate(Species =gsub("_", " ", Species)) %>% filter(Species %in% baltic_ciliates$species)
nrow(classified_4_Cil_spec)
nrow(baltic_4_Cil_spec)
################################################################################
################################################################################

### Comparisons

#Overview over all databases (species level) for all sequences assigned to 
#ciliates
overview_1 <- as.data.frame.table(table(classified_1_Cil$Species))
overview_2 <- as.data.frame.table(table(classified_2_Cil$Species))
overview_3 <- as.data.frame.table(table(classified_3_Cil$Species))
overview_4 <- as.data.frame.table(table(classified_4_Cil$Species))

#comparison of numbers of ASVs for each assignment
comparison_num_seq <- full_join(overview_1, overview_2, by = "Var1")
comparison_tax <- full_join(overview_2, overview_3, by = "Var1") 
comparison_seq <- full_join(overview_3, overview_4, by = "Var1")

#Comparison of assignment for each ASV (all, species level, genus level)
comparison_ASVs <- full_join(classified_1_Cil, classified_2_Cil, by = "V1", suffix = c(".v1", ".v2")) %>%
  full_join(full_join(classified_3_Cil, classified_4_Cil, by = "V1", suffix = c(".v3", ".v4")), by = "V1") 

comparison_ASV_Spec <- full_join(classified_1_Cil_spec, classified_2_Cil_spec, by = "V1", suffix = c(".v1", ".v2")) %>%
  full_join(full_join(classified_3_Cil_spec, classified_4_Cil_spec, by = "V1", suffix = c(".v3", ".v4")), by = "V1") 

comparison_ASV_Gen <- full_join(classified_1_Cil_gen, classified_2_Cil_gen, by = "V1", suffix = c(".v1", ".v2")) %>%
  full_join(full_join(classified_3_Cil_gen, classified_4_Cil_gen, by = "V1", suffix = c(".v3", ".v4")), by = "V1") 


############ Compare  the different genus assignments
comparison_genus <-classified_1_Cil_gen %>% mutate(database="v1")%>%
                      bind_rows(classified_2_Cil_gen %>% mutate(database="v2"))%>%
                      bind_rows(classified_3_Cil_gen %>% mutate(database="v3"))%>%
                      bind_rows(classified_4_Cil_gen %>% mutate(database="v4"))%>%
                      group_by(database, Genus)%>%
                      summarise(ASV= n(),
                                reads=sum(total_reads))%>%
                      ungroup()%>%
                      arrange(Genus, database)%>%
                      mutate(Baltic = case_when(Genus %in% baltic_ciliates$genus ~ "TRUE",
                                                TRUE ~ "FALSE"))

# only keep those genera in which the read number changed
comparison_genus_changes <- comparison_genus%>%
  group_by(Genus) %>%
  mutate(v1_reads = ifelse(database == "v1", reads, NA),
         v4_reads = ifelse(database == "v4", reads, NA)) %>%
  summarize(v1_reads = first(na.omit(v1_reads)),
            v4_reads = first(na.omit(v4_reads)),
            .groups = "drop") %>%
  filter(v1_reads != v4_reads | is.na(v1_reads) | is.na(v4_reads)) %>%
  select(Genus) %>%
  inner_join(comparison_genus, by = "Genus")

comparison_genus_reads <- comparison_genus %>% #filter(database %in% c("v2", "v1"))%>%
  select(-ASV)%>%
  pivot_wider(names_from = database, values_from = reads, names_prefix = "reads_")%>%
  mutate(
    reads_v1 = ifelse(is.na(reads_v1), 0, reads_v1),
    reads_v2 = ifelse(is.na(reads_v2), 0, reads_v2),
    reads_v3 = ifelse(is.na(reads_v3), 0, reads_v3),
    reads_v4 = ifelse(is.na(reads_v4), 0, reads_v4)
  ) %>%
  mutate(reads_diff_v1.2 = reads_v2 - reads_v1,
         reads_diff_v1.3 = reads_v3 - reads_v1,
         reads_diff_v2.3 = reads_v3 - reads_v2,
         reads_diff_v1.4 = reads_v4 - reads_v1,
         reads_diff_v3.4 = reads_v4 - reads_v3)


comparison_genus_ASV <- comparison_genus %>% #filter(database %in% c("v2", "v1"))%>%
  select(-reads)%>%
  pivot_wider(names_from = database, values_from = ASV, names_prefix = "ASV_")%>%
  mutate(
    ASV_v1 = ifelse(is.na(ASV_v1), 0, ASV_v1),
    ASV_v2 = ifelse(is.na(ASV_v2), 0, ASV_v2),
    ASV_v3 = ifelse(is.na(ASV_v3), 0, ASV_v3),
    ASV_v4 = ifelse(is.na(ASV_v4), 0, ASV_v4)
  ) %>%
  mutate(ASV_diff_v1.2 = ASV_v2 - ASV_v1,
         ASV_diff_v1.3 = ASV_v3 - ASV_v1,
         ASV_diff_v2.3 = ASV_v3 - ASV_v2,
         ASV_diff_v1.4 = ASV_v4 - ASV_v1,
         ASV_diff_v3.4 = ASV_v4 - ASV_v3)

comparison_genus_ASV_reads <- comparison_genus_ASV %>% full_join(comparison_genus_reads, by =c("Baltic", "Genus"))

write.table(comparison_genus_ASV_reads, file=file.path(output_path,"13_Ciliate_RefDBcomparison__genus_readsNasvs.tsv"), sep="\t", row.names = F)

rm(comparison_genus_ASV,comparison_genus_reads)


#Number classified as ciliates (ASVs)
number_classified <- c(sum(!is.na(comparison_ASVs$Species.x)), sum(!is.na(comparison_ASVs$Species.y)), sum(!is.na(comparison_ASVs$Species.x.x)), sum(!is.na(comparison_ASVs$Spec_perc.y.y)))

#Comparison of ASVs that were assigned identical with all 4 versions
comparison_ASVs_id <- comparison_ASVs %>% 
  filter(Species.v1 == Species.v2 & Species.v2 == Species.v3 & Species.v3 == Species.v4) %>% 
  mutate(Spec_perc.v1 = as.numeric(Spec_perc.v1), Spec_perc.v2 = as.numeric(Spec_perc.v2), 
         Spec_perc.v3 = as.numeric(Spec_perc.v3), Spec_perc.v4 = as.numeric(Spec_perc.v4))

comparison_ASV_Spec_id <- comparison_ASV_Spec %>% 
  filter(Species.v1 == Species.v2 & Species.v2 == Species.v3 & Species.v3 == Species.v4) %>% 
  mutate(Spec_perc.v1 = as.numeric(Spec_perc.v1), Spec_perc.v2 = as.numeric(Spec_perc.v2), 
         Spec_perc.v3 = as.numeric(Spec_perc.v3), Spec_perc.v4 = as.numeric(Spec_perc.v4))

comparison_ASV_Gen_id <- comparison_ASV_Gen %>% 
  filter(Genus.v1 == Genus.v2 & Genus.v2 == Genus.v3 & Genus.v3 == Genus.v4) %>%
  mutate(Genu_perc.v1 = as.numeric(Genu_perc.v1), Genu_perc.v2 = as.numeric(Genu_perc.v2), 
         Genu_perc.v3 = as.numeric(Genu_perc.v3), Genu_perc.v4 = as.numeric(Genu_perc.v4))

#Changes in bootstrap values for identical assigned ASVs
boxplot(comparison_ASV_Gen_id$Genu_perc.v2 - comparison_ASV_Gen_id$Genu_perc.v1, 
        comparison_ASV_Gen_id$Genu_perc.v3 - comparison_ASV_Gen_id$Genu_perc.v2, 
        comparison_ASV_Gen_id$Genu_perc.v4 - comparison_ASV_Gen_id$Genu_perc.v3, 
        names = c("Remove Pr2 sequences", "Change taxonomy", "Add NCBI sequences"), main = "Changes in bootstrap values")

data.frame(v1_v2 = comparison_ASV_Gen_id$Genu_perc.v2 - comparison_ASV_Gen_id$Genu_perc.v1,
           v2_v3 = comparison_ASV_Gen_id$Genu_perc.v3 - comparison_ASV_Gen_id$Genu_perc.v2,
           v3_v4 = comparison_ASV_Gen_id$Genu_perc.v4 - comparison_ASV_Gen_id$Genu_perc.v3)%>%
  pivot_longer(everything(),names_to = "comparison", values_to = "difference")%>%
  mutate(Comparison= case_when(
    comparison=="v1_v2" ~ "Remove PR2 sequences",
    comparison=="v2_v3" ~ "Change taxonomy",
    comparison=="v3_v4" ~ "Add NCBI sequences"  ))%>%
  mutate(Comparison= factor(Comparison, levels=c("Remove PR2 sequences", "Change taxonomy", "Add NCBI sequences")))%>%
  ggplot(aes(x=Comparison, y= difference))+
    geom_boxplot(fill="#7AC5CD")+theme_light()+
    labs(title="Changes in bootstrap values on the genus level", y="Difference", x="Changes in the reference database")+
    theme(axis.title.x = element_text(vjust = -2), axis.title.y = element_text(vjust = 2), title = element_text(vjust = 2)) 
ggsave(file.path(output_path, "13_Figure_bootstrap_changes.jpg"), dpi = 300, width=8, height = 6)


################################################################################
################################################################################

# print supplementary tables

write.table(assignment_stats, file=file.path(output_path,"13_SUMMARY_Read_ASV_statistics.tsv"), sep="\t", row.names = F)
write.table(comparison_ASVs_id, file=file.path(output_path,"13_Ciliate_Analysis_RefDB_comparison__assigned.tsv"), sep="\t", row.names = F)
write.table(comparison_ASV_Spec_id, file=file.path(output_path,"13_Ciliate_Analysis_RefDB_comparison__species.tsv"), sep="\t", row.names = F)
write.table(comparison_ASV_Gen_id, file=file.path(output_path,"13_Ciliate_Analysis_RefDB_comparison__genus.tsv"), sep="\t", row.names = F)

################################################################################
################################################################################


#"real" species level assignments and additonal! genus level assignments (read 
#that are assigned on species level are not included on genus level)
spec_4 <- classified_4_Cil_spec %>% filter(str_detect(Species, "_sp.", negate = TRUE))
gen_4 <- classified_4_Cil_gen %>% filter(!V1 %in% spec_4$V1)

spec_3 <- classified_3_Cil_spec %>% filter(str_detect(Species, "_sp.", negate = TRUE))
gen_3 <- classified_3_Cil_gen %>% filter(!V1 %in% spec_3$V1)

spec_2 <- classified_2_Cil_spec %>% filter(str_detect(Species, "_sp.", negate = TRUE))
gen_2 <- classified_2_Cil_gen %>% filter(!V1 %in% spec_2$V1)

spec_1 <- classified_1_Cil_spec %>% filter(str_detect(Species, "_sp.", negate = TRUE))
gen_1 <- classified_1_Cil_gen %>% filter(!V1 %in% spec_1$V1)


############# Plot results
## Number of classified ASVs
tax_unit <- rep(c("species", "genus"), 4)
database <- sort(rep(c("PR2\n(V1)", "PR2 (WoRMS exists)\n(V2)", "PR2 + WoRMS\n(V3)", "PR2 + WoRMS + NCBI Sequences\n(V4)"), 2))
num_classified <- c(nrow(spec_1), nrow(gen_1), nrow(spec_2), nrow(gen_2), nrow(spec_3), nrow(gen_3), nrow(spec_4), nrow(gen_4))

stats_classified <- data.frame(tax_unit, database, num_classified) 

ggplot(stats_classified, aes(fill = tax_unit, y = num_classified, x = database)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = "Database version", y = "Total ASV", fill = "Taxonomic level") +
  theme_light() +
  scale_fill_manual(values=c("darkseagreen3", "#00688B"), 
                    labels=c("genus"="Genus", "species"="Species"))+
  ggtitle("Number of classified ASVs")+
  theme(axis.title.x = element_text(vjust = -2), axis.title.y = element_text(vjust = 2), title = element_text(vjust = 2)) 
ggsave(file.path(output_path, "13_Figure_stats_classified__identifiedASVs.jpg"), dpi = 300, width=9, height = 5)

## Number Ciliates classified (Taxa)
num_classified_tax <- c(n_distinct(spec_1$Species), n_distinct(classified_1_Cil_gen$Genus), n_distinct(spec_2$Species), n_distinct(classified_2_Cil_gen$Genus), n_distinct(spec_3$Species), n_distinct(classified_3_Cil_gen$Genus), n_distinct(spec_4$Species), n_distinct(classified_4_Cil_gen$Genus))

stats_classified_tax <- data.frame(tax_unit, database, num_classified_tax)

ggplot(stats_classified, aes(fill = tax_unit, y = num_classified_tax, x = database)) +
  geom_bar(position = "dodge", stat = "identity") +
  #ggtitle("Number of species and genera identified") +
  labs(x = "Database version", y = "Number of identified taxa", fill = "Taxonomic level") +
  theme_light() +
  scale_fill_manual(values=c("darkseagreen3", "#00688B"), 
                    labels=c("genus"="Genus", "species"="Species"))+
  theme(axis.title.x = element_text(vjust = -2), axis.title.y = element_text(vjust = 2), title = element_text(vjust = 2)) 

ggsave(file.path(output_path, "13_Figure_stats_classified__identifiedTaxa.jpg"), dpi = 300, width=9, height = 5)
  
## Number Ciliates classified (Reads)
num_classified_reads <- c(sum(spec_1$total_reads), sum(gen_1$total_reads), sum(spec_2$total_reads), sum(gen_2$total_reads), sum(spec_3$total_reads), sum(gen_3$total_reads), sum(spec_4$total_reads), sum(gen_4$total_reads))
stats_classified_reads <- data.frame(tax_unit, database, num_classified_reads)

ggplot(stats_classified, aes(fill = tax_unit, y = num_classified_reads, x = database)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = "Database version", y = "Total read number", fill = "Taxonomic level") +
  theme_light() +
  scale_fill_manual(values=c("darkseagreen3", "#00688B"), 
                    labels=c("genus"="Genus", "species"="Species"))+
  ggtitle("Number of reads identified")+scale_y_continuous(labels = scales::comma) +
  theme(axis.title.x = element_text(vjust = -2), axis.title.y = element_text(vjust = 2), title = element_text(vjust = 2)) 
ggsave(file.path(output_path, "13_Figure_stats_classified__identifiedReads.jpg"), dpi = 300, width=9, height = 5)

################################################################################
################################################################################

#Save cleaned results for Treemap script
saveRDS(results_1_clean, file.path(output_path,"13__Results_DB1"))
saveRDS(results_2_clean, file.path(output_path,"13__Results_DB2"))
saveRDS(results_3_clean, file.path(output_path,"13__Results_DB3"))
saveRDS(results_4_clean, file.path(output_path,"13__Results_DB4"))
