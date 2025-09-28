################################################################################
#This script can be used to calculate the summary and the 
# comparison of the different databases
################################################################################
require(tidyr)
require(dplyr)
require(stringr)



#setwd("/PATH/TO/02_Scripts_folder") # to path above the Script folder
setwd("/PATH/TO/11_Benchmark_Scripts")
source("00_Function_Library.R") # read file with functions
read_simple_ini("00_login_data.ini")

# path of the databases
path <- "11_Benchmark_Database_Versions"
crabs_database4 <- "11_Euka02_Taxonomy_DB4.tax"

##
overview <- read.table("10_FINAL_results/9.6_Overview_Sequences_FINAL.csv", header = TRUE)

#original pr2 database - taxonomy file for mothur
path_pr2 <- "../"
pr2_tax<- read.table(file.path(path_pr2,"pr2_version_5.1.0_SSU_mothur.tax"))


#output
output_path <- "13_Analyses_2025"


## list with Baltic Species
baltic_list <- file.path(input_path, input_table)
species_list <- read.table(baltic_list, header = TRUE, sep = ",")

#######################################################################################

###########
#Remove sequences that do not have sufficient taxonomic information
pr2_tax <- pr2_tax%>%
  mutate(Kingdom = str_split_i(V2, ";", 1)) %>%
  mutate(Phylum = str_split_i(V2, ";", 2)) %>%
  mutate(Subphylum = str_split_i(V2, ";", 3)) %>%
  mutate(Class = str_split_i(V2, ";", 4)) %>%
  mutate(Subclass = str_split_i(V2, ";", 5)) %>%
  mutate(Order = str_split_i(V2, ";", 6)) %>%
  mutate(Family = str_split_i(V2, ";", 7)) %>%
  mutate(Genus = str_split_i(V2, ";", 8)) %>%
  mutate(Species = str_split_i(V2, ";", 9)) %>%
  select(-V2)

pr2_badTax <- pr2_tax %>%
  #Remove sequences without taxonomic information at least on genus level
  filter(str_detect(Species, "_X")|  str_detect(Species, "_group_")|
           str_detect(Species, "^[[:upper:]][[:upper:]]")) 
#deleted bad taxonomy
nrow(pr2_badTax)
  
pr2_badTax_Cil <-pr2_badTax %>% filter(Class=="Ciliophora")

#deleted bad taxonomy Ciliophora
nrow(pr2_badTax_Cil)


pr2_badTax_Cil_sum <- pr2_badTax_Cil %>% group_by(Subclass, Order, Family,Species)%>%
      summarise(seq= n())%>% 
      group_by(Subclass, Order, Family)%>%
      summarise(Species= n(),
                seq= sum(seq))

write.table(pr2_badTax_Cil_sum, file = file.path(output_path, "15_Ciliate_stats_deletedSeq_noTax_familyLevel.tsv"), row.names = F, sep="\t")



#split database into origin
ncbi <- overview %>% filter(Database == "NCBI") %>%
  rowwise() %>%
  mutate(dummy = Acc_Name) %>%
  mutate(dummy = ifelse(is.na(Acc_Name), Clean_Name, Acc_Name))
  
pr2 <- overview %>% filter(Database == "pr2")

# sequences iwth accepted and clean name
overview_acc <- overview %>%
  filter(!is.na(Acc_Name)) %>%
  mutate(dummy = Acc_Name)
overview_cle <- overview %>%
  filter(is.na((Acc_Name))) %>%
  mutate(dummy = Clean_Name)
print(nrow(overview))
overview <- bind_rows(overview_acc, overview_cle)
print(nrow(overview))

overview_pr2 <- overview %>%
  filter(Database == "pr2")

print("New Species")
length(unique(overview$dummy)) - length(unique(overview_pr2$dummy)) 
print("New Sequences")
nrow(ncbi)

stats_df <- data.frame(database= "V4", variable= "PR2  Species No.", value= length(unique(overview_pr2$dummy)) )
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="PR2 Sequences No.", value= length(overview$dummy)))
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Species No.", value= length(unique(overview$dummy))))
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="NCBI Species No.", value= length(unique(overview$dummy)) - length(unique(overview_pr2$dummy)) ))
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="NCBI Sequences No.", value= nrow(ncbi)))


#sequences without taxonomy
overview_tax <- overview %>%
  filter(Taxonomy != 5)
length(unique(overview_tax$dummy))

stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Species No. - no taxonomy", value= length(unique(overview_tax$dummy))))
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Sequences No. - no taxonomy", value= length(overview_tax$dummy)))

# stats for the different taxonomy types
one <- overview %>%
  filter(Taxonomy == 1)
two <- overview %>%
  filter(Taxonomy == 2)
three <- overview %>%
  filter(Taxonomy == 3)
four <- overview %>%
  filter(Taxonomy == 4)
six <- overview %>%
  filter(Taxonomy == 6)

length(unique(overview$dummy))

stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Species No. - taxonomy=1", value= length(unique(one$dummy))))
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Sequences No. - taxonomy=1", value= length(one$dummy)))
#
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Species No. - taxonomy=2", value= length(unique(two$dummy))))
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Sequences No. - taxonomy=2", value= length(two$dummy)))
#
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Species No. - taxonomy=3", value= length(unique(three$dummy))))
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Sequences No. - taxonomy=3", value= length(three$dummy)))
#
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Species No. - taxonomy=4", value= length(unique(four$dummy))))
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Sequences No. - taxonomy=4", value= length(four$dummy)))
#
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Species No. - taxonomy=6", value= length(unique(six$dummy))))
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Sequences No. - taxonomy=6", value= length(six$dummy)))


tax <- read.table("10_FINAL_results/10.2_Taxonomy_FINAL.tax", sep = ";")
genus <- tax %>%
  filter(str_detect(V9, "sp."))
species <- tax %>%
  filter(str_detect(V9, "sp.", negate = TRUE))

length(unique(tax$V8)) # no unique genera
length(unique(tax$V7)) # no unique families

#
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Genus No.", value= length(unique(tax$V8)) ))
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Family No.", value= length(unique(tax$V7))))


#######################################################

## extract baltic genera from database taxonomy
baltic <- tax %>% filter(V8 %in% species_list$genus)

#####
# extract how many of the Baltic listed genera are in the final database
included <- species_list %>% filter(genus %in% tax$V8)

length(unique(species_list$genus))
stats_df <- rbind(stats_df, data.frame(database= "Species list", variable ="Genus No.", value= length(unique(species_list$genus))))


# no. of Baltic genera in DB
length(unique(baltic$V8))
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Baltic Genus No.", value= length(unique(baltic$V8))))

#
length(unique(included$genus))


#perecentage
length(unique(included$genus))/length(unique(species_list$genus)) * 100
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Represented Baltic Genus (%)", value= length(unique(included$genus))/length(unique(species_list$genus)) * 100))

#####
# extract how many of the Baltic listed genera are in the final database
ciliate_list <- species_list %>% filter(taxon == "Ciliophora")
included_cil <- ciliate_list %>% filter(genus %in% tax$V8)

length(unique(ciliate_list$genus))
stats_df <- rbind(stats_df, data.frame(database= "Species list", variable ="Ciliate Genus No.", value= length(unique(ciliate_list$genus))))

length(unique(included_cil$genus))
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Baltic Ciliate Genus No.", value= length(unique(included_cil$genus))))

length(unique(included_cil$genus))/length(unique(ciliate_list$genus)) * 100
stats_df <- rbind(stats_df, data.frame(database= "V4", variable ="Represented Baltic Ciliate Genus (%)", value= length(unique(included_cil$genus))/length(unique(ciliate_list$genus)) * 100))


##stats compared to orginal PR2 database
pr2_tax <- read.table("00_input/pr2_version_5.0.0_SSU_mothur.tax", sep="\t") %>%
  mutate(genus = str_split_i(V2, ";", 8))

included_pr2 <- species_list %>% filter(genus %in% pr2_tax$genus)
included_cil_pr2 <- ciliate_list %>% filter(genus %in% pr2_tax$genus)

# extract how many of the Baltic listed genera are in the original PR2 database
length(unique(species_list$genus))
length(unique(included_pr2$genus))
stats_df <- rbind(stats_df, data.frame(database= "V1", variable ="Baltic Genus No.", value= length(unique(included_pr2$genus))))

length(unique(included_pr2$genus))/length(unique(species_list$genus)) *100
stats_df <- rbind(stats_df, data.frame(database= "V1", variable ="Represented Baltic Genus (%)", value= length(unique(included_pr2$genus))/length(unique(species_list$genus)) *100))

# extract how many of the Baltic listed ciliate genera are in the original PR2 database
length((unique(included_cil_pr2$genus)))
length((unique(included_cil_pr2$genus)))/length(unique(ciliate_list$genus)) *100

stats_df <- rbind(stats_df, data.frame(database= "V1", variable ="Baltic Ciliate Genus No.", value= length((unique(included_cil_pr2$genus)))))
stats_df <- rbind(stats_df, data.frame(database= "V1", variable ="Represented Baltic Ciliate Genus (%)", value= length((unique(included_cil_pr2$genus)))/length(unique(ciliate_list$genus)) *100))
########################################
#extract stats about removed sequences
no_tax <- overview %>% 
  filter(Taxonomy == 5)

nrow(no_tax)


deleted_pr2_tax <- pr2_tax %>% filter(V1 %in% no_tax$pr2_accession)

deleted_pr2_tax <- deleted_pr2_tax%>%
  mutate(supergroup = str_split_i(V2, ";", 2))%>%
  mutate(division = str_split_i(V2, ";", 3))%>%
  mutate(subdivision = str_split_i(V2, ";", 4))%>%
  mutate(class = str_split_i(V2, ";", 5))%>%
  mutate(order = str_split_i(V2, ";", 6))%>%
  mutate(family = str_split_i(V2, ";", 7))%>%
  mutate(species = str_split_i(V2, ";", 9))%>%
  relocate(genus, .after= family)%>%
  select(-V2)

deleted_pr2_tax %>% group_by(supergroup, division)%>%
  summarise(no =n())%>%
  arrange(desc(no))

deleted_pr2_tax %>% group_by(supergroup, division,subdivision)%>%
  summarise(no =n())%>%
  arrange(desc(no))

deleted_pr2_tax %>% group_by(supergroup, division,subdivision,class)%>%
  summarise(no =n())%>%
  arrange(desc(no))

deleted_pr2_tax %>% group_by(supergroup, division,subdivision,class, order)%>%
  summarise(no =n())%>%
  arrange(desc(no))

deleted_pr2_tax %>% group_by(supergroup, division,subdivision,class, order,family)%>%
  summarise(no =n())%>%
  arrange(desc(no))


## save

noTax_stats_subdvision <- deleted_pr2_tax %>% group_by(supergroup, division,subdivision)%>%
  summarise(no_seq =n())%>%
  arrange(desc(no_seq))%>%as.data.frame()

noTax_stats_order <- deleted_pr2_tax %>% group_by(supergroup, division,subdivision,class, order)%>%
  summarise(no_seq =n())%>%
  arrange(desc(no_seq))%>%as.data.frame()


noTax_stats_family <- deleted_pr2_tax %>% group_by(supergroup, division,subdivision,class, order,family)%>%
  filter(subdivision != "Fungi")%>%
  summarise(no =n())%>%
  arrange(desc(no))


noTax_stats_plants <- deleted_pr2_tax %>% group_by(supergroup, division,subdivision,class, order,family,genus)%>%
  filter(order == "Embryophyceae_X")%>%
  summarise(no =n())%>%
  arrange(desc(no))


write.table(noTax_stats_subdvision, file = file.path(output_path, "15__stats_deletedSeq_noTax_subdivisionLevel.tsv"), row.names = F, sep="\t")
write.table(noTax_stats_order, file = file.path(output_path, "15__stats_deletedSeq_noTax_orderLevel.tsv"), row.names = F, sep="\t")

rownames(noTax_stats_subdvision) <- noTax_stats_subdvision$subdivision
rownames(noTax_stats_order) <- noTax_stats_order$order

noTax_stats_subdvision["Streptophyta_X","no_seq"]/ nrow(deleted_pr2_tax) *100
noTax_stats_subdvision["Fungi","no_seq"]/ nrow(deleted_pr2_tax) *100
noTax_stats_subdvision["Metazoa","no_seq"]/ nrow(deleted_pr2_tax) *100
noTax_stats_order["Hexapoda","no_seq"]/ nrow(deleted_pr2_tax) *100
noTax_stats_order["Chelicerata","no_seq"]/ nrow(deleted_pr2_tax) *100


stats_df <- rbind(stats_df, data.frame(database= "V2", variable ="Deleted Streptophyta_X Sequences No.", value= noTax_stats_subdvision["Streptophyta_X","no_seq"]))
stats_df <- rbind(stats_df, data.frame(database= "V2", variable ="Deleted Streptophyta_X Sequences compared to all deleted (%)", 
                                       value= noTax_stats_subdvision["Streptophyta_X","no_seq"]/ nrow(deleted_pr2_tax) *100))

stats_df <- rbind(stats_df, data.frame(database= "V2", variable ="Deleted Fungi Sequences No.", value= noTax_stats_subdvision["Fungi","no_seq"]))
stats_df <- rbind(stats_df, data.frame(database= "V2", variable ="Deleted Fungi Sequences compared to all deleted (%)", 
                                       value= noTax_stats_subdvision["Fungi","no_seq"]/ nrow(deleted_pr2_tax) *100))

stats_df <- rbind(stats_df, data.frame(database= "V2", variable ="Deleted Metazoa Sequences No.", value= noTax_stats_subdvision["Metazoa","no_seq"]))
stats_df <- rbind(stats_df, data.frame(database= "V2", variable ="Deleted Metazoa Sequences compared to all deleted (%)", 
                                       value= noTax_stats_subdvision["Metazoa","no_seq"]/ nrow(deleted_pr2_tax) *100))

stats_df <- rbind(stats_df, data.frame(database= "V2", variable ="Deleted Hexapoda Sequences No.", value= noTax_stats_order["Hexapoda","no_seq"]))
stats_df <- rbind(stats_df, data.frame(database= "V2", variable ="Deleted Hexapoda Sequences compared to all deleted (%)", 
                                       value= noTax_stats_order["Hexapoda","no_seq"]/ nrow(deleted_pr2_tax) *100))

stats_df <- rbind(stats_df, data.frame(database= "V2", variable ="Deleted Chelicerata Sequences No.", value= noTax_stats_order["Chelicerata","no_seq"]))
stats_df <- rbind(stats_df, data.frame(database= "V2", variable ="Deleted Chelicerata Sequences compared to all deleted (%)", 
                                       value= noTax_stats_order["Chelicerata","no_seq"]/ nrow(deleted_pr2_tax) *100))

## stats databases:

### database 1 
db1 <-  read.table(file.path(path,"DB_V1_Taxonomy.tax"), sep = " ")%>%
  mutate(order = str_split_i(V2, ";", 6))%>%
  mutate(family = str_split_i(V2, ";", 7))%>%
  mutate(genus = str_split_i(V2, ";",8))%>%
  mutate(species = str_split_i(V2, ";", 9))%>%
  mutate(Family= ifelse(grepl("_X", family), NA, family))%>%
  mutate(Genus= ifelse(grepl("_X", genus), NA, genus))%>%
  mutate(Species= ifelse(grepl("_sp.", species), NA, species))


nrow(db1)
length(unique(db1$order))
length(unique(db1$Family)%>%na.omit())
length(unique(db1$Genus)%>%na.omit())
length(unique(db1$Species)%>%na.omit())

stats_df_DB <- data.frame()

stats <- calculate_refstats(db1, "V1")
stats_df_DB <- rbind(stats_df_DB, stats)

### database 2 
db2<-  read.table(file.path(path,"DB_V2_Taxonomy.tax"), sep = " ")%>%
  mutate(order = str_split_i(V2, ";", 6))%>%
  mutate(family = str_split_i(V2, ";", 7))%>%
  mutate(genus = str_split_i(V2, ";",8))%>%
  mutate(species = str_split_i(V2, ";", 9))%>%
  mutate(Family= ifelse(grepl("_X", family), NA, family))%>%
  mutate(Genus= ifelse(grepl("_X", genus), NA, genus))%>%
  mutate(Species= ifelse(grepl("_sp.", species), NA, species))


nrow(db2)
length(unique(db2$order))
length(unique(db2$Family)%>%na.omit())
length(unique(db2$Genus)%>%na.omit())
length(unique(db2$Species)%>%na.omit())


stats <- calculate_refstats(db2, "V2")
stats_df_DB <- rbind(stats_df_DB, stats)

### database 3
db3<-  read.table(file.path(path,"DB_V3_Taxonomy.tax"), sep = " ")%>%
  mutate(order = str_split_i(V2, ";", 6))%>%
  mutate(family = str_split_i(V2, ";", 7))%>%
  mutate(genus = str_split_i(V2, ";",8))%>%
  mutate(species = str_split_i(V2, ";", 9))%>%
  mutate(Family= ifelse(grepl("_X", family), NA, family))%>%
  mutate(Genus= ifelse(grepl("_X", genus), NA, genus))%>%
  mutate(Species= ifelse(grepl("_sp.", species), NA, species))


nrow(db3)
length(unique(db3$order))
length(unique(db3$Family)%>%na.omit())
length(unique(db3$Genus)%>%na.omit())
length(unique(db3$Species)%>%na.omit())



stats <- calculate_refstats(db3, "V3")
stats_df_DB <- rbind(stats_df_DB, stats)

### database 4
db4<-  read.table(file.path(path,"DB_V4_Taxonomy.tax"), sep = " ")%>%
  mutate(order = str_split_i(V2, ";", 6))%>%
  mutate(family = str_split_i(V2, ";", 7))%>%
  mutate(genus = str_split_i(V2, ";",8))%>%
  mutate(species = str_split_i(V2, ";", 9))%>%
  mutate(Family= ifelse(grepl("_X", family), NA, family))%>%
  mutate(Genus= ifelse(grepl("_X", genus), NA, genus))%>%
  mutate(Species= ifelse(grepl("_sp.", species), NA, species))

#records
nrow(db4)
length(unique(db4$order))
length(unique(db4$Family)%>%na.omit())
length(unique(db4$Genus)%>%na.omit())
length(unique(db4$Species)%>%na.omit())
#sequences
length(db4 %>% filter(is.na(Species)) %>% pull(Genus)%>%na.omit()) #genus level assignment
length(db4$Species%>%na.omit())

stats <- calculate_refstats(db4, "V4")
stats_df_DB <- rbind(stats_df_DB, stats)

db4_amp<-  read.table(file.path(path,crabs_database4), sep = "\t")%>%
  mutate(order = str_split_i(V2, ";", 6))%>%
  mutate(family = str_split_i(V2, ";", 7))%>%
  mutate(genus = str_split_i(V2, ";",8))%>%
  mutate(species = str_split_i(V2, ";", 9))%>%
  mutate(Family= ifelse(grepl("_X", family), NA, family))%>%
  mutate(Genus= ifelse(grepl("_X", genus), NA, genus))%>%
  mutate(Species= ifelse(grepl("_sp.", species), NA, species))


nrow(db4_amp)
length(unique(db4_amp$order))
length(unique(db4_amp$Family)%>%na.omit())
length(unique(db4_amp$Genus)%>%na.omit())
length(unique(db4_amp$Species)%>%na.omit())

stats <- calculate_refstats(db4_amp, "V4 CRABS")
stats_df_DB <- rbind(stats_df_DB, stats)





#deleted sequences
nrow(db1)- nrow(db2)
nrow(db4)- nrow(db4_amp)


stats_deletedSeq <- stats_df%>% unique() %>% 
      filter(grepl("Deleted", variable))
write.table(stats_deletedSeq, file = file.path(output_path, "15_Summary_stats_deletedSeq.tsv"), row.names = F, sep="\t")


stats_taxo <- stats_df%>% unique() %>% 
  filter(grepl("taxonomy=", variable))
write.table(stats_taxo, file = file.path(output_path, "15_Summary_stats_V4_taxonomyStrategy.tsv"), row.names = F, sep="\t")


stats_many <- stats_df%>% unique() %>% 
  filter(!grepl("taxonomy=", variable) & !grepl("Deleted", variable))
write.table(stats_many, file = file.path(output_path, "15_Summary_stats_joined.tsv"), row.names = F, sep="\t")


stats_V4 <- stats_df%>% unique() %>% 
  filter(grepl("V4", database)) %>% 
  filter(!grepl("taxonomy=", variable))
write.table(stats_V4, file = file.path(output_path, "15_Summary_stats_DBV4.tsv"), row.names = F, sep="\t")


stats_rest <- stats_df%>% unique() %>% 
  filter(!grepl("V4", database)) 
write.table(stats_rest, file = file.path(output_path, "15_Summary_stats_NON_DBV4.tsv"), row.names = F, sep="\t")


stats_df_DB <- stats_df_DB %>% pivot_wider(names_from = "database", values_from = "value")%>% arrange(variable)
write.table(stats_df_DB, file = file.path(output_path, "15_Summary_stats_general.tsv"), row.names = F, sep="\t")


