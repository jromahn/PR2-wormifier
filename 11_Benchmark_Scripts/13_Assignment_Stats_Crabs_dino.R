################################################################################
#This script can be used to asses how the different database versions performed
#regarding taxonomic identification
################################################################################
#Run after mothur or 12_Convert_header_Crabs_part2.R
#TAxon: Dinophyceae
################################################################################
rm(list = ls())

require(tidyverse)
require(stringr)
require(ggplot2)
require(ggpubr)
require(tibble)
library(gapminder) 

options(scipen = 999) # stop scientific notation
#setwd("/PATH/TO/02_Scripts_folder") # to path above the Script folder
source("00_Function_Library.R") # read file with functions
source("11_Benchmark_Scripts/00_Functions_Benchmark.R") # read file with functions
read_simple_ini("00_login_data.ini")

## create output folder
input_path <- "11_Benchmark_Database_Versions_Euka02"
output_path <- "13_Analyses_2026_Euka02"

output_path_tab <- file.path(output_path, "tables")
if (!dir.exists(output_path)){ dir.create(output_path) }
if (!dir.exists(output_path_tab)){ dir.create(output_path_tab) }

cleaned_community_file <- "00_input/Supplementary_CommunityMatrix_Cleaned_Euka02.csv"
mothur_results_db1 <- file.path(input_path, "9_PhytoArk_Euka02_MUC__Taxonomy_amplicons__DB1.wang.taxonomy")
mothur_results_db2 <- file.path(input_path,"9_PhytoArk_Euka02_MUC__Taxonomy_amplicons__DB2.wang.taxonomy")
mothur_results_db3 <- file.path(input_path,"9_PhytoArk_Euka02_MUC__Taxonomy_amplicons__DB3.wang.taxonomy")
mothur_results_db4 <- file.path(input_path,"9_PhytoArk_Euka02_MUC__Taxonomy_amplicons__DB4.wang.taxonomy")

color_ciliate <- c("#669bbc", "#073b4c")
color_dino <- c("darkseagreen3", "#52796f")



## list of Baltic Species
baltic_list <- file.path(input_path, input_table)
baltic_dinophyceae <- read.table(baltic_list, header = T, sep=",") %>% filter(taxon=="Dinophyceae")

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

################################## (1) Assignment done with original PR2-database and PR2 taxonomy
#ASVs and Taxonomic Assignment
results_1_mothur <- read.table(mothur_results_db1)
results_1_clean <- clean_mothur_results( results_1_mothur)

#Classified Ciliates (all, species level, genus level)
identif <- identify_classified_reads(results_1_clean, counttable, Subclass, "Dinophyceae", "v1")
stats <- identif[[1]]
classified_1_Dino_fam <- identif[[2]]
classified_1_Dino_gen <- identif[[3]]
classified_1_Dino_spec <- identif[[4]]
classified_1_Dino <- identif[[5]]
assignment_stats <- rbind(assignment_stats,stats)

################################## (2) Assignment with original PR2 taxonomy but only sequences where a WoRMS entry exist
#ASVs and Taxonomic Assignment
results_2_mothur <- read.table(mothur_results_db2)
results_2_clean <- clean_mothur_results( results_2_mothur)


#Classified Ciliates (all, species level, genus level)
identif <- identify_classified_reads(results_2_clean,counttable, Subclass, "Dinophyceae", "v2")
stats <- identif[[1]]
classified_2_Dino_fam <- identif[[2]]
classified_2_Dino_gen <- identif[[3]]
classified_2_Dino_spec <- identif[[4]]
classified_2_Dino <- identif[[5]]
assignment_stats <- rbind(assignment_stats,stats)


################################## (3) Assignment with cleaned PR2 database and WoRMS Taxonomy
#ASVs and Taxonomic Assignment
results_3_mothur <- read.table(mothur_results_db3)
results_3_clean <- clean_mothur_results( results_3_mothur)

#Classified Ciliates (all, species level, genus level)
identif <- identify_classified_reads(results_3_clean,counttable, Class, "Dinophyceae", "v3")
stats <- identif[[1]]
classified_3_Dino_fam <- identif[[2]]
classified_3_Dino_gen <- identif[[3]]
classified_3_Dino_spec <- identif[[4]]
classified_3_Dino <- identif[[5]]
assignment_stats <- rbind(assignment_stats,stats)

################################## (4) Assignment done with updated PR2 database and WoRMS taxonomy
#ASVs and Taxonomic Assignment
results_4_mothur <- read.table(mothur_results_db4)
results_4_clean <- clean_mothur_results( results_4_mothur)


#Classified Ciliates (all, species level, genus level)
identif <- identify_classified_reads(results_4_clean, counttable,Class, "Dinophyceae", "v4")
stats <- identif[[1]]
classified_4_Dino_fam <- identif[[2]]
classified_4_Dino_gen <- identif[[3]]
classified_4_Dino_spec <- identif[[4]]
classified_4_Dino <- identif[[5]]
assignment_stats <- rbind(assignment_stats,stats)


##compare families
stats_family1 <- classified_1_Dino_fam %>% full_join(classified_2_Dino_fam, by = "V1", suffix = c(".db1",".db2"))
stats_family2 <- classified_3_Dino_fam %>% full_join(classified_4_Dino_fam, by = "V1", suffix = c(".db3",".db4"))
stats_family <- stats_family1 %>% full_join(stats_family2, by = "V1")
write.table(stats_family, file=file.path(output_path_tab,"13_Dinophyceae_family_list_RefDBcomparison.tsv"), sep="\t", row.names = F)
rm(stats_family1, stats_family2)

#investigating
new_family_ASV_db2 <-stats_family %>% filter( is.na(Family.db1)& !is.na(Family.db2))%>% pull(V1)
results_1_clean %>% filter(V1 %in%new_family_ASV_db2 )%>% select(Class, Subclass, Order, Family)%>%unique

new_family_ASV_db3 <-stats_family %>% filter( is.na(Family.db2)& !is.na(Family.db3))%>% pull(V1)
results_2_clean %>% filter(V1 %in%new_family_ASV_db3 )%>% select(Class, Subclass, Order, Family)%>%unique()


### baltic genera
baltic_4_Dino_gen <- classified_4_Dino_gen %>% filter(Genus %in% baltic_dinophyceae$genus | grepl("Biecheleria|Apocalathium|Tripos", Genus) )
balticNon_4_Dino_gen <- classified_4_Dino_gen %>% filter(!Genus %in% baltic_dinophyceae$genus & !grepl("Biecheleria|Apocalathium|Tripos", Genus))
nrow(classified_4_Dino_gen)
nrow(baltic_4_Dino_gen)
baltic <-nrow(baltic_4_Dino_gen)/nrow(classified_4_Dino_gen)*100
print(baltic)


baltic_4_Dino_spec <- classified_4_Dino_spec %>% mutate(Species =gsub("_", " ", Species)) %>% filter(Species %in% baltic_dinophyceae$species)
nrow(classified_4_Dino_spec)
nrow(baltic_4_Dino_spec)
################################################################################
################################################################################

### Comparisons

#Overview over all databases (species level) for all sequences assigned to 
#ciliates
overview_1 <- as.data.frame.table(table(classified_1_Dino$Species))
overview_2 <- as.data.frame.table(table(classified_2_Dino$Species))
overview_3 <- as.data.frame.table(table(classified_3_Dino$Species))
overview_4 <- as.data.frame.table(table(classified_4_Dino$Species))

#comparison of numbers of ASVs for each assignment
comparison_num_seq <- full_join(overview_1, overview_2, by = "Var1")
comparison_tax <- full_join(overview_2, overview_3, by = "Var1") 
comparison_seq <- full_join(overview_3, overview_4, by = "Var1")

#Comparison of assignment for each ASV (all, species level, genus level)
comparison_ASVs <- full_join(classified_1_Dino, classified_2_Dino, by = "V1", suffix = c(".v1", ".v2")) %>%
  full_join(full_join(classified_3_Dino, classified_4_Dino, by = "V1", suffix = c(".v3", ".v4")), by = "V1") 

comparison_ASV_Spec <- full_join(classified_1_Dino_spec, classified_2_Dino_spec, by = "V1", suffix = c(".v1", ".v2")) %>%
  full_join(full_join(classified_3_Dino_spec, classified_4_Dino_spec, by = "V1", suffix = c(".v3", ".v4")), by = "V1") 

comparison_ASV_Gen <- full_join(classified_1_Dino_gen, classified_2_Dino_gen, by = "V1", suffix = c(".v1", ".v2")) %>%
  full_join(full_join(classified_3_Dino_gen, classified_4_Dino_gen, by = "V1", suffix = c(".v3", ".v4")), by = "V1") 


############ Compare  the different genus assignments
compare_genus_level(classified_1_Dino_gen, classified_2_Dino_gen, classified_3_Dino_gen, classified_4_Dino_gen,baltic_dinophyceae, "Dinophyceae", output_path_tab)


#Number classified as ciliates (ASVs)
number_classified <- c(sum(!is.na(comparison_ASVs$Species.x)), sum(!is.na(comparison_ASVs$Species.y)), sum(!is.na(comparison_ASVs$Species.x.x)), sum(!is.na(comparison_ASVs$Spec_perc.y.y)))

#Comparison of ASVs that were assigned identical with all 4 versions
comparison_ASVs_id <- comparison_ASVs %>% 
  filter(Species.v1 == Species.v2 & Species.v2 == Species.v3 & Species.v3 == Species.v4) %>% 
  mutate(Spec_perc.v1 = as.numeric(Spec_perc.v1), 
         Spec_perc.v2 = as.numeric(Spec_perc.v2), 
         Spec_perc.v3 = as.numeric(Spec_perc.v3), 
         Spec_perc.v4 = as.numeric(Spec_perc.v4))

comparison_ASV_Spec_id <- comparison_ASV_Spec %>% 
  filter(Species.v1 == Species.v2 & Species.v2 == Species.v3 & Species.v3 == Species.v4) %>% 
  mutate(Spec_perc.v1 = as.numeric(Spec_perc.v1), 
         Spec_perc.v2 = as.numeric(Spec_perc.v2), 
         Spec_perc.v3 = as.numeric(Spec_perc.v3), 
         Spec_perc.v4 = as.numeric(Spec_perc.v4))

comparison_ASV_Gen_id <- comparison_ASV_Gen %>% 
  filter(Genus.v1 == Genus.v2 & Genus.v2 == Genus.v3 & Genus.v3 == Genus.v4) %>%
  mutate(Genu_perc.v1 = as.numeric(Genu_perc.v1), 
         Genu_perc.v2 = as.numeric(Genu_perc.v2), 
         Genu_perc.v3 = as.numeric(Genu_perc.v3), 
         Genu_perc.v4 = as.numeric(Genu_perc.v4))

#Changes in bootstrap values for identical assigned ASVs
plot_data <- data.frame(v1_v2 = comparison_ASV_Gen_id$Genu_perc.v2 - comparison_ASV_Gen_id$Genu_perc.v1,
           v2_v3 = comparison_ASV_Gen_id$Genu_perc.v3 - comparison_ASV_Gen_id$Genu_perc.v2,
           v3_v4 = comparison_ASV_Gen_id$Genu_perc.v4 - comparison_ASV_Gen_id$Genu_perc.v3)%>%
  pivot_longer(everything(),names_to = "comparison", values_to = "difference")%>%
  mutate(Comparison= case_when(
    comparison=="v1_v2" ~ "Remove PR2 sequences",
    comparison=="v2_v3" ~ "Change taxonomy",
    comparison=="v3_v4" ~ "Add NCBI sequences"  ))%>%
  mutate(Comparison= factor(Comparison, levels=c("Remove PR2 sequences", "Change taxonomy", "Add NCBI sequences")))

ggplot(plot_data, aes(x=Comparison, y= difference))+
    geom_boxplot(fill="#7AC5CD")+theme_light()+
    labs(title="Changes in bootstrap values on the genus level", y="Difference", x="Changes in the reference database")+
    theme(axis.title.x = element_text(vjust = -2), axis.title.y = element_text(vjust = 2), title = element_text(vjust = 2)) 
ggsave(file.path(output_path, "13_Figure_bootstrap_changes_Dinophyceae.jpg"), dpi = 300, width=8, height = 6)
saveRDS(plot_data, file.path(output_path,"13_plotData_bootstrap_stats_Euka02_Din.RDS"))


################################################################################
################################################################################

# print supplementary tables

write.table(assignment_stats, file=file.path(output_path_tab,"13_SUMMARY_Read_ASV_statistics_Dinophyceae.tsv"), sep="\t", row.names = F)
write.table(comparison_ASVs_id, file=file.path(output_path_tab,"13_Dinophyceae_Analysis_RefDB_comparison__assigned.tsv"), sep="\t", row.names = F)
write.table(comparison_ASV_Spec_id, file=file.path(output_path_tab,"13_Dinophyceae_Analysis_RefDB_comparison__species.tsv"), sep="\t", row.names = F)
write.table(comparison_ASV_Gen_id, file=file.path(output_path_tab,"13_Dinophyceae_Analysis_RefDB_comparison__genus.tsv"), sep="\t", row.names = F)

################################################################################
################################################################################


#"real" species level assignments and additonal! genus level assignments (read 
#that are assigned on species level are not included on genus level)
spec_4 <- classified_4_Dino_spec %>% filter(str_detect(Species, "_sp.", negate = TRUE))
gen_4 <- classified_4_Dino_gen %>% filter(!V1 %in% spec_4$V1)

spec_3 <- classified_3_Dino_spec %>% filter(str_detect(Species, "_sp.", negate = TRUE))
gen_3 <- classified_3_Dino_gen %>% filter(!V1 %in% spec_3$V1)

spec_2 <- classified_2_Dino_spec %>% filter(str_detect(Species, "_sp.", negate = TRUE))
gen_2 <- classified_2_Dino_gen %>% filter(!V1 %in% spec_2$V1)

spec_1 <- classified_1_Dino_spec %>% filter(str_detect(Species, "_sp.", negate = TRUE))
gen_1 <- classified_1_Dino_gen %>% filter(!V1 %in% spec_1$V1)


############# Plot results
## Number of classified ASVs
tax_unit <- rep(c("species", "genus"), 4)
database <- sort(rep(c("PR2\n(V1)", "PR2 (WoRMS exists)\n(V2)", "PR2 + WoRMS\n(V3)", "PR2 + WoRMS + NCBI Sequences\n(V4)"), 2))
num_classified <- c(nrow(spec_1), nrow(gen_1), nrow(spec_2), nrow(gen_2), nrow(spec_3), nrow(gen_3), nrow(spec_4), nrow(gen_4))

stats_classified <- data.frame(tax_unit, database, num_classified) 

plot1 <- ggplot(stats_classified, aes(fill = tax_unit, y = num_classified, x = database)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = "Database version", y = "Total ASV", fill = "Taxonomic level",
       title="Number of classified ASVs",subtitle = paste("V4 covers ~", round(baltic,2),"% of Baltic genera")) +
  theme_light() +
  scale_fill_manual(values=color_dino, 
                    labels=c("genus"="Genus", "species"="Species"))+
  theme(axis.title.x = element_text(vjust = -2), axis.title.y = element_text(vjust = 2), title = element_text(vjust = 2)) 
print(plot1)
ggsave(file.path(output_path, "13_Figure_stats_classified__identifiedASVs_Dinophyceae.jpg"), dpi = 300, width=9, height = 5)

## Number Ciliates classified (Taxa)
num_classified_tax <- c(n_distinct(spec_1$Species), n_distinct(classified_1_Dino_gen$Genus), n_distinct(spec_2$Species), n_distinct(classified_2_Dino_gen$Genus), n_distinct(spec_3$Species), n_distinct(classified_3_Dino_gen$Genus), n_distinct(spec_4$Species), n_distinct(classified_4_Dino_gen$Genus))
stats_classified_tax <- data.frame(tax_unit, database, num_classified_tax)

plot2 <- ggplot(stats_classified_tax, aes(fill = tax_unit, y = num_classified_tax, x = database)) +
  geom_bar(position = "dodge", stat = "identity") +
  #ggtitle("Number of species and genera identified") +
  labs(x = "Database version", y = "Number of identified taxa", fill = "Taxonomic level") +
  theme_light() +
  scale_fill_manual(values=color_dino, 
                    labels=c("genus"="Genus", "species"="Species"))+
  theme(axis.title.x = element_text(vjust = -2), axis.title.y = element_text(vjust = 2), title = element_text(vjust = 2)) 
print(plot2)
ggsave(file.path(output_path, "13_Figure_stats_classified__identifiedTaxa_Dinophyceae.jpg"), dpi = 300, width=9, height = 5)
  
## Number Ciliates classified (Reads)
num_classified_reads <- c(sum(spec_1$total_reads), sum(gen_1$total_reads), sum(spec_2$total_reads), sum(gen_2$total_reads), sum(spec_3$total_reads), sum(gen_3$total_reads), sum(spec_4$total_reads), sum(gen_4$total_reads))
stats_classified_reads <- data.frame(tax_unit, database, num_classified_reads)

plot3 <- ggplot(stats_classified_reads, aes(fill = tax_unit, y = num_classified_reads, x = database)) +
  geom_bar(position = "stack", stat = "identity") +
  labs(x = "Database version", y = "Total read number", fill = "Taxonomic level") +
  theme_light() +
  scale_fill_manual(values=color_dino, 
                    labels=c("genus"="Genus", "species"="Species"))+
  ggtitle("Number of reads identified")+scale_y_continuous(labels = scales::comma) +
  theme(axis.title.x = element_text(vjust = -2), axis.title.y = element_text(vjust = 2), title = element_text(vjust = 2)) 
print(plot3)
ggsave(file.path(output_path, "13_Figure_stats_classified__identifiedReads_Dinophyceae.jpg"), dpi = 300, width=9, height = 5)

plot_arrange <- theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 25))
ggarrange(plot1 + plot_arrange, # Adjust left margin (l = 50)
          plot2  + plot_arrange ,
          plot3 + plot_arrange , 
          nrow = 3, labels = c("A","B", "C"),
          font.label = list(size = 18), hjust = -0.7)
ggsave(file.path(output_path, "13_Figure_stats_classified__all_Dinophyceae.jpg"), dpi = 300, width=9, height = 15)
################################################################################
################################################################################

#Save cleaned results for Treemap script
saveRDS(results_1_clean, file.path(output_path,"13__Results_DB1_Dino.RDS"))
saveRDS(results_2_clean, file.path(output_path,"13__Results_DB2_Dino.RDS"))
saveRDS(results_3_clean, file.path(output_path,"13__Results_DB3_Dino.RDS"))
saveRDS(results_4_clean, file.path(output_path,"13__Results_DB4_Dino.RDS"))

################################################################################
#keep stats for plotting

stats_classified <- stats_classified %>% 
  mutate(organism= "Dinophyceae", dataset= "TAREuk", method="ASV")%>%
  rename(values= num_classified)

stats_classified_tax<- stats_classified_tax %>% 
  mutate(organism= "Dinophyceae", dataset= "TAREuk", method="Taxa")%>%
  rename(values= num_classified_tax)

stats_classified_reads <- stats_classified_reads%>% 
  mutate(organism= "Dinophyceae", dataset= "TAREuk", method="Reads")%>%
  rename(values= num_classified_reads)

stats_classified_combined <- stats_classified %>% bind_rows(stats_classified_tax) %>% bind_rows(stats_classified_reads)
saveRDS(stats_classified_combined, file.path(output_path,"13__combined_stats_Euka02_Dino.RDS"))
