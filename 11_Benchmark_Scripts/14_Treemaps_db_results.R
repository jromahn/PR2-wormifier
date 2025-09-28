################################################################################
#This script can be used to create Treemaps aout of the databases or the species
#classification results
################################################################################
#Run after '13_Assignment_Stats_CRABS.R'
################################################################################
rm(list = ls())

#setwd("/PATH/TO/02_Scripts_folder") # to path above the Script folder
setwd("/PATH/TO/11_Benchmark_Scripts")
source("00_Function_Library.R") # read file with functions

require(dplyr)
require(ggplot2)
require(treemapify)
require(stringr)
require(tidyr)
require(tibble)
require(ggpubr)
require(ghibli)

output_path <- "13_Analyses_2025"
cleaned_community_file <- "00_input/Supplementary_Table_5_CommunityMatrix_Cleaned.csv"
crabs_database1 <- "11_Benchmark_Database_Versions/11_Euka02_Taxonomy_DB1.tax"
crabs_database2 <- "11_Benchmark_Database_Versions/11_Euka02_Taxonomy_DB2.tax"
crabs_database3 <- "11_Benchmark_Database_Versions/11_Euka02_Taxonomy_DB3.tax"
crabs_database4 <- "11_Benchmark_Database_Versions/11_Euka02_Taxonomy_DB4.tax"


#######################################################################################
#Import the CRABS versions of all 4 database versions
db_1 <- read.table(crabs_database1, sep = "\t", header = FALSE) %>%
  mutate(domain = str_split_i(V2, ";", 1)) %>%
  mutate(supergroup = str_split_i(V2, ";", 2)) %>%
  mutate(division = str_split_i(V2, ";", 3)) %>%
  mutate(subdivision = str_split_i(V2, ";", 4)) %>%
  mutate(class = str_split_i(V2, ";", 5)) %>%
  mutate(order = str_split_i(V2, ";", 6)) %>%
  mutate(family = str_split_i(V2, ";", 7)) %>%
  mutate(genus = str_split_i(V2, ";", 8)) %>%
  mutate(species = str_split_i(V2, ";", 9)) %>%
  select(-V2) %>%
  filter

db_2 <- read.table(crabs_database2, sep = "\t", header = FALSE) %>%
  mutate(domain = str_split_i(V2, ";", 1)) %>%
  mutate(supergroup = str_split_i(V2, ";", 2)) %>%
  mutate(division = str_split_i(V2, ";", 3)) %>%
  mutate(subdivision = str_split_i(V2, ";", 4)) %>%
  mutate(class = str_split_i(V2, ";", 5)) %>%
  mutate(order = str_split_i(V2, ";", 6)) %>%
  mutate(family = str_split_i(V2, ";", 7)) %>%
  mutate(genus = str_split_i(V2, ";", 8)) %>%
  mutate(species = str_split_i(V2, ";", 9)) %>%
  select(-V2)

db_3 <- read.table(crabs_database3, sep = "\t", header = FALSE) %>%
  mutate(kingdom = str_split_i(V2, ";", 1)) %>%
  mutate(phylum = str_split_i(V2, ";", 2)) %>%
  mutate(subphylum = str_split_i(V2, ";", 3)) %>%
  mutate(class = str_split_i(V2, ";", 4)) %>%
  mutate(subclass = str_split_i(V2, ";", 5)) %>%
  mutate(order = str_split_i(V2, ";", 6)) %>%
  mutate(family = str_split_i(V2, ";", 7)) %>%
  mutate(genus = str_split_i(V2, ";", 8)) %>%
  mutate(species = str_split_i(V2, ";", 9)) %>%
  select(-V2)

db_4 <- read.table(crabs_database4, sep = "\t", header = FALSE) %>%
  mutate(kingdom = str_split_i(V2, ";", 1)) %>%
  mutate(phylum = str_split_i(V2, ";", 2)) %>%
  mutate(subphylum = str_split_i(V2, ";", 3)) %>%
  mutate(class = str_split_i(V2, ";", 4)) %>%
  mutate(subclass = str_split_i(V2, ";", 5)) %>%
  mutate(order = str_split_i(V2, ";", 6)) %>%
  mutate(family = str_split_i(V2, ";", 7)) %>%
  mutate(genus = str_split_i(V2, ";", 8)) %>%
  mutate(species = str_split_i(V2, ";", 9)) %>%
  select(-V2)

#Import results of the species classification for cilitate
results_db1 <- readRDS(file.path(output_path,"13__Results_DB1"))
results_db2 <- readRDS(file.path(output_path,"13__Results_DB2"))
results_db3 <- readRDS(file.path(output_path,"13__Results_DB3"))
results_db4 <- readRDS(file.path(output_path,"13__Results_DB4"))

#Import counttable
counttable <- read.table(cleaned_community_file, header = TRUE, sep = ",") %>%
  select(-c(sample_id, tag, station, depth, replicate, total_reads)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  rowwise() %>%
  mutate(total_reads = sum(across(starts_with("V")))) %>%
  select(ASV, total_reads) %>%
  dplyr::rename(V1 = ASV)


#Prepare results so that they only contain ciliates, assigned to at least Family
#level, not removed by cleaning. Number of rows indicate number of reads
results_cil_db1_fam <- results_db1 %>%
  filter(Class == "Ciliophora") %>%
  filter(str_detect(Family, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Family, "_X", negate = TRUE)) %>%
  left_join(counttable, by = "V1") %>%
  filter(!is.na(total_reads)) %>%
  select(-V1, -ends_with("perc")) %>%
  uncount(total_reads)

results_cil_db4_fam <- results_db4 %>%
  filter(Phylum == "Ciliophora") %>%
  filter(str_detect(Family, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Family, "_X", negate = TRUE)) %>%
  left_join(counttable, by = "V1") %>%
  filter(!is.na(total_reads)) %>%
  select(-V1, -ends_with("perc")) %>%
  uncount(total_reads)


#Prepare results so that they only contain ciliates, assigned to at least genus
#level, not removed by cleaning. Number of rows indicate number of reads
results_cil_db1 <- results_db1 %>%
  filter(Class == "Ciliophora") %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE)) %>%
  left_join(counttable, by = "V1") %>%
  filter(!is.na(total_reads)) %>%
  select(-V1, -ends_with("perc")) %>%
  uncount(total_reads)

results_cil_db2 <- results_db2 %>%
  filter(Class == "Ciliophora") %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE)) %>%
  left_join(counttable, by = "V1") %>%
  filter(!is.na(total_reads)) %>%
  select(-V1, -ends_with("perc")) %>%
  uncount(total_reads)

results_cil_db3 <- results_db3 %>%
  filter(Phylum == "Ciliophora") %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE)) %>%
  left_join(counttable, by = "V1") %>%
  filter(!is.na(total_reads)) %>%
  select(-V1, -ends_with("perc")) %>%
  uncount(total_reads)

results_cil_db4 <- results_db4 %>%
  filter(Phylum == "Ciliophora") %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE)) %>%
  left_join(counttable, by = "V1") %>%
  filter(!is.na(total_reads)) %>%
  select(-V1, -ends_with("perc")) %>%
  uncount(total_reads)


#pr2_treemap function from PR2 website, with some adaptations for maps of the 
#reference database
pr2_treemap_db <- function(pr2, level1, level2) {
  length_col <- pr2 %>% count({{level1}}) %>% nrow()
  # Group
  pr2_class <- pr2 %>%
    count({{level1}},{{level2}}) %>% 
    ungroup()
  # Do a treemap
  ggplot(pr2_class, aes(area = n, 
                        fill = {{level2}}, 
                        subgroup = {{level1}}, 
                        label = {{level2}})) +
    treemapify::geom_treemap()
  
  plot <- ggplot(pr2_class, aes(area = n, 
                                fill= {{level1}}, 
                                subgroup = {{level1}}, 
                                label = {{level2}})) +
    treemapify::geom_treemap(size = 1) +
    treemapify::geom_treemap_text(colour = "white", place = "centre", grow = TRUE) +
    treemapify::geom_treemap_subgroup_border(size = 3) +
    treemapify::geom_treemap_subgroup_text(place = "centre", grow = FALSE,
                                           alpha = 0.8, colour = "grey85", 
                                           min.size = 0) +
    theme_bw() +
    # ghibli stuff
    scale_fill_manual(values = colorRampPalette(ghibli_palettes$MarnieMedium2)(length_col))+
    guides(fill = FALSE)
  print(plot)
  return(plot)
}

#same function adapted for maps of the genus level results, single reads are 
#removed, allows to add single rows of other results, to achieve consistent 
#colors
pr2_treemap_results <- function(pr2, level1, level2) {
  length_col <- pr2 %>% count({{level1}}) %>% nrow()
  # Group
  pr2_class <- pr2 %>%
    count({{level1}},{{level2}}) %>% 
    mutate(n = replace(n, n == 1, 0)) %>%
    ungroup()

  # Do a treemap
  ggplot(pr2_class, aes(area = n, 
                        fill = {{level2}}, 
                        subgroup = {{level1}}, 
                        label = {{level2}})) +
    treemapify::geom_treemap()
  
  plot <- ggplot(pr2_class, aes(area = n, 
                        fill= {{level1}}, 
                        subgroup = {{level1}}, 
                        label = {{level2}})) +
    treemapify::geom_treemap(size = 0) +
    treemapify::geom_treemap_subgroup_border(size = 2) +
    treemapify::geom_treemap_subgroup_text(place = "centre", grow = FALSE,
                                           alpha = 1, colour = "white", 
                                           min.size = 1) +
    theme_bw() +
    # ghibli stuff
    scale_fill_manual(values = colorRampPalette(ghibli_palettes$YesterdayMedium[-1])(length_col))+
    guides(fill = FALSE)
  return(plot)
}

pr2_treemap_stats<- function(pr2, level1, level2){
  pr2_class <- pr2 %>%
    count({{level1}},{{level2}}) %>% 
    mutate(n = replace(n, n == 1, 0)) %>%
    ungroup()%>%
    arrange(desc(n))
  
  return(pr2_class)
}

## treemap stats --> readcount count
A_stats <- pr2_treemap_stats(db_4 %>% filter(phylum == "Ciliophora"), class, genus)
#to family level
B_stats_fam <- pr2_treemap_stats(results_cil_db1_fam, Order, Family)
C_stats_fam <- pr2_treemap_stats(results_cil_db4_fam, Order, Family)
#to genus level
B_stats_gen <- pr2_treemap_stats(results_cil_db1, Family, Genus)
C_stats_gen <- pr2_treemap_stats(results_cil_db4, Family, Genus)
# to species level
B_stats_spec <- pr2_treemap_stats(results_cil_db1, Genus, Species)
C_stats_spec <- pr2_treemap_stats(results_cil_db4, Genus, Species)

##ASV
## to family 
ASV_family_db1 <- results_db1 %>%
  filter(Class == "Ciliophora") %>%
  filter(str_detect(Family, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Family, "_X", negate = TRUE))%>%nrow()

ASV_family_db4 <-results_db4 %>%
  filter(Phylum == "Ciliophora") %>%
  filter(str_detect(Family, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Family, "_X", negate = TRUE))%>%nrow()

## to genus 
ASV_genus_db1 <-results_db1 %>%
  filter(Class == "Ciliophora") %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE))%>%nrow()

ASV_genus_db4 <-results_db4 %>%
  filter(Phylum == "Ciliophora") %>%
  filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Genus, "_X", negate = TRUE))%>%nrow()

## to species 
ASV_species_db1 <-results_db1 %>%
  filter(Class == "Ciliophora") %>%
  filter(str_detect(Species, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Species, "_X", negate = TRUE))%>%
  filter(str_detect(Species, "_sp", negate = TRUE))%>%nrow()

ASV_species_db4 <- results_db4 %>%
  filter(Phylum == "Ciliophora") %>%
  filter(str_detect(Species, "unclassified", negate = TRUE)) %>%
  filter(str_detect(Species, "_X", negate = TRUE))%>%
  filter(str_detect(Species, "_sp", negate = TRUE))%>%nrow()


stats_assignment <- data.frame(db_version= c("v1","v4"))%>%
          mutate(family_level= c(ASV_family_db1, ASV_family_db4),
                 genus_level= c(ASV_genus_db1, ASV_genus_db4),
                 species_level= c(ASV_species_db1, ASV_species_db4),
                 measurement="ASV")
rm(ASV_family_db1, ASV_family_db4, ASV_genus_db1, ASV_genus_db4,ASV_species_db1, ASV_species_db4)

# reads assigned to family
read_family_db1<- sum(B_stats_fam$n)
read_family_db4<- sum(C_stats_fam$n)

# reads assigned to genus
read_genus_db1<- sum(B_stats_gen$n)
read_genus_db4<- sum(C_stats_gen$n)

# reads assigned to species
read_species_db1<- sum(B_stats_spec$n)
read_species_db4<-sum(C_stats_spec$n)

stats_assignment <- stats_assignment %>%
  bind_rows(data.frame(db_version= c("v1","v4"))%>%
  mutate(family_level= c(read_family_db1, read_family_db4),
         genus_level= c(read_genus_db1, read_genus_db4),
         species_level= c(read_species_db1, read_species_db4),
         measurement="reads"))

write.table(stats_assignment, file=file.path(output_path,"14_Ciliate_RefDB_summary.tsv"), sep="\t", row.names = F)
rm(read_family_db1,read_family_db4, read_genus_db1, read_genus_db4,read_species_db1, read_species_db4)


stats_family <- B_stats_fam %>% select(Family, n)%>%
                  full_join(C_stats_fam %>% select(Family, n), by ="Family", suffix = c(".db1", ".db4"))

stats_genus <- B_stats_gen %>% select(Genus, n)%>%
  full_join(C_stats_gen %>% select(Genus, n), by ="Genus", suffix = c(".db1", ".db4"))

stats_species <- B_stats_spec %>% select(Species, n)%>%
  full_join(C_stats_spec %>% select(Species, n), by ="Species", suffix = c(".db1", ".db4"))

rm(B_stats_fam, C_stats_fam,B_stats_gen, C_stats_gen, B_stats_spec, C_stats_spec )



write.table(stats_family, file=file.path(output_path,"14_Ciliate_RefDB_comparison__reads_family.tsv"), sep="\t", row.names = F)
write.table(stats_genus, file=file.path(output_path,"14_Ciliate_RefDB_comparison__reads_genus.tsv"), sep="\t", row.names = F)
write.table(stats_species, file=file.path(output_path,"14_Ciliate_RefDB_comparison__reads_species.tsv"), sep="\t", row.names = F)


####################################### Plot ####################################### 
#Treemaps of the database versions
pr2_treemap_db(db_1 %>% filter(subdivision == "Ciliophora"), class, genus)
pr2_treemap_db(db_2, subdivision, class)
pr2_treemap_db(db_3, kingdom, class)
A <- pr2_treemap_db(db_4 %>% filter(phylum == "Ciliophora"), class, genus)




length_col <- 12





#Treemaps of the ciliate results
B <- pr2_treemap_results(bind_rows(results_cil_db1, distinct(results_cil_db4)), Genus, Species) 
pr2_treemap_results(results_cil_db2, Genus, Species)
pr2_treemap_results(results_cil_db3, Genus, Species)
C <- pr2_treemap_results(bind_rows(results_cil_db4, distinct(results_cil_db1)), Genus, Species)
#To B and C distinct row of the other dataframe are added so that both will have
#the same amount of levels. Consequently same names will have same colors. These
#"fake" reads are removed by the function above

distinct(bind_rows(results_cil_db4, results_cil_db1))


plot_arrange <- theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 25))

ggarrange(A + plot_arrange, # Adjust left margin (l = 50)
          ggarrange(B  + plot_arrange , C + plot_arrange , 
                    ncol = 2, labels = c("B","C"), 
                    font.label = list(size = 18), hjust = -0.7), 
          nrow = 2, labels = c("A",""), heights = c(2, 1),
          font.label = list(size = 18), hjust = -0.7)
ggsave(file.path(output_path, "14_Figure_stats_treemap.jpg"), dpi = 300, width=12, height = 8)

