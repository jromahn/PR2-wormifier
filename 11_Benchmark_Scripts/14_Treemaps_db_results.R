################################################################################
#This script can be used to create Treemaps aout of the databases or the species
#classification results
################################################################################
#Run after '13_Assignment_Stats_CRABS.R'
#TAxon: CIliophora
################################################################################
rm(list = ls())

#setwd("/PATH/TO/02_Scripts_folder") # to path above the Script folder
source("00_Function_Library.R") # read file with functions
source("11_Benchmark_Scripts/00_Functions_Benchmark.R") # read file with functions

require(dplyr)
require(ggplot2)
require(treemapify)
require(stringr)
require(tidyr)
require(tibble)
require(ggpubr)
require(ghibli)

input_path <- "13_Analyses_2026_Euka02"
output_path <- "13_Analyses_2026_Euka02"
crabs_path <- "11_Benchmark_Database_Versions_Euka02"
primer <- "Euka02"


cleaned_community_file <- "00_input/Supplementary_CommunityMatrix_Cleaned_Euka02.csv"
crabs_database1 <- file.path(crabs_path,"11_Euka02_Taxonomy_DB1.tax")
crabs_database2 <- file.path(crabs_path,"11_Euka02_Taxonomy_DB2.tax")
crabs_database3 <- file.path(crabs_path,"11_Euka02_Taxonomy_DB3.tax")
crabs_database4 <- file.path(crabs_path,"11_Euka02_Taxonomy_DB4.tax")

#######################################################################################
output_path_tab <- file.path(output_path,"tables")
path_creation(output_path )
path_creation(output_path_tab )
#######################################################################################
#Import the CRABS versions of all 4 database versions
db_1 <- transform_DBs(crabs_database1) 
db_2 <- transform_DBs(crabs_database2) 
db_3 <- transform_DBs(crabs_database3) 
db_4 <- transform_DBs(crabs_database4) 


#Import results of the species classification for cilitate
results_db1 <- readRDS(file.path(input_path,"13__Results_DB1_Cil.RDS"))
results_db2 <- readRDS(file.path(input_path,"13__Results_DB2_Cil.RDS"))
results_db3 <- readRDS(file.path(input_path,"13__Results_DB3_Cil.RDS"))
results_db4 <- readRDS(file.path(input_path,"13__Results_DB4_Cil.RDS"))

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
results_cil_db1_fam <- filter_taxon(results_db1, counttable, Family, Class,"Ciliophora")
results_cil_db4_fam <- filter_taxon(results_db4, counttable, Family, Phylum,"Ciliophora")

#Prepare results so that they only contain ciliates, assigned to at least genus
#level, not removed by cleaning. Number of rows indicate number of reads
results_cil_db1 <- filter_taxon(results_db1, counttable, Genus, Class,"Ciliophora")
results_cil_db2 <- filter_taxon(results_db2, counttable, Genus, Class,"Ciliophora")
results_cil_db3 <- filter_taxon(results_db3, counttable, Genus, Phylum,"Ciliophora")
results_cil_db4 <- filter_taxon(results_db4, counttable, Genus, Phylum,"Ciliophora")

## treemap stats --> readcount count
A_stats <- pr2_treemap_stats(db_4 %>% filter(supergroup == "Ciliophora"), class, genus)

#to family level
B_stats_fam <- pr2_treemap_stats(results_cil_db1_fam, Order, Family)
C_stats_fam <- pr2_treemap_stats(results_cil_db4_fam, Order, Family)

#to genus level
B_stats_gen <- pr2_treemap_stats(results_cil_db1, Family, Genus)
C_stats_gen <- pr2_treemap_stats(results_cil_db4, Family, Genus)
# to species level
B_stats_spec <- pr2_treemap_stats(results_cil_db1, Genus, Species)
C_stats_spec <- pr2_treemap_stats(results_cil_db4, Genus, Species)

######------------- treemap stats --> ASV count
## to family 
ASV_family_db1 <- summarize_stats(results_db1, Family, Class, "Ciliophora" )
ASV_family_db4 <- summarize_stats(results_db4, Family, Phylum, "Ciliophora" )

## to genus 
ASV_genus_db1 <- summarize_stats(results_db1, Genus, Class, "Ciliophora" )
ASV_genus_db4 <- summarize_stats(results_db4, Genus, Phylum, "Ciliophora" )

## to species 
ASV_species_db1 <- summarize_stats(results_db1, Species, Class, "Ciliophora" )
ASV_species_db4 <- summarize_stats(results_db4, Species, Phylum, "Ciliophora" )


stats_assignment <- data.frame(db_version= c("v1","v4"))%>%
          mutate(family_level= c(ASV_family_db1, ASV_family_db4),
                 genus_level= c(ASV_genus_db1, ASV_genus_db4),
                 species_level= c(ASV_species_db1, ASV_species_db4),
                 measurement="ASV")
rm(ASV_family_db1, ASV_family_db4, ASV_genus_db1, ASV_genus_db4,ASV_species_db1, ASV_species_db4)

######-------------
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

write.table(stats_assignment, file=file.path(output_path_tab,"14_RefDB_summary_Ciliate.tsv"), sep="\t", row.names = F)
rm(read_family_db1,read_family_db4, read_genus_db1, read_genus_db4,read_species_db1, read_species_db4)

######-------------
stats_family <- B_stats_fam %>% select(Family, n)%>%
                  full_join(C_stats_fam %>% select(Family, n), by ="Family", suffix = c(".db1", ".db4"))

stats_genus <- B_stats_gen %>% select(Genus, n)%>%
  full_join(C_stats_gen %>% select(Genus, n), by ="Genus", suffix = c(".db1", ".db4"))

stats_species <- B_stats_spec %>% select(Species, n)%>%
  full_join(C_stats_spec %>% select(Species, n), by ="Species", suffix = c(".db1", ".db4"))


write.table(stats_family, file=file.path(output_path_tab,"14_RefDB_comparison__reads_family_Ciliate.tsv"), sep="\t", row.names = F)
write.table(stats_genus, file=file.path(output_path_tab,"14_RefDB_comparison__reads_genus_Ciliate.tsv"), sep="\t", row.names = F)
write.table(stats_species, file=file.path(output_path_tab,"14_RefDB_comparison__reads_species_Ciliate.tsv"), sep="\t", row.names = F)

rm(B_stats_fam, C_stats_fam,B_stats_gen, C_stats_gen, B_stats_spec, C_stats_spec )

####################################### Plot ####################################### 
#Treemaps of the database versions
pr2_treemap_db(db_1 %>% filter(subdivision == "Ciliophora"), class, genus)
pr2_treemap_db(db_2, subdivision, class)
pr2_treemap_db(db_3, domain, class)
A <- pr2_treemap_db(db_4 %>% filter(supergroup == "Ciliophora"), class, genus)+
      labs(title=paste(primer,"Ciliophora", sep=" - "))+
      theme(plot.title= element_text(face="bold", size=28))
  



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
ggsave(file.path(output_path, "14_Figure_stats_treemap_Ciliophora.jpg"), dpi = 300, width=12, height = 8)

