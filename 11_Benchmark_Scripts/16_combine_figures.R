################################################################################
#This script can be used to combine the data and the plot of both datasets and both taxonomic groups
################################################################################
#Run after you run '13_Assignment_Stats_CRABS.R' and 
### '13_Assignment_Stats_CRABS_dino.R' for both (TAREuk & Euka02) dataset
################################################################################
rm(list = ls())

require(tidyverse)
require(stringr)
require(ggplot2)
require(ggpubr)
require(tibble)
library(gapminder) 
library(legendry) ## group legend
library(ggh4x) ## facet_grid2()

options(scipen = 999) # stop scientific notation
#setwd("/PATH/TO/02_Scripts_folder") # to path above the Script folder
setwd("/Users/jromahn/Documents/XX_PAPERS/2025sub_stefanie_wormifier/PR2-wormifier_april26")
source("00_Function_Library.R") # read file with functions
source("11_Benchmark_Scripts/00_Functions_Benchmark.R") # read file with functions


output_path_TAREuk <- "13_Analyses_2026_TAREuk"
output_path_Euka02 <- "13_Analyses_2026_Euka02"
output_path <- "13_Analyses_2026_TAREuk"


###### input data
stats_Cil_Euka02 <- readRDS( file.path(output_path_Euka02,"13__combined_stats_Euka02_Cil.RDS"))
stats_Dino_Euka02 <- readRDS( file.path(output_path_Euka02,"13__combined_stats_Euka02_Dino.RDS"))
stats_Cil_TAREuk <- readRDS( file.path(output_path_TAREuk,"13__combined_stats_TAREuk_Cil.RDS"))
stats_Dino_TAREuk <- readRDS( file.path(output_path_TAREuk,"13__combined_stats_TAREuk_Dino.RDS"))

plot_Cil_Euka02 <- readRDS( file.path(output_path_Euka02,"13_plotData_bootstrap_stats_Euka02_Cil.RDS"))
plot_Dino_Euka02 <- readRDS( file.path(output_path_Euka02,"13_plotData_bootstrap_stats_Euka02_Din.RDS"))
plot_Cil_TAREuk <- readRDS( file.path(output_path_TAREuk,"13_plotData_bootstrap_stats_TAREuk_Cil.RDS"))
plot_Dino_TAREuk <- readRDS( file.path(output_path_TAREuk,"13_plotData_bootstrap_stats_TAREuk_Din.RDS"))



##############################################################################
# Plot 1-3 Assignment sucess to genus and species level based on assigned taxa, ASVs and reads
##############################################################################
stats_all <- stats_Cil_Euka02 %>%bind_rows(stats_Dino_Euka02) %>% bind_rows(stats_Cil_TAREuk) %>% bind_rows(stats_Dino_TAREuk)%>%
    mutate(tax_unit = case_when( tax_unit=="genus" ~ "Genus",
                                 tax_unit =="species" ~"Species"))%>%
    mutate(database = gsub(" NCBI", "\nNCBI", database))

stats_all$database %>% unique

figure_theme <- theme(title = element_text(vjust = 2, size = 20),
                      axis.title.x = element_text(vjust = -2, size = 18,  hjust = 1), 
                      axis.title.y = element_text(vjust = 2, size = 18), 
                      axis.text.x = element_text( size = 12),
                      axis.text.y = element_text( size = 16),
                      strip.background = element_rect(fill="#636363"),
                      strip.text = element_text( size = 16),
                      legend.title = element_text( hjust = 0.5, size= 18),
                      legend.text = element_text( size = 16),
                      legend.position = "bottom",
                      # Adds space below the plot area (pushes legend down)
                      plot.margin = margin(t = 5, r = 5, b = 20, l = 5, unit = "pt"),
                      
                      # Adds space above the legend box itself
                      legend.box.margin = margin(t = 10, unit = "pt"),
                      legend.title.position = "top") 


color_ciliate <- c("#669bbc", "#073b4c")
color_dino <- c("darkseagreen3", "#52796f")


plot1 <- ggplot(stats_all %>% filter(method=="ASV"), 
                aes(fill = paste(organism,tax_unit), y = values, x = database)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_light()+
  facet_grid(organism ~ dataset, axes = "all_y", axis.labels = "all_x",  scales = "free_y") +
  scale_fill_manual(values=c(color_ciliate,color_dino) )+
  labs(x = "Database version", y = "Total ASV", fill = "Taxonomic level") +
  guides(fill = "legend_group")+
  figure_theme +
  theme(legendry.group.spacing = unit(1, "cm"))
print(plot1)
ggsave(file.path(output_path, "17_Figure_stats_classified__identifiedASVs.jpg"), dpi = 300, width=14, height = 10)


plot2 <- ggplot(stats_all %>% filter(method=="Taxa"), 
                aes(fill = paste(organism,tax_unit), y = values, x = database)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_light()+
  facet_grid(organism ~ dataset, axes = "all_y", axis.labels = "all_x",  scales = "free_y") +
  scale_fill_manual(values=c(color_ciliate,color_dino) )+
  labs(x = "Database version", y = "Total Taxa", fill = "Taxonomic level") +
  guides(fill = "legend_group")+
  figure_theme +
  theme(legendry.group.spacing = unit(1, "cm"))+
  scale_y_continuous(labels = scales::comma)
print(plot2)
ggsave(file.path(output_path, "17_Figure_stats_classified__identifiedTaxa.jpg"), dpi = 300, width=14, height = 10)
ggsave(file.path(output_path, "17_Figure_stats_classified__identifiedTaxa.pdf"), dpi = 300, width=14, height = 10)

plot3 <- ggplot(stats_all %>% filter(method=="Reads"), 
                aes(fill = paste(organism,tax_unit), y = values, x = database)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_light()+
  facet_grid2(organism ~dataset, scales = "free_y", independent = "y") +
  scale_fill_manual(values=c(color_ciliate,color_dino) )+
  labs(x = "Database version", y = "Total Reads", fill = "Taxonomic level") +
  guides(fill = "legend_group")+
  figure_theme +
  theme(legendry.group.spacing = unit(1, "cm"))+
  scale_y_continuous(labels = scales::comma)
print(plot3)
ggsave(file.path(output_path, "17_Figure_stats_classified__identifiedReads.jpg"), dpi = 300, width=14, height = 10)


##############################################################################
# Plot 4 changes in mothur bootstrap value depending on the changes in the database
##############################################################################

plot_data <- plot_Cil_Euka02 %>% mutate(primer="Euka02", taxon = "Ciliophora")%>%
                  bind_rows(plot_Dino_Euka02 %>% mutate(primer="Euka02", taxon = "Dinophyceae"))%>%
                  bind_rows(plot_Cil_TAREuk %>% mutate(primer="TAREuk", taxon = "Ciliophora"))%>%
                  bind_rows(plot_Dino_TAREuk %>% mutate(primer="TAREuk", taxon = "Dinophyceae"))

ggplot(plot_data, aes(x=Comparison, y= difference, fill= taxon))+
  geom_boxplot()+theme_light()+
  labs(title="Changes in bootstrap values on the genus level", 
       y="Difference", 
       x="Changes in the reference database", fill= "Taxonomic group")+
  facet_grid2(taxon ~primer, scales = "free_y", independent = "y") +
  theme(axis.title.x = element_text(vjust = -2), 
        axis.title.y = element_text(vjust = 2), title = element_text(vjust = 2)) +
  figure_theme +
  scale_fill_manual(values=c("#669bbc","darkseagreen3"), 
                    labels=c("Ciliophora", "Dinophyceae"))
ggsave(file.path(output_path, "17_Figure_bootstrap_changes_ALL.jpg"), dpi = 300, width=14, height = 9)
