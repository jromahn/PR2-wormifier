
####### clean mothur results
clean_mothur_results <- function(data){
  data_clean <- data %>%
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
  
  return(data_clean)
}

identify_classified_reads <- function( data, counttable, rank_filter, taxon_filter, version){
  classified_1_Cil <- data %>%
    filter({{rank_filter}} == taxon_filter) %>%
    select(V1, Species, Spec_perc) %>%
    inner_join(counttable, by = "V1") %>%
    filter(total_reads > 0)
  print(nrow(classified_1_Cil))
  
  classified_1_Cil_spec <- data %>%
    filter({{rank_filter}} == taxon_filter) %>%
    filter(str_detect(Species, "unclassified", negate = TRUE)) %>%
    filter(str_detect(Species, "sp.", negate = TRUE)) %>%
    filter(str_detect(Species, "sp", negate = TRUE)) %>%
    filter(str_detect(Species, "_X", negate = TRUE)) %>%
    select(V1, Species, Spec_perc) %>%
    inner_join(counttable, by = "V1") %>%
    filter(total_reads > 0)
  
  classified_1_Cil_gen <- data %>%
    filter({{rank_filter}} == taxon_filter) %>%
    filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
    filter(str_detect(Genus, "_X", negate = TRUE)) %>%
    select(V1, Genus, Genu_perc) %>%
    inner_join(counttable, by = "V1") %>%
    filter(total_reads > 0)
  
  classified_1_Cil_fam <- data %>%
    filter({{rank_filter}} == taxon_filter) %>%
    filter(str_detect(Family, "unclassified", negate = TRUE)) %>%
    filter(str_detect(Family, "_X", negate = TRUE)) %>%
    select(V1, Family, Fami_perc) %>%
    inner_join(counttable, by = "V1") %>%
    filter(total_reads > 0)
  
  #Classified (all ASVs classified to species/genus level)
  classified_1_spec <- data %>%
    filter(str_detect(Species, "unclassified", negate = TRUE)) %>%
    filter(str_detect(Species, "sp.", negate = TRUE)) %>%
    filter(str_detect(Species, "sp", negate = TRUE)) %>%
    filter(str_detect(Species, "_X", negate = TRUE)) %>%
    select(V1, Species, Spec_perc) %>%
    inner_join(counttable, by = "V1") %>%
    filter(total_reads > 0)
  
  classified_1_gen <- data %>%
    filter(str_detect(Genus, "unclassified", negate = TRUE)) %>%
    filter(str_detect(Genus, "_X", negate = TRUE)) %>%
    select(V1, Genus, Genu_perc) %>%
    inner_join(counttable, by = "V1") %>%
    filter(total_reads > 0)
  summed_reads <- sum(counttable$total_reads)
  summed_reads_target <- data %>%
    filter({{rank_filter}} == taxon_filter) %>%
    select(V1) %>%
    inner_join(counttable, by = "V1") %>%
    pull(total_reads) %>% sum()
  
  stats <- data.frame(database=version) %>% 
    mutate(genus= nrow(classified_1_gen),
           species= nrow(classified_1_spec),
           ciliate= nrow(classified_1_Cil),
           cil_fam = nrow(classified_1_Cil_fam),
           cil_gen_no = length(unique(classified_1_Cil_gen$Genus)),
           cil_gen = nrow(classified_1_Cil_gen),
           cil_spec = length(unique(classified_1_Cil_spec$Species)),
           cil_ASV = nrow(classified_1_Cil_spec),
           measurement= "ASV")%>%
    bind_rows( data.frame(database=version, measurement="Reads") %>%
                 mutate(genus= sum(classified_1_gen$total_reads),
                        species= sum(classified_1_spec$total_reads),
                        ciliate= sum(classified_1_Cil$total_reads),
                        cil_fam = sum(classified_1_Cil_fam$total_reads),
                        cil_gen = sum(classified_1_Cil_gen$total_reads),
                        cil_spec = sum(classified_1_Cil_spec$total_reads)))%>%
    bind_rows( data.frame(database=version, measurement="ReadsProportion") %>%
                 mutate(genus= round(sum(classified_1_gen$total_reads)/summed_reads *100,1),
                        species= round(sum(classified_1_spec$total_reads)/summed_reads *100,1),
                        ciliate= round(sum(classified_1_Cil$total_reads)/summed_reads *100,1),
                        cil_fam = round(sum(classified_1_Cil_fam$total_reads)/summed_reads_target *100,1),
                        cil_gen = round(sum(classified_1_Cil_gen$total_reads)/summed_reads_target *100,1),
                        cil_spec = round(sum(classified_1_Cil_spec$total_reads)/summed_reads_target *100,1)))
  
  data_new <- list(stats, classified_1_Cil_fam, classified_1_Cil_gen, classified_1_Cil_spec, classified_1_Cil)
  return(data_new)
  
}

## comparison asv and read number on genus level
compare_genus_level <- function(db1_gen, db2_gen, db3_gen, db4_gen, baltic_list, taxon ,output_path){
  comparison_genus <-db1_gen %>% mutate(database="v1")%>%
    bind_rows(db2_gen %>% mutate(database="v2"))%>%
    bind_rows(db3_gen %>% mutate(database="v3"))%>%
    bind_rows(db4_gen %>% mutate(database="v4"))%>%
    group_by(database, Genus)%>%
    summarise(ASV= n(),
              reads=sum(total_reads))%>%
    ungroup()%>%
    arrange(Genus, database)%>%
    mutate(Baltic = case_when(Genus %in% baltic_list$genus ~ "TRUE",
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
      reads_v4 = ifelse(is.na(reads_v4), 0, reads_v4)) %>%
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
      ASV_v4 = ifelse(is.na(ASV_v4), 0, ASV_v4)) %>%
    mutate(ASV_diff_v1.2 = ASV_v2 - ASV_v1,
           ASV_diff_v1.3 = ASV_v3 - ASV_v1,
           ASV_diff_v2.3 = ASV_v3 - ASV_v2,
           ASV_diff_v1.4 = ASV_v4 - ASV_v1,
           ASV_diff_v3.4 = ASV_v4 - ASV_v3)
  
  comparison_genus_ASV_reads <- comparison_genus_ASV %>% full_join(comparison_genus_reads, by =c("Baltic", "Genus"))
  write.table(comparison_genus_ASV_reads, file=file.path(output_path,paste("13",taxon,"RefDBcomparison__genus_readsNasvs.tsv", sep="_")), sep="\t", row.names = F)
  write.table(comparison_genus_changes, file=file.path(output_path,paste("13",taxon,"RefDBcomparison__genus_changedReadNo.tsv", sep="_")), sep="\t", row.names = F)
  
}

## transform the DBs

transform_DBs <- function(data){
  data_db <- read.table(data, sep = "\t", header = FALSE) %>%
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
  return(data_db)
}

#Prepare results so that they only contain of a certian class taxon, assigned to at least a given taxonomic level
filter_taxon <- function(data, counttable, rank_col, filter_rank, filter_taxon) {
  data2 <- data %>%
    filter({{filter_rank}} == filter_taxon) %>%
    filter(str_detect({{ rank_col }}, "unclassified", negate = TRUE)) %>%
    filter(str_detect({{ rank_col }}, "_X", negate = TRUE)) %>%
    left_join(counttable, by = "V1") %>%
    filter(!is.na(total_reads)) %>%
    select(-V1, -ends_with("perc")) %>%
    uncount(total_reads)
  
  return(data2)
}

# count ASVs of a certain level for a specific taxon
summarize_stats <-function(data, rank, filter_rank, filter_taxon){
  data2 <- data %>%
    filter({{filter_rank}} == filter_taxon) %>%
    filter(str_detect({{rank}}, "unclassified", negate = TRUE)) %>%
    filter(str_detect({{rank}}, "_X", negate = TRUE))%>%nrow()
}

#######---------------------- Treemap functions

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