#################
# create one or several dirs
path_creation <- function(path){
  if(length(path)==1){
    if (!dir.exists(path)){ dir.create(path) }
  }else{
    for (p in path){
      if (!dir.exists(p)){ dir.create(p) }
    }
  }
}

##################
# Filter

#filter out wrongly matched taxa based on pr2 taxonomy
## dependency : dplyr

clean_pr2_taxonomy <- function(data){
  
  data2 <- data %>%
    mutate(Clean_Name = species) %>%
    relocate(.before = domain) %>%
    
    #Remove sequences from cell organelles
    filter(!str_detect(Clean_Name, ":mito")) %>%                                                        
    filter(!str_detect(Clean_Name, ":nucl")) %>%
    filter(!str_detect(Clean_Name, ":plas")) %>%
    filter(!str_detect(Clean_Name, ":chro")) %>%
    filter(!str_detect(Clean_Name, ":apic")) %>%
    
    #Remove sequences that do not have sufficient taxonomic information
    filter(!str_detect(Clean_Name, "_X")) %>%                          #Remove sequences without taxonomic information at least on genus level
    filter(!str_detect(Clean_Name, "_group_")) %>%                     #Remove sequences with "_group_" in their name, indicating there is no proper taxonomic information     
    filter(!str_detect(Clean_Name, "^[[:upper:]][[:upper:]]")) %>%     #Removing sequences where the species name start with two upper letters, indicating it is not a proper name but a fabricated construct
    
    #Clean up species names
    mutate(Clean_Name = str_remove_all(Clean_Name, "-?\\d")) %>%                                     #Remove all digits from species names, including a possible "-" leading up to them
    mutate(Clean_Name = str_remove_all(Clean_Name, "Candidatus_")) %>%                               #Removing the label "Candidatus_" leading up to some bacterial genera
    mutate(Clean_Name = str_replace_all(Clean_Name, "_[[:upper:]]+[[:lower:]]?_", "_")) %>%          #Removing blocks that contain upper letters
    mutate(Clean_Name = str_replace_all(Clean_Name, "_._", "_")) %>%                                 #Removing artefacts
    mutate(Clean_Name = str_replace_all(Clean_Name, "__", "_")) %>%                                  #Replacing double "_"
    mutate(Clean_Name = str_replace_all(Clean_Name, "_\\([[:upper:]][[:lower:]]+\\)", "_sp.")) %>%   #Removing Subgenus information found in brackets
    
    #Correct names for certain species with obvious problems
    mutate(Clean_Name = str_replace(Clean_Name, "Chlorella_Chlorella", "Chlorella_sp.")) %>%
    mutate(Clean_Name = str_replace(Clean_Name, "Micractinium_Micractinium", "Micractinium_sp.")) %>%
    mutate(Clean_Name = str_replace(Clean_Name, "Oncidium_Gower", "Oncidium_sp.")) %>%
    mutate(Clean_Name = str_replace(Clean_Name, "Tergestiellaadriatica", "Tergestiella_adriatica")) %>%
    mutate(Clean_Name = str_replace(Clean_Name, "Isochrysislitoralis", "Isochrysis_litoralis")) %>%
    mutate(Clean_Name = str_replace(Clean_Name, "Cardiosp.oridium_cionae", "Cardiosporidium_cionae")) %>%
    
    filter(species != "Candidatus_Chlorothrix_sp.") %>%                                              #Remove to avoid confusion with Chlorothrix_sp.
    
    #Filter out everything remaining that still does not fit the normal structure of a species name
    filter(str_detect(Clean_Name, "^[[:upper:]][[:lower:]]+-?[[:lower:]]+_[[:lower:]]+\\.?$")) 
  
  return(data2)
}

filter_for_right_taxa <- function(data){
  data2 <- data %>%
  filter(!(order == "Trematoda" & Worms_record_class != "Trematoda")) %>%  
  filter(!(subdivision == "Fungi" & Worms_record_kingdom != "Fungi")) %>%                       
  filter(!(subdivision == "Metazoa" & Worms_record_kingdom != "Animalia")) %>%
  filter(!(subdivision == "Euglenozoa" & Worms_record_phylum != "Euglenozoa")) %>%
  filter(!(subdivision == "Ciliophora" & Worms_record_phylum != "Ciliophora")) %>%
  filter(!(division == "Chlorophyta" & Worms_record_phylum != "Chlorophyta")) %>%
  filter(!(domain == "Bacteria" & Worms_record_kingdom != "Bacteria")) %>%
  filter(!(supergroup == "Archaeplastida" & Worms_record_kingdom != "Plantae")) %>%
  filter(!(supergroup == "TSAR" & !Worms_record_kingdom %in% c("Protozoa","Chromista"))) %>%
  filter(!(supergroup == "Amoebozoa" & Worms_record_phylum != "Amoebozoa")) %>%
  #animalia
    filter(!(family == "Amphibia" & Worms_record_class != "Amphibia")) %>%
    filter(!(class == "Annelida" & Worms_record_phylum != "Annelida")) %>%
    filter(!(class == "Arthropoda" & Worms_record_phylum != "Arthropoda")) %>%
    filter(!(class == "Bryozoa" & Worms_record_phylum != "Bryozoa")) %>%
    filter(!(order == "Ceramiales" & Worms_record_order != "Ceramiales")) %>%
    filter(!(class == "Ctenophora" & Worms_record_phylum != "Ctenophora")) %>%# new
    filter(!(class == "Cnidaria" & Worms_record_phylum != "Cnidaria")) %>%
    filter(!(class == "Craniata" & Worms_record_phylum != "Chordata")) %>%  # new
    filter(!(order == "Hexapoda" & Worms_record_class != "Hexapoda")) %>%
    filter(!(class == "Mollusca" & Worms_record_phylum != "Mollusca")) %>%
    filter(!(class == "Nematoda" & Worms_record_phylum != "Nematoda")) %>%
    filter(!(class == "Porifera" & Worms_record_phylum != "Porifera")) %>%
    filter(!(class == "Rotifera" & Worms_record_phylum != "Rotifera")) %>%# new
    filter(!(family == "Rhabdocoela" & Worms_record_order != "Rhabdocoela")) %>%
    filter(!(order == "Trematoda" & Worms_record_class != "Trematoda")) 
  
  # just to douple check if maybe we lost some
  data3 <- data %>%
    filter(!(  Worms_record_phylum != division | 
               Worms_record_phylum != subdivision | 
               Worms_record_phylum  != class | 
               Worms_record_class   != order |
               Worms_record_class   != class |
               Worms_record_class   != family |
               Worms_record_order   != order  |
               Worms_record_order   != family  ))%>%
    bind_rows(data2)%>%
    unique()
  
  print(paste( "data2:", nrow(data2), "; data3:", nrow(data3)))
  return(data3)
}

filter_for_right_taxa_algaebase <- function(data){
  data2 <- data %>%filter(mapply( grepl, pr2_division , paste(`dwc:phylum`, `dwc:class`, `dwc:order`, `dwc:family`, sep=";"))|
                            mapply( grepl, pr2_subdivision , paste(`dwc:phylum`, `dwc:class`, `dwc:order`, `dwc:family`, sep=";"))|
                            mapply( grepl, pr2_class , paste(`dwc:phylum`, `dwc:class`, `dwc:order`, `dwc:family`, sep=";"))|
                            mapply( grepl, substr(pr2_class,1,8) , paste(`dwc:phylum`, `dwc:class`, `dwc:order`, `dwc:family`, sep=";"))|
                            (pr2_subdivision =="Euglenozoa"& `dwc:phylum`=="Euglenophyta")|
                            (pr2_division =="Streptophyta"& `dwc:kingdom`=="Plantae")|
                            is.na(pr2_subdivision))
  return(data2)
}

filter_algaebase_for_algae <- function(data){
  data2 <- data %>% 
    filter(`dwc:phylum` =="Rhodophyta"|
             `dwc:phylum` =="Heterokontophyta"|
             `dwc:phylum` =="Haptophyta"|
             `dwc:phylum` =="Euglenophyta"|
             `dwc:phylum` =="Dinoflagellata"|
             `dwc:phylum` =="Chromeridophyta"|
             (`dwc:phylum` =="Chlorophyta" & `dwc:class`  =="Ulvophyceae" )|
             `dwc:phylum` =="Charophyta")
  
  return(data2)
}

filter_pr2_for_algae <- function(data){
  data2 <- data %>% 
    filter(subdivision =="Dinoflagellata" |
             division =="Chlorophyta" |
             division =="Rhodophyta" |
             supergroup == "Cryptista" |
             division =="Cyanobacteria" |
             division =="Haptophyta"|
             division =="Stramenopiles" |
             division =="Streptophyta"|
             subdivision =="Euglenozoa" )
  
  return(data2)
}

##################
#Worms functions  - single entry no vectors or list
# dependency: "worrms"

try_worms <- function(name, marine_only = FALSE) {
  tryCatch(wm_records_taxamatch(name, marine_only = FALSE), error = function(err){NA})
}

try_synonyms <- function(AphiaID) {
  tryCatch(wm_synonyms(AphiaID), error = function(err){NA})
}

try_record <- function(AphiaID) {
  tryCatch(wm_record(AphiaID), error = function(err){NA})
}

try_taxonomy <- function(AphiaID) {
  tryCatch(wm_classification(AphiaID), error = function(err){NA})
}

try_names <- function(name, marine_only = FALSE) {
  tryCatch(wm_records_name(name, marine_only = FALSE, fuzzy = FALSE), error = function(err){NA})
}


#################
# Algaebase functions
### dependencies jsonlite & curl

# input: data frame with at least three columns: genus, species epithel, taxonomic group
check_species_against_Algaebase <- function(data, key, password){
  colnames(data)[c(1:3)] <- c("genus", "spec", "taxon")
  data <- unique(data)
  data <- data %>% mutate(spec =gsub("(cf. |\n|\\[|\\]|\\(.*\\)|∗)","", spec))
  
  algaebase.species <- data.frame()
  non.algaebase.df <- data.frame()
  failed <- data.frame()
  
  pb <- txtProgressBar(min = 1, max = nrow(data), style = 3)
  for (i in c(1:nrow(data))){

    date_today <- Sys.Date() 
    genus <- data$genus[i]; spec <- data$spec[i]; Taxon <-data$taxon[i]
    
    # api adresss
    quest_spec <- paste( "https://api.algaebase.org/v1.3/species?scientificname=", genus,"%20",spec, "&taxonrank=species", sep="")
    
    #curl command to extract
    curl_command <- c("curl",  "-s",  "-X", "GET",  "-H", paste(key,": " ,password, sep=""),  # HTTP header with your key and password
                      quest_spec ) #"https://api.algaebase.org/v1.3/species?scientificname=Apocalathium%20malmogiense" # URL
    
    # Execute the curl command
    result <- system2(curl_command, stdout = TRUE)
    #check for bad requests
    if (validate(result)){
        objext <- fromJSON(result)
        if("result" %in% names(objext)) {
          data.spec  <- fromJSON(result)$result %>% 
            mutate(request = paste(genus, spec),
                   taxon = Taxon, 
                   date.spec = as.character(date_today))
          algaebase.species <- plyr::rbind.fill(algaebase.species,data.spec )
        }else{
          non.algaebase.df <-rbind(non.algaebase.df,data[i,] )
        }
      }else{
        failed <- rbind(failed, data[i,])
    }
    
    setTxtProgressBar(pb, i)
  }; close(pb); rm(pb, i, data.spec)
  new_data <- list(algaebase.species,non.algaebase.df ,failed )
  return(new_data)
}

algaebase_species_clean_up <- function( data){
  #data <- tax_6E
  # keep  perfect match or matches sharing the first three characters of genus and species to get results for species with spelling mistakes
  data <- data %>%
    mutate(result= paste(`dwc:genus`, `dwc:specificEpithet`)) %>%
    mutate(request_spec = gsub("([\\w\\-ë]+) ([\\w\\-ë\\.]+).*","\\2",request, perl=T)) %>%
    filter( request == result |
              ( ! grepl('sp\\.', request, perl =T) & 
                  substr(request, 1, 3) == substr(`dwc:genus`, 1, 3) &
                  substr(request_spec, 1, 3) ==substr(`dwc:specificEpithet`, 1, 3))) %>%
    mutate(request_spec =NULL)
  
  
  # & keep just one entry per request!!!! so far several may exist because of similar names
  step2 <- data %>%
    group_by(request, taxon, `dwc:genus`,`dwc:specificEpithet`) %>%
    filter(row_number(plyr::desc(`dcterms:modified`)) == 1)
  
  # keep perfect fit
  step2.1 <- step2 %>% group_by(request, taxon)%>%
    filter(request == paste(`dwc:genus`,`dwc:specificEpithet`)) 
  # win back which we loose because of spelling mistake
  step2.2 <- step2[!step2$request %in% step2.1$request,]
  # Combine the results
  data <- bind_rows(step2.1, step2.2) %>%   distinct()
  
  # make competible with alternative names
  data.accepted <- data %>% 
    filter(`dwc:taxonomicStatus` =="currently accepted taxonomically") %>%
    mutate(algaebase.species_name = paste(`dwc:genus`, `dwc:specificEpithet`),
           algaebase.scientificName = `dwc:scientificName`)%>%
    mutate(algaebase.species_name = as.character(algaebase.species_name))
  
  #if Neotype & Holotype exists keep accepted
  data.accepted <- data.accepted %>%
    group_by(taxon, request)%>%
    mutate(n_records =n())
  
  data.accepted <- data.accepted%>%
    filter(n_records >1)%>%
    group_by(taxon, request)%>%
    filter(`dwc:nomenclaturalStatus`=="nom. et typ. cons.")%>%
    bind_rows(data.accepted %>% filter(n_records == 1))%>%
    mutate(n_records=NULL)
  
  ###  if name is not accepted  use accepted name and if status is unclear use search termin and  save as uncertian
  data.alternativ <- data %>% 
    filter(`dwc:taxonomicStatus` !="currently accepted taxonomically") %>% 
    filter(! request %in% data.accepted$request) %>%
    mutate(algaebase.species_name = ifelse(is.na(`dwc:acceptedNameUsage`), paste(`dwc:genus`, `dwc:specificEpithet`),   
                                           gsub("([\\w\\-ë]+) ([\\w\\-ë]+).*","\\1 \\2",iconv(`dwc:acceptedNameUsage`, to='UTF-8'), perl=T)),
           algaebase.scientificName = ifelse(is.na(`dwc:acceptedNameUsage`), "Uncertain",  `dwc:acceptedNameUsage`))%>%
    mutate(result=NULL)

  
  
  ####################################################################
  ##### merge & clean both
  if(nrow(data.alternativ)>0){
    data_FINAL <- rbind(data.accepted, data.alternativ)
  }else{
    data_FINAL <- data.accepted
  }
  
  
  data_FINAL <- data_FINAL %>% 
    mutate(`dwc:scientificName` =NULL,
           `dwc:scientificName` =NULL, 
           `dwc:genus` =NULL, 
           `dwc:specificEpithet`=NULL,
           `dwc:taxonomicStatus`=NULL, 
           `dwc:acceptedNameUsage` = NULL, 
           `dwc:acceptedNameUsageID`=NULL)%>%
    dplyr::rename("Algaebase_search"= request,
                  "acc_Name"= algaebase.species_name)%>%
    mutate(Algaebase_entry= rep("Species", length(Algaebase_search))) %>%
    mutate(acc_Genus= gsub("([\\w\\-ë]+) ([\\w\\-ë]+).*","\\1",acc_Name, perl=T))
  return(data_FINAL)
}

check_genus_n_taxonomy_algaebase <- function(data, key, password ){
  
  colnames(data)[c(1:2)] <- c("taxon", "genus")
  
  data_known <- data.frame()
  data_unknown <- data.frame()
  failed <- data.frame()
  # query
  pb <- txtProgressBar(min = 1, max = nrow(data), style = 3)
  for(i in c(1:nrow(data))){
    date_today <-Sys.Date()
    genus <- data$genus[i]; taxon <- data$taxon[i]
    
    #api link
    #genus <- "Chondria"
    quest_genus <- paste( "https://api.algaebase.org/v1.3/genus?scientificname=", genus, sep="")
    
    # -g enable globbing (if needed) -s silent "-X", "GET",       # HTTP request method (GET)
    curl_command <- c("curl",     "-s", "-X", "GET", "-H", paste(key,": " ,password, sep=""),  # HTTP header with your key and password
                      quest_genus ) #"https://api.algaebase.org/v1.3/genus?scientificname=Apocalathium"  # URL
    
    # Execute the curl command
    result <- system2(curl_command, stdout = TRUE)
    
    if (validate(result)){
      objext <- fromJSON(result)
      result
      
      if("result" %in% names(objext)) {
        data.spec  <- fromJSON(result)$result %>% 
          mutate(request_genus = genus,
                 request_taxon = taxon,
                 date.genus = as.character(date_today))
        data_known <- plyr::rbind.fill(data_known,data.spec )
      }else{
        data_unknown <- rbind(data_unknown, data[i,])
      }
    }else{
      failed <- rbind(failed, data[i,])
    }
    setTxtProgressBar(pb, i)
  }; close(pb); rm(pb, i, data.spec)
  new_data <- list(data_known,data_unknown, failed )
  return(new_data)
}

####### Stats
calculcate_uniq_entries <- function(data, column){
  number <- data %>%pull({{column}}) %>% unique() %>%na.omit() %>%length()
  return(number)
}
calculcate_entries <- function(data, column){
  number <- data %>%pull({{column}})%>%na.omit() %>%length()
  return(number)
}

# Read and parse key-value pairs from ini-style file
read_simple_ini <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[!grepl("^\\s*#", lines)]  # remove comments
  lines <- lines[nzchar(lines)]           # remove empty lines
  
  kv_pairs <- strsplit(lines, "=", fixed = TRUE)
  kv_pairs <- lapply(kv_pairs, function(x) setNames(trimws(gsub('"', '', x[2])), trimws(x[1])))
  values <- do.call(c, kv_pairs)
  
  list2env(as.list(values), envir = .GlobalEnv)
}


calculate_refstats <- function(data, db_type){
  stats_df_DB <- data.frame()
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Sequence No.", value= nrow(data)))
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Order No.", value= calculcate_uniq_entries(data, "order")))
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Family No.", value=calculcate_uniq_entries(data, "Family")))
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Genus No.", value=calculcate_uniq_entries(data, "Genus")))
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Species No.", value=calculcate_uniq_entries(data, "Species")))
  
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Ciliate Sequence No.", value= nrow(data %>% filter(grepl("Ciliophora",V2)))))
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Ciliate Order No.", value= calculcate_uniq_entries(data %>% filter(grepl("Ciliophora",V2)), "order")))
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Ciliate Family No.", value=calculcate_uniq_entries(data %>% filter(grepl("Ciliophora",V2)), "Family")))
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Ciliate Genus No.", value=calculcate_uniq_entries(data %>% filter(grepl("Ciliophora",V2)), "Genus")))
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Ciliate Species No.", value=calculcate_uniq_entries(data %>% filter(grepl("Ciliophora",V2)), "Species")))
  
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Order Seq. No.", value= calculcate_entries(data, "order")))
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Family Seq. No.", value=calculcate_entries(data, "Family")))
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Genus Seq. No.", value=calculcate_entries(data, "Genus")))
  stats_df_DB <- rbind(stats_df_DB, data.frame(database= db_type, variable ="Species Seq. No.", value=calculcate_entries(data, "Species")))
  
  return(stats_df_DB)
  
}