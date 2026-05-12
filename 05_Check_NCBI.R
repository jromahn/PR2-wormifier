rm(list=ls())

source("00_Function_Library.R") # read file with functions
read_simple_ini("00_login_data.ini")
#####

library(rentrez)
library(stringr)

##################################################################################
# Aim: Download species from NCBI to create a Reference Database.
# Function:  Searches NCBI for nuclear 18S sequences of missing or related species, filters results, and downloads relevant sequences. 
##### The created folder and metadata will contain the data when the downloading was started.
######
# 1.) Species-level NCBI search
# 2.) Genus-level fallback search 
##################################################################################


 
### flags  ────────────────────────────────────────────────────────────────

FLAG_search_spec <- TRUE   # Step 1 – search by species
FLAG_search_gen  <- TRUE   # Step 2 – fall back to genus

### NCBI search parameters  ────────────────────────────────────────────────────────────────

search_suffix <- paste0(
  " AND (18S OR small ribosomal)",
  " NOT (chloroplast[filter] OR mitochondrion[filter])",
  " AND (350[SLEN] : 5000[SLEN])"
)
ncbi_db <- "nucleotide"

## input file paths ────────────────────────────────────────────────────────────────
Input_folder  <- file.path(path_to_output, "01_intermediate_results/")

file_missing            <- paste0(Input_folder, "4.1_Missing_Species.csv")
file_missing_but_genus  <- paste0(Input_folder, "4.2_Present_Species.csv")
file_missing_but_genus_new <- paste0(Input_folder,"4.2_Present_Species_completed_withNCBI.csv")
file.copy(file_missing_but_genus, file_missing_but_genus_new, overwrite = TRUE)

# ── entrez configuration ─────────────────────────────────────────────────────────────
TIMEOUT      <- 30          # seconds (used as httr / curl option where possible)

# Read email from config
set_entrez_key(NULL)        # clear any stale API key
options(ENTREZ_KEY = NULL)
entrez_email <- Entrez_email   # rentrez uses this automatically

# ── create new file paths ────────────────────────────────────────────────────────────────

#date <- format(Sys.Date(), "%d_%b_%Y")
date <- "20_Apr_2026"


NCBI_folder         <- file.path(path_to_output,paste0("02_NCBI_", date, "/"))
NCBI_results_folder <- file.path(path_to_output,paste0("02_NCBI_", date, "_results/"))
output_prefix       <- "05_Missing_PR2"

path_creation(c(NCBI_folder,NCBI_results_folder ))

file_still_missing   <- paste0(NCBI_results_folder, output_prefix,"_Species_and_missing_NCBI__",  date, ".csv")
file_missing_genera  <- paste0(NCBI_results_folder, output_prefix,"_incomplete_genera_despite_NCBI__", date, ".csv")
file_missing_genera2 <- paste0(NCBI_results_folder, output_prefix,"_MISSING_genera_despite_NCBI__", date, ".csv")
file_accession_info  <- paste0(NCBI_results_folder, output_prefix,"_Downloaded_species_NCBI_info__", date, ".tsv")
file_LOG             <- paste0(NCBI_results_folder, output_prefix,"LOG__", date, ".txt")
file_LOG2            <- paste0(NCBI_results_folder, output_prefix,"LOG__", date, "_successfull_genera.txt")

##############################################################################
# STEP 1 – Search and download by species
##############################################################################

if (FLAG_search_spec) {
  
  message("Start downloading species")
  
  con_missing   <- file(file_missing,            "r")
  con_still     <- file(file_still_missing,      "w")
  con_found     <- file(file_missing_but_genus_new, "a")
  con_log       <- file(file_LOG,                "a")
  con_assec     <- file(file_accession_info,     "w")
  
  lines <- readLines(con_missing)
  lines <- lines[-1]                    # drop header
  
  #line <- lines[2]
  for (line in lines) {
    line <- trimws(line)
    line <- gsub('"', '', line)
    parts   <- strsplit(line, ",")[[1]]
    genus   <- parts[1]
    species <- parts[2]
    
    message(species)
    new_search <- paste0(species, "[ORGANISM]", search_suffix)
    writeLines(new_search, con_log)
    
    success <- FALSE
    for (attempt in seq_len(3)) {
      message("  attempt ", attempt)
      result <- tryCatch({
        
        #tax <- entrez_search(db = "taxonomy", term = "Chironomus aprilinus")
        #tax$ids
        
        #new_search <- "Chironomus aprilinus[ORGANISM:noexp] AND (18S OR small ribosomal) NOT (chloroplast[filter] OR mitochondrion[filter]) AND (350[SLEN] : 5000[SLEN])"
        # --- search ---
        search_res <- entrez_search(db = ncbi_db, term = new_search,retmax = 10, use_history = FALSE)
        cat(search_res$QueryTranslation)
        
        search_res <- entrez_search(db = ncbi_db, term = new_search, retmax = 10, use_history = FALSE)
        summaries <- entrez_summary(db = ncbi_db, id = search_res$ids)
        
        # extract organism info
        sapply(summaries, function(x) x$organism)
        id_list <- search_res$ids
        
        if (length(id_list) == 0) {
          writeLines(line, con_still)
          success <- TRUE
          break
        }
        
        message("  downloading...")
        out_file <- paste0(NCBI_folder,gsub(" ", "_", species), "__NCBI.fasta")
        
        # --- fetch FASTA ---
        fasta_text <- entrez_fetch(db = ncbi_db, id = id_list,rettype = "fasta", retmode = "text")
        writeLines(fasta_text, out_file)
        
        # --- fetch GenBank ---
        if (file.exists(out_file) && file.info(out_file)$size > 0) {
          writeLines(paste0(genus, ",", species), con_found)
          
          gb_text <- entrez_fetch(db = ncbi_db, id = id_list,rettype = "genbank", retmode = "text")
          
          info_df <- parse_genbank_info(gb_text, species)
          info_df$file <- out_file
          write.table(info_df, con_assec,sep = "\t", col.names = FALSE, row.names = FALSE,quote = FALSE, append = TRUE)
        }
        
        success <- TRUE
        "ok"
        
      }, error = function(e) {
        msg <- paste0(line, "\tproblem: ", class(e)[1], " - ", conditionMessage(e))
        message("  ERROR: ", msg)
        writeLines(msg, con_log)
        Sys.sleep(2)
        e
      })
      
      if (isTRUE(success)) break
    }
  }
  
  close(con_missing)
  close(con_still)
  close(con_found)
  close(con_log)
  close(con_assec)
}



##############################################################################
# STEP 2 – Fall back to genus-level search
##############################################################################

if (FLAG_search_gen) {
  message("Start downloading species from the same genus")
  
  # -- 2a) Build genera_dict: genus -> vector of already-known species --------
  genera_dict <- list()
  lines_genus <- readLines(file_missing_but_genus_new)
  lines_genus <- lines_genus[-1]         # drop header
  
  for (line in lines_genus) {
    line  <- trimws(gsub('"', '', line))
    parts <- strsplit(line, ",")[[1]]
    g <- parts[1]; sp <- parts[2]
    if (g %in% names(genera_dict)) {
      genera_dict[[g]] <- c(genera_dict[[g]], sp)
    } else {
      genera_dict[[g]] <- sp
    }
  }
  
  # -- 2b) Build list of genera that still need searching --------------------
  genera_search <- character(0)
  if (file.exists(file_still_missing)) {
    lines_still <- readLines(file_still_missing)
    for (line in lines_still) {
      line  <- trimws(gsub('"', '', line))
      parts <- strsplit(line, ",")[[1]]
      g <- parts[1]
      if (!(g %in% genera_search)) genera_search <- c(genera_search, g)
    }
  }
  
  # -- 2c) Search & download by genus ----------------------------------------
  con_incomp <- file(file_missing_genera,  "w")
  con_miss   <- file(file_missing_genera2, "w")
  con_log    <- file(file_LOG,             "a")
  con_assec  <- file(file_accession_info,  "a")
  con_log2   <- file(file_LOG2,            "w")
  
  genera_search[genera_search %in% names(genera_dict)]
  genus <- "Acartia"
  for (genus in genera_search) {
    message("\n", genus)
    
    
    # Build search term excluding already-downloaded species
    if (genus %in% names(genera_dict)) {
      exclusions <- paste0(genera_dict[[genus]], "[ORGANISM]", collapse = " OR ")
      new_search <- paste0(genus, "[ORGANISM] NOT (", exclusions, ") ", search_suffix)
    } else {
      new_search <- paste0(genus, "[ORGANISM] ", search_suffix)
    }
    
    success <- FALSE
    for (attempt in seq_len(3)) {
      result <- tryCatch({
        
        # --- search ---
        search_res <- entrez_search(db = ncbi_db, term = new_search,
                                    retmax = 50, use_history = FALSE)
        id_list <- search_res$ids
        
        if (length(id_list) == 0) {
          message("Genus unknown: ", genus)
          if (genus %in% names(genera_dict)) {
            writeLines(paste0(genus, ",",paste(genera_dict[[genus]], collapse = "&"),new_search), con_miss)
          } else {
            writeLines(paste0(genus, ",", new_search), con_incomp)
          }
        } else {
          writeLines(paste0(genus, ",", new_search), con_log2)
          
          out_file <- paste0(NCBI_folder, genus, "__NCBI.fasta")
          
          # --- fetch FASTA ---
          fasta_text <- entrez_fetch(db = ncbi_db, id = id_list, rettype = "fasta", retmode = "text")
          writeLines(fasta_text, out_file)
          
          # --- fetch GenBank and write accession info ---
          if (file.exists(out_file) && file.info(out_file)$size > 0) {
            gb_text <- entrez_fetch(db = ncbi_db, id = id_list,rettype = "genbank", retmode = "text")
            
            info_df <- parse_genbank_info(gb_text, genus)
            info_df$file <- out_file
            write.table(info_df, con_assec,
                        sep = "\t", col.names = FALSE, row.names = FALSE,
                        quote = FALSE, append = TRUE)
          }
        }
        
        success <- TRUE
        "ok"
        
      }, error = function(e) {
        message("Something went wrong: ", conditionMessage(e))
        writeLines(paste0(genus, ",", new_search), con_log)
        Sys.sleep(2)
        e
      })
      
      if (isTRUE(success)) break
    }
  }
  
  close(con_incomp)
  close(con_miss)
  close(con_log)
  close(con_assec)
  close(con_log2)
}

message("Done.")
