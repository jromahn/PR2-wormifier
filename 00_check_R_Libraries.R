################################################################################
## Script: 00_Package_Versions.R
## Aim:    Load all packages used across the pipeline and report their versions.
## Usage:  Run this script independently to document the software environment.
################################################################################

# ── All packages used across scripts 1–10 ─────────────────────────────────────

packages <- c(
  # Data wrangling
  "dplyr",
  "tidyr",
  "stringr",
  "stringi",
  "plyr",
  "purrr",
  
  # Bioinformatics / sequence handling
  "Biostrings",
  "seqinr",
  "phylotools",
  "pr2database",
  
  # Taxonomy databases
  "worrms",
  "taxonomizr",
  
  # NCBI access
  "rentrez",
  
  # Web / API for Algaebase
  "jsonlite",
  "curl",
  
  
  ##vizualize & output creation
  "ggplot2",
  "tibble",
  "gapminder",
  "treemapify",
  "ggpubr",
  "ghibli",
  "openxlsx", 
  "legendry",
  "ggh4x"
)

# ── Load each package and collect version info ─────────────────────────────────

results <- data.frame(
  Package = character(),
  Version = character(),
  stringsAsFactors = FALSE
)

for (pkg in packages) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  ver <- as.character(packageVersion(pkg))
  results <- rbind(results, data.frame(
    Package = pkg,
    Version = ver,
    stringsAsFactors = FALSE
  ))
}

# ── R version ─────────────────────────────────────────────────────────────────

cat("\n")
cat("=======================================================\n")
cat("  R version\n")
cat("=======================================================\n")
print(R.version.string)

# ── Print results table ────────────────────────────────────────────────────────

cat("\n")
cat("=======================================================\n")
cat("  Package versions\n")
cat("=======================================================\n")
print(results, row.names = FALSE)

# ── Save citations into BibLatex ────────────────────────────────────────────────────────

sink("00_Package_Citations.bib")
for (pkg in packages) {
  tryCatch({
    cit <- citation(pkg)
    cat(toBibtex(cit), "\n\n")
  }, error = function(e) {
    message("No citation available for: ", pkg)
  })
}
sink()


# ── Optionally write to file ───────────────────────────────────────────────────

write.table(
  results,
  file      = "00_Package_Versions.tsv",
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

cat("\nVersion table saved to: 00_Package_Versions.tsv\n")
