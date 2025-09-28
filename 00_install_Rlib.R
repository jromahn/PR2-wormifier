options(repos = c(CRAN = "https://cloud.r-project.org/"))

install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("Biostrings")
devtools::install_github("pr2database/pr2database")
install.packages("stringi")
