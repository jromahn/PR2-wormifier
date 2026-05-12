options(repos = c(CRAN = "https://cloud.r-project.org/"))

install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("Biostrings")

# data handling
install.packages("stringi")
install.packages("stringr")
install.packages("tidyr")
install.packages("dplyr")


## install blaster from archive (dependency for pr2database) since it is not supported anymore
url <- "https://cran.r-project.org/src/contrib/Archive/blaster/blaster_1.0.7.tar.gz"
pkgFile <- "blaster_1.0.7.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
file.remove(pkgFile)

devtools::install_github("pr2database/pr2database")

## for downloading from algaebase
install.packages("jsonlite")
install.packages("curl")

#worms taxonomy
install.packages("worrms")

# download from ncbi
install.packages("rentrez")

#ncbi assembly to taxonomy
install.packages("taxonomizr")



