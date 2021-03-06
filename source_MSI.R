## S2MET source for MSI
## 
## A script that automatically loads the data relevant for the S2MET project


# Load packages
packages <- c("dplyr", "tidyr", "tibble", "stringr", "readxl", "readr", "parallel",
              "PopVar", "purrr", "pbsim", "Matrix")

## Determine the package directory by assessing the version of R
vers <- paste0(strsplit(x = paste0(version$major, ".", version$minor), split = "\\.")[[1]][1:2], collapse = ".")
# Test directories
dir1 <- file.path("/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-unknown-linux-gnu-library", vers)
dir2 <- file.path("/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/", vers)
dir_list <- c(dir1, dir2)

package_dir <- dir_list[which(sapply(dir_list, dir.exists))]
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "/path/to/project/on/supercomputer"

geno_dir <-  "/path/to/genotype/data/on/supercomputer"
pheno_dir <- "/path/to/phenotypic/data/on/supercomputer"

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")


## Source the 'source.R' script from a define starting point
source_lines <- readLines(file.path(repo_dir, "source.R"))
source_lines_discard <- seq(which(grepl(pattern = "^# MSI", x = source_lines)))
source_lines_run <- source(textConnection(source_lines[-source_lines_discard]))

