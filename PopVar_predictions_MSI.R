# Use PopVar to simulate crosses between the S2C1 and predict population mean,
## population V_G, and superior progeny mean

# Are we working on MSI?
MSI = TRUE

if (MSI) {
  geno.file <- "/panfs/roc/groups/6/smithkp/neyhartj/Genomic_Selection/Data/GBS_Genos/S2_genos_hmp.RData"
  pheno.file <- "/panfs/roc/groups/6/smithkp/neyhartj/Genomic_Selection/Data/Phenos/PopVarVal_BLUE.RData"
  
  setwd("/panfs/roc/groups/6/smithkp/neyhartj/Genomic_Selection/PopVarVal")
  
  package.dir <- "/panfs/roc/groups/6/smithkp/neyhartj/R/x86_64-pc-linux-gnu-library/3.3/"
  
  # Load packages
  library(tidyverse, quietly = T, package.dir)
  library(stringr, quietly = T, package.dir)
  library(rrBLUP, quietly = T, package.dir)
  library(PopVar, quietly = T, package.dir)
  library(parallel, quietly = T, package.dir)
  
} else {
  geno.file <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/Genotypic Data/GBS Genotype Data/S2_genos_hmp.RData"
  pheno.file <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/Phenotypic Data/Final/Master Phenotypes/PopVarVal_BLUE.RData"
  
  # Load packages
  library(tidyverse)
  library(stringr)
  library(rrBLUP)
  library(PopVar)
  library(parallel)
  
}

# Read in data
load(geno.file)
load(pheno.file)

# Detect cores
n.cores <- detectCores()

# Training and prediction line names
tp <- names(s2_imputed_genos) %>%
  str_subset("^[0-9]{2}")

vp <- names(s2_imputed_genos) %>%
  str_subset("^2MS")

# Generate an ordered crossing table
crossing.table <- combn(vp, 2) %>% t()
colnames(crossing.table) <- c("Par1", "Par2")

# Trim the phenotype data.frame
y.in <- PopVarVal.tidy.BLUE %>%
  filter(line_name %in% tp) %>%
  spread(trait, value) %>%
  as.data.frame()
  
# Edit the snp info for the map
map.in <- s2_imputed_genos %>%
  select(`rs#`, chrom, cM_pos) %>%
  rename(mkr = `rs#`) %>%
  as.data.frame()
  
# Pre-allocate the G.in data.frame
G.in <- matrix(NA, nrow = length(tp) + length(vp) + 1, ncol = nrow(map.in) + 1) %>%
  as.data.frame()

# Add entry names
G.in[,1] <- c("mkr", tp, vp)
# Add marker names
G.in[1,-1] <- map.in$mkr
# Add genotypic data
G.in[-1, -1] <- t(s2_imputed_genos[,-1:-5]) %>% as.data.frame()

# Split the crossing table into parallelizable chunks
crossing.table.list <- seq(nrow(crossing.table)) %>%
  split(f = cut(seq(nrow(crossing.table)), breaks = n.cores))

# Execute the simulations in parallel

popvar_out <- mclapply(X = crossing.table.list, FUN = function(cross.table) {
  
  cross.table <- crossing.table[cross.table,]
  
  pop.predict(G.in = G.in, y.in = y.in, map.in = map.in, 
              crossing.table = cross.table, tail.p = 0.1, nInd = 150, 
              min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, 
              impute = "pass", nSim = 25, nCV.iter = 0, models = "rrBLUP")
  
}, mc.cores = n.cores)


# Save the output
save("popvar_out", file = "Output/PopVar_predictions_20161206.RData")
