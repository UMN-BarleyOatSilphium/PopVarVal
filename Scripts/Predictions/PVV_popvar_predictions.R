## PopVar Predictions
## 
## These are predictions of genetic variance using all available C1 lines (to be
## listed below).
## 
## This script can be run locally or on MSI
## 

# List of packages
packages <- c("dplyr", "tidyr", "tibble", "purrr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "PopVar", "purrrlyr", "boot")

# Set the directory of the R packages
package_dir <- NULL
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/PopVarVal/" 
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/PopVarVal/" 

## Geno and pheno directories
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/"
pheno_dir <- file.path(proj_dir, "Phenotype_Data/")

geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/GBS_Genos"
pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"

# Load genotypic and phenotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))
load(file.path(pheno_dir, "PVV_train_BLUEs.RData"))

# Other directories
result_dir <- file.path(proj_dir, "Results")
entry_dir <- file.path(proj_dir, "Plant_Material")


# Load project entries
entry_file <- file.path(entry_dir, "PopVarVal_project_entries.xlsx")
entry_list <- read_excel(entry_file)


# Cross data
cross_list <- read_excel(file.path(entry_dir, "PopVarVal_cross_info.xlsx")) %>%
  mutate(Cross_no = str_c("4", F1_Pot_No))

## Crosses that we made
post_hoc_cross <- cross_list %>%
  filter(!is.na(Summer16_Trial)) %>%
  select(cross = Cross_no, parent1 = Parent1, parent2 = Parent2, experiment = Experiment, 
         exp_no_lines = Harvested_No_lines)

# Extract the TP
tp <- entry_list %>%
  filter(Purpose == "S2TP") %>%
  pull(Line)

tp_geno <- intersect(tp, row.names(s2_imputed_mat))

# Extract the potential parents
potential_parents <- entry_list %>% 
  filter(Purpose == "Potential_parent") %>%
  pull(Line)

# Find the potential parents with genotypic information
pp_geno <- intersect(potential_parents, row.names(s2_imputed_mat))

# Find the used parents
pp_used <- post_hoc_cross %>% 
  select(parent1, parent2) %>% 
  unlist() %>% 
  unique()

# Intersect with available genos
pp_used_geno <- intersect(pp_used, row.names(s2_imputed_mat))

# Extract phenotypic data for the genotyped TP
tp_geno_BLUEs <- PVV_train_BLUEs_summ %>% 
  ungroup() %>% 
  filter(line_name %in% tp_geno)


# Subset the genos
genos_use <- s2_imputed_mat[c(tp_geno, pp_used_geno),]

# Reformat the genotypic data
# Create data.frame for output
G_in <- as.data.frame(cbind( c("", row.names(genos_use)), rbind(colnames(genos_use), genos_use)) )

# Format the phenos
phenos_use <- tp_geno_BLUEs %>% 
  spread(trait, value) %>%
  as.data.frame()

# Format the map
map_use <- snp_info %>% 
  select(marker = `rs#`, chrom, cM_pos) %>%
  as.data.frame()
  

# Format the crossing block
crossing_block <- post_hoc_cross %>% 
  filter(parent1 %in% pp_used_geno,
         parent2 %in% pp_used_geno) %>%
  select(parent1, parent2) %>%
  as.data.frame()

## Predict
pred_out <- pop.predict(G.in = G_in, y.in = phenos_use, map.in = map_use, 
                        crossing.table = crossing_block, tail.p = 0.1, nInd = 150,
                        min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
                        nSim = 25, nCV.iter = 1, models = "rrBLUP", impute = "pass")

# Define a function to convert the output list from PopVar to a combined data.frame
# with all traits
tidy_popvar <- function(x) {
  x$predictions %>% 
    map(as_data_frame) %>% 
    map(~mutate_all(., unlist)) %>%
    list(., names(.)) %>%
    pmap(.f = function(df, tr) {
      mutate(df, trait = tr) }) %>%
    bind_rows() %>%
    select(trait, names(.)) %>%
    mutate(trait = str_replace(trait, "_param.df", ""))
}

# Convert to DF
PVV_used_family_pred <- tidy_popvar(pred_out) %>%
  left_join(post_hoc_cross, ., by = c("parent1" = "Par1", "parent2" = "Par2"))


# Save the results
save_file <- file.path(result_dir, "PVV_post_hoc_prediction_results.RData")
save("PVV_used_family_pred", file = save_file)


# Determine the number of cores
n_cores <- detectCores()


## Predict all possible bp families
genos_use <- s2_imputed_mat[c(tp_geno, pp_geno),]
# Create data.frame for output
G_in <- as.data.frame(cbind( c("", row.names(genos_use)), rbind(colnames(genos_use), genos_use)) )

# Create the crossing block of all pp_geno individuals
# Remove the crosses that were already predicted
pp_crossing_block <- combn(x = pp_geno, m = 2) %>%
  t() %>%
  as.data.frame() %>%
  rename(parent1 = V1, parent2 = V2) %>%
  setdiff(., crossing_block)


# Add cores to the crossing block and split by core
pred_out_list <- pp_crossing_block %>% 
  mutate(core = sort(rep(seq(n_cores), length.out = nrow(.)))) %>% 
  split(.$core) %>%
  # Apply the function by core
  mclapply(X = ., FUN = function(core_cb) {
    
    # Return predictions
    pop.predict(G.in = G_in, y.in = phenos_use, map.in = map_use, 
                crossing.table = core_cb, tail.p = 0.1, nInd = 150,
                min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
                nSim = 25, nCV.iter = 1, models = "rrBLUP", impute = "pass")

  }, mc.cores = n_cores)


PVV_family_pred <- pred_out_list %>% 
  map_df(tidy_popvar) %>%
  rename(parent1 = Par1, parent2 = Par2)

# Combine data
PVV_all_pred <- bind_rows(PVV_used_family_pred, PVV_family_pred)

## Save
save_file <- file.path(result_dir, "PVV_prediction_results.RData")
save("PVV_all_pred", file = save_file)



