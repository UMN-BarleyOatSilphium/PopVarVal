## PopVar Predictions
## 
## These are predictions of genetic variance using all available C1 lines (to be
## listed below).
## 
## This script can be run locally or on MSI
## 

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/PopVarVal/" 
source(file.path(repo_dir, "source_MSI.R"))

# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))

# Load genotypic and phenotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))
load(file.path(pheno_dir, "PVV_BLUE.RData"))


# Cross data
cross_list <- entry_list %>%
  filter(Group == "Experimental") %>%
  separate(Pedigree, c("parent1", "parent2"), sep = "/") %>%
  rename_all(str_to_lower)


# Subset the genos
genos_use <- s2_imputed_mat[c(tp_geno, pot_pars_geno),]

# Reformat the genotypic data
# Create data.frame for output
G_in <- as.data.frame(cbind( c("", row.names(genos_use)), rbind(colnames(genos_use), genos_use)) )

# Format the phenos
tp_prediction_phenos <- tp_prediction_BLUE %>% 
# phenos_use <- tp_relevant_BLUE %>% 
  spread(trait, value) %>%
  as.data.frame()

# Format the map
map_use <- snp_info %>% 
  select(marker = `rs#`, chrom, cM_pos) %>%
  as.data.frame()
  

# # Format the crossing block
# crossing_block <- cross_list %>%
#   distinct(family, parent1, parent2) %>%
#   as.data.frame() %>%
#   column_to_rownames("family")
# 
# ## Predict
# pred_out <- pop.predict(G.in = G_in, y.in = tp_prediction_phenos, map.in = map_use, 
#                         crossing.table = crossing_block, tail.p = 0.1, nInd = 150,
#                         min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
#                         nSim = 25, nCV.iter = 1, models = "rrBLUP", impute = "pass")
# 
# 
# # Convert to DF
# posthoc_family_pred <- tidy.popvar(pred_out) %>%
#   left_join(rownames_to_column(crossing_block, "family"), ., by = c("parent1" = "Par1", "parent2" = "Par2"))
# 
# 
# # Save the results
# save_file <- file.path(result_dir, "post_hoc_prediction_results.RData")
# save("posthoc_family_pred", file = save_file)


# Determine the number of cores
n_cores <- detectCores()


## Predict all possible bp families

# Create the crossing block of all pp_geno individuals
# Remove the crosses that were already predicted
all_par_crossing_block <- combn(x = pot_pars_geno, m = 2) %>%
  t() %>%
  as.data.frame() %>%
  rename(parent1 = V1, parent2 = V2) %>%
  setdiff(., crossing_block)


# Add cores to the crossing block and split by core
all_par_crossing_block_split <- all_par_crossing_block %>% 
  assign_cores(df = ., n_core = n_cores) %>%
  split(.$core)



# Apply the function by core
all_par_pred_out <- mclapply(X = all_par_crossing_block_split, FUN = function(core_df) {
    
    # Return predictions
    out <- pop.predict(G.in = G_in, y.in = phenos_use, map.in = map_use, 
                crossing.table = core_df, tail.p = 0.1, nInd = 150,
                min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
                nSim = 25, nCV.iter = 1, models = "rrBLUP", impute = "pass")
    
    tidy.popvar(out)

  }, mc.cores = n_cores)

all_family_pred <- bind_rows(all_par_pred_out)


## Save
save_file <- file.path(result_dir, "all_family_prediction_results.RData")
save("all_family_pred", file = save_file)



