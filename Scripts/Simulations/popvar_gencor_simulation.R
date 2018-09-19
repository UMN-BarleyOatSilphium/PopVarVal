## PopVarVal simulations of prediction accuracy for the genetic
## correlation between two traits
## 
## 
## Author: Jeff Neyhart
## Last modified: September 2, 2018
## 

# Run the source script
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/PopVarVal/"
source(file.path(repo_dir, "source_MSI.R"))

# Load the two-row simulation genotypes
load(file.path(geno_dir, "s2_cap_simulation_data.RData"))




# # Run on a local machine
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))
# # Additional libraries
# library(pbsim)
# library(PopVar)
# 
# # Load the two-row simulation genotypes
# load(file.path(gdrive_dir, "BarleyLab/Projects/SideProjects/Resources/s2_cap_simulation_data.RData"))


# Remove monomorphic SNPs
s2_cap_genos <- s2_cap_genos[,!colMeans(s2_cap_genos) %in% c(0, 2)]
s2_snp_info <- subset(s2_snp_info, rs %in% colnames(s2_cap_genos))





# Number cores
n_cores <- 16 # Local machine for demo
n_cores <- detectCores()


#### Simulation to test different genetic architectures and training population parameters


## Fixed parameters
sim_pop_size <- 150
n_iter <- 25
n_env <- 3
n_rep <- 1
n_crosses <- 50
k_sp <- 1.76

## Outline the parameters to perturb
trait1_h2_list <- trait2_h2_list <- c(0.3, 0.6, 1)
nQTL_list <- c(30, 100)
tp_size_list <- seq(150, 600, by = 150)
gencor_list <- c(-0.5, 0, 0.5)
probcor_list <- data_frame(arch = c("pleio", "close_link", "loose_link"),
                           input = list(cbind(0, 1), cbind(5, 1), cbind(30, 1) ))

# Create a data.frame of parameters
param_df <- crossing(trait1_h2 = trait1_h2_list, trait2_h2 = trait2_h2_list, nQTL = nQTL_list, tp_size = tp_size_list,
                     gencor = gencor_list, probcor = probcor_list, iter = seq(n_iter))

map_sim <- s2_snp_info %>%
  split(.$chrom) %>%
  map(~setNames(.$cM_pos, .$rs)) %>%
  map(~structure(., class = "A"))  %>%
  structure(., class = "map") %>%
  # Jitter
  qtl::jittermap(.) %>%
  `names<-`(., seq_along(.))

# Create the base genome - this is fixed for all simulations
genome <- sim_genome(map = map_sim)


# Split the parameter df
param_df_split <- param_df %>%
  assign_cores(n_cores) %>%
  split(.$core)

# Parallelize
simulation_out <- mclapply(X = param_df_split, FUN = function(core_df) {

  # ## For local machine
  # i <- 1
  # core_df <- param_df_split[[i]]
  # i <- 65
  # ##

  # Create a results list
  results_out <- vector("list", nrow(core_df))

  # Iterate over the rows of the param_df
  for (i in seq_along(results_out)) {

    trait1_h2 <- core_df$trait1_h2[i]
    trait2_h2 <- core_df$trait2_h2[i]
    L <- core_df$nQTL[i]
    maxL <- max(param_df$nQTL)
    tp_size <- core_df$tp_size[i]
    gencor <- core_df$gencor[i]
    probcor <- core_df$input[[i]]

    # Simulate QTL
    qtl_model <- replicate(n = 2, matrix(NA, ncol = 4, nrow = L), simplify = FALSE)
    genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "geometric", max.qtl = maxL,
                                   corr = gencor, prob.corr = probcor)

    ## Adjust the genetic architecture - only if pleiotropy is not present
    if (probcor[1] != 0) {
      genome1 <- adj_multi_gen_model(genome = genome1, geno = s2_cap_genos, gencor = gencor)
    }


    ## Create the TP by selecting on the larger base population
    ## Hopefully this will enforce the desired genetic correlation
    # Create the larger base population
    tp1 <- create_pop(genome = genome1, geno = s2_cap_genos) %>% 
      # Phenotype the base population
      sim_phenoval(pop = ., h2 = c(trait1_h2, trait2_h2), n.env = n_env, n.rep = n_rep) %>%
      # Select using equal weight
      select_pop(pop = ., intensity = tp_size, index = c(1, 1), type = "phenotypic")
    
    # Measure the genetic correlation in the TP
    tp_cor <- cor(tp1$geno_val[,-1])[1,2]
    # Measure the phenotypic correlation in the TP
    tp_pheno_cor <- cor(tp1$pheno_val$pheno_mean[,-1])[1,2]
    par_pop <- tp1

    # ##### Create a set of cycle 1 progeny from which to predict crosses
    # ## Select the best tp individuals for both traits
    # tp_select <- select_pop(pop = tp1, intensity = 0.1, index = c(1, 1), type = "phenotypic")
    # # Randomly create crosses from the TP individuals
    # crossing_block <- sim_crossing_block(parents = indnames(tp_select), n.crosses = 40)
    # # Pedigree to accompany the crosses
    # ped <- sim_pedigree(n.ind = 25, n.selfgen = Inf)
    # 
    # # Make theses crosses
    # par_pop <- sim_family_cb(genome = genome1, pedigree = ped, founder.pop = tp_select, crossing.block = crossing_block)
    # #####

    # Randomly create crosses from the cycle1 individuals
    crossing_block <- sim_crossing_block(parents = indnames(par_pop), n.crosses = n_crosses)
    # Pedigree to accompany the crosses
    ped <- sim_pedigree(n.ind = sim_pop_size, n.selfgen = Inf)


    ## Predict genetic variance and correlation
    pred_out <- pred_genvar(genome = genome1, pedigree = ped, training.pop = tp1, founder.pop = par_pop,
                            crossing.block = crossing_block) %>%
      mutate(pred_musp = pred_mu + (k_sp * sqrt(pred_varG)))



    # ## Check the predictions of correlation versus PopVar
    # # Convert items for PopVar
    # # Pheno
    # pheno_use <- tp1$pheno_val$pheno_mean
    # geno_use <- genotype(genome1, combine_pop(list(tp1, par_pop)))
    # map_use <- genome2$gen_model$trait1 %>%
    #   select(qtl_name, chr, pos)
    #
    #
    # # # Pass to PopVar
    # # Convert genotypes into something useable for PopVar
    # geno_use <- as.data.frame(cbind( c("", row.names(geno_use)), rbind(colnames(geno_use), geno_use)) )
    #
    # pred_out_pv <- pop.predict(G.in = geno_use, y.in = pheno_use, map.in = map_use,
    #                           crossing.table = crossing_block, tail.p = 0.1, nInd = sim_pop_size,
    #                           min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
    #                           nSim = n_pops, nCV.iter = 1, models = "rrBLUP", impute = "pass")
    #
    # # Tidy
    # tidy_pred_out <- pred_out_pv$predictions %>%
    #   map(as.data.frame) %>%
    #   map(~mutate_all(., unlist) %>% rename_at(vars(contains("cor")), ~"pred.corG")) %>%
    #   list(., names(.)) %>%
    #   pmap_df(~mutate(.x, trait = str_extract(.y, "trait[0-9]{1}"))) %>%
    #   select(parent1 = Par1, parent2 = Par2, trait, pred.varG, pred.corG)
    #
    # ## Correlate predictions of genetic variance and correlation
    # pred_out %>%
    #   left_join(., tidy_pred_out) %>%
    #   distinct(parent1, parent2, pred_corG, pred.corG) %>%
    #   summarize(acc = cor(pred_corG, pred_corG))
    #
    # pred_out %>%
    #   left_join(., tidy_pred_out) %>%
    #   group_by(trait) %>%
    #   summarize(acc = cor(pred_varG, pred_varG))
    #
    # ## Completely accurate!




    ## Calculate the expected genetic variance in these populations
    expected_var <- calc_exp_genvar(genome = genome1, pedigree = ped, founder.pop = par_pop, crossing.block = crossing_block) %>%
      mutate(exp_musp = exp_mu + (k_sp * sqrt(exp_varG)))
    
    ## Combine the expected and predicted results
    expected_predicted <- full_join(
      x = pred_out %>% rename_all(~str_replace(., "pred_", "")) %>% gather(parameter, prediction, -parent1:-trait),
      y = expected_var %>% rename_all(~str_replace(., "exp_", "")) %>% gather(parameter, expectation, -parent1:-trait),
      by = c("parent1", "parent2", "trait", "parameter")
    ) %>% 
      # If looking at the correlation, change the trait to simply trait1
      filter(!(trait == "trait2" & parameter == "corG"))

    # Summarize
    results_summ <- expected_predicted %>% 
      group_by(trait, parameter) %>% 
      summarize(accuracy = cor(prediction, expectation), 
                bias = mean(prediction - expectation) / mean(expectation))

    ## Add the accuracy results to the results list
    results_out[[i]] <- list(
      summary = results_summ,
      other = data.frame(variable = c("tp_gencor", "tp_phencor"), value = c(tp_cor, tp_pheno_cor))
    )

  }

  # Add the results to the core_df, remove core
  core_df %>%
    mutate(results = results_out) %>%
    select(-core)

}, mc.cores = n_cores)

# Bind and save
popvar_prediction_simulation_out <- bind_rows(simulation_out)

# Save
save_file <- file.path(result_dir, "popvar_gencor_simulation_prediction_results.RData")
save("popvar_prediction_simulation_out", file = save_file)







# #### Simulation to test different degrees of linkage and pleiotropy
# 
# 
# ## Fixed parameters
# sim_pop_size <- 150
# n_iter <- 10
# n_env <- 3
# n_rep <- 1
# n_crosses <- 50
# k_sp <- 1.76
# 
# ## Outline the parameters to perturb
# trait1_h2_list <- trait2_h2_list <- 0.5
# nQTL_list <- c(100)
# tp_size_list <- 300
# gencor_list <- c(-0.5, 0, 0.5)
# 
# ## Percentage of pleiotropy versus degree of linkage
# pPleio <- seq(0, 0.9, by = 0.1)
# dLinkage <- seq(5, 50, by = 5)
# 
# probcor_list <- crossing(pPleio, dLinkage) %>%
#   add_row(pPleio = 1, dLinkage = 0) %>%
#   mutate(pLinkage = 1 - pPleio) %>%
#   pmap(~rbind(cbind(0, ..1), cbind(..2, ..3))) %>%
#   # Remove any rows with 0 probability
#   map(~.[.[,2] != 0,,drop = FALSE])
# 
# 
# probcor_df <- data_frame(probor = probcor_list)
# 
# 
# # Create a data.frame of parameters
# param_df <- crossing(trait1_h2 = trait1_h2_list, trait2_h2 = trait2_h2_list, nQTL = nQTL_list, tp_size = tp_size_list,
#                      gencor = gencor_list, iter = seq(n_iter), probcor = probcor_df)
# 
# map_sim <- s2_snp_info %>%
#   split(.$chrom) %>%
#   map(~setNames(.$cM_pos, .$rs)) %>%
#   map(~structure(., class = "A"))  %>%
#   structure(., class = "map") %>%
#   # Jitter
#   qtl::jittermap(.) %>%
#   `names<-`(., seq_along(.))
# 
# # Create the base genome - this is fixed for all simulations
# genome <- sim_genome(map = map_sim)
# 
# 
# # Split the parameter df
# param_df_split <- param_df %>%
#   assign_cores(n_cores) %>%
#   split(.$core)
# 
# # Parallelize
# simulation_out <- mclapply(X = param_df_split, FUN = function(core_df) {
# 
#   # ## For local machine
#   # i <- 1
#   # core_df <- param_df_split[[i]]
#   # ##
# 
#   # Create a results list
#   results_out <- vector("list", nrow(core_df))
# 
#   # Iterate over the rows of the param_df
#   for (i in seq_along(results_out)) {
# 
#     trait1_h2 <- core_df$trait1_h2[i]
#     trait2_h2 <- core_df$trait2_h2[i]
#     L <- core_df$nQTL[i]
#     maxL <- max(param_df$nQTL)
#     tp_size <- core_df$tp_size[i]
#     gencor <- core_df$gencor[i]
#     probcor <- core_df$probor[[i]]
# 
# 
#     # Simulate QTL
#     qtl_model <- replicate(n = 2, matrix(NA, ncol = 4, nrow = L), simplify = FALSE)
#     genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "geometric", max.qtl = maxL,
#                                    corr = gencor, prob.corr = probcor)
# 
#     ## Adjust the genetic architecture - only if pleiotropy is not present
#     if (all(probcor[,1] == 0)) {
#       genome1 <- adj_multi_gen_model(genome = genome1, geno = s2_cap_genos, gencor = gencor)
#     }
# 
#     ## Create the TP by selecting on the larger base population
#     ## Hopefully this will enforce the desired genetic correlation
#     # Create the larger base population
#     tp1 <- create_pop(genome = genome1, geno = s2_cap_genos) %>% 
#       # Phenotype the base population
#       sim_phenoval(pop = ., h2 = c(trait1_h2, trait2_h2), n.env = n_env, n.rep = n_rep) %>%
#       # Select using equal weight
#       select_pop(pop = ., intensity = tp_size, index = c(1, 1), type = "phenotypic")
#     
#     # Measure the genetic correlation in the TP
#     tp_cor <- cor(tp1$geno_val[,-1])[1,2]
#     # Measure the phenotypic correlation in the TP
#     tp_pheno_cor <- cor(tp1$pheno_val$pheno_mean[,-1])[1,2]
#     par_pop <- tp1
# 
#     # ##### Create a set of cycle 1 progeny from which to predict crosses
#     # ## Select the best tp individuals for both traits
#     # tp_select <- select_pop(pop = tp1, intensity = 0.1, index = c(1, 1), type = "phenotypic")
#     # # Randomly create crosses from the TP individuals
#     # crossing_block <- sim_crossing_block(parents = indnames(tp_select), n.crosses = 40)
#     # # Pedigree to accompany the crosses
#     # ped <- sim_pedigree(n.ind = 25, n.selfgen = Inf)
#     #
#     # # Make theses crosses
#     # par_pop <- sim_family_cb(genome = genome1, pedigree = ped, founder.pop = tp_select, crossing.block = crossing_block)
#     # #####
# 
#     # Randomly create crosses from the cycle1 individuals
#     crossing_block <- sim_crossing_block(parents = indnames(par_pop), n.crosses = n_crosses)
#     # Pedigree to accompany the crosses
#     ped <- sim_pedigree(n.ind = sim_pop_size, n.selfgen = Inf)
# 
# 
# 
#     ## Predict genetic variance and correlation
#     pred_out <- pred_genvar(genome = genome1, pedigree = ped, training.pop = tp1, founder.pop = par_pop,
#                             crossing.block = crossing_block) %>%
#       mutate(pred_musp = pred_mu + (k_sp * sqrt(pred_varG)))
# 
# 
# 
# 
#     ## Calculate the expected genetic variance in these populations
#     expected_var <- calc_exp_genvar(genome = genome1, pedigree = ped, founder.pop = par_pop, crossing.block = crossing_block) %>%
#       mutate(exp_musp = exp_mu + (k_sp * sqrt(exp_varG)))
# 
# 
#     ## Summarize the expected and predicted values
#     expected_predicted <- full_join(pred_out, expected_var, by = c("parent1", "parent2", "trait"))
# 
#     # Summarize variance and means
#     expected_predicted_mean_var <- expected_predicted %>%
#       select(-contains("corG")) %>%
#       gather(parameter, value, pred_mu:exp_musp) %>%
#       separate(parameter, c("type", "parameter"), sep = "_") %>%
#       spread(type, value)
# 
#     # Summarize correlations
#     expected_predicted_cor <- expected_predicted %>%
#       select(parent1, parent2, contains("corG")) %>%
#       distinct() %>%
#       gather(parameter, value, pred_corG, exp_corG) %>%
#       separate(parameter, c("type", "parameter"), sep = "_") %>%
#       spread(type, value)
# 
#     # Combine
#     expected_predicted1 <- bind_rows(expected_predicted_mean_var, expected_predicted_cor)
# 
#     # Summarize
#     results_summ <- expected_predicted1 %>%
#       group_by(trait, parameter) %>%
#       nest(exp, pred) %>%
#       mutate(accuracy = map_dbl(data, ~cor(.$exp, .$pred, use = "complete.obs")),
#              bias = map_dbl(data, ~mean(.$pred - .$exp, na.rm = T) / mean(.$exp, na.rm = T)),
#              nMissing = map_dbl(data, ~sum(is.na(.$exp))))
# 
#     ## Add the accuracy results to the results list
#     results_out[[i]] <- list(
#       summary = results_summ,
#       other = data.frame(variable = c("tp_gencor", "tp_phencor"), value = c(tp_cor, tp_pheno_cor))
#     )
# 
#   }
# 
#   # Add the results to the core_df, remove core
#   core_df %>%
#     mutate(results = results_out) %>%
#     select(-core)
# 
# }, mc.cores = n_cores)
# 
# # Bind and save
# popvar_corG_space_simulation_out <- bind_rows(simulation_out)
# 
# # Save
# save_file <- file.path(result_dir, "popvar_gencor_space_simulation_results.RData")
# save("popvar_corG_space_simulation_out", file = save_file)









