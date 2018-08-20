## This is the script to run the simulations for 'PVV_data_analysis_simulations'

# Load packages
library(tidyverse)
library(stringr)
library(PopVar)
library(qtl)
library(pbsim)
library(lme4)
library(gws)
library(sommer)
library(car)

# Load some data
data("s6_cap_haploid")
data("s6_snp_info")

proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/PopVarVal"

sim_dir <- file.path(proj_dir, "Simulations")


# Set some parameters
n_qtl <- 100
h2 <- 0.5 
n_crosses <- 15
n_ind_per_cross <- 20
n_env <- 4
n_rep <- 1
n_check_rep <- 3

n_iter <- 25

# Readjust the map
map_true <- s6_snp_info %>%
  dplyr::select(-alleles) %>%
  as.data.frame() %>%
  column_to_rownames("rs") %>%
  table_to_map()

# Create the genome
genome <- sim_genome(map = map_true, type = "hypred")

# Pedigree
ped <- sim_pedigree(n.ind = n_ind_per_cross, n.selfgen = 3)


## Loop starts here
simulation_results <- replicate(n = n_iter, expr = {


  # Simulate genetic architecture
  qtl.model <- matrix(NA, nrow = n_qtl, ncol = 4)
  genome <- sim_gen_model(genome = genome, qtl.model = qtl.model, add.dist = "geometric")
  
  # Extract the map (minus QTL)
  map_use <- genome %>% 
    find_markerpos(marker = markernames(.)) %>%
    rownames_to_column("rs")
  
  # Create the training population and phenotype
  training_pop <- create_pop(genome = genome, geno = s6_cap_haploid) %>%
    sim_phenoval(h2 = h2, n.env = n_env, n.rep = n_rep)
  
  # Grab the variance components
  V_E <- training_pop$pheno_val$var_comp$V_E
  V_R <- training_pop$pheno_val$var_comp$V_R
  
  # Randomly select 10 individuals to become checks
  check_pop <- training_pop %>%
    subset_pop(individual = sample(indnames(training_pop), 10))
  
  
  # Randomly select parents
  cb <- sim_crossing_block(parents = indnames(training_pop), n.crosses = n_crosses)
  
  
  ## Predict genetic variance
  G <- genotype(genome = genome, pop = training_pop)
  y_in <- training_pop$pheno_val$pheno_mean
  
  # Expectation
  sim_fam <- pop_predict_quick(G.in = G, y.in = y_in, map.in = map_use, crossing.table = cb)
  
  # # Via PopVar
  # G_in <- geno_to_popvar(genome = genome, geno = training_pop$geno)
  # pred_out <- pop.predict(G.in = G_in, y.in = y_in, map.in = map_use, crossing.table = cb,
  #                         nInd = 50, min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, nSim = 3,
  #                         nCV.iter = 1, models = "rrBLUP")
  # 
  # # Extract data
  # sim_fam_PV <- pred_out$predictions %>% 
  #   as.data.frame() %>% 
  #   mutate_all(unlist)
  # 
  # # Sanity checks
  # cor(sim_fam$pred_mu, sim_fam_PV$pred.mu)
  # cor(sim_fam$pred_varG, sim_fam_PV$pred.varG)
  # cor(sim_fam$pred_mu_sp_high, sim_fam_PV$mu.sp_high)
  
  
  # Create those populations
  real_fam <- sim_family_cb(genome = genome, pedigree = ped, founder.pop = training_pop, crossing.block = cb)
  
  # Simulate phenotypes
  real_fam1 <- sim_trial(pop = real_fam, h2 = h2, n.env = 2, n.rep = n_rep,
                         check.pop = check_pop, check.rep = n_check_rep, V_E = V_E,
                         V_R = V_R)
  
  # Separate by family and calculate the genetic variance
  rf_summary <- real_fam$geno_val %>% 
    mutate(fam = str_extract(ind, "[0-9]{4}")) %>% 
    group_by(fam) %>% 
    summarize(mu = mean(trait1),
              V_G = var(trait1),
              mu_sp = mean(sort(trait1, decreasing = T)[1:10]))
  
  # Correlation with the true values
  mu_acc <- cor(sim_fam$pred_mu, rf_summary$mu)
  varG_acc <- cor(sim_fam$pred_varG, rf_summary$V_G)
  mu_sp_acc <- cor(sim_fam$pred_mu_sp_high, rf_summary$mu_sp)
  
  # Extract the phenotypes and assign families
  phenos <- real_fam1$pheno_val$pheno_obs %>% 
    mutate(family = ifelse(ind %in% indnames(check_pop), "check", str_extract(ind, "[0-9]{4}")))
  
  
  # Model
  # Fit genotype as fixed and env + gxe as random
  rf_mod <- mmer2(fixed = phenoval ~ 1 + env, random = ~ ind + at(family):ind + ind:env, 
                  data = phenos, method = "NR",
                  constraint = TRUE)
                  
  # Summary
  summ_rf_mod <- summary(rf_mod)
  summ_rf_mod1 <- summary(rf_mod) 
        
  # Estimates of genetic variance
  varG_hat <- summ_rf_mod1$w %>% 
    rownames_to_column("comp") %>% 
    filter(str_detect(comp, "[0-9]{4}")) %>%
    pull(VarComp)
    
  

  
  # Create a DF and return
  data.frame(mu_acc = mu_acc, varG_acc = varG_acc, mu_sp_acc = mu_sp_acc)
  
}, simplify = FALSE)
 

# Interpretation
out_df <- simulation_results[,,] %>% 
  t() %>% 
  as.data.frame() %>%
  mutate_all(unlist) %>%
  mutate(iter = seq(nrow(.)))

# Tidy
results_tidy <- out_df %>%
  gather(measure, value, -iter)

# Summarize
results_tidy %>%
  group_by(measure) %>%
  summarize(mean = mean(value),
            sd = sd(value))

# Model
mod <- lm(mu_sp_acc ~ mu_acc, data = out_df)
mod1 <- lm(mu_sp_acc ~ mu_acc + varG_acc, data = out_df)

# Compare models
anova(mod, mod1)
