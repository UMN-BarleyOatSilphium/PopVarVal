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



# Number cores
n_cores <- detectCores()



## Fixed parameters
sim_pop_size <- 150
n_pops <- 25
n_iter <- 25
n_env <- 2
n_rep <- 1
n_crosses <- 50

## Outline the parameters to perturb
trait1_h2_list <- trait2_h2_list <- c(0.2, 0.5, 0.8)
nQTL_list <- c(30, 100)
tp_size_list <- seq(150, 600, by = 150)
gencor_list <- c(0.2, 0.5, 0.8)

# Create a data.frame of parameters
param_df <- crossing(trait1_h2 = trait1_h2_list, trait2_h2 = trait2_h2_list, nQTL = nQTL_list, tp_size = tp_size_list, 
                     gencor = gencor_list, iter = seq(n_iter))

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
    
    # Simulate QTL
    qtl_model <- replicate(n = 2, matrix(NA, ncol = 4, nrow = L), simplify = FALSE)
    genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "geometric", max.qtl = maxL,
                                   corr = gencor, prob.corr = cbind(0, 1))
    
    ## Create the TP
    # Sample from the genotypes
    sample_genos <- s2_cap_genos[sort(sample(nrow(s2_cap_genos), size = tp_size)),]
    tp <- create_pop(genome = genome1, geno = sample_genos, ignore.gen.model = FALSE)
    # Measure the genetic correlation in the TP
    tp_cor <- cor(tp$geno_val[,-1])[1,2]
    
    
    # Phenotype
    tp1 <- sim_phenoval(pop = tp, h2 = c(trait1_h2, trait2_h2), n.env = n_env, n.rep = n_rep)
    # Measure the phenotypic correlation in the TP
    tp_pheno_cor <- cor(tp1$pheno_val$pheno_mean[,-1])[1,2]
    par_pop <- tp1
    
    ##### Create a set of cycle 1 progeny from which to predict crosses
    ## Select the best tp individuals for both traits
    tp_select <- select_pop(pop = tp1, intensity = 0.1, index = c(1, 1), type = "phenotypic")
    # Randomly create crosses from the TP individuals
    crossing_block <- sim_crossing_block(parents = indnames(tp_select), n.crosses = 40)
    # Pedigree to accompany the crosses
    ped <- sim_pedigree(n.ind = 25, n.selfgen = Inf)
    
    # Make theses crosses
    par_pop <- sim_family_cb(genome = genome1, pedigree = ped, founder.pop = tp_select, crossing.block = crossing_block)
    #####
    
    # Randomly create crosses from the cycle1 individuals
    crossing_block <- sim_crossing_block(parents = indnames(par_pop), n.crosses = n_crosses)
    # Pedigree to accompany the crosses
    ped <- sim_pedigree(n.ind = sim_pop_size, n.selfgen = Inf)
    
    
    # # # Pass to PopVar
    # # Convert genotypes into something useable for PopVar
    # geno_use <- as.data.frame(cbind( c("", row.names(pv_geno)), rbind(colnames(pv_geno), pv_geno)) )
    # 
    # pred_out <- pop.predict(G.in = geno_use, y.in = pheno_use, map.in = map_use,
    #                         crossing.table = crossing_block, tail.p = 0.1, nInd = sim_pop_size,
    #                         min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
    #                         nSim = n_pops, nCV.iter = 1, models = "rrBLUP", impute = "pass")
    # 
    # # Tidy
    # tidy_pred_out <- pred_out$predictions %>%
    #   mutate_all(unlist) %>%
    #   as_data_frame() %>%
    #   select(parent1 = Par1, parent2 = Par2, mpv = midPar.Pheno, pred_mu = pred.mu,
    #          pred_varG = pred.varG, pred_musp = mu.sp_high)
    
    
    ## Use the expected variance calculation as a sanity check - Actually it's much quicker than popvar
    # First predict marker effects, then predict the genotypic value of the TP
    tp_mar_eff <- pred_mar_eff(genome = genome1, training.pop = tp1)
    par_pop <- pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = par_pop)
    par_pop$pred_val <- par_pop$pred_val %>% 
      mutate(ind = as.character(ind)) %>%
      gather(trait, prediction, -ind)
    
    # Create a new genome
    genome2 <- genome1
    # Tidy the predicted marker effects
    marker_effs <- tp_mar_eff$mar_eff %>%
      gather(trait, add_eff, -marker)
    
    # Assign the qtl model with every marker as a QTL
    marker_qtl_model <- find_markerpos(genome = genome1, marker = tp_mar_eff$mar_eff$marker) %>%
      mutate(qtl_name = row.names(.)) %>%
      left_join(., marker_effs, by = c("qtl_name" = "marker")) %>%
      mutate(dom_eff = 0, qtl1_pair = NA) %>%
      select(trait, chr, pos, add_eff, dom_eff, qtl_name, qtl1_pair) %>%
      split(.$trait) %>%
      map(select, -trait)
    
    ## Edit the qtl pair for trait2
    marker_qtl_model$trait2 <- marker_qtl_model$trait2 %>% 
      mutate(qtl1_pair = qtl_name)
    
    # Add to the genome
    genome2$gen_model <- marker_qtl_model
    
    
    # Use the genetic variance prediction function
    pred_out <- calc_exp_genvar(genome = genome2, pedigree = ped, founder.pop = par_pop, crossing.block = crossing_block) %>%
      left_join(x = ., y = par_pop$pred_val, by = c("parent1" = "ind", "trait")) %>% 
      left_join(x = ., y = par_pop$pred_val, by = c("parent2" = "ind", "trait")) %>%
      mutate(exp_mu = (prediction.x + prediction.y) / 2,
             exp_musp = exp_mu + (1.76 * sqrt(exp_varG))) %>%
      select(-contains("prediction")) %>%
      rename_all(~str_replace(., "exp", "pred"))
    
    ##
    
    
    ## Calculate the expected genetic variance in these populations
    expected_var <- calc_exp_genvar(genome = genome1, pedigree = ped, founder.pop = par_pop, crossing.block = crossing_block) %>%
      mutate(exp_musp = exp_mu + (1.76 * sqrt(exp_varG)))
    
    
    ## Summarize the expected and predicted values
    expected_predicted <- full_join(pred_out, expected_var, by = c("parent1", "parent2", "trait"))
    
    # Summarize variance and means
    expected_predicted_mean_var <- expected_predicted %>%
      select(-contains("corG")) %>%
      gather(parameter, value, pred_mu:exp_musp) %>% 
      separate(parameter, c("type", "parameter"), sep = "_") %>% 
      spread(type, value)
    
    # Summarize correlations
    expected_predicted_cor <- expected_predicted %>%
      select(parent1, parent2, contains("corG")) %>%
      distinct() %>%
      gather(parameter, value, pred_corG, exp_corG) %>% 
      separate(parameter, c("type", "parameter"), sep = "_") %>% 
      spread(type, value)
    
    # Combine
    expected_predicted1 <- bind_rows(expected_predicted_mean_var, expected_predicted_cor)
    
    # Summarize
    results_summ <- expected_predicted1 %>%
      group_by(trait, parameter) %>% 
      nest(exp, pred) %>%
      mutate(accuracy = map_dbl(data, ~cor(.$exp, .$pred)),
             bias = map_dbl(data, ~mean(.$pred - .$exp) / mean(.$exp)))
    
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
popvar_simulation_out <- bind_rows(simulation_out)

# Save
save_file <- file.path(result_dir, "popvar_gencor_simulation_results.RData")
save("popvar_simulation_out", file = save_file)















