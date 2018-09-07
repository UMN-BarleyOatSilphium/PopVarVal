## PopVarVal: Simulation to assess conditions where PopVar succeeds
## 
## This script will run simulations to perturb different parameters of a
## quantitative trait or a genetic map in order to find the conditions that
## might explain the success or failure of PopVar
## 
## Author: Jeff Neyhart
## Last modified: August 30, 2018
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
# library(Matrix)
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
h2_list <- c(0.2, 0.5, 0.8)
nQTL_list <- c(30, 100)
map_error_list <- c(0, 0.01, 1)
tp_size_list <- seq(150, 600, by = 150)

# Create a data.frame of parameters
param_df <- crossing(h2 = h2_list, nQTL = nQTL_list, map_error = map_error_list, tp_size = tp_size_list, iter = seq(n_iter))

map_sim <- s2_snp_info %>% 
  split(.$chrom) %>% 
  map(~setNames(.$cM_pos, .$rs)) %>%
  map(~structure(., class = "A"))  %>%
  structure(., class = "map") %>%
  # Jitter
  qtl::jittermap(.)

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
    
    h2 <- core_df$h2[i]
    L <- core_df$nQTL[i]
    maxL <- max(param_df$nQTL)
    map_error <- core_df$map_error[i]
    tp_size <- core_df$tp_size[i]
  
    # Simulate QTL
    qtl_model <- matrix(NA, ncol = 4, nrow = L)
    genome1 <- sim_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "geometric", max.qtl = maxL)
    
    ## Create the TP
    # Sample from the genotypes
    sample_genos <- s2_cap_genos[sort(sample(nrow(s2_cap_genos), size = tp_size)),]
    tp <- create_pop(genome = genome1, geno = sample_genos, ignore.gen.model = FALSE)
    # # Genotype
    # pv_geno <- genotype(genome = genome1, pop = tp)
    
    # Phenotype
    tp1 <- sim_phenoval(pop = tp, h2 = h2, n.env = n_env, n.rep = n_rep)
    par_pop <- tp1
    
    ##### Create a set of cycle 1 progeny from which to predict crosses
    ## Select the best tp individuals
    tp_select <- select_pop(pop = tp1, intensity = 0.1, index = 1, type = "phenotypic")
    # Randomly create crosses from the TP individuals
    crossing_block <- sim_crossing_block(parents = indnames(tp_select), n.crosses = 40)
    # Pedigree to accompany the crosses
    ped <- sim_pedigree(n.ind = 25, n.selfgen = Inf)
    
    # Make theses crosses
    par_pop <- sim_family_cb(genome = genome1, pedigree = ped, founder.pop = tp_select, crossing.block = crossing_block)
    #####

    # Randomly create crosses from the TP individuals
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
      mutate(ind = as.character(ind))
    
    # Create a new genome
    genome2 <- genome1
    # Assign the qtl model with every marker as a QTL
    marker_qtl_model <- find_markerpos(genome = genome1, marker = tp_mar_eff$mar_eff$marker) %>%
      mutate(qtl_name = row.names(.)) %>%
      left_join(., tp_mar_eff$mar_eff, by = c("qtl_name" = "marker")) %>%
      mutate(dom_eff = 0, qtl1_pair = NA) %>%
      select(chr, pos, add_eff = trait1, dom_eff, qtl_name, qtl1_pair)
    
    # Add to the genome
    genome2$gen_model <- list(marker_qtl_model)
    # Add the map to the genome
    genome2$map <- genome2$map %>%
      # Remove the QTL
      map(~. + rnorm(n = length(.), mean = 0, sd = sqrt(map_error))) %>% 
      map(sort)
    
    
    # Use the genetic variance prediction function
    exp_gen_var <- calc_exp_genvar1(genome = genome2, pedigree = ped, founder.pop = par_pop, crossing.block = crossing_block)
    
    pred_out <- exp_gen_var$crossing_block %>%
      left_join(x = ., y = par_pop$pred_val, by = c("parent1" = "ind")) %>% left_join(x = ., y = par_pop$pred_val, by = c("parent2" = "ind")) %>%
      mutate(exp_mu = (trait1.x + trait1.y) / 2,
             exp_musp = exp_mu + (1.76 * sqrt(exp_varG))) %>%
      select(-contains("trait")) %>%
      rename_all(~str_replace(., "exp", "pred"))

    ## What is the average QTL variance and covariance?
    pred_var_covar_summ <- map(exp_gen_var$var_covar, ~rbind(.[1], .[2] / 2)) %>% 
      do.call("cbind", .) %>% 
      rowMeans()

    
    ## Calculate the expected genetic variance in these populations
    expected_var_out <- calc_exp_genvar1(genome = genome1, pedigree = ped, founder.pop = par_pop, crossing.block = crossing_block)
    expected_var <- expected_var_out$crossing_block %>%
      mutate(exp_musp = exp_mu + (1.76 * sqrt(exp_varG))) %>%
      select(-trait)
    
    exp_var_covar_summ <- map(expected_var_out$var_covar, ~rbind(.[1], .[2] / 2)) %>% 
      do.call("cbind", .) %>% 
      rowMeans()
    
    
    
    ## Summarize the expected and predicted values
    expected_predicted <- full_join(pred_out, expected_var, by = c("parent1", "parent2")) %>% 
      gather(parameter, value, pred_mu:exp_musp) %>% 
      separate(parameter, c("type", "parameter"), sep = "_") %>% 
      spread(type, value)
    
    # Summarize
    results_summ <- expected_predicted %>%
      group_by(parameter) %>% 
      nest(exp, pred) %>%
      mutate(accuracy = map_dbl(data, ~cor(.$exp, .$pred)),
             bias = map_dbl(data, ~mean(.$pred - .$exp) / mean(.$exp)))
    
    ## Add the accuracy results to the results list
    results_out[[i]] <- list(
      results = results_summ,
      var_covar_summ = list(pred = pred_var_covar_summ, exp = exp_var_covar_summ)
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
save_file <- file.path(result_dir, "popvar_suitability_simulation_results.RData")
save("popvar_simulation_out", file = save_file)















