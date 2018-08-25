## PopVarVal: Simulation to assess conditions where PopVar succeeds
## 
## This script will run simulations to perturb different parameters of a
## quantitative trait or a genetic map in order to find the conditions that
## might explain the success or failure of PopVar
## 
## Author: Jeff Neyhart
## Last modified: August 20, 2018
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
n_iter <- 10
n_env <- 2
n_rep <- 1
n_crosses <- 50

## Outline the parameters to perturb
h2_list <- c(0.2, 0.5, 0.8)
nQTL_list <- c(30, 100)
map_error_list <- c(0.01, 0.05, 0.10)

# Create a data.frame of parameters
param_df <- crossing(h2 = h2_list, nQTL = nQTL_list, map_error = map_error_list, iter = seq(n_iter))

map_sim <- s2_snp_info %>% 
  split(.$chrom) %>% 
  map(~setNames(.$cM_pos, .$rs))

# Create the base genome - this is fixed for all simulations
genome <- sim_genome(map = map_sim)


# Split the parameter df
param_df_split <- param_df %>%
  assign_cores(n_cores) %>%
  split(.$core)

# Parallelize
simulation_out <- mclapply(X = param_df_split, FUN = function(core_df) {
  
  # Create a results list
  results_out <- vector("list", nrow(core_df))
  
  # Iterate over the rows of the param_df
  for (i in seq_along(results_out)) {
    
    h2 <- param_df$h2[i]
    L <- param_df$nQTL[i]
    maxL <- max(param_df$nQTL)
    map_error <- param_df$map_error[i]
  
    # Simulate QTL
    qtl_model <- matrix(NA, ncol = 4, nrow = L)
    genome1 <- sim_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "geometric", max.qtl = maxL)
    
    # Create the TP
    tp <- create_pop(genome = genome1, geno = s2_cap_genos, ignore.gen.model = FALSE)
    
    # Phenotype
    tp1 <- sim_phenoval(pop = tp, h2 = h2, n.env = n_env, n.rep = n_rep)
    # Convert the phenotypes into something PopVar will handle
    pheno_use <- tp1$pheno_val$pheno_mean
    
    ## Select the best among the TP
    tp_best <- select_pop(pop = tp1, intensity = 0.1, index = 1,  type = "phenotypic")
    
    # Randomly generate 50 crosses to simulate, then actually create
    crossing_block <- sim_crossing_block(parents = indnames(tp_best), n.crosses = n_crosses)
    ped <- sim_pedigree(n.ind = 40, n.selfgen = Inf)
    
    # Create these crosses
    cycle1 <- sim_family_cb(genome = genome1, pedigree = ped, founder.pop = tp_best, crossing.block = crossing_block)
    
    # Combine the tp with cycle1
    pv_pop <- combine_pop(pop_list = list(tp1, cycle1))
  
    ## Add error to the genetic map
    map_use <- map_sim %>% 
      # Remove the QTL
      map(~.[names(.) %in% markernames(genome1, include.qtl = FALSE)]) %>%
      map(~. + rnorm(n = length(.), mean = 0, sd = sqrt(map_error))) %>% 
      map(sort) %>%
      map(~data_frame(marker = names(.), pos = .)) %>% 
      map2_df(.x = ., .y = names(.), ~mutate(.x, chrom = .y)) %>%
      select(marker, chrom, pos) %>%
      as.data.frame()
      
    
    ## Genotype and convert for PopVar
    geno_data <- genotype(genome = genome1, pop = pv_pop)
    geno_use <- as.data.frame(cbind( c("", row.names(geno_data)), rbind(colnames(geno_data), geno_data)) )
    
    # Randomly create crosses from the cycle1 individuals
    crossing_block <- sim_crossing_block(parents = indnames(cycle1), n.crosses = n_crosses)

    # Pass to PopVar
    pred_out <- pop.predict(G.in = geno_use, y.in = pheno_use, map.in = map_use,
                            crossing.table = crossing_block, tail.p = 0.1, nInd = sim_pop_size,
                            min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
                            nSim = n_pops, nCV.iter = 1, models = "rrBLUP", impute = "pass")
    
    # Tidy
    tidy_pred_out <- pred_out$predictions %>% 
      mutate_all(unlist) %>% 
      as_data_frame() %>%
      select(parent1 = Par1, parent2 = Par2, mpv = midPar.Pheno, pred_mu = pred.mu,
             pred_varG = pred.varG, mu_sp = mu.sp_high)
    
    
    
    
    ## Calculate the expected genetic variance in these populations
    
    # Get the QTL information
    qtl_info <- pull_qtl(genome1) %>%
      # Pull only effective QTL
      filter(add_eff != 0)
    
    # Get the map and genotypes of only the QTL
    qtl_geno <- pull_genotype(genome = genome1, geno = pv_pop$geno, loci = qtl_info$qtl_name) - 1
    
    
    ## Calculate the expected genetic variance at each locus, assuming p = q = 0.5
    qtl_info1 <- qtl_info %>% 
      mutate(exp_var = add_eff^2) %>%
      split(.$chr)
    
    
    ## Second calculate the expected covariance between QTL
    # Create a distance matrix
    qtl_dist <- qtl_info1 %>% 
      map(~{
        temp <- .
        as.matrix(dist(temp$pos)) %>% 
          `dimnames<-`(., list(temp$qtl_name, temp$qtl_name))
      })
    
    # Calculate pairwise D
    pairwise_D <- qtl_dist %>%
      map(~qtl::mf.h(.)) %>%
      map(~((1 - (2 * .)) / (1 + (2 * .))) )
    
    # Calculate the pairwise product of all marker effects
    pairwise_prod <- qtl_info1 %>%
      map(~{
        temp <- .
        crossprod(t(temp$add_eff)) %>% `dimnames<-`(., list(temp$qtl_name, temp$qtl_name))
      })
    
    # The covariance is the QTL effect product multiplied by the expected D
    exp_covar <- map2(pairwise_prod, pairwise_D, `*`) %>% 
      Matrix::.bdiag() %>%
      as.matrix() %>%
      `dimnames<-`(., list(qtl_info$qtl_name, qtl_info$qtl_name))
    exp_var <- qtl_info1 %>% 
      bind_rows() %>% 
      column_to_rownames("qtl_name") %>% 
      select(exp_var) %>% 
      as.matrix()
    
    # List of expected variance, mean, etc
    expected_out <- vector("list", nrow(crossing_block))
    
    # Iterate over the crossing block
    for (j in seq(nrow(crossing_block))) {
      
      pars <- as.character(crossing_block[j,])
      
      # Get the QTL genotypes and find those that are polymorphic
      poly_qtl <- colMeans(qtl_geno[pars,]) %>% {.[.==0]} %>% names()
      exp_var1 <- exp_var[poly_qtl,,drop = F]
      exp_covar1 <- exp_covar[poly_qtl, poly_qtl]
      
      # The expected variance is the sum of the variances at the polymorphic QTL, plus 2 times
      # the expected covariance between all polymorphic QTL
      exp_var_j <- sum(exp_var1) + (2 * sum(exp_covar1[upper.tri(exp_covar1)]))
      
      # The expected mu is simply the mean of the genotypic values of the two parents
      exp_mu_j <- mean(subset(pv_pop$geno_val, ind %in% pars, trait1, drop = T))
      
      # The expected mu_sp is a function of the expected variance and the expected mean
      exp_musp_j <- exp_mu_j + (1.76 * sqrt(exp_var_j))
      
      # Create a df and add to the list
      expected_out[[j]] <- crossing_block[j,] %>% mutate(exp_mu = exp_mu_j, exp_var = exp_var_j, exp_musp = exp_musp_j)
      
    }
    
    
    expected_out_df <- bind_rows(expected_out)
    
    ## Caculate the accuracy of mean, varG, and mu_sp
    h2 <- data_frame(
      term = c("mean", "varG", "mu_sp"),
      accuracy = c(cor(expected_out_df$exp_mu, tidy_pred_out$pred_mu),
                   cor(expected_out_df$exp_var, tidy_pred_out$pred_varG),
                   cor(expected_out_df$exp_musp, tidy_pred_out$mu_sp))
    )
    
    ## Add the accuracy results to the results list
    results_out[[i]] <- list(accuracy = accuracy)
    
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















