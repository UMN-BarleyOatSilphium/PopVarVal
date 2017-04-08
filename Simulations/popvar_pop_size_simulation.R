## PopVar Size Simulations

library(tidyverse)
library(stringr)
library(PopVar)
library(qtl)
library(pbsim)
library(simcross)
library(lme4)
library(parallel)


proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/PopVarVal"

sim_dir <- file.path(proj_dir, "Simulations")


data("s2_genos")
data("s2_snp_info")


## Immutable parameter 
# Number of simulation iterations
n_iter <- 50

# Number of replicates in environments
n_rep <- 2

# Number of environments to phenotype
n_env <- 4

# Error variance for genetic map positions
map_error <- 0.002

# Number of bi-parental populations
n_cross <- 20


## Mutable parameters
# Heritability
herit <- c(0.2, 0.5, 0.8)

# Number of QTL
nQTL <- c(30, 100)

# Number of individuals in the bi-parental family
n_ind <- c(25, 75, 125, 200, 275)

# Create a data.frame of the mutable parameters to iterate over
param_df <- expand.grid(herit = herit, nQTL = nQTL, n_ind = n_ind, iter = seq(n_iter))

# Detect cores
n_cores <- detectCores()

# Split the data.frame
df_split <- split(x = param_df, f = sort(rep_len(x = seq(n_cores), length.out = nrow(param_df))))





# Map length based on Munoz-Amatriain et al 2011
map_len <- c(143.3, 172.93, 180.13, 146.5, 189.9, 142.3, 162.5)

# Number of markers per chromosome
n_mar <- table(s2_snp_info$chrom)

# Map list
s2_snp_map <- s2_snp_info %>% 
  split(.$chrom) %>% 
  map(function(chrom) {
    chrom_pos <- chrom$cM_pos
    names(chrom_pos) <- chrom$rs
    class(chrom_pos) <- "A"
    return(chrom_pos) })

class(s2_snp_map) <- "map"

# Create the genome
genome <- sim_genome(len = map_len, n.mar = n_mar, map = s2_snp_map)

# Vector of base population line names
base_pop <- row.names(s2_genos)

# All combinations of base lines to form the populations
par_combn <- as.data.frame(t(combn(x = base_pop, m = 2)))
names(par_combn) <- c("par1", "par2")




# Iterate over sections of the data.frame
sim_results <- mclapply(X = df_split, FUN = function(df) {
  
  # Empty list to store results
  results <- vector("list", nrow(df))
  
  # Iterate over mutable simulation parameters
  for (p in seq(nrow(df))) {
    
    params <- param_df[p,]
    
    h2 <- unlist(params[1])
    n_QTL <- unlist(params[2])
    nind <- unlist(params[3])
    
    # Define genetic architecture
    # Blank QTL model
    qtl.model <- matrix(NA, nrow = n_QTL, ncol = 4)
    genome <- sim_gen_model(genome = genome, qtl.model = qtl.model, add.dist = "geometric")
    
    
    # Recode the genotypes
    s2_genos_recode <- s2_genos + 1
    # Phenotype the training population
    tp_pheno <- sim_pheno(genome = genome, geno = s2_genos_recode, h2 = h2, n.env = n_env, n.rep = 3)
    
    # Extract variance components
    V_G_base <- tp_pheno$var_comp$V_G
    V_E_base <- tp_pheno$var_comp$V_E
    V_R_base <- tp_pheno$var_comp$V_R
    
    # Select crosses
    crosses_sel <- par_combn %>% 
      sample_n(size = n_cross)
    
    # Create a RIL pedigree
    ped <- sim_pedigree(n.ind = 5000)
    
    # Iterate over the crosses and calculate the "true" genetic variance
    true_varG <- crosses_sel %>% 
      rowwise() %>%
      do(V_G_true = {
        
        # Gather the founder genotypes
        pars_genos <- s2_genos_recode[c(.$par1, .$par2),]
        # Split the genos based on the map
        pars_genos <- map(genome$map, function(chr) pars_genos[,names(chr)])
        
        # Simulate a family
        fam <- sim_family(genome = genome, pedigree = ped, founder_geno = pars_genos)
        
        # Get the genotypic values
        geno_val <- sim_pheno(genome = genome, geno = fam, h2 = h2)$geno_val
        geno_val %>% 
          summarize(varG = var(geno_val))
        
      })
    
    # Recombine data.frame
    crosses_sel_true <- true_varG %>% 
      mutate(V_G_true = unlist(V_G_true)) %>%
      ungroup() %>% 
      mutate(par1 = crosses_sel$par1, par2 = crosses_sel$par2) %>% 
      select(par1, par2, V_G_true)
    
    
    
    ## Run PopVar
    
    # Save the genetic map for use in PopVar
    popvar_map <- s2_snp_info %>% 
      select(-alleles) %>% 
      mutate(chrom = as.numeric(as.factor(chrom))) %>%
      filter(!rs %in% genome$gen_model$snp) %>%
      as.data.frame()
    
    # Subset the genos for the non-QTL markers
    s2_genos_noQTL <- s2_genos[,popvar_map$rs]
    
    # Save the training genotype data for use in PopVar
    popvar_Gin <- as.data.frame(cbind( c("", base_pop), 
                                       rbind(colnames(s2_genos_noQTL), s2_genos_noQTL)) )
    
    # Save the training phenotype data for use in PopVar
    popvar_pheno <- tp_pheno$pheno_mean
    
    # Save the crossing block for use in PopVar
    popvar_cb <- crosses_sel
    
    # Run
    bp_predict <- pop.predict(G.in = popvar_Gin, y.in = popvar_pheno, map.in = popvar_map, 
                              crossing.table = popvar_cb, tail.p = 0.1, nInd = 150, min.maf = 0,
                              mkr.cutoff = 1, entry.cutoff = 1, impute = "pass", nCV.iter = 1, 
                              models = "rrBLUP")
    
    # Extract the results
    columns <- c("Par1", "Par2", "pred.mu", "pred.varG", "pred.varG_sd", "mu.sp_high")
    bp_predictions <- bp_predict$predictions[, columns] %>%
      apply(MARGIN = 2, FUN = unlist) %>% 
      as.data.frame() %>% 
      mutate_at(-1:-2, .funs = funs(as.numeric(as.matrix(.))))
    
    
    ## Simulate making the actual crosses
    # Assume the families will be inbred to the F_4
    f4_ped <- sim_pedigree(n.ind = nind, n.selfgen = 3)
    
    # Create the families
    families_sel <- crosses_sel %>% 
      # Add family number
      mutate(fam_num = seq(n())) %>%
      rowwise() %>%
      do(family_geno = {
        
        # Gather the founder genotypes
        pars_genos <- s2_genos_recode[c(.$par1, .$par2),]
        # Split the genos based on the map
        pars_genos <- map(genome$map, function(chr) pars_genos[,names(chr)])
        
        # Simulate a family
        sim_family(genome = genome, pedigree = f4_ped, founder_geno = pars_genos, 
                   family.num = .$fam_num)
        
      })
    
    # Combine the genotypes
    pop_geno <- families_sel %>% 
      unlist(recursive = F) %>% 
      do.call("rbind", .)
    
    ## Simulate field testing
    pop_phenos <- sim_pheno(genome = genome, geno = pop_geno, h2 = h2, n.env = n_env, 
                            n.rep = n_rep, V_E = V_E_base, V_R = V_R_base)
    
    ## Estimate the genetic variance by family
    # Assign family name
    pop_pheno_vals <- pop_phenos$pheno_val %>%
      mutate(family = str_sub(entry, 1, 7))
    
    # Fit a mixed model to estimate genetic variance
    # Also calculate the top 10% of BLUPs
    crosses_sel_hat <- pop_pheno_vals %>% 
      group_by(family) %>% 
      do({
        mod <- lmer(pheno_val ~ (1|entry) + env + (1|entry:env), data = .)
        # Find variance of random effects
        vcov_fit <- VarCorr(mod)
        # Extract genetic variance estimate
        varG_hat <- as.numeric(subset(as.data.frame(vcov_fit), grp == "entry", vcov))
        # Find random effects of genotypes
        random_eff <- as.matrix(ranef(mod)$entry)
        # Calculate mean and mu_sp
        mu <- mean(random_eff)
        mu_sp <- mean(head(sort(random_eff, decreasing = T), n = 0.1 * length(random_eff)))
        # Export data.frame
        data.frame(varG_hat = varG_hat, mu_sp = mu_sp, mu = mu) })
    
    # Combine predictions, truth, and observations
    crosses_combined <- bind_cols(crosses_sel, crosses_sel_hat) %>%
      select(-family) %>% 
      full_join(., crosses_sel_true) %>% 
      full_join(., bp_predictions, by = c("par1" = "Par1", "par2" = "Par2"))
    
    # Add the df to the list
    results[[p]] <- crosses_combined
    
  }
  
  # Return the results
  return(results) })


# Save the output
save("sim_results", file = "simulation_results.RData")



