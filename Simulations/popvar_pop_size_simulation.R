## PopVar Size Simulations
## 
## This script will run some simulations to determine the effect of various 
## perturbations on the ability to statistically differentiate genetic variance 
## in different populations
## 

packages <- c("dplyr", "tidyr", "tibble", "purrr", "readr", "stringr", "modelr", 
              "pbsim", "pbr", "purrrlyr", "lme4", "parallel", "gws")

# Set the directory of the R packages
package_dir <- NULL
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/PopVarVal"
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/PopVarVal"

sim_dir <- file.path(proj_dir, "Simulations")

# Load the S6 data
data(s6_cap_haploid)
data(s6_snp_info)




## Immutable parameter 
# Number of simulation iterations
n_iter <- 5
# Number of replicates in environments
n_rep <- 3
# Number of environments to phenotype
n_env <- 4
# Number of bi-parental populations
n_cross <- 20
# Number of selfing generations
self_gen <- 2
# Selection intensity for superior progeny calculation
i_sp <- 0.10 

# Define the biologically "real" map
real_map <- s6_snp_info %>% 
  split(.$chrom) %>% 
  map(~ set_names(x = .$cM_pos, nm = .$rs)) %>%
  structure(class = "map")

# Create the real genome
genome <- sim_genome(map = real_map, type = "hypred")


## Mutable parameters
# Error variance for genetic map positions
map_error <- c(1, 5, 10)
# Heritability
herit <- c(0.2, 0.5, 0.8)
# Number of QTL
nQTL <- c(30, 100)
# Number of individuals in the bi-parental family
n_ind <- c(25, 100, 500)

# Create a data.frame of the mutable parameters to iterate over
param_df <- expand.grid(iter = seq(n_iter), map_error = map_error, herit = herit, 
                        nQTL = nQTL, n_ind = n_ind)

# The maximum number of qtl is the max of the nQTL vector
max_qtl <- max(nQTL)


# Detect cores
n_cores <- detectCores()

# Add a column of cores and split the data.frame
param_df %>% 
  mutate(core = sort(rep(seq(n_cores), length.out = nrow(.)))) %>%
  split(.$core) %>%
  mclapply(function(core_df) {
    
    # Apply by row
    core_df %>%
      by_row(function(i) {
        
        ## Variables
        map_error <- i$map_error
        herit <- i$herit
        n_qTL <- i$nQTL
        n_ind <- i$n_ind
        
        ## Create the genomes with genetic models
        
        # Add error to the genetic map
        useable_map <- real_map %>% 
          map(~ . + rnorm(n = length(.), mean = 0, sd = sqrt(map_error))) %>%
          # Correct for negative map positions
          map(~ . + abs(min(.))) %>%
          map(sort) %>%
          structure(class = "map")
        
        # Create the usable genome
        genome_usable <- sim_genome(map = useable_map, type = "hypred")
        
        # Genetic model matrix
        qtl.model <- matrix(nrow = n_qtl, ncol = 4)
        
        # Add genetic architecture
        genome <- sim_gen_model(genome, qtl.model, add.dist = "geometric", max.qtl = max_qtl)
        # Copy the genetic model to the useable genome
        genome_usable$gen_model <- genome$gen_model
        
        # Define a pedigree for the families of selection candidates
        pedigree <- sim_pedigree(n.ind = n_ind, n.selfgen = self_gen)
        
        
        
        ## PopVar prep
        # Pull out marker names
        marker_names <- markernames(genome)
        
        # Popvar prep for the usable genome
        popvar_map <- map_to_popvar(genome_usable) %>%
          filter(!marker %in% qtlnames(genome_usable))
        
        # Create the base population
        training_pop <- create_pop(genome = genome, geno = s6_cap_haploid)
        # Phenotypic observations for the training population
        training_pop <- sim_phenoval(pop = training_pop, h2 = herit, n.env = n_env, n.rep = n_rep)
        
        # Select the best 50 MN and best 50 ND
        NDnames <- indnames(training_pop) %>% 
          str_subset("^ND")
        MNnames <- setdiff(indnames(training_pop), NDnames)
        
        NDbest <- subset_pop(pop = training_pop, individual = NDnames) %>%
          select_pop(intensity = 50, index = 1)
        MNbest <- subset_pop(pop = training_pop, individual = MNnames) %>%
          select_pop(intensity = 50, index = 1)
        
        potential_parents <- combine_pop(list(NDbest, MNbest)) %>% indnames()

        
        
        # Extract the variance components
        trait_var_comp <- training_pop$pheno_val$var_comp
        
        # Designate the training population
        ## The 'popvar_geno_pop' is for use in PopVar (before parental selection)
        ## The 'selection_training_pop' is for predicting the selection candidates
        popvar_geno_pop <- training_pop
        
        # Use the phenotypic means as the predicted genotypic values
        popvar_geno_pop$pred_val <- popvar_geno_pop$pheno_val$pheno_mean
        
        
        ### Parental selection
        # The maximum number of crosses given the potential parents
        max_crosses <- choose(length(potential_parents), 2)
        
        # Simualte a crossing block
        crossing_block <- sim_crossing_block(parents = potential_parents, n.crosses = max_crosses)
        
        # Extract the genotype data for use in the function
        G <- genotype(genome = genome, pop = popvar_geno_pop)
        # Phenotype data
        phenos <- popvar_geno_pop$pheno_val$pheno_mean
        
        
        sim_fam <- pop_predict_quick(G.in = G, y.in = phenos, map.in = popvar_map, 
                                     crossing.table = crossing_block, tail.p = i_sp)
        
        
        # Subset the families by those with predicted genetic variances that are
        # less than or greater than the 10% upper and lower quantiles
        # Then sort on the mean
        sim_fam_subset <- sim_fam %>% 
          filter(pred_vaxrG < quantile(pred_varG, 0.10) | pred_varG > quantile(pred_varG, 0.90)) %>%
          arrange(pred_mu)
        
        
        
        
        
        
        ## Data.frame of information to collect for each cycle and iteration
        tracked_data <- vector("list", n_cycles)
        




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



