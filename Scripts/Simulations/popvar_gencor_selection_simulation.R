## PopVarVal simulations of correlated response to selection
## by taking advantage of the predicted genetic correlation
## 
## 
## Author: Jeff Neyhart
## Last modified: September 7, 2018
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



# Number of cores
n_cores <- 8 # Local machine
n_cores <- detectCores()


## Selection on the genetic correlation versus random selection or selection on the mean

## Fixed parameters
sim_pop_size <- 100
n_iter <- 10
n_env <- 3
n_rep <- 1
n_tp_select <- 50
n_crosses <- 50
k_sp <- 1.76

## Outline the parameters to perturb
trait1_h2_list <- trait2_h2_list <- seq(0.25, 1, by = 0.25)
nQTL_list <- 100
tp_size_list <- 300
gencor_list <- c(-0.5, 0.5)

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
  # ##
  
  # Create a results list
  results_out <- vector("list", nrow(core_df))
  
  # Iterate over the rows of the param_df
  for (i in seq_along(results_out)) {
    
    trait1_h2 <- core_df$trait1_h2[i]
    trait2_h2 <- core_df$trait2_h2[i]
    nQTL <- core_df$nQTL[i]
    maxQTL <- max(param_df$nQTL)
    tp_size <- core_df$tp_size[i]
    gencor <- core_df$gencor[i]
    probcor <- core_df$input[[i]]
    
    # Simulate QTL
    qtl_model <- replicate(n = 2, matrix(NA, ncol = 4, nrow = nQTL), simplify = FALSE)
    genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "geometric", max.qtl = maxQTL,
                                   corr = gencor, prob.corr = probcor)
    
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
    
    ## Predict the genotypic values for the TP, then select the best from the tp
    par_pop <- pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = tp1) %>%
      # Use indepdent culling levels
      select_pop(pop = ., intensity = sqrt(n_tp_select/nind(tp1)), index = c(1, 0), type = "phenotypic") %>%
      select_pop(pop = ., intensity = sqrt(n_tp_select/nind(tp1)), index = c(0, 1), type = "phenotypic") %>%
      pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = .)
    
    # Get the PGVs
    par_pop_pgv <- par_pop$pred_val %>%
      gather(trait, pgv, -ind) %>%
      mutate(ind = as.character(ind))
    
    # ## Plot the GVs of the tp and selections
    # tp1$geno_val %>% 
    #   mutate(selected = ind %in% indnames(par_pop)) %>% 
    #   qplot(x = trait1, y = trait2, data = ., color = selected)
    
    
    
    # Measure the genetic correlation in the TP
    tp_select_cor <- cor(par_pop$geno_val[,-1])[1,2]
    
    ## Create a crossing block with all of the possible non-reciprocal combinations of the TP
    crossing_block <- sim_crossing_block(parents = indnames(par_pop), n.crosses = choose(nind(par_pop), 2))
    # Calculate the mean PGVs for the combinations
    crossing_block_pgv <- crossing_block %>% 
      left_join(., par_pop_pgv, by = c("parent1" = "ind")) %>% 
      left_join(., par_pop_pgv, by = c("parent2" = "ind", "trait")) %>% 
      mutate(mean_pgv = (pgv.x + pgv.y) / 2) %>% 
      select(-pgv.x, -pgv.y) %>%
      spread(trait, mean_pgv) %>%
      ## Calculate an index
      mutate(index = trait1 + trait2)
    
    # Select the best
    # crossing_block_use <- crossing_block_pgv[order(crossing_block_pgv$index, decreasing = T)[1:cross_preselect],]
    crossing_block_use <- crossing_block_pgv
    # Pedigree for later family development
    ped <- sim_pedigree(n.ind = sim_pop_size, n.selfgen = Inf)
    
    
    ## Predict genetic variance and correlation
    pred_out <- pred_genvar(genome = genome1, pedigree = ped, training.pop = tp1, founder.pop = par_pop, 
                            crossing.block = select(crossing_block_use, contains("parent")))
    
    
    ## Select crosses that maximize the predicted correlation
    pred_out_select <- pred_out %>% distinct(parent1, parent2, pred_corG)
    pred_out_select <- pred_out_select[order(pred_out_select$pred_corG, decreasing = TRUE)[1:n_crosses],]
    cb_select_corG <- select(pred_out_select, contains("parent"))
    
    
    ## Now select random crosses
    cb_select_rand <- pred_out %>% 
      select(contains("parent")) %>% 
      sample_n(n_crosses)
    
    
    ## Combine crossing blocks
    cb_list <- list(corG = cb_select_corG, rand = cb_select_rand)
    
    ### Create all of the crosses
    cycle1_list <- cb_list %>%
      map(~sim_family_cb(genome = genome1, pedigree = ped, founder.pop = par_pop, crossing.block = .) %>%
            pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = .))
    
    ## Correlations of the whole population
    cycle1_cor <- map_df(cycle1_list, ~ cor(.$geno_val[,-1])[1,2]) %>%
      gather(type, corG, corG, rand)
    
    
    ## Select progeny on a range of selection intensities
    selections <- data_frame(intensity = seq(0.05, 0.95, by = 0.05))
    selections$summ <- selections$intensity %>% map(~{
      int <- .
      cycle1_list %>%
        map(~select_pop(pop = ., intensity = sqrt(int), index = c(1,0), type = "genomic") %>%
              select_pop(pop = ., intensity = sqrt(int), index = c(0, 1), type = "genomic")) %>%
        map("geno_val")
    })
    
    ## Summarize the mean and cor of each/both traits
    response <- selections %>%
      mutate(corG = map(summ, "corG"), rand = map(summ, "rand")) %>% 
      select(-summ) %>% gather(selection, out, corG, rand) %>% 
      unnest() %>% 
      group_by(intensity, selection) %>% 
      summarize(trait1_varG = var(trait1), trait2_varG = var(trait2), 
                trait1_mean = mean(trait1), trait2_mean = mean(trait2), 
                corG = cor(trait1, trait2)) %>%
      ungroup()
    
    
    ## How many of the crosses are the same?
    common_crosses <- combn(x = cb_list, m = 2, simplify = FALSE) %>% 
      map_df(~data.frame(selection1 = names(.)[1], selection2 = names(.)[2], common_cross = nrow(reduce(., intersect)),
                         stringsAsFactors = FALSE))
    
    # Other data
    cors <- data.frame(type = c("tp_base_cor", "tp_pheno_cor", "tp_select_cor"), corG = c(tp_cor, tp_pheno_cor, tp_select_cor),
                       stringsAsFactors = FALSE, row.names = NULL) %>%
      bind_rows(., cycle1_cor)
    
    ## Add the response and other results
    results_out[[i]] <- list(
      response = response,
      common_cross = common_crosses,
      correlations = cor
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
save_file <- file.path(result_dir, "popvar_gencor_selection_simulation_results.RData")
save("popvar_simulation_out", file = save_file)




# 
# 
# 
# 
# 
# ## Fixed parameters
# tp_size <- 300
# tp_select <- 100
# cross_preselect <- 150
# cross_select <- 25
# n_progeny <- 100
# k_sp <- 1.76
# k_sp1 <- 1.12
# 
# L <- 100
# n_iter <- 100
# n_env <- 3
# n_rep <- 1
# 
# 
# ## Outline the parameters to perturb
# trait1_h2_list <- c(0.5, 0.8)
# trait2_h2_list <- c(0.2, 0.5)
# gencor_list <- c(-0.75, -0.25, 0.25, 0.75)
# probcor_list <- data_frame(arch = c("pleio", "close_link", "loose_link"),
#                            input = list(cbind(0, 1), cbind(5, 1), cbind(30, 1) ))
# 
# # Create a data.frame of parameters
# param_df <- crossing(trait1_h2 = trait1_h2_list, trait2_h2 = trait2_h2_list, gencor = gencor_list, 
#                      probcor = probcor_list, iter = seq(n_iter))
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
#     gencor <- core_df$gencor[i]
#     probcor <- core_df$input[[i]]
#     
#     # Simulate QTL
#     qtl_model <- replicate(n = 2, matrix(NA, ncol = 4, nrow = L), simplify = FALSE)
#     genome1 <- sim_multi_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "geometric", max.qtl = L,
#                                    corr = gencor, prob.corr = probcor)
#     
#     ## Create the TP
#     # Sample from the genotypes
#     sample_genos <- s2_cap_genos[sort(sample(nrow(s2_cap_genos), size = tp_size)),]
#     tp <- create_pop(genome = genome1, geno = sample_genos, ignore.gen.model = FALSE)
#     # Measure the genetic correlation in the TP
#     tp_cor <- cor(tp$geno_val[,-1])[1,2]
#     
#     
#     # Phenotype
#     tp1 <- sim_phenoval(pop = tp, h2 = c(trait1_h2, trait2_h2), n.env = n_env, n.rep = n_rep)
#     # Measure the phenotypic correlation in the TP
#     tp_pheno_cor <- cor(tp1$pheno_val$pheno_mean[,-1])[1,2]
#     
#     ## Predict the genotypic values for the TP, then select the best from the tp
#     par_pop <- pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = tp1) %>%
#       # Use indepdent culling levels
#       select_pop(pop = ., intensity = sqrt(tp_select/nind(tp1)), index = c(1, 0), type = "genomic") %>%
#       select_pop(pop = ., intensity = sqrt(tp_select/nind(tp1)), index = c(0, 1), type = "genomic")
#     # Get the PGVs
#     par_pop_pgv <- par_pop$pred_val %>%
#       gather(trait, pgv, -ind) %>%
#       mutate(ind = as.character(ind))
#   
#     # ## Plot the GVs of the tp and selections
#     # tp1$geno_val %>% 
#     #   mutate(selected = ind %in% indnames(par_pop)) %>% 
#     #   qplot(x = trait1, y = trait2, data = ., color = selected)
#     
#     
#     
#     # Measure the genetic correlation in the TP
#     tp_select_cor <- cor(par_pop$geno_val[,-1])[1,2]
#     
#     ## Create a crossing block with all of the possible non-reciprocal combinations of the TP
#     crossing_block <- sim_crossing_block(parents = indnames(par_pop), n.crosses = choose(nind(par_pop), 2))
#     # Calculate the mean PGVs for the combinations
#     crossing_block_pgv <- crossing_block %>% 
#       left_join(., par_pop_pgv, by = c("parent1" = "ind")) %>% 
#       left_join(., par_pop_pgv, by = c("parent2" = "ind", "trait")) %>% 
#       mutate(mean_pgv = (pgv.x + pgv.y) / 2) %>% 
#       select(-pgv.x, -pgv.y) %>%
#       spread(trait, mean_pgv) %>%
#       ## Calculate an index
#       mutate(index = trait1 + trait2)
#     
#     # Select the best
#     crossing_block_use <- crossing_block_pgv[order(crossing_block_pgv$index, decreasing = T)[1:cross_preselect],]
#     # Pedigree for later family development
#     ped <- sim_pedigree(n.ind = n_progeny, n.selfgen = Inf)
#     
#     
#     ## Predict genetic variance and correlation
#     pred_out <- pred_genvar(genome = genome1, pedigree = ped, training.pop = tp1, founder.pop = par_pop, 
#                             crossing.block = select(crossing_block_use, contains("parent"))) %>%
#       mutate(pred_musp = pred_mu + (k_sp * sqrt(pred_varG))) %>%
#       group_by(parent1, parent2) %>%
#       # mutate(pred_muspC = rev(pred_mu + (k_sp * pred_corG * pred_corG))) %>%
#       mutate(pred_muspC = rev(pred_mu + (k_sp1 * pred_corG * pred_corG))) %>%
#       ungroup()
#     
#     
#     # ## Check the predictions of correlation versus PopVar
#     # # Convert items for PopVar
#     # # Pheno
#     # pheno_use <- tp1$pheno_val$pheno_mean
#     # geno_use <- genotype(genome1, tp1)
#     # map_use <- genome2$gen_model$trait1 %>%
#     #   select(qtl_name, chr, pos)
#     # 
#     # 
#     # # # Pass to PopVar
#     # # Convert genotypes into something useable for PopVar
#     # geno_use <- as.data.frame(cbind( c("", row.names(geno_use)), rbind(colnames(geno_use), geno_use)) )
#     # 
#     # pred_out_pv <- pop.predict(G.in = geno_use, y.in = pheno_use, map.in = map_use,
#     #                           crossing.table = select(crossing_block_use, contains("parent")), tail.p = 0.1, nInd = 150,
#     #                           min.maf = 0, mkr.cutoff = 1, entry.cutoff = 1, remove.dups = FALSE,
#     #                           nSim = 25, nCV.iter = 1, models = "rrBLUP", impute = "pass")
#     # 
#     # # Tidy
#     # tidy_pred_out <- pred_out_pv$predictions %>%
#     #   map(as.data.frame) %>%
#     #   map(~mutate_all(., unlist) %>% rename_at(vars(contains("cor")), ~"pred.corG")) %>%
#     #   list(., names(.)) %>%
#     #   pmap_df(~mutate(.x, trait = str_extract(.y, "trait[0-9]{1}"))) %>%
#     #   select(parent1 = Par1, parent2 = Par2, trait, pred.mu, mu.sp_high, pred.varG, contains("cor"), contains("high.resp")) %>%
#     #   unite(high_resp, c("high.resp_trait1", "high.resp_trait2")) %>% mutate(high_resp = str_replace(high_resp, "_NA", "") %>% parse_number())
#     # 
#     # ## Correlate predictions of genetic variance and correlation
#     # pred_out %>%
#     #   left_join(., tidy_pred_out) %>% 
#     #   group_by(trait) %>%
#     #   summarize(corG_acc = cor(pred_corG, pred.corG), 
#     #             varG_acc = cor(pred_varG, pred.varG),
#     #             muspC_acc = cor(pred_muspC, high_resp),
#     #             mu_acc = cor(pred_mu, pred.mu),
#     #             musp_acc = cor(pred_musp, mu.sp_high))
#     # 
#     # ## Very accurate!
#     
# 
#     # ## Get the expected correlation
#     # expected_out <- calc_exp_genvar(genome = genome1, pedigree = ped, founder.pop = par_pop, 
#     #                                 crossing.block = select(crossing_block_use, contains("parent"))) %>%
#     #   mutate(exp_musp = exp_mu + (pred_ksp * sqrt(exp_varG)))
#     
#       
#     
#     ## Select the crosses that maximize the sum of the correlated responses
#     ## Alternatively, select the cross that gives the largest sum of the trait1 mu_sp and 
#     ## predicted correlated reponse in trait2
#     pred_out_select <- pred_out %>%
#       mutate(pred_musp_muspC = pred_musp + pred_muspC) %>%
#       filter(trait == "trait1") %>%
#       distinct(parent1, parent2, pred_muspC, pred_musp_muspC)
#     # pred_out_select <- pred_out_select[order(pred_out_select$pred_muspC_index, decreasing = TRUE)[1:cross_select],]
#     pred_out_select <- pred_out_select[order(pred_out_select$pred_musp_muspC, decreasing = TRUE)[1:cross_select],]
#     cb_select_corG <- select(pred_out_select, contains("parent"))
#     
#     ## Now select crosses based on an index of the predicted mean PGV
#     pred_out_select_mean <- pred_out %>% 
#       group_by(parent1, parent2) %>%
#       summarize(index = sum(pred_mu)) %>%
#       ungroup()
#     pred_out_select_mean <- pred_out_select_mean[order(pred_out_select_mean$index, decreasing = T)[1:cross_select],]
#     cb_select_mean <- select(pred_out_select_mean, contains("parent"))
#     
#     ## Now select on an index of the superior progeny means (this effectively ignores the correlation)
#     ## Select the crosses that maximize the correlated superior progeny mean of trait2
#     pred_out_select_musp <- pred_out %>% 
#       group_by(parent1, parent2) %>%
#       summarize(index = sum(pred_musp)) %>%
#       ungroup()
#     pred_out_select_musp <- pred_out_select_musp[order(pred_out_select_musp$index, decreasing = TRUE)[1:cross_select],]
#     cb_select_musp <- select(pred_out_select_musp, contains("parent"))
#     
#     ## Now select random crosses
#     cb_select_rand <- pred_out %>% 
#       select(contains("parent")) %>% 
#       sample_n(cross_select)
#     
#     
#     ## Combine crossing blocks
#     cb_list <- list(corG = cb_select_corG, mean = cb_select_mean, musp = cb_select_musp, rand = cb_select_rand)
#     
#     ### Create all of the crosses
#     cycle1_list <- cb_list %>%
#       map(~sim_family_cb(genome = genome1, pedigree = ped, founder.pop = par_pop, crossing.block = .) %>%
#             pred_geno_val(genome = genome1, training.pop = tp1, candidate.pop = .))
#     
#     
#     ## Create a list of selections
#     ## Use independent culling levels
#     selections <- cycle1_list %>%
#       map(~select_pop(pop = ., intensity = sqrt(0.1), index = c(1,0), type = "genomic") %>%
#             select_pop(pop = ., intensity = sqrt(0.1), index = c(0, 1), type = "genomic")) %>%
#       map("geno_val")
#     
#     # ## Select progeny on a range of selection intensities
#     # selections <- data_frame(k_sp = seq(0.05, 1, by = 0.05))
#     # selections$summ <- selections$k_sp %>% map(~{
#     #   sel <- select_pop(pop = cycle1, intensity = ., index = c(1,1), type = "genomic")$geno_val[,-1]
#     #   # Return the mean of each trait and the correlation
#     #   as.data.frame(t(c(colMeans(sel), corG = cor(sel)[1,2])))
#     # })
#     
#     response <- selections %>%
#       map(~as.data.frame(t(c(colMeans(.[,-1]), corG = cor(.[,-1])[1,2])))) %>% 
#       list(., names(.)) %>% 
#       pmap_df(~mutate(.x, selection = .y)) %>%
#       # Add the base SD
#       gather(trait, mean_gv, contains("trait")) %>% 
#       left_join(., gather(sqrt(tp1$pheno_val$var_comp$V_G), trait, base_sdG), by = "trait") %>%
#       left_join(., data.frame(trait = names(par_pop$geno_val)[-1], par_mean_gv = colMeans(par_pop$geno_val[,-1]),
#                               stringsAsFactors = FALSE), by = "trait")
#     
#     
#     ## Summarize the response
#     response_summary <- response %>% 
#       left_join(., subset(response, selection == "mean", c(corG, mean_gv, trait)), by = c("trait")) %>% 
#       select(trait, selection, corG = corG.x, mean_gv = mean_gv.x, corG_mean = corG.y, mean_gv_mean = mean_gv.y, par_mean_gv, base_sdG) %>%
#       mutate(relative_response = (mean_gv - mean_gv_mean) / base_sdG,
#              stand_response = (mean_gv - par_mean_gv) / base_sdG)
#     
#     ## How many of the crosses are the same?
#     common_crosses <- combn(x = cb_list, m = 2, simplify = FALSE) %>% 
#       map_df(~data.frame(selection1 = names(.)[1], selection2 = names(.)[2], common_cross = nrow(reduce(., intersect)),
#                          stringsAsFactors = FALSE))
#     
#     # Other data
#     cors <- data.frame(type = c("tp_base_cor", "tp_pheno_cor", "tp_select_cor"), corG = c(tp_cor, tp_pheno_cor, tp_select_cor),
#                        stringsAsFactors = FALSE, row.names = NULL)
#     
#     ## Add the response and other results
#     results_out[[i]] <- list(
#       response = response_summary,
#       common_cross = common_crosses,
#       correlations = cors
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
# popvar_simulation_out <- bind_rows(simulation_out)
# 
# # Save
# save_file <- file.path(result_dir, "popvar_gencor_selection_simulation_results.RData")
# save("popvar_simulation_out", file = save_file)















