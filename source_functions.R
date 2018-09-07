## S2MET Functions
## 
## A script with useful functions used in the prediction analysis of the S2MET
## 



## Other/utility functions
# A function to assign cores to a data.frame
assign_cores <- function(df, n_core) {
  df$core <- sort(rep(seq(n_core), length.out = nrow(df)))
  return(df)
}


# A function to calculate heritability, BLUEs, and variance components from tidy phenotype data
summarize_pheno <- function(data, blue.model = c("lmer", "sommer")) {
  
  # Make sure necessary columns are in 'data'
  needed_cols <- c("trial", "environment", "location", "year", "line_name", "value", "std.error")
  stopifnot(all(needed_cols %in% names(data)))
  
  blue.model <- match.arg(blue.model)
  
  # If the number of trials/environment is > 1, fit a model to get the genotype mean
  # for a trait-environment combination
  mto_trial <- group_by(data, environment) %>% 
    summarize(n_trial = n_distinct(trial)) %>%
    filter(n_trial > 1)
  
  if (nrow(mto_trial) > 0) {
    
    env_mean <- data %>%
      filter(environment %in% mto_trial$environment) %>%
      group_by(environment) %>%
      do({
        data1 <- .
        fit1 <- lm(value ~ -1 + line_name + trial, data = data1)
        
        # Tidy
        tidy(fit1) %>% 
          select(term, estimate, std.error) %>% 
          filter(str_detect(term, "line_name")) %>% 
          mutate(term = str_replace(term, "line_name", "")) %>% 
          rename(line_name = term, value = estimate)
        
      })
    
    # Combine these results with the original data
    data1 <- bind_rows(
      data %>% filter(environment %in% mto_trial$environment) %>% distinct(environment, location, year, trait) %>% left_join(., env_mean, by = "environment"),
      data %>% filter(!environment %in% mto_trial$environment)
    ) %>%
      select(trial, names(.)) %>%
      arrange(environment)
    
  } else {
    data1 <- data
    
  }
  
  control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  wts <- data1$std.error^2
  
  # If the number of environments is < 2, drop relevant random effects
  if (n_distinct(data1$environment) < 2) {
    formula <- value ~ (1|line_name)
    exp <- "line_name / (line_name + (Residual / (n_e * n_r)))"
    
    ## Drop terms
    fit_noge <- lmer(formula = formula, data = data1, control = control, weights = wts)
    fit_nog <- lm(formula = value ~ 1, data = data1)
    
  } else {
    formula <- value ~ (1|line_name) + environment + (1|line_name:environment)
    exp <- "line_name / (line_name + (line_name:environment / n_e) + (Residual / (n_e * n_r)))"
    
    ## Drop terms
    fit_noge <- lmer(formula = value ~ (1|line_name) + environment, data = data1, control = control, weights = wts)
    fit_nog <- lmer(formula = value ~  environment + (1|line_name:environment), data = data1, control = control, weights = wts)
    
    
  }

  fit <- lmer(formula = formula, data = data1, control = control, weights = wts)
  
  plot_table <- xtabs(~line_name + environment, data1)
  
  # Get the harmonic mean of the number of environments / reps
  n_e <- plot_table %>% 
    ifelse(test = . > 1, 1, .) %>% 
    rowSums() %>% 
    harm_mean()
  
  n_r <- plot_table %>% 
    harm_mean()
  
  # Estimate heritability
  h2 <- herit(object = fit, exp = exp, n_r = n_r, n_e = n_e)
  
  

  # Calculate significance
  ge_sig <- lr_test(fit, fit_noge)
  g_sig <- lr_test(fit, fit_nog)
  sig_test <- bind_rows(g_sig, ge_sig) %>% 
    mutate(term = c("g", "ge")) %>% 
    select(term, names(.), -full_model)
  
  
  ## Split on whether to use lmer or sommer
  if (blue.model == "lmer") {
  
    ## Modify formula so line_name is fixed, then fit the model
    new_form <- tail(as.character(formula), 1) %>% 
      str_replace(string = ., pattern = "\\(1 \\| line_name\\)", "line_name") %>% str_c("value ~ -1 + ", .) %>% 
      as.formula()
    
    if (any(str_detect(new_form, "\\("))) {
      ## Now refit the model, but change genotype from random to fixed
      fit_blue <- lmer(formula = new_form, data = data1, control = control, weights = wts)
      
    } else {
      ## Now refit the model, but change genotype from random to fixed
      fit_blue <- lm(formula = new_form, data = data1)
      
    }
    
  
    
    # Tidy
    tidy_blue <- tidy(fit_blue) %>% 
      filter(str_detect(term, "line_name"), !str_detect(term, "sd")) %>%
      mutate(line_name = str_replace(term, "line_name", "")) %>% 
      select(line_name, value = estimate)
    
  } else if (blue.model == "sommer") {
    
    
    stopifnot(n_distinct(data$environment) > 1)
    
    ## Use sommer to calculate the genotype BLUEs
    mf <- model.frame(value ~ line_name + environment, data1)
    y <- model.response(mf)
    X <- model.matrix(~ -1 + line_name + environment, mf)
    Z <- model.matrix(~ -1 + line_name:environment, mf)
    K <- diag(ncol(Z))
    R <- diag(wts)
    
    fit_blue <- sommer::mmer(Y = y, X = X, Z = list(ge = list(Z = Z, K = diag(ncol(Z)))), R = list(unit = R))
    
    # Tidy
    tidy_blue <- fit_blue$beta.hat %>% 
      as.data.frame() %>% 
      rownames_to_column("term") %>% 
      rename(estimate = T1) %>%
      filter(str_detect(term, "line_name")) %>% 
      mutate(line_name = str_replace(term, "line_name", "")) %>% 
      select(line_name, value = estimate)
    
    
  }
  
  # Return all this nonsense
  data_frame(BLUE = list(tidy_blue), n_e = n_distinct(data$environment), h2 = list(h2), sig_test = list(sig_test))
  
}


# A function to calculate genetic variance
calc_varG <- function(data, method = c("lmer", "sommer")) {
  
  # Check the data input
  data <- droplevels(as.data.frame(data))
  method <- match.arg(method)
  
  # Check column names for the required columns
  needed_cols <- c("environment", "location", "year", "line_name", "value", "std.error", "family")
  stopifnot(all(needed_cols %in% names(data)))
  
  
  # Number of lines in the family
  n_lines <- n_distinct(data$line_name)
  n_env <- n_distinct(data$environment)
  
  plot_table <- xtabs(~line_name + environment, data)
  
  # Split based on the number of environments
  if (n_env > 1) {
  
    # Get the harmonic mean of the number of environments / reps
    n_e <- plot_table %>% 
      ifelse(test = . > 1, 1, .) %>% 
      rowSums() %>% 
      harm_mean()
    
    n_r <- plot_table %>% 
      harm_mean()
    
    wts <- data$std.error^2
    
    
    # Split flow based on method
    if (method == "lmer") {
      
      control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
      formula <- value ~ (1|line_name) + environment + (1|line_name:environment)
      
      fit <- lmer(formula = formula, data = data, control = control, weights = wts, contrasts = list(environment = "contr.sum"))
      
      # Estimate heritability
      h2 <- herit(object = fit, exp = "line_name / (line_name + (line_name:environment / n_e) + (Residual / (n_e * n_r)))",
                  n_e = n_e, n_r = n_r)
      
      ## Drop terms
      fit_noge <- lmer(formula = value ~ (1|line_name) + environment, data = data, control = control, weights = wts)
      fit_nog <- lmer(formula = value ~  environment + (1|line_name:environment), data = data, control = control, weights = wts)
      
      # Calculate significance
      ge_sig <- lr_test(fit, fit_noge)
      fam_sig <- lr_test(fit, fit_nog)
      sig_test <- bind_rows(fam_sig, ge_sig) %>% 
        mutate(full_model = c("family", "ge")) %>% 
        rename(term_red = full_model)
      
      family_mean <- fixef(fit)[[1]]
      
    } else if (method == "sommer") {
      
      # Create the model matrices
      mf <- model.frame(value ~ line_name + environment, data)
      y <- model.response(mf)
      X <- model.matrix(~ 1 + environment, mf, contrasts.arg = list(environment = "contr.sum"))
      
      Zg <- model.matrix(~ -1 + line_name, mf)
      Kg <- diag(ncol(Zg))
      Zge <- model.matrix(~ -1 + line_name:environment, mf)
      Kge <- diag(ncol(Zge))
      
      R <- solve(diag(wts))
      
      # Fit the model
      fit <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = Kg), ge = list(Z = Zge, K = Kge)))
      
      varG <- fit$var.comp$g[1]
      varGE <- fit$var.comp$ge[1]
      varR <- fit$var.comp$units[1]
      
      h2 <- varG / (varG + (varGE / n_e) + (varR / (n_e + n_r)))
      var_comp <- data_frame(source = c("line_name:environment", "line_name", "Residual"),
                             variance = c(varGE, varG, varR))
      
      h2 <- list(heritability = h2, var_comp = var_comp)
      
      
      
      ## Drop terms
      fit_noge <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = Kg)))
      fit_nog <- sommer::mmer(Y = y, X = X, Z = list(ge = list(Z = Zge, K = Kge)))
      
      # Calculate significance
      ge_sig <- data_frame(term_red = "ge", statistic = -2 * (fit_noge$LL - fit$LL)) %>%
        mutate(df = 1, p_value = pchisq(q = statistic, df = df, lower.tail = FALSE))
      fam_sig <- data_frame(term_red = "family", statistic = -2 * (fit_nog$LL - fit$LL)) %>%
        mutate(df = 1, p_value = pchisq(q = statistic, df = df, lower.tail = FALSE))
      sig_test <- bind_rows(fam_sig, ge_sig)
      
      family_mean <- fit$beta.hat[1]
      
    }
    
  } else {
    
    n_r <- plot_table %>% 
      harm_mean()
    
    wts <- data$std.error^2
    
    
    # Split flow based on method
    if (method == "lmer") {
      
      control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
      formula <- value ~ (1|line_name)
      
      fit <- lmer(formula = formula, data = data, control = control, weights = wts)
      
      # Estimate heritability
      h2 <- herit(object = fit, exp = "line_name / (line_name + (Residual / (n_r)))",
                  n_r = n_r)
      
      ## Drop terms
      fit_noge <- fit
      fit_nog <- lm(formula = value ~ 1, data = data)
      
      # Calculate significance
      ge_sig <- data_frame(term_red = "ge", statistic = -2 * (as.numeric(logLik(fit_noge)) - as.numeric(logLik(fit)))) %>%
        mutate(df = 1, p_value = pchisq(q = statistic, df = df, lower.tail = FALSE))
      fam_sig <- data_frame(term_red = "family", statistic = -2 * (as.numeric(logLik(fit_nog)) - as.numeric(logLik(fit)))) %>%
        mutate(df = 1, p_value = pchisq(q = statistic, df = df, lower.tail = FALSE))
      sig_test <- bind_rows(fam_sig, ge_sig)
      
      family_mean <- fixef(fit)[[1]]
      
    } else if (method == "sommer") {
      
      # Create the model matrices
      mf <- model.frame(value ~ line_name, data)
      y <- model.response(mf)
      X <- model.matrix(~ 1, mf, contrasts.arg = list(environment = "contr.sum"))
      
      Zg <- model.matrix(~ -1 + line_name, mf)
      Kg <- diag(ncol(Zg))
      
      R <- solve(diag(wts))
      
      # Fit the model
      fit <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Zg, K = Kg)))
      
      varG <- fit$var.comp$g[1]
      varR <- fit$var.comp$units[1]
      
      h2 <- varG / (varG + (varR / (n_r)))
      var_comp <- data_frame(source = c("line_name", "Residual"),
                             variance = c(varG, varR))
      
      h2 <- list(heritability = h2, var_comp = var_comp)
      
      
      
      ## Drop terms
      fit_nog <- sommer::mmer(Y = y, X = X, Z = list(g = list(Z = diag(length(y)), K = diag(length(y)))))

      # Calculate significance
      ge_sig <- data_frame(term_red = "ge", statistic = -2 * (fit$LL - fit$LL)) %>%
        mutate(df = 1, p_value = pchisq(q = statistic, df = df, lower.tail = FALSE))
      fam_sig <- data_frame(term_red = "family", statistic = -2 * (fit_nog$LL - fit$LL)) %>%
        mutate(df = 1, p_value = pchisq(q = statistic, df = df, lower.tail = FALSE))
      sig_test <- bind_rows(fam_sig, ge_sig)
      
      family_mean <- fit$beta.hat[1]
      
    }
    
    
  }
  
  # Return all this nonsense
  data_frame(family_mean = family_mean, fit = list(fit), h2 = list(h2), sig_test = list(sig_test))
  
}



# A function to calculate family mean and superior progeny mean
calc_mean <- function(data, i = 0.1) {
  
  # Check the data input
  data <- droplevels(as.data.frame(data))
  stopifnot(between(i, 0, 1))
  
  # Check column names for the required columns
  needed_cols <- c("environment", "location", "year", "line_name", "value", "std.error", "family")
  stopifnot(all(needed_cols %in% names(data)))
  
  # Convert some variables to factors
  data1 <- mutate_at(data, vars(line_name, environment), as.factor)
  
  
  # Number of lines in the family
  n_lines <- n_distinct(data1$line_name)
  n_env <- n_distinct(data1$environment)
  
  plot_table <- xtabs(~line_name + environment, data1)
  
  # Split based on the number of environments
  if (n_env > 1) {
    
    wts <- data1$std.error^2

    control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
    formula <- value ~ 1 + line_name + environment + (1|line_name:environment)
    
    fit <- lmer(formula = formula, data = data1, control = control, weights = wts, 
                contrasts = list(environment = "contr.sum", line_name = "contr.sum"))
    
    # Get the fixed effect estimates
    coefs <- tidy(fit) %>% 
      filter(str_detect(term, "line_name|Intercept"), group == "fixed") %>%
      select(term, estimate) %>%
      add_row(term = "last_line", estimate = -sum(tail(.$estimate, -1))) %>%
      mutate(term = c("family_mean", levels(data1$line_name)),
             mean = c(estimate[1], estimate[-1] + estimate[1]))

    
  } else {
    
    formula <- value ~ 1 + line_name
    
    fit <- lm(formula = formula, data = data1, contrasts = list(line_name = "contr.sum"))
    
    # Get the fixed effect estimates
    coefs <- tidy(fit) %>% 
      filter(str_detect(term, "line_name|Intercept")) %>%
      select(term, estimate) %>%
      add_row(term = "last_line", estimate = -sum(tail(.$estimate, -1))) %>%
      mutate(term = c("family_mean", levels(data1$line_name)),
             mean = c(estimate[1], estimate[-1] + estimate[1]))
    
  }
  
  geno_means <- subset(coefs, term != "family_mean", mean, drop = T)
  
  # Calculate the superior progeny mean
  mu_sp <- mean(geno_means[geno_means <= quantile(geno_means, i)])
    
  # Return all this nonsense
  data_frame(means = list(coefs), family_mean = coefs$mean[1], mu_sp = mu_sp)
  
}




# A function to return a tidy output from PopVar
tidy.popvar <- function(x) {
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
  
  

## Modification of the the pbsim function for calculating the expected genetic variance.
calc_exp_genvar1 <- function(genome, pedigree, founder.pop, crossing.block) {
  
  # Error handling
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")

  # Check the crossing block
  if (ncol(crossing.block) != 2) {
    stop("The crossing block should have two columns.")
  } else {
    crossing.block <- as.data.frame(crossing.block)
  }
  
  # founder.pop needs to be a pop object
  if (!inherits(founder.pop, "pop"))
    stop("The input 'founder.pop' must be of class 'pop'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = founder.pop$geno))
    stop("The geno did not pass. See warning for reason.")
  
  ## How many traits
  n_traits <- length(genome$gen_model)
  
  # If it is more than 2, error out
  stopifnot(n_traits <= 2)
  
  
  ## Calculate the expected genetic variance
  
  
  ## What are the expected allele frequencies in the population?
  ## Is there any backcrossing?
  mom_ped <- pedigree[pedigree$mom == 1,]
  dad_ped <- pedigree[pedigree$mom == 2,]
  
  mom_dist_gen <- length(unique(mom_ped$gen))
  dad_dist_gen <- length(unique(dad_ped$gen))
  
  max_bc_gen <- pmax(mom_dist_gen, dad_dist_gen) - 1
  
  # The expected frequency of the minor allele is 0.5 ^ n_bc_gen + 1
  exp_q <- 0.5^(max_bc_gen + 1)
  exp_p <- 1 - exp_q
  
  # Get the QTL information - drop unused levels
  qtl_info <- pull_qtl(genome, unique = FALSE)
  # Filter out QTL with no additive effect
  qtl_info <- droplevels(qtl_info[qtl_info$add_eff != 0,,drop = FALSE])
  # Split by trait
  qtl_info_split <- split(qtl_info, qtl_info$trait)
  
  
  ## Iterate over traits
  qtl_covariance <- lapply(X = qtl_info_split, FUN = function(trait_qtl) {
    
    # Get the map and genotypes of only the QTL
    qtl_geno <- pull_genotype(genome = genome, geno = founder.pop$geno, loci = trait_qtl$qtl_name) - 1
    
    ## Calculate the expected genetic variance and covariance of QTL
    trait_qtl_split <- split(trait_qtl, trait_qtl$chr)
    
    
    ## Calculate the expected covariance between QTL
    qtl_covar <- lapply(X = trait_qtl_split, FUN = function(qtl_chr) {
      d <- as.matrix(dist(qtl_chr$pos))
      
      # Calculate pairwise D (see Zhong and Jannink, 2007)
      # First convert cM to recombination fraction
      c <- qtl:::mf.h(d)
      D <- ((1 - (2 * c)) / (1 + (2 * c)))
      # # The diagonals are 0
      # diag(D) <- 0
      
      # Calculate the pairwise product of all QTL effects
      qtl_crossprod <- tcrossprod(qtl_chr$add_eff)
      dimnames(qtl_crossprod) <- list(qtl_chr$qtl_name, qtl_chr$qtl_name)
      
      # The covariance is the QTL effect product multiplied by the expected D
      qtl_crossprod * D
      
    })
    
    # Combine into a block diagonal, since the covariance between QTL on different chromosomes
    # is expected to be 0
    Cov <- .bdiag(qtl_covar)
    dimnames(Cov) <- list(trait_qtl$qtl_name, trait_qtl$qtl_name)
    # Return
    return(Cov)
    
  })
  
  ## Calculate the genetic covariance between QTL for different traits
  # Split by chromosome
  qtl_chr_split <- split(qtl_info, qtl_info$chr)
  
  
  if (n_traits > 1) {
    
    ## Iterate over chromosomes
    qtl_trait_covariance <- lapply(X = qtl_chr_split, FUN = function(chr_qtl) {
      
      # Split by trait
      trait_split <- split(chr_qtl, chr_qtl$trait)
      
      ## QTL names for each trait
      qtl_names <- lapply(X = trait_split, FUN = "[[", "qtl_name")
      qtl_pos <- lapply(X = trait_split, FUN = "[[", "pos")
      qtl_eff <- lapply(X = trait_split, FUN = function(q) as.matrix(q$add_eff))
      
      ## Calculate the pairwise distance
      d <- abs(outer(X = qtl_pos[[1]], Y = qtl_pos[[2]], FUN = `-`))
      # Calculate pairwise D (see Zhong and Jannink, 2007)
      # First convert cM to recombination fraction
      c <- qtl:::mf.h(d)
      D <- ((1 - (2 * c)) / (1 + (2 * c)))
      
      # Product of QTL effects
      qtl_crossprod <- tcrossprod(qtl_eff[[1]], qtl_eff[[2]])
      dimnames(qtl_crossprod) <- qtl_names
      
      # The covariance is the QTL effect product multiplied by the expected D
      qtl_crossprod * D
      
    })
    
    # Combine into a block diagonal, since the covariance between QTL on different chromosomes
    # is expected to be 0
    traitCov <- .bdiag(qtl_trait_covariance)
    dimnames(traitCov) <- list(qtl_info_split[[1]]$qtl_name, qtl_info_split[[2]]$qtl_name)
    
  } else {
    traitCov <- NULL
    
  }
  
  
  ## Now we iterate over the parent pairs to determine the QTL that are segregating
  
  # Replicate the crossing block 
  
  ## Add columns to the crossing.block for exp mu and exp varG
  crossing_block <- crossing(crossing.block, trait = paste0("trait", seq(length(genome$gen_model))))
  exp_mu <- exp_varG <- exp_corG <- list()
  var_and_covar <- list()
  
  # Iterate over the crossing block
  for (j in seq(nrow(crossing.block))) {
    
    pars <- as.character(crossing.block[j,1:2])
    
    ## Pull out the qtl genotypes for each trait
    qtl_names <- lapply(X = qtl_info_split, FUN = "[[", "qtl_name")
    qtl_geno <- lapply(X = qtl_names, pull_genotype, genome = genome, geno = founder.pop$geno)
    
    ## Get a list of the polymorphic QTL
    poly_qtl_list <- lapply(X = qtl_geno, FUN = function(tr_qtl) {
      
      # Subset the parents
      par_qtl_geno <- tr_qtl[pars,,drop = FALSE] - 1
      qtl_means <- colMeans(par_qtl_geno)
      poly_qtl <- names(qtl_means)[qtl_means == 0]
      # Get the marker states for parent 1 - use this to determine phase
      par_qtl_geno[1, poly_qtl, drop = F]
      
    })
    
    # Iterate over the traits
    trait_var1 <- mapply(poly_qtl_list, qtl_covariance, FUN = function(poly_mat, qtl_cov) {
      
      poly_qtl <- colnames(poly_mat)
      exp_covar1 <- as.matrix(qtl_cov)[poly_qtl, poly_qtl] * crossprod(poly_mat)
      
      # The expected variance is the sum of the variances at the polymorphic QTL, plus 2 times
      # the expected covariance between all polymorphic QTL
      
      # The sum of the variance of polymorphic QTL is the sum of the diagonal elements
      qtl_var <- sum(diag(exp_covar1))
      # The sum of the off-diagonal elements is 2 times the sum of the covariance
      qtl_covar <- sum(exp_covar1) - qtl_var
      
      # Return these two elements
      c(var = qtl_var, covar = qtl_covar)
    })
    
    # Sum
    trait_var <- colSums(trait_var1)
    
    if (!is.null(traitCov)) {
      
      poly_qtl_names <- lapply(poly_qtl_list, colnames)
      ## Calculate the expected covariance
      trait_cov <- as.matrix(traitCov)[poly_qtl_names[[1]], poly_qtl_names[[2]]] * crossprod(poly_qtl_list[[1]], poly_qtl_list[[2]])
      # The covariance is 2 times the sum of this matrix (it is not diagonal)
      trait_cov1 <- sum(trait_cov)
      
      # The expected correlation is calculated using the expected sd and expected cov
      exp_corG_j <- trait_cov1 / prod(sqrt(trait_var)) 
      exp_corG[[j]] <- rep(exp_corG_j, 2)
      
    }
    
    # The expected mu is simply the mean of the genotypic values of the two parents
    exp_mu_j <- colMeans(founder.pop$geno_val[founder.pop$geno_val$ind %in% pars,-1,drop = F])
    
    ## Add to the lists
    exp_mu[[j]] <- exp_mu_j
    exp_varG[[j]] <- trait_var
    var_and_covar[[j]] <- trait_var1
    
  }
  
  ## Add the variances and means to the crossing block
  crossing_block$exp_mu <- unlist(exp_mu)
  crossing_block$exp_varG <- unlist(exp_varG)  
  crossing_block$exp_corG <- unlist(exp_corG)
  
  # Return the crossing block
  return(list(crossing_block = crossing_block, var_covar = var_and_covar))
  
}










## A function to calculate the expected genetic variance in a bi-parental population
calc_exp_genvar2 <- function(G.in = G_in, y.in = phenos_use, map.in = map_use, crossing.table = core_df, tail.p = 0.1,
                             models = "rrBLUP", impute = "pass") {
  
  # Error handling
  if (!inherits(genome, "genome"))
    stop("The input 'genome' must be of class 'genome.'")
  
  # Check the pedigree
  if (!check_pedigree(pedigree, ignore_sex = TRUE))
    stop("The pedigree is not formatted correctly.")
  
  # Check the crossing block
  if (ncol(crossing.block) != 2) {
    stop("The crossing block should have two columns.")
  } else {
    crossing.block <- as.data.frame(crossing.block)
  }
  
  # founder.pop needs to be a pop object
  if (!inherits(founder.pop, "pop"))
    stop("The input 'founder.pop' must be of class 'pop'")
  
  # Check the genome and geno
  if (!check_geno(genome = genome, geno = founder.pop$geno))
    stop("The geno did not pass. See warning for reason.")
  
  ## How many traits
  n_traits <- length(genome$gen_model)
  
  # If it is more than 2, error out
  stopifnot(n_traits <= 2)
  
  
  ## Calculate the expected genetic variance
  
  
  ## What are the expected allele frequencies in the population?
  ## Is there any backcrossing?
  mom_ped <- pedigree[pedigree$mom == 1,]
  dad_ped <- pedigree[pedigree$mom == 2,]
  
  mom_dist_gen <- length(unique(mom_ped$gen))
  dad_dist_gen <- length(unique(dad_ped$gen))
  
  max_bc_gen <- pmax(mom_dist_gen, dad_dist_gen) - 1
  
  # The expected frequency of the minor allele is 0.5 ^ n_bc_gen + 1
  exp_q <- 0.5^(max_bc_gen + 1)
  exp_p <- 1 - exp_q
  
  # Get the QTL information - drop unused levels
  qtl_info <- pull_qtl(genome, unique = FALSE)
  # Filter out QTL with no additive effect
  qtl_info <- droplevels(qtl_info[qtl_info$add_eff != 0,,drop = FALSE])
  # Split by trait
  qtl_info_split <- split(qtl_info, qtl_info$trait)
  
  
  ## Iterate over traits
  qtl_covariance <- lapply(X = qtl_info_split, FUN = function(trait_qtl) {
    
    # Get the map and genotypes of only the QTL
    qtl_geno <- pull_genotype(genome = genome, geno = founder.pop$geno, loci = trait_qtl$qtl_name) - 1
    
    ## Calculate the expected genetic variance and covariance of QTL
    trait_qtl_split <- split(trait_qtl, trait_qtl$chr)
    
    
    ## Calculate the expected covariance between QTL
    qtl_covar <- lapply(X = trait_qtl_split, FUN = function(qtl_chr) {
      d <- as.matrix(dist(qtl_chr$pos))
      
      # Calculate pairwise D (see Zhong and Jannink, 2007)
      # First convert cM to recombination fraction
      c <- qtl:::mf.h(d)
      D <- ((1 - (2 * c)) / (1 + (2 * c)))
      # # The diagonals are 0
      # diag(D) <- 0
      
      # Calculate the pairwise product of all QTL effects
      qtl_crossprod <- tcrossprod(qtl_chr$add_eff)
      dimnames(qtl_crossprod) <- list(qtl_chr$qtl_name, qtl_chr$qtl_name)
      
      # The covariance is the QTL effect product multiplied by the expected D
      qtl_crossprod * D
      
    })
    
    # Combine into a block diagonal, since the covariance between QTL on different chromosomes
    # is expected to be 0
    Cov <- .bdiag(qtl_covar)
    dimnames(Cov) <- list(trait_qtl$qtl_name, trait_qtl$qtl_name)
    # Return
    return(Cov)
    
  })
  
  ## Calculate the genetic covariance between QTL for different traits
  # Split by chromosome
  qtl_chr_split <- split(qtl_info, qtl_info$chr)
  
  
  if (n_traits > 1) {
    
    ## Iterate over chromosomes
    qtl_trait_covariance <- lapply(X = qtl_chr_split, FUN = function(chr_qtl) {
      
      # Split by trait
      trait_split <- split(chr_qtl, chr_qtl$trait)
      
      ## QTL names for each trait
      qtl_names <- lapply(X = trait_split, FUN = "[[", "qtl_name")
      qtl_pos <- lapply(X = trait_split, FUN = "[[", "pos")
      qtl_eff <- lapply(X = trait_split, FUN = function(q) as.matrix(q$add_eff))
      
      ## Calculate the pairwise distance
      d <- abs(outer(X = qtl_pos[[1]], Y = qtl_pos[[2]], FUN = `-`))
      # Calculate pairwise D (see Zhong and Jannink, 2007)
      # First convert cM to recombination fraction
      c <- qtl:::mf.h(d)
      D <- ((1 - (2 * c)) / (1 + (2 * c)))
      
      # Product of QTL effects
      qtl_crossprod <- tcrossprod(qtl_eff[[1]], qtl_eff[[2]])
      dimnames(qtl_crossprod) <- qtl_names
      
      # The covariance is the QTL effect product multiplied by the expected D
      qtl_crossprod * D
      
    })
    
    # Combine into a block diagonal, since the covariance between QTL on different chromosomes
    # is expected to be 0
    traitCov <- .bdiag(qtl_trait_covariance)
    dimnames(traitCov) <- list(qtl_info_split[[1]]$qtl_name, qtl_info_split[[2]]$qtl_name)
    
  } else {
    traitCov <- NULL
    
  }
  
  
  ## Now we iterate over the parent pairs to determine the QTL that are segregating
  
  # Replicate the crossing block 
  
  ## Add columns to the crossing.block for exp mu and exp varG
  crossing_block <- crossing(crossing.block, trait = paste0("trait", seq(length(genome$gen_model))))
  exp_mu <- list()
  exp_varG <- list()
  exp_corG <- list()
  
  # Iterate over the crossing block
  for (j in seq(nrow(crossing.block))) {
    
    pars <- as.character(crossing.block[j,1:2])
    
    ## Pull out the qtl genotypes for each trait
    qtl_names <- lapply(X = qtl_info_split, FUN = "[[", "qtl_name")
    qtl_geno <- lapply(X = qtl_names, pull_genotype, genome = genome, geno = founder.pop$geno)
    
    ## Get a list of the polymorphic QTL
    poly_qtl_list <- lapply(X = qtl_geno, FUN = function(tr_qtl) {
      
      # Subset the parents
      par_qtl_geno <- tr_qtl[pars,,drop = FALSE] - 1
      qtl_means <- colMeans(par_qtl_geno)
      poly_qtl <- names(qtl_means)[qtl_means == 0]
      # Get the marker states for parent 1 - use this to determine phase
      par_qtl_geno[1, poly_qtl, drop = F]
      
    })
    
    # Iterate over the traits
    trait_var <- mapply(poly_qtl_list, qtl_covariance, FUN = function(poly_mat, qtl_cov) {
      
      poly_qtl <- colnames(poly_mat)
      exp_covar1 <- as.matrix(qtl_cov)[poly_qtl, poly_qtl] * crossprod(poly_mat)
      
      # The expected variance is the sum of the variances at the polymorphic QTL, plus 2 times
      # the expected covariance between all polymorphic QTL
      # One can simply sum up the elements in the covariance matrix
      sum(exp_covar1)
      
    })
    
    if (!is.null(traitCov)) {
      
      poly_qtl_names <- lapply(poly_qtl_list, colnames)
      ## Calculate the expected covariance
      trait_cov <- as.matrix(traitCov)[poly_qtl_names[[1]], poly_qtl_names[[2]]] * crossprod(poly_qtl_list[[1]], poly_qtl_list[[2]])
      # The covariance is 2 times the sum of this matrix (it is not diagonal)
      trait_cov1 <- sum(trait_cov)
      
      # The expected correlation is calculated using the expected sd and expected cov
      exp_corG_j <- trait_cov1 / prod(sqrt(trait_var)) 
      exp_corG[[j]] <- rep(exp_corG_j, 2)
      
    }
    
    # The expected mu is simply the mean of the genotypic values of the two parents
    exp_mu_j <- colMeans(founder.pop$geno_val[founder.pop$geno_val$ind %in% pars,-1,drop = F])
    
    ## Add to the lists
    exp_mu[[j]] <- exp_mu_j
    exp_varG[[j]] <- trait_var
    
  }
  
  ## Add the variances and means to the crossing block
  crossing_block$exp_mu <- unlist(exp_mu)
  crossing_block$exp_varG <- unlist(exp_varG)  
  crossing_block$exp_corG <- unlist(exp_corG)
  
  # Return the crossing block
  return(crossing_block)
  
}


    
 