## PopVar simulation to test the size of populations necessary to validate PopVar
# April 8, 2016

# Load packages
library(hypred)
library(PopVar)
library(rrBLUP)
library(boot)

# Set working directory
setwd("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/PopVar/POPSEQ_Pred_VarG/PopVar_validation_simulations/")

# Load some other functions
source("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Side Projects/Simulations/BarleySimGS-TPUpdate/Code/hypred_simulation_FUNCTIONS.R")

boot.cor <- function(data, boot.reps, CI = 0.95) {
  
  # Define a function for the correlation
  boot.cor <- function(input.data, i) {
    rep.data <- input.data[i,]
    return(cor(rep.data[,1], rep.data[,2]))
  }
  
  # Perform the bootstrapping
  boot.results <- boot(data = data, statistic = boot.cor, R = boot.reps)
  
  # Calculate a confidence interval
  CI.ind <- boot.reps * ((1 - CI) / 2)
  
  CI.upper <- sort(boot.results$t, decreasing = T)[CI.ind]
  CI.lower <- sort(boot.results$t, decreasing = F)[CI.ind]
  
  # Assemble list
  results.out <- list(r = boot.results$t0, r.sd.hat = sd(boot.results$t), CI = c(CI.lower, CI.upper))
  
  # Parse the results and return
  return(results.out)
}

# Load the CAP data
load("../../../../Side Projects/Simulations/BarleySimGS-TPUpdate/Files/Barley_CAP_simuation_starting_material.RData")

# Experimental design
##### User-controlled parameters #####
# Number of iterations
n.iterations = 500

mutation.rate.snp = 7e-8
mutation.rate.qtl = 7e-8

# Define the experimental design environments and reps
## In other words, how many environments and reps will be used to estimate the variance and superior progeny mean?
n.exp.env = 2
n.exp.reps = 2

# Define whether the tails will be selected prior to estimating variance
select.tails = TRUE
# Define the selection intensity
sel.intensity = 0.1

# Define the trait parameters
h2.vec <- c(0.3, 0.6)
n.QTL <- 100

# Define the population sizes to test
pop.sizes <- c(50, 100, 500, 1000)

# Define the error variance from which to draw error for the genetic map
marker.V_e = 0.002

# Select the appropriate haploid genome matrix
CAP.gametes <- CAP.gametes.0.03
sampled.markers <- sampled.markers.0.03

chr.len <- as.numeric(tapply(X = sampled.markers$pos, INDEX = sampled.markers$chrom, FUN = max))

# Define the initial genome
hv.genome <- make.genome( n.chr = 7, 
                          chr.len = chr.len, 
                          n.chr.snps = 100,
                          genetic.map = sampled.markers$pos)

# For loop over heritabilities
for (h2 in h2.vec) {
  
  # Sapply over iterations
  sapply(X = 1:n.iterations, FUN = function(i) {

    #### Define trait parameters ####
    hv.genome <- trait.architecture(genome = hv.genome,
                                    n.QTL = n.QTL, 
                                    qtl.ids = NULL, 
                                    qtl.dom.ids = NULL, 
                                    qtl.per.ids = NULL, 
                                    qtl.add.eff = "geometric", 
                                    qtl.dom.eff = NULL,
                                    keep.all.snps = FALSE)
    
    # Since more loci have been added to the genome, we need to add columns to the gamete matricies to make up for this
    TP.gametes <- CAP.gametes
    # Conver the gametes to genotypes
    TP.genos <- genotype.markers(gametes = TP.gametes, genome = hv.genome, include.QTL = F)
    
    ## Add error to the genetic map
    # Pull out the Morgan position of the marker loci
    pos.snps <- slot(hv.genome, "pos.snp")
    pos.qtl <- slot(hv.genome, "pos.add.qtl")$ID
    pos.markers <- data.frame(chr = rep(1:slot(hv.genome, "num.chr"), each = slot(hv.genome, "num.snp.chr")), pos = pos.snps[-pos.qtl])
    # Take out the marker names
    marker.names <- paste("M", 1:nrow(pos.markers), sep = "")
    
    # To add error, we will add a deviation from the known position of the markers based on the given error variance
    pos.markers.new <- pos.markers$pos + rnorm(n = nrow(pos.markers), 0, sqrt(marker.V_e))
    
    # Make sure no negative positions were given
    while(any(pos.markers.new < 0)) {
      pos.markers.new <- pos.markers$pos + rnorm(n = nrow(pos.markers), 0, sqrt(marker.V_e))
    }
    
    # Save the genetic map for use in PopVar
    popvar.map <- data.frame(marker.names = marker.names,
                             chr = pos.markers$chr,
                             cM.pos = pos.markers.new * 100)
    
    # Calculate the error variance between the known and estimated marker positions
    marker.ev <- var(pos.markers.new - pos.markers$pos)
    
    # Measure geno and pheno values
    TP.values <- evaluate.population(genome = hv.genome,
                                     gametes = TP.gametes,
                                     h2 = h2,
                                     n.env = n.exp.env,
                                     n.rep = n.exp.env)
    # Subset the phenotypic values
    TP.phenos <- TP.values$mean.pheno.values
  
    # Find the top 100 for use in crossing
    parent.lines <- select.population(TP.phenos, sel.intensity = 100, selection = "best")$lines.sel
    # Subset gametes
    parent.gametes <- subset.gametes(gametes = TP.gametes, line.names = parent.lines)
    
    # Reformat the data for PopVar
    TP.genos.popvar <- as.data.frame(cbind( c("", row.names(TP.genos)), rbind(marker.names,
                                                                              TP.genos)) )
    
    TP.phenos.popvar <- data.frame(genos = row.names(TP.phenos), sim.trait = TP.phenos)
  
    
    ##### Predict the genetic variance through PopVar
    # Create a random crossing block of 30 crosses
    C0.crossing.block <- make.crossing.block(parent1.lines = parent.lines,
                                             parent2.lines = parent.lines,
                                             n.crosses = 30,
                                             method = "random")
    
    # Use popvar to predict genetic variance
    C1.pop.predict <- pop.predict(G.in = TP.genos.popvar,
                                  y.in = TP.phenos.popvar,
                                  map.in = popvar.map,
                                  crossing.table = C0.crossing.block, 
                                  tail.p = sel.intensity,
                                  min.maf = 0,
                                  mkr.cutoff = 1,
                                  entry.cutoff = 1,
                                  impute = "pass",
                                  nSim = 25,
                                  nInd = 150,
                                  models = "rrBLUP",
                                  nCV.iter = 1,
                                  remove.dups = F)
    
    # Create a matrix of the Cross, predicted variance and mu_sp, and sd of those estimates
    columns <- c("Par1", "Par2", "pred.mu", "pred.varG", "pred.varG_sd", "mu.sp_high")
    C1.predictions <- C1.pop.predict$predictions[, columns]
    # Subset the genetic variance
    C1.predicted.V_g <- unlist(C1.predictions$pred.varG)
    
    # Create a list of comparisons and hypothesis testing
    # This will involve creating 30 populations of various sizes and correlating the observed
    ## and predicted genetic variance
    correlation.list <- list()
    
    # Iterate over the vector of population sizes
    sapply(X = pop.sizes, FUN = function(pop) {
    
      # Make the crosses and different population sizes
      ## Progeny will be inbred to the F_3
      C1.population <- make.population(hv.genome,
                                       named.parental.gametes = parent.gametes,
                                       crossing.block = C0.crossing.block,
                                       N = pop,
                                       cycle.number = 1,
                                       generations = 2,
                                       pop.type = "inbred",
                                       mutation.rate.snp = mutation.rate.snp,
                                       mutation.rate.qtl = mutation.rate.qtl)
      
      # Genotype the markers
      C1.genos <- genotype.markers(gametes = C1.population, genome = hv.genome, include.QTL = F)
      
      # First compare using the whole populations
      # Measure genotypic values in the actual population
      C1.values <- evaluate.population(genome = hv.genome,
                                       gametes = C1.population,
                                       h2 = h2, 
                                       V_e.scale = 8, 
                                       n.env = n.exp.env, 
                                       n.rep = n.exp.reps)
                                            
      
      # Calculate genetic variance within familes
      C1.observed.V_g <- tapply(X = C1.known.geno.values, INDEX = substring(text = row.names(C1.known.geno.values), first = 6, last = 7), FUN = var)
      
      # Correlate - use bootstrapping for CI
      boot.results <- boot.cor(data = cbind(C1.observed.V_g, C1.predicted.V_g), boot.reps = 1000, CI = 0.95)
      # KS test
      ks.test.results <- ks.test(x = C1.observed.V_g, y = C1.predicted.V_g)
      
      
      ### Test predictions after phenotyping trials
      # Phenotype
      
      
      
      
      
      
      # Add to the list
      correlation.list[as.character(pop)] <- whole.pop.cor
      
    } # Close the per-pop.size loop
    
    
    
    ### Now correlate based on using genomic prediction to measure GEBVs,
    # find the tails, grow those out, and measure genetic variance
    
    # Get GEBVs - SOMETHING MAY BE WRONG HERE
    C1.GEBVs <- make.predictions(pheno.train = TP.phenos, 
                                 geno.train = TP.genos, 
                                 geno.pred = C1.genos, 
                                 model = "RRBLUP")$GEBV
    
    # Reformat to data.frame
    C1.GEBVs <- data.frame(family = substring(text = row.names(C1.known.values), first = 6, last = 7),
                           GEBV = C1.GEBVs)
    
    # Find the top and bottom 10% in each family
    C1.selection <- tapply(X = C1.GEBVs$GEBV, INDEX = C1.GEBVs$family, FUN = function(family) {
      top <- sort
      bottom <- select.population(pheno.mat = family, sel.intensity = sel.intensity, selection = "worst")$lines.sel
      c(top,bottom)
    })
    
    C1.top <- select.population(pheno.mat = C1.GEBVs,
                                sel.intensity = sel.intensity,
                                selection = "best")
    C1.top.gametes <- subset.gametes(gametes = C1.population, line.names = C1.top$lines.sel)
    
    C1.bottom <- select.population(pheno.mat = C1.GEBVs,
                                   sel.intensity = sel.intensity,
                                   selection = "worst")
    C1.bottom.gametes <- subset.gametes(gametes = C1.population, line.names = C1.bottom$lines.sel)
    
    ## Evaluate the top and bottom lines in the field
    C1.selection.gametes <- rbind(C1.top.gametes, C1.bottom.gametes)
    C1.top.phenos <- measure.values(genome = hv.genome,
                                    gametes = )
    
    
    
    # If the tails are selected, follow this procedure:
    ## 1. Calculate GEBVs using the TP
    ## 2. Select the top 'sel.intensity' and bottom 'sel.intensity'
    ## 3. Measure the genotypic and phenotypic values of these sets
    ## 4. Calculate the mean using the tails
    ## 5. Measure the mean of the top tail
    ## 6. Use the superior progeny mean equation to determine the genetic variance
      
    # Apply the function over each family
    C1.estimation <- lapply(X = C1.population, FUN = function(family) {
      # Determine the true genetic variance
      family.values <- measure.values(hv.genome, family, h2, just.geno.values = T)
      family.V_g <- family.values$var.components$V_g
      # Determine the superior progeny genotypic values
      mu_sp.geno <- mean(sort(family.values$geno.values, decreasing = T)[1:(sel.intensity * pop)])
      
      # Genotype the progeny in the family
      family.genos <- genotype.markers(gametes = family, genome = hv.genome)
      # Make predictions
      family.predictions <- make.predictions(pheno.train = TP.phenos, geno.train = TP.genos, geno.pred = family.genos, model = "RRBLUP")

      
      # Select the top and bottom
      bottom.tail <- select.population(family.predictions$GEBV, sel.intensity = sel.intensity, selection = "worst")
      top.tail <- select.population(family.predictions$GEBV, sel.intensity = sel.intensity, selection = "best")
      
      # Subset the gametes for the top and bottom tails
      tails.gametes <- subset.gametes(gametes = family, line.names = c(bottom.tail$lines.sel, top.tail$lines.sel))
      
      # Perform the field experiment
      tails.values <- measure.values(hv.genome, gametes = tails.gametes, h2 = h2, n.env = n.exp.env, n.rep = n.exp.reps, model = "anova")
      
      # Calculate the superior and inferior progeny mean based on the GEBVs
      s_m1 <- mean(tails.values$pheno.mu[top.tail$lines.sel,])
      i_m1 <- mean(tails.values$pheno.mu[bottom.tail$lines.sel,])
      mu_m = s_m1 - i_m1
      
      # Calculate V_g using the standardized selection coefficient
      V_g.hat1 = ( (s_m1 - mu_m ) / sel.intensity )^2
      
      # Return the estimates of the mean, s_m, and V_g
      output.list <- list(true = list(V_g = family.V_g, mu_sp.geno = mu_sp.geno),
                          estimated = list(V_g.hat = V_g.hat1, mu_sp.hat = s_m1) )
      return(output.list)
    })
    
    # Pull out the true values for each cross
    C1.true <- t(sapply(X = C1.estimation, FUN = function(cross) data.frame(V_g = cross$true$V_g, mu_sp = cross$true$mu_sp.geno)))
    C1.estimated <- t(sapply(X = C1.estimation, FUN = function(cross) data.frame(V_g.hat = cross$estimated$V_g.hat, mu_sp.hat = cross$estimated$mu_sp.hat)))
    
    # Add to the comparison list
    comparison.list[[as.character(pop)]] <- list(true = C1.true, estimated = C1.estimated)
    
    # notify
    print(paste("Population size:", pop, "complete."))
    
    # Close the pop.size loop
  }
  
  C1.predictions
  
  
  
  
  
                           








