## PopVarValidation
## Calculation and analysis of genetic correlation
##
## 
## Author: Jeff Neyhart
## Last modified: August 29, 2018
## 

# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load the pbr library
library(pbr)
library(broom)
library(ggridges)
library(cowplot)
library(modelr)


# Load the PVV BLUEs
load(file.path(data_dir, "PVV_BLUE.RData"))
# Load the S2 BLUEs and filter
load(file.path(gdrive_dir, "BarleyLab/Breeding/PhenotypicData/Final/MasterPhenotypes/S2_tidy_BLUE.RData"))

vp_family_tomodel <- s2_tidy_BLUE %>%
  filter(trait %in% traits, line_name %in% c(pot_pars, exper), str_detect(trial, "PVV"))


boot_rep <- 1000
alpha <- 0.05
i <- 0.1


## First if the same line is measured in two different trials, calculate an environment mean
vp_family_tomodel1 <- vp_family_tomodel %>%
  filter(line_name %in% exper) %>%
  group_by(trait, environment) %>%
  do({
    df <- .
    
    if (n_distinct(df$trial) > 1) {
      
      fit <- lm(value ~ -1 + line_name + trial, data = df)
      fit_tidy <- tidy(fit) %>% filter(str_detect(term, "line_name")) %>% 
        mutate(line_name = str_replace(term, "line_name", "")) %>% 
        select(line_name, value = estimate, std.error)
      
      df %>% 
        mutate(trial = NA) %>% 
        distinct(trial, environment, location, year, trait) %>% 
        cbind(., fit_tidy)
      
    } else {
      
      df
    }
    
  }) %>% ungroup()


## Calculate genetic correlation via REML by fitting a model with G and GE

## Use the stage-one BLUEs to calculate correlation
vp_family_tomodel1 <- vp_family_tomodel1 %>% 
  mutate(family = str_extract(line_name, "4[0-9]{3}"))

## Iterate over pairs of traits
trait_pairs <- combn(x = traits, m = 2, simplify = F)

# What families were measured for both traits in the pair
vp_family_tomodel2 <- trait_pairs %>%
  map(~{
    trs <- .
    families <- vp_family_tomodel1 %>% 
      distinct(trait, family) %>% 
      group_by(family) %>%
      filter(trait %in% trs) %>% 
      filter(n() == length(trs))
    
    left_join(families, vp_family_tomodel1, by = c("trait", "family")) %>% 
      ungroup()
  })

## Fit a model per family

vp_family_corG1 <- vp_family_tomodel2 %>% 
  map_df(~{
    df <- .
    
    df %>%
      select(family, environment, line_name, trait, value) %>%
      spread(trait, value) %>%
      group_by(family) %>%
      do({
        df2 <- .
        # Bind the traits
        Y <- as.matrix(df2[-1:-3])
        fit2 <- sommer::mmer2(fixed = Y ~ environment, 
                              random = ~us(trait):line_name + us(trait):line_name:environment, 
                              data = df2, rcov = ~ us(traits):units, silent = T)
        
        ## Return variance components
        varcomp <- fit2$var.comp$line_name
        data.frame(trait1 = colnames(Y)[1], trait2 = colnames(Y)[2], correlation = varcomp[1,2] / prod(sqrt(diag(varcomp))),
                   row.names = NULL, stringsAsFactors = FALSE)
        
      })
    
  })

# Distributions
vp_family_corG1 %>%
  ggplot(aes(x = correlation)) +
  geom_histogram() +
  facet_grid(trait1 ~ trait2)


    
v## Use the stage-two BLUEs to calculate correlation
vp_family_tomodel3 <- vp_BLUE %>% 
  filter(line_name %in% exper) %>% 
  mutate(family = str_extract(line_name, "4[0-9]{3}"))


vp_family_tomodel3 <- trait_pairs %>%
  map(~{
    trs <- .
    families <- vp_family_tomodel3 %>% 
      distinct(trait, family) %>% 
      group_by(family) %>%
      filter(trait %in% trs) %>% 
      filter(n() == length(trs))
    
    left_join(families, vp_family_tomodel3, by = c("trait", "family")) %>% 
      ungroup()
  })

# Fit models per family
vp_family_corG2 <- vp_family_tomodel3 %>% 
  map_df(~{
    df <- .
    
    df %>%
      group_by(family) %>%
      do({
        df2 <- spread(., trait, value)
        
        # Bind the traits
        Y <- as.matrix(df2[-1:-2])
        fit2 <- sommer::mmer2(fixed = Y ~ 1, 
                              random = ~us(trait):line_name, 
                              data = df2, rcov = ~ us(traits):units, silent = T)
        
        ## Return variance components
        varcomp <- fit2$var.comp$line_name
        data.frame(trait1 = colnames(Y)[1], trait2 = colnames(Y)[2], correlation = varcomp[1,2] / prod(sqrt(diag(varcomp))),
                   row.names = NULL, stringsAsFactors = FALSE)
        
      })
    
  })
  




















