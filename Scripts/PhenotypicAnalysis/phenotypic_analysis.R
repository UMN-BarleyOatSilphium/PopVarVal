## PopVarValidation
## Phenotypic analysis
##
## This script will run statistical analyses of the phenotypic data from the PopVarVal project. This
## will include:
## 1. Phenotypic analysis of training population data
## 2. Phenotypic analysis of validation family data
## 
## Author: Jeff Neyhart
## Last modified: August 19, 2018
## 

# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load the pbr library
library(pbr)
library(broom)
library(ggridges)
library(cowplot)


## Load the S2 BLUEs
load(file.path(gdrive_dir, "BarleyLab/Breeding/PhenotypicData/Final/MasterPhenotypes/S2_tidy_BLUE.RData"))




### Phenotypic analysis of training population data

# First gather data that would have been used to make the predictions
tp_prediction_tomodel <- s2_tidy_BLUE %>% 
  filter(trait %in% traits, line_name %in% tp, year %in% 2014:2015,
         location %in% c("STP", "CRM")) %>%
  arrange(year, location, trial, trait, line_name)

# Now gather data on the training population that would be relevant to the predictions, regardless of year
tp_relevant_tomodel <- s2_tidy_BLUE %>%
  filter(trait %in% traits, line_name %in% tp, year %in% 2014:2017, 
         location %in% c("STP", "CRM", "FND", "BCW"))


## Run models
tp_analysis <- tp_prediction_tomodel %>%
  group_by(trait) %>%
  do({
    df <- .
    print(unique(df$trait))
    summarize_pheno(data = df)
  })


## Look at variance components and heritability
(g_tp_prediction_h2 <- tp_analysis %>% 
  mutate(h2 = map_dbl(h2, "heritability")) %>% 
  qplot(x = trait, y = h2, geom = "col", fill = "blue", data = .) +
  geom_text(aes(label = str_c("Envs: ", n_e)), vjust = 2) +
  geom_text(aes(label = str_c("h2: ", round(h2, 3))), vjust = 4) + 
  scale_fill_discrete(guide = FALSE))

ggsave(filename = "tp_prediction_h2.jpg", plot = g_tp_prediction_h2, path = fig_dir, height = 5, width = 5, dpi = 1000)
  

tp_var_prop <- tp_analysis %>% 
  mutate(varcomp = map(h2, "var_comp")) %>%
  unnest(varcomp) %>% 
  group_by(trait) %>% 
  mutate(var_prop = variance / sum(variance)) 

# trait       source                variance  var_prop
# 1 FHBSeverity line_name:environment 49.2     0.825    
# 2 FHBSeverity line_name             10.4     0.175    
# 3 FHBSeverity Residual               0.00456 0.0000764
# 4 HeadingDate line_name:environment  0.725   0.0665   
# 5 HeadingDate line_name              9.20    0.844    
# 6 HeadingDate Residual               0.972   0.0892   
# 7 PlantHeight line_name:environment 19.5     0.642    
# 8 PlantHeight line_name             10.9     0.358    
# 9 PlantHeight Residual               0.00245 0.0000808

g_tp_prediction_varprop <- tp_var_prop %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col(position = "dodge")

ggsave(filename = "tp_prediction_varprop.jpg", plot = g_tp_prediction_varprop, path = fig_dir, height = 5, width = 5, dpi = 1000)



tp_analysis %>%
  unnest(sig_test) %>%
  mutate(annotate = case_when(p_value <= 0.01 ~ "***", p_value <= 0.05 ~ "**", p_value <= 0.1 ~ "*", TRUE ~ ""))

# trait         n_e term     df statistic   p_value
# 1 FHBSeverity     4 g         1     26.5  2.59e-  7 ***     
# 2 FHBSeverity     4 ge        1     75.8  3.13e- 18 ***     
# 3 HeadingDate     4 g         1    721.   8.07e-159 ***     
# 4 HeadingDate     4 ge        1      5.87 1.54e-  2 **      
# 5 PlantHeight     2 g         1     24.9  6.00e-  7 ***     
# 6 PlantHeight     2 ge        1      4.67 3.07e-  2 ** 

## G and GxE are significant for all traits, though more so for FHB severity.

# Unnest the blues
tp_prediction_BLUE <- tp_analysis %>%
  unnest(BLUE) %>%
  ungroup() %>%
  select(-n_e)



## Do the same thing for the relevant tp data

## Run models
tp_analysis <- tp_relevant_tomodel %>%
  group_by(trait) %>%
  do({
    df <- .
    print(unique(df$trait))
    summarize_pheno(data = df)
  })


## Look at variance components and heritability
tp_analysis %>% select(trait, h2, n_e) %>% mutate(h2 = map_dbl(h2, "heritability"))

# trait          h2   n_e
# 1 FHBSeverity 0.446     4
# 2 HeadingDate 0.972    10
# 3 PlantHeight 0.746     9

(g_tp_relevant_h2 <- tp_analysis %>% 
    mutate(h2 = map_dbl(h2, "heritability")) %>% 
    qplot(x = trait, y = h2, geom = "col", fill = "blue", data = .) +
    geom_text(aes(label = str_c("Envs: ", n_e)), vjust = 2) +
    geom_text(aes(label = str_c("h2: ", round(h2, 3))), vjust = 4) + 
    scale_fill_discrete(guide = FALSE))

ggsave(filename = "tp_relevant_h2.jpg", plot = g_tp_relevant_h2, path = fig_dir, height = 5, width = 5, dpi = 1000)

tp_var_prop <- tp_analysis %>% 
  mutate(varcomp = map(h2, "var_comp")) %>%
  unnest(varcomp) %>% 
  group_by(trait) %>% 
  mutate(var_prop = variance / sum(variance)) 

# trait         n_e source                   variance    var_prop
# 1 FHBSeverity     4 line_name:environment 49.2        0.825      
# 2 FHBSeverity     4 line_name             10.4        0.175      
# 3 FHBSeverity     4 Residual               0.00456    0.0000764  
# 4 HeadingDate    10 line_name:environment  2.98       0.222      
# 5 HeadingDate    10 line_name             10.4        0.778      
# 6 HeadingDate    10 Residual               0.00000856 0.000000638
# 7 PlantHeight     9 line_name:environment 22.2        0.752      
# 8 PlantHeight     9 line_name              7.33       0.248      
# 9 PlantHeight     9 Residual               0.000507   0.0000172


g_tp_relevant_varprop <- tp_var_prop %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col(position = "dodge")

ggsave(filename = "tp_relevant_varprop.jpg", plot = g_tp_relevant_varprop, path = fig_dir, height = 5, width = 5, dpi = 1000)



tp_analysis %>%
  unnest(sig_test)  %>%
  mutate(annotate = case_when(p_value <= 0.01 ~ "***", p_value <= 0.05 ~ "**", p_value <= 0.1 ~ "*", TRUE ~ ""))

## G and GxE are significant for all traits, though more so for FHB severity.

# trait         n_e term     df statistic   p_value annotate
# 1 FHBSeverity     4 g         1      26.5 2.59e-  7 ***     
# 2 FHBSeverity     4 ge        1      75.8 3.13e- 18 ***     
# 3 HeadingDate    10 g         1    2066.  0.        ***     
# 4 HeadingDate    10 ge        1     553.  3.27e-122 ***     
# 5 PlantHeight     9 g         1     215.  1.30e- 48 ***     
# 6 PlantHeight     9 ge        1     538.  5.42e-119 *** 



# Unnest the blues
tp_relevant_BLUE <- tp_analysis %>%
  unnest(BLUE) %>%
  ungroup() %>%
  select(-n_e)





### Phenotypic analysis of validation family data
# Gather data

vp_family_tomodel <- s2_tidy_BLUE %>%
  filter(trait %in% traits, line_name %in% c(pot_pars, exper), str_detect(trial, "PVV"))


## Run models
vp_analysis <- vp_family_tomodel %>%
  group_by(trait) %>%
  do({
    df <- .
    print(unique(df$trait))
    summarize_pheno(data = df, blue.model = "sommer")
  })


## Look at variance components and heritability
vp_analysis %>% 
  select(trait, h2, n_e) %>% 
  mutate(h2 = map_dbl(h2, "heritability"))

# trait             h2   n_e
# 1 FHBSeverity 0.114     4
# 2 HeadingDate 0.784     4
# 3 PlantHeight 0.740     4

# Heritability is low for FHB, but adequate for HD and PH.

(g_vp_h2 <- vp_analysis %>% 
    mutate(h2 = map_dbl(h2, "heritability")) %>% 
    qplot(x = trait, y = h2, geom = "col", fill = "blue", data = .) +
    geom_text(aes(label = str_c("Envs: ", n_e)), vjust = 2) +
    geom_text(aes(label = str_c("h2: ", round(h2, 3))), vjust = 4) + 
    scale_fill_discrete(guide = FALSE))

ggsave(filename = "vp_h2.jpg", plot = g_vp_h2, path = fig_dir, height = 5, width = 5, dpi = 1000)


vp_var_prop <- vp_analysis %>% 
  mutate(varcomp = map(h2, "var_comp")) %>%
  unnest(varcomp) %>% 
  group_by(trait) %>% 
  mutate(var_prop = variance / sum(variance)) 

# trait         n_e source                 variance  var_prop
# 1 FHBSeverity     4 line_name:environment 137.      0.959    
# 2 FHBSeverity     4 line_name               5.83    0.0408   
# 3 FHBSeverity     4 Residual                0.00805 0.0000562
# 4 HeadingDate     4 line_name:environment   4.43    0.309    
# 5 HeadingDate     4 line_name               8.38    0.584    
# 6 HeadingDate     4 Residual                1.53    0.107    
# 7 PlantHeight     4 line_name:environment  19.5     0.465    
# 8 PlantHeight     4 line_name              22.5     0.535    
# 9 PlantHeight     4 Residual                0.00377 0.0000896


g_vp_varprop <- vp_var_prop %>% 
  ggplot(aes(x = trait, y = var_prop, fill = source)) + 
  geom_col(position = "dodge")

ggsave(filename = "vp_varprop.jpg", plot = g_vp_varprop, path = fig_dir, height = 5, width = 5, dpi = 1000)


vp_analysis %>%
  unnest(sig_test)  %>%
  mutate(annotate = case_when(p_value <= 0.01 ~ "***", p_value <= 0.05 ~ "**", p_value <= 0.1 ~ "*", TRUE ~ ""))

## G and GxE are significant for all traits

# trait         n_e term     df statistic   p_value annotate
# 1 FHBSeverity     4 g         1      10.4 1.28e-  3 ***
# 2 FHBSeverity     4 ge        1    1775.  0.        ***
# 3 HeadingDate     4 g         1    3499.  0.        ***
# 4 HeadingDate     4 ge        1      41.9 9.38e- 11 ***
# 5 PlantHeight     4 g         1    1698.  0.        ***
# 6 PlantHeight     4 ge        1     558.  2.40e-123 ***
# 


# Unnest the blues
vp_BLUE <- vp_analysis %>%
  unnest(BLUE) %>%
  ungroup() %>%
  select(-n_e)

## For FHB severity, try estimating variance components separately by location
vp_analysis_FHB <- vp_family_tomodel %>%
  filter(trait == "FHBSeverity") %>%
  group_by(location) %>%
  do(summarize_pheno(data = ., blue.model = "sommer"))

# Look at heritability
vp_analysis_FHB %>% 
  mutate(h2 = map_dbl(h2, "heritability")) %>%
  select(location, n_e, h2)

# location   n_e     h2
# 1 CRM          2 0.267 
# 2 STP          2 0.0417


# Is GE still significant when analyzing locations independently?
vp_analysis_FHB %>% unnest(sig_test)

# location   n_e term     df statistic   p_value
# 1 CRM          2 g         1    33.2   8.46e-  9
# 2 CRM          2 ge        1   332.    4.26e- 74
# 3 STP          2 g         1     0.755 3.85e-  1
# 4 STP          2 ge        1   850.    7.16e-187

## Yes


# Extract the BLUEs
vp_FHB_BLUE <- vp_analysis_FHB %>%
  unnest(BLUE) %>% 
  select(-n_e)


## Save the BLUEs
save("tp_prediction_BLUE", "tp_relevant_BLUE", "vp_BLUE", "vp_FHB_BLUE", file = file.path(data_dir, "PVV_BLUE.RData"))

# Load
load(file.path(data_dir, "PVV_BLUE.RData"))


## TP mean and range
tp_prediction_BLUE %>% 
  group_by(trait) %>% 
  summarize_at(vars(value), funs(min, max, mean))

# trait         min   max  mean
# 1 FHBSeverity  5.08  38.6  16.1
# 2 HeadingDate 43.9   56.5  50.2
# 3 PlantHeight 60.1   87.4  73.9


## Other plots


# Plot the distributions of each trait and for each family
vp_family_BLUE <- vp_BLUE %>%
  filter(line_name %in% exper) %>%
  mutate(family = str_extract(line_name, "4[0-9]{3}")) %>%
  group_by(family, trait) %>%
  mutate(family_mean_est = mean(value)) %>%
  ungroup()


## Calculate the mean, min, and max for all traits for the VP
vp_family_BLUE %>% 
  group_by(trait) %>% 
  summarize(mean = mean(value), min = min(value), max = max(value))

# trait        mean   min   max
# 1 FHBSeverity  21.1  4.19  56.3
# 2 HeadingDate  53.5 44.2   66.6
# 3 PlantHeight 101.  82.2  117.




# Combine
vp_family_BLUE_toplot <- bind_rows(
  mutate(vp_family_BLUE, type = "family"), 
  mutate(vp_family_BLUE, family = "All", type = "all"))

# Plot per trait
g_vp_family_density <- vp_family_BLUE_toplot %>%
  split(.$trait) %>%
  map(~{
    temp <- .
    # Order the family based on the family mean
    family_order <- distinct(temp, family, family_mean_est) %>% 
      filter(family != "All") %>% 
      arrange(family_mean_est) %>% 
      pull(family)
    
    # Convert family to factor
    temp$family <- factor(temp$family, levels = c(family_order, "All"))
    
    # Create a color scheme
    family_color <- setNames(c(all_colors(n = nlevels(temp$family) - 1), "grey"), levels(temp$family))
    
    # Plot
    temp %>%
      ggplot(aes(x = value, y = type, fill = family)) +
      geom_density_ridges(alpha = 0.2) +
      facet_wrap(~ trait, ncol = 1, scale = "free") +
      scale_fill_manual(values = family_color, name = "Family", guide = FALSE) +
      theme_acs() +
      theme(axis.title = element_blank())
    
  })

# Cowplot
g_density_plot <- plot_grid(plotlist = g_vp_family_density, ncol = 1, align = "hv")
ggsave(filename = "vp_pheno_mean_density.jpg", plot = g_density_plot, path = fig_dir,
       height = 7, width = 4, dpi = 1000)









### Calculate genetic variance within families
### 

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

## Two methods to be used:
## 1. Calculate variance by fitting a model with G and GE

## Use the stage-one BLUEs to calculate variance
vp_family_tomodel1 <- vp_family_tomodel1 %>% 
  mutate(family = str_extract(line_name, "4[0-9]{3}"))


# Fit a model per family
vp_family_varG1 <- vp_family_tomodel1 %>% 
  group_by(trait, family) %>%
  do(calc_varG(data = ., method = "lmer")) %>%
  ungroup()
    

## Plot heritability for each trait/family
vp_family_varG1 %>% 
  mutate(h2 = map_dbl(h2, "heritability")) %>% 
  qplot(x = family, y = h2, data = ., fill = trait, geom = "col") + 
  facet_grid(trait ~ .)


# Plot mean versus varG
vp_family_varG_method1 <- vp_family_varG1 %>% 
  ungroup() %>%
  mutate(varG = map(h2, "var_comp")) %>% 
  unnest(varG) %>% 
  filter(source == "line_name") %>% 
  select(-source)

vp_family_varG_method1 %>% 
  qplot(x = family_mean, y = variance, data = .) + 
  facet_wrap(~ trait, ncol = 2, scales = "free")


## Calculate the superior progeny mean
vp_family_musp <- vp_family_tomodel1 %>% 
  group_by(trait, family) %>%
  do(calc_mean(data = .)) %>%
  ungroup()




## Calculate mean and genetic variance for FHB Severity separately for
## STP and CRM
vp_family_varG1_FHB <- vp_family_tomodel1 %>% 
  filter(trait == "FHBSeverity") %>%
  group_by(location, family) %>%
  do(calc_varG(data = ., method = "lmer")) %>%
  ungroup()

## Calculate the superior progeny mean
vp_family_FHB_musp <- vp_family_tomodel1 %>% 
  filter(trait == "FHBSeverity") %>%
  group_by(location, family) %>%
  do(calc_mean(data = .)) %>%
  ungroup()


## Plot heritability for each location/family
vp_family_varG1_FHB %>% 
  mutate(h2 = map_dbl(h2, "heritability")) %>% 
  qplot(x = family, y = h2, data = ., geom = "col") + 
  facet_grid(location ~ .)

# How many families
n_distinct(vp_family_varG1_FHB$family)


# Plot mean versus varG
vp_family_varG_FHB_method1 <- vp_family_varG1_FHB %>% 
  ungroup() %>%
  mutate(varG = map(h2, "var_comp")) %>% 
  unnest(varG) %>% 
  filter(source == "line_name") %>% 
  select(-source)

vp_family_varG_FHB_method1 %>% 
  qplot(x = family_mean, y = variance, data = .) + 
  facet_wrap(~ location, ncol = 2, scales = "free")



## Save the variance calculations
save("vp_family_varG_method1", "vp_family_varG_FHB_method1", "vp_family_musp", "vp_family_FHB_musp",
     file = file.path(result_dir, "vp_family_analysis.RData"))


load(file.path(result_dir, "vp_family_analysis.RData"))





### Family analysis 

## What is the coefficient of variation for the mean and variance for each trait
vp_family_varG_method1 %>% 
  gather(parameter, estimate, family_mean, variance) %>% 
  group_by(trait, parameter) %>% 
  summarize(cv = sd(estimate) / mean(estimate))



## Using the pedigree information and the BLUEs, correlate the MPV with the family mean
vp_pedigree <- entry_list %>% 
  filter(Group == "Experimental") %>% 
  separate(Pedigree, c("Par1", "Par2"), "/") %>%
  distinct(Par1, Par2, Family)


# Add the BLUEs for each parent and calculate MPV
family_MPV <- vp_pedigree %>%
  rename_all(str_to_lower) %>%
  left_join(., vp_BLUE, by = c("par1" = "line_name")) %>% 
  left_join(., vp_BLUE, by = c("par2" = "line_name", "trait")) %>% 
  select(family, par1, par2, trait, par1_value = value.x, par2_value = value.y) %>% 
  arrange(trait, family) %>%
  mutate(mpv = (par1_value + par2_value) / 2,
         family = str_c("4", family))

family_MPV_mean <- family_MPV %>% 
  left_join(., vp_family_varG_method1)

# Plot
g_mean_mpv <- family_MPV_mean %>% 
  # filter(!(trait == "FHBSeverity" & family_mean > 22.5)) %>% # Remove outlier 
  qplot(x = mpv, y = family_mean, data = .) + 
  facet_wrap(~ trait, ncol = 2, scales = "free")

ggsave(filename = "family_mean_mpv.jpg", plot = g_mean_mpv, path = fig_dir, height = 5, width = 5, dpi = 1000)


# Correlate
family_MPV_mean %>%
  group_by(trait) %>%
  # filter(!(trait == "FHBSeverity" & family_mean > 22.5)) %>% # Remove outlier
  summarize(corr = cor(mpv, family_mean))

# trait         corr
# 1 FHBSeverity 0.243
# 2 HeadingDate 0.723
# 3 PlantHeight 0.613



## What is the relationship between family mean and superior progeny mean?
vp_family_estimates <- vp_family_musp %>% 
  select(-means) %>% 
  left_join(., select(vp_family_varG_method1, -family_mean))

# Mean and range of means and variances
vp_family_estimates %>% 
  gather(parameter, estimate, family_mean:variance) %>% 
  group_by(trait, parameter) %>% 
  summarize_at(vars(estimate), funs(min, max, mean))

# trait       parameter        min   max  mean
# 1 FHBSeverity family_mean 2.14e+ 1 30.8  25.1 
# 2 FHBSeverity mu_sp       1.06e+ 1 19.4  14.8 
# 3 FHBSeverity variance    4.87e-14 15.4   4.95
# 4 HeadingDate family_mean 4.79e+ 1 55.8  52.5 
# 5 HeadingDate mu_sp       4.38e+ 1 53.4  49.5 
# 6 HeadingDate variance    4.33e- 1  8.57  2.27
# 7 PlantHeight family_mean 7.30e+ 1 90.8  81.8 
# 8 PlantHeight mu_sp       6.79e+ 1 84.0  75.3 
# 9 PlantHeight variance    0.       18.4   6.37



vp_family_estimates %>% 
  group_by(trait) %>% 
  summarize(mean_cor = cor(family_mean, mu_sp))

# trait       mean_cor
# 1 FHBSeverity    0.755
# 2 HeadingDate    0.963
# 3 PlantHeight    0.956

# Plot
g_vp_means <- vp_family_estimates %>% 
  ggplot(aes(x = family_mean, y = mu_sp)) +
  geom_point() +
  facet_wrap(~trait, nrow = 1, scales = "free")

ggsave(filename = "vp_mean_corr.jpg", plot = g_vp_means, path = fig_dir, height = 4, width = 8, dpi = 1000)

# Relationship between variance and mu_sp
vp_family_estimates %>% 
  group_by(trait) %>% 
  summarize(cor = cor(variance, mu_sp))

# trait          cor
# 1 FHBSeverity -0.629
# 2 HeadingDate -0.599
# 3 PlantHeight -0.190

## A negative relationship is expected, since a larger variance should lead to lower mu_sp

# Plot
g_vp_musp_var <- vp_family_estimates %>% 
  ggplot(aes(x = variance, y = mu_sp)) +
  geom_point() +
  facet_wrap(~trait, nrow = 1, scales = "free")

ggsave(filename = "vp_musp_var_corr.jpg", plot = g_vp_means, path = fig_dir, height = 4, width = 8, dpi = 1000)


## Relationship of mean and variance
g_vp_mean_var <- vp_family_estimates %>% 
  ggplot(aes(x = family_mean, y = variance)) +
  geom_point() +
  facet_wrap(~trait, nrow = 1, scales = "free")

ggsave(filename = "vp_mean_var_corr.jpg", plot = g_vp_mean_var, path = fig_dir, height = 4, width = 8, dpi = 1000)



## Model mu_sp as a contribution of family_mean and variance
vp_musp_models <- vp_family_estimates %>%
  group_by(trait) %>%
  do({
    data <- .
    
    # fit a linear model
    fit <- lm(mu_sp ~ family_mean + variance, data = data)
    # Tidy
    fit_tidy <- tidy(fit)
    
    # Get the r squared
    r_squared <- summary(fit)$r.squared
    
    # Fit just variance or just the family_mean
    fit_no_var <- lm(mu_sp ~ family_mean, data = data)
    family_mean_rsquared <- summary(fit_no_var)$r.squared
    
    fit_no_mean <- lm(mu_sp ~ variance, data = data)
    var_rsquared <- summary(fit_no_mean)$r.squared
    
    # Return a df
    data_frame(summary = list(fit_tidy), R2 = r_squared, R2_mean = family_mean_rsquared, R2_var = var_rsquared)
      
  })

vp_musp_models %>%
  unnest() %>% 
  filter(term != "(Intercept)")


# trait          R2 R2_mean R2_var term        estimate std.error statistic  p.value
# 1 FHBSeverity 0.733   0.570 0.396  family_mean    0.540    0.139       3.89 2.15e- 3
# 2 FHBSeverity 0.733   0.570 0.396  variance      -0.196    0.0724     -2.70 1.92e- 2
# 3 HeadingDate 0.977   0.927 0.359  family_mean    0.953    0.0374     25.5  6.81e-19
# 4 HeadingDate 0.977   0.927 0.359  variance      -0.369    0.0506     -7.29 1.56e- 7
# 5 PlantHeight 0.944   0.915 0.0363 family_mean    0.911    0.0460     19.8  2.22e-16
# 6 PlantHeight 0.944   0.915 0.0363 variance      -0.169    0.0471     -3.58 1.52e- 3




## Find the number of transgressive segregants per family
vp_and_parent_BLUE <- family_MPV %>% 
  left_join(., vp_family_BLUE) %>% 
  select(trait, family:par2_value, line_name, value)

vp_trans_seg <- vp_and_parent_BLUE %>% 
  group_by(trait, family) %>% 
  mutate(min_parent = min(par1_value, par2_value)) %>% 
  summarize(p_ts = mean(value < min_parent),
            par1_value = unique(par1_value),
            par2_value = unique(par2_value),
            min_parent = unique(min_parent)) %>%
  ungroup()

# Add the family variance
vp_trans_seg1 <- vp_trans_seg %>% 
  left_join(vp_family_varG_method1) %>% 
  select(p_ts, variance)










### Appendix

# ## Use the genotype BLUEs to calculate variance (i.e. already accounting for E and GE)
# 
# # Subset the BLUEs and assign family designators
# vp_family_tomodel2 <- vp_BLUE %>%
#   filter(line_name %in% exper) %>%
#   mutate(family = str_extract(line_name, "4[0-9]{3}"))
# 
# # Fit a model per family
# vp_family_varG2 <- vp_family_tomodel2 %>% 
#   group_by(trait, family) %>%
#   do({
#     df <- droplevels(.)
#     df1 <- df
#     
#     print(unique(df$family))
#     print(unique(df$trait))
#     
#     # Number of lines in the family
#     n_lines <- n_distinct(df1$line_name)
#     
#     control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
#     formula <- value ~ (1|line_name)
#     
#     fit <- lmer(formula = formula, data = df1, control = control)
#     
#     plot_table <- xtabs(~line_name, model.frame(fit))
#     
#     # Get the harmonic mean of the number of environments / reps
#     n_r <- plot_table %>% 
#       harm_mean()
#     
#     # Estimate heritability
#     h2 <- herit(object = fit, exp = "line_name / (line_name + (Residual / (n_r)))", n_r = n_r)
#     
#     fit_null <- lm(value ~ 1, df1)
#     lr <- -2 * as.numeric(logLik(fit_null) - logLik(fit))
#     p_value <- pchisq(q = lr, df = 1, lower.tail = FALSE)
#     
#     sig_test <- data_frame(term = "family", statistic = lr, df = 1, pvalue = p_value)
#     
#     # The intercept of the model is the grand/family mean
#     family_mean <- fixef(fit)[[1]]
#     
#     # Return
#     data_frame(family_mean = family_mean, h2 = list(h2), sig_test = list(sig_test))
#     
#   })
# 
# 
# ## Plot heritability for each trait/family
# vp_family_varG2 %>% 
#   mutate(h2 = map_dbl(h2, "heritability")) %>% 
#   qplot(x = family, y = h2, data = ., fill = trait, geom = "col") + 
#   facet_grid(trait ~ .)
# 
# 
# # Plot mean versus varG
# vp_family_varG_method2 <- vp_family_varG2 %>% 
#   ungroup() %>%
#   mutate(varG = map(h2, "var_comp")) %>% 
#   unnest(sig_test) %>%
#   unnest(varG) %>% 
#   filter(source == "line_name") %>% 
#   select(-source, -term)
# 
# vp_family_varG_method2 %>% 
#   qplot(x = family_mean, y = variance, data = .) + 
#   facet_wrap(~ trait, ncol = 2, scales = "free")
# 
