## PopVarValidation - analysis of predictions
##
## This script will look at the prediction output from PopVar and examine the distribution of predictions and
## the relationship of predicted mean and variance
## 
## Author: Jeff Neyhart
## Last modified: August 20, 2018
## 
## 


# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

library(cowplot)

# Load the predictions
load(file.path(result_dir, "prediction_results_realistic.RData"))

## First subset the relevant columns
popvar_pred <- PVV_all_pred %>%
  filter(trait %in% traits) %>%
  select(family = cross, parent1, parent2, trait, pred_mu = pred.mu, pred_varG = pred.varG, musp_high = mu.sp_high,
         musp_low = mu.sp_low, cor_HeadingDate = `cor_w/_HeadingDate`, cor_PlantHeight = `cor_w/_PlantHeight`,
         cor_FHBSeverity = `cor_w/_FHBSeverity`)


### Plot distributions
### 

# Distribution of predicted mu, varG, and mu_sp
popvar_pred_toplot <- popvar_pred %>% 
  select(family:trait, family_mean = pred_mu, variance = pred_varG, mu_sp = musp_low) %>%
  gather(parameter, prediction, family_mean:mu_sp)

## Map over trait-parameter combinations
g_plot_list <- popvar_pred_toplot %>%
  split(list(.$trait, .$parameter)) %>%
  map(~{
    ggplot(., aes(x = prediction)) +
      geom_histogram() + 
      facet_grid(trait ~ parameter, scale = "free", switch = "y") +
      theme_acs() +
      theme(axis.title = element_blank(), strip.placement = "outside", 
            axis.text.y = element_blank(), axis.ticks.y = element_blank())
  })

# Reorder
g_plot_list1 <- g_plot_list[c(1,4,7,2,5,8,3,6,9)]

## Edit some plots
# Remove all strips
g_plot_list1[c(5,6,8,9)] <- g_plot_list1[c(5,6,8,9)] %>%
  map(~. + theme(strip.text = element_blank()))
# Remove the x strip
g_plot_list1[c(4,7)] <- g_plot_list1[c(4,7)] %>%
  map(~. + theme(strip.text.x = element_blank()))
# Remove the y strip
g_plot_list1[c(2,3)] <- g_plot_list1[c(2,3)] %>%
  map(~. + theme(strip.text.y = element_blank()))


## Create a plot grid
g_pred_hist <- plot_grid(plotlist = g_plot_list1, ncol = 3, align = "hv")

# Save
ggsave(filename = "predictions_histogram.jpg", plot = g_pred_hist, path = fig_dir,
       height = 6, width = 8, dpi = 1000)









