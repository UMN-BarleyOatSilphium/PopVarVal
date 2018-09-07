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
load(file.path(result_dir, "prediction_results.RData"))

## First subset the relevant columns
popvar_pred <- list(pred_results_realistic, pred_results_relevant) %>% 
  setNames(c("realistic", "relevant")) %>%
  map(~{
    filter(., trait %in% traits) %>%
      left_join(., cross_list, by = c("Par1" = "parent1", "Par2" = "parent2")) %>%
      select(parent1 = Par1, parent2 = Par2, family, trait, note, pred_mu = pred.mu, pred_varG = pred.varG, musp_high = mu.sp_high,
             musp_low = mu.sp_low, cor_HeadingDate = `cor_w/_HeadingDate`, cor_PlantHeight = `cor_w/_PlantHeight`,
             cor_FHBSeverity = `cor_w/_FHBSeverity`)
  })



### Plot distributions
### 

# Distribution of predicted mu, varG, and mu_sp
popvar_pred_toplot <- popvar_pred %>% 
  map(~select(., parent1:trait, family_mean = pred_mu, variance = pred_varG, mu_sp = musp_low) %>%
        gather(parameter, prediction, family_mean:mu_sp) %>%
        mutate(parameter = str_replace_all(parameter, param_replace),
               parameter = factor(parameter, levels = param_replace)) )

# Mean and range of predictions
popvar_pred_summ <- popvar_pred_toplot %>% 
  map(~group_by(., trait, parameter) %>% 
        summarize_at(vars(prediction), funs(min, max, mean)))

# trait       parameter        min   max   mean
# 1 FHBSeverity family_mean 11.3     19.8  15.3  
# 2 FHBSeverity mu_sp       10.4     19.2  14.2  
# 3 FHBSeverity variance     0.00979  1.77  0.417
# 4 HeadingDate family_mean 47.3     55.6  51.2  
# 5 HeadingDate mu_sp       46.4     55.1  50.0  
# 6 HeadingDate variance     0.0111   2.05  0.452
# 7 PlantHeight family_mean 70.3     77.8  74.0  
# 8 PlantHeight mu_sp       69.8     77.4  73.2  
# 9 PlantHeight variance     0.00590  1.09  0.265
# 
# $relevant
# trait       parameter        min   max   mean
# 1 FHBSeverity family_mean 11.3     19.8  15.3  
# 2 FHBSeverity mu_sp       10.3     19.2  14.2  
# 3 FHBSeverity variance     0.00974  1.76  0.418
# 4 HeadingDate family_mean 46.9     55.6  51.4  
# 5 HeadingDate mu_sp       46.2     55.1  50.2  
# 6 HeadingDate variance     0.0116   2.09  0.481
# 7 PlantHeight family_mean 84.9     92.5  88.8  
# 8 PlantHeight mu_sp       84.0     92.0  87.8  
# 9 PlantHeight variance     0.00869  1.43  0.364

# Save as a table
write_csv(x = popvar_pred_summ$realistic, path = file.path(fig_dir, "prediction_summary.csv"))

## Ranges and mean for selected crosses
popvar_pred_summ_select <- popvar_pred_toplot %>%   
  map(~filter(., !is.na(family)) %>% 
        group_by(trait, parameter) %>% 
        summarize_at(vars(prediction), funs(min, max, mean)))
# 
# $`realistic`
# trait       parameter       min    max   mean
# 1 FHBSeverity family_mean 12.6    19.4   16.4  
# 2 FHBSeverity mu_sp       11.6    18.6   15.4  
# 3 FHBSeverity variance     0.0522  1.28   0.388
# 4 HeadingDate family_mean 47.8    54.1   50.2  
# 5 HeadingDate mu_sp       47.0    53.1   49.3  
# 6 HeadingDate variance     0.0899  1.29   0.370
# 7 PlantHeight family_mean 70.9    76.0   73.7  
# 8 PlantHeight mu_sp       70.3    75.4   72.9  
# 9 PlantHeight variance     0.0552  0.697  0.252
# 
# $relevant
# trait       parameter       min    max   mean
# 1 FHBSeverity family_mean 12.6    19.4   16.4  
# 2 FHBSeverity mu_sp       11.6    18.6   15.4  
# 3 FHBSeverity variance     0.0508  1.28   0.393
# 4 HeadingDate family_mean 47.4    53.5   50.2  
# 5 HeadingDate mu_sp       46.4    52.6   49.2  
# 6 HeadingDate variance     0.109   1.52   0.425
# 7 PlantHeight family_mean 84.9    90.9   88.5  
# 8 PlantHeight mu_sp       84.2    90.1   87.5  
# 9 PlantHeight variance     0.0774  0.738  0.353


## Fold difference in parameters for all crosses
popvar_pred_summ$realistic %>% mutate(fold = max / min)


## What are the correlations between the realistic versus relevant predictions?
popvar_pred_compare <- popvar_pred_toplot %>% 
  list(., names(.)) %>% 
  pmap(~{names(.x)[ncol(.x)] <- .y; .x}) %>% 
  reduce(left_join)

## Correlation
popvar_pred_compare %>% 
  group_by(trait, parameter) %>% 
  summarize(pred_cor = cor(realistic, relevant))

# trait       parameter   pred_cor
# 1 FHBSeverity family_mean    1.000
# 2 FHBSeverity mu_sp          1.000
# 3 FHBSeverity variance       0.998
# 4 HeadingDate family_mean    0.973
# 5 HeadingDate mu_sp          0.976
# 6 HeadingDate variance       0.958
# 7 PlantHeight family_mean    0.877
# 8 PlantHeight mu_sp          0.867
# 9 PlantHeight variance       0.773

# Plot
g_pred_compare <- popvar_pred_compare %>% 
  # group_by(trait, parameter) %>%
  # sample_n(1000) %>%
  ggplot(aes(x = realistic, y = relevant)) +
  geom_point() +
  facet_wrap(~ trait + parameter, scales = "free", labeller = labeller(parameter = label_parsed)) +
  theme_acs()

ggsave(filename = "prediction_compare.jpg", plot = g_pred_compare, path = fig_dir, height = 6, width = 8, dpi = 1000)



# ## Compare variance of full-sib crosses versus non full-sib
# popvar_pred_toplot1 <- popvar_pred_toplot %>% 
#   map(~mutate_at(., vars(contains("parent")), funs(family = str_extract(., "[0-9]{4}"))) %>%
#         mutate(full_sib = parent1_family == parent2_family))
# 
# popvar_pred_toplot1 %>% 
#   map(~filter(., parameter == "variance") %>% 
#         group_by(trait, full_sib) %>% 
#         summarize(mean = mean(prediction)) )
# 
# # Plot
# sib_var_list <- popvar_pred_toplot1 %>% 
#   map(~filter(., parameter == "variance") %>% 
#         qplot(x = trait, y = prediction, data = ., geom = "boxplot", fill = full_sib) +
#         theme_acs() )







## Map over trait-parameter combinations
g_plot_list <- popvar_pred_toplot %>%
  map(~{
    df <- .
    plot_list <- df %>%
      split(.$trait) %>%
      map(~{
        ggplot(., aes(x = prediction)) +
          geom_histogram() + 
          geom_point(data = subset(df, !is.na(family) & trait == unique(.$trait)), 
                     aes(y = 10000, color = "selected"), size = 1) +
          # geom_segment()
          facet_grid(trait ~ parameter, scale = "free_x", switch = "y", labeller = labeller(parameter = label_parsed)) +
          scale_color_manual(guide = FALSE, values = umn_palette(n = 3)[3]) +
          # ylim(c(0, 55000)) +
          theme_acs() +
          theme(axis.title = element_blank(), strip.placement = "outside", 
                axis.text.y = element_blank(), axis.ticks.y = element_blank())
        
      })
        
      plot_grid(
        plot_list$FHBSeverity,
        plot_list$HeadingDate + theme(strip.text.x = element_blank()),
        plot_list$PlantHeight + theme(strip.text.x = element_blank()),
        nrow = 3, align = "v", rel_heights = c(1, 0.9, 0.9))

  })

# Save the plots
for (i in seq_along(g_plot_list)) {
  filename <- str_c(names(g_plot_list)[i], "_predictions_histogram.jpg")
  ggsave(filename = filename, plot = g_plot_list[[i]], path = fig_dir, height = 4, width = 5, dpi = 1000)
}





### Plot the mean versus the variance
g_mean_var_list <- popvar_pred_toplot %>%
  map(~{
    spread(., parameter, prediction) %>%
      # group_by(trait) %>% sample_n(1000) %>%
      ggplot(aes(x = mu, y = `V[G]`)) +
      geom_point(size = 1) +
      ylab(expression(V[G])) +
      xlab(expression(mu)) +
      facet_wrap(~ trait, ncol = 3, scales = "free") + 
      theme_acs()
  })
  
# Save the plots
for (i in seq_along(g_mean_var_list)) {
  filename <- str_c(names(g_plot_list)[i], "_prediction_mean_variance.jpg")
  ggsave(filename = filename, plot = g_mean_var_list[[i]], path = fig_dir, height = 3, width = 8, dpi = 1000)
}



## There is this interesting section in the mean v variance plot for FHB severity.
## Is this a particular individual?
# popvar_pred_toplot_FHB <- popvar_pred_toplot %>%
#   filter(parameter != "mu_sp", trait == "FHBSeverity") %>%
#   spread(parameter, prediction)

noted_parent <- "2MS14_3323-011"
label <- str_c(noted_parent, "/*")

g_mean_var_noted_list <- popvar_pred_toplot %>% 
  map(~{
    df <- spread(., parameter, prediction)
    
    df %>%
      # group_by(trait) %>% sample_n(1000) %>%
      ggplot(aes(x = mu, y = `V[G]`)) +
      geom_point(size = 1) +
      geom_point(data = subset(df, parent1 == noted_parent | parent2 == noted_parent),
                 aes(color = "selected"), size = 1) +
      facet_wrap(~ trait, ncol = 3, scales = "free") + 
      scale_color_manual(name = NULL, values = c("selected" = umn_palette(n = 4)[4]), labels = label) +
      ylab(expression(V[G])) +
      xlab(expression(mu)) +
      theme_acs() +
      theme(legend.position = c(0.07, 0.95), legend.margin = margin(), legend.key.width = unit(0, "lines"),
            legend.background = element_rect(fill = alpha("white", 0)))
  })

# Save the plots
for (i in seq_along(g_mean_var_noted_list)) {
  filename <- str_c(names(g_mean_var_noted_list)[i], "_prediction_mean_variance_noted.jpg")
  ggsave(filename = filename, plot = g_mean_var_noted_list[[i]], path = fig_dir, height = 3, width = 8, dpi = 1000)
}




# ## Predict the value of all vp individuals using the TP
# load(file.path(pheno_dir, "PVV_BLUE.RData"))
# 
# K <- A.mat(s2_imputed_mat_use, min.MAF = 0, max.missing = 1)
# 
# tp_prediction_BLUE1 <- tp_prediction_BLUE %>%
#   filter(trait == "FHBSeverity", line_name %in% tp_geno) %>%
#   select(-trait)
# 
# preds <- kin.blup(data = as.data.frame(tp_prediction_BLUE1), 
#                   geno = "line_name", pheno = "value", GAUSS = FALSE, K = K)
# 
# vp_pred <- preds$pred %>% 
#   data_frame(line_name = names(.), pred = .) %>%
#   filter(line_name %in% pot_pars_geno)
# 
# ## Look for a marker that might explain the high predicted value for this individual
# marker_eff <- mixed.solve(y = pull(tp_prediction_BLUE1, value), Z = s2_imputed_mat_use[tp_geno,])$u



## What is the relationship between family mean and mu_sp?
## 

# Plot
g_pred_mean_musp <- popvar_pred %>% 
  ggplot(aes(x = pred_mu, y = musp_low)) + 
  geom_point(size = 1) + 
  geom_abline(slope = 1) + 
  facet_wrap(~ trait, ncol = 3, scales = "free") + 
  theme_acs()

ggsave(filename = "prediction_mean_musp.jpg", plot = g_pred_mean_musp, path = fig_dir, 
       height = 3, width = 7.5, dpi = 1000)


## Fit models
models <- popvar_pred %>%
  map(~group_by(., trait) %>%
        do(data_frame(family_mean = list(lm(musp_low ~ pred_mu, data = .)),
                      variance = list(lm(musp_low ~ pred_varG, data = .)),
                      both = list(lm(musp_low ~ pred_mu + pred_varG, data = .)))) %>%
        ungroup() %>%
        gather(model, fit, -trait) %>%
        mutate(R2 = map_dbl(fit, ~summary(.)$r.squared)) )
  
# $`realistic`
# trait       model       fit          R2
# 1 FHBSeverity family_mean <S3: lm> 0.947 
# 2 HeadingDate family_mean <S3: lm> 0.961 
# 3 PlantHeight family_mean <S3: lm> 0.970 
# 4 FHBSeverity variance    <S3: lm> 0.0691
# 5 HeadingDate variance    <S3: lm> 0.0815
# 6 PlantHeight variance    <S3: lm> 0.0608
# 7 FHBSeverity both        <S3: lm> 0.998 
# 8 HeadingDate both        <S3: lm> 0.999 
# 9 PlantHeight both        <S3: lm> 0.999 
# 
# $relevant
# trait       model       fit          R2
# 1 FHBSeverity family_mean <S3: lm> 0.947 
# 2 HeadingDate family_mean <S3: lm> 0.963 
# 3 PlantHeight family_mean <S3: lm> 0.966 
# 4 FHBSeverity variance    <S3: lm> 0.0669
# 5 HeadingDate variance    <S3: lm> 0.116 
# 6 PlantHeight variance    <S3: lm> 0.0667
# 7 FHBSeverity both        <S3: lm> 0.998 
# 8 HeadingDate both        <S3: lm> 0.999 
# 9 PlantHeight both        <S3: lm> 0.999 


## Select the best families based on mean, then fit the same models
# models_select <- popvar_pred %>%
#   group_by(trait) %>%
#   top_n(n = -10000, wt = pred_mu) %>%
#   do(data_frame(family_mean = list(lm(musp_low ~ pred_mu, data = .)),
#                 variance = list(lm(musp_low ~ pred_varG, data = .)),
#                 both = list(lm(musp_low ~ pred_mu + pred_varG, data = .)))) %>%
#   ungroup() %>%
#   gather(model, fit, -trait) %>%
#   mutate(R2 = map_dbl(fit, ~summary(.)$r.squared))


## What is the ratio of the variance of means to the variance of variances?
param_var <- popvar_pred %>%
  map(~group_by(., trait) %>% 
        summarize_at(vars(pred_mu, pred_varG, musp_low), var) %>%
        mutate(ratio = pred_mu / pred_varG) )

# $`realistic`
# trait       pred_mu pred_varG musp_low ratio
# 1 FHBSeverity    1.15    0.0415     1.22  27.7
# 2 HeadingDate    1.33    0.0384     1.44  34.5
# 3 PlantHeight    1.02    0.0123     1.08  83.0
# 
# $relevant
# trait       pred_mu pred_varG musp_low ratio
# 1 FHBSeverity    1.15    0.0418     1.21  27.5
# 2 HeadingDate    1.49    0.0457     1.66  32.7
# 3 PlantHeight    1.12    0.0212     1.20  53.0



