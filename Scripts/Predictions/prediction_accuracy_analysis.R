## PopVarValidation - prediction accuracy and analysis
##
## This script will look at the prediction accuracy of family means, variance, and superior
## progeny mean.
## 
## Author: Jeff Neyhart
## Last modified: August 20, 2018
## 

# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

library(cowplot)
library(gridExtra)

# Boot reps
boot_reps <- 1000
alpha <- 0.05

## Load the predictions and the results
# Predictions
load(file.path(result_dir, "prediction_results.RData"))
# Results
load(file.path(result_dir, "vp_family_analysis.RData"))

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

# Filter for the empirical crosses
# Then select only the relevant parameters and tidy it up
popvar_pred_cross <- popvar_pred %>%
  map(~mutate(., family = str_c(4, family)) %>%
        filter(family %in% vp_family_musp$family) %>%
        select(family, parent1, parent2, trait, note, family_mean = pred_mu, variance = pred_varG, mu_sp = musp_low) %>%
        gather(parameter, prediction, family_mean:mu_sp) ) %>%
  list(., names(.)) %>%
  pmap(~{names(.x)[ncol(.x)] <- .y; .x}) %>% 
  reduce(left_join) %>%
  gather(tp_set, prediction, realistic, relevant)

## Edit the bootstrapping results
vp_family_results <- vp_family_boot %>% 
  mutate(parameter = case_when(parameter == "mu" ~ "family_mean", parameter == "varG" ~ "variance", TRUE ~ "mu_sp"),
         expectation = base + bias) %>%
  select(trait, family, parameter, estimate = base, expectation)

# Combine the predictions with the estimates - remove NAs
popvar_pred_obs <- left_join(popvar_pred_cross, vp_family_results) %>%
  filter(!is.na(estimate))






### Measure prediction accuracy
# Calculate the correlation between predictions and observations
set.seed(242)
pred_acc <- popvar_pred_obs %>% 
  group_by(trait, parameter, tp_set) %>% 
  do(cbind(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = boot_reps, alpha = alpha), n_fam = length(.$prediction))) %>%
  rowwise() %>%
  mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", "")) %>%
  ungroup()

# ## Do predictions of the variance improve when using the unbiased expectation?
# pred_acc_exp <- popvar_pred_obs %>%
#   group_by(trait, parameter, tp_set) %>%
#   do(bootstrap(x = .$prediction, y = .$expectation, fun = "cor", boot.reps = boot_reps, alpha = alpha)) %>%
#   rowwise() %>%
#   mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", "")) %>%
#   ungroup()


# trait       parameter   tp_set    statistic    base     se     bias ci_lower ci_upper n_fam annotation
# 1 FHBSeverity family_mean realistic cor       0.462   0.203   0.00728   0.0714    0.846    14 *         
# 2 FHBSeverity family_mean relevant  cor       0.460   0.208   0.0125    0.0446    0.850    14 *         
# 3 FHBSeverity mu_sp       realistic cor       0.693   0.161  -0.0171    0.277     0.891    14 *         
# 4 FHBSeverity mu_sp       relevant  cor       0.690   0.153  -0.0123    0.279     0.884    14 *         
# 5 FHBSeverity variance    realistic cor       0.00693 0.222   0.0353   -0.357     0.556    14 ""        
# 6 FHBSeverity variance    relevant  cor       0.00454 0.237   0.0284   -0.367     0.580    14 ""        
# 7 HeadingDate family_mean realistic cor       0.618   0.0829  0.00523   0.447     0.764    26 *         
# 8 HeadingDate family_mean relevant  cor       0.672   0.0807  0.00411   0.502     0.819    26 *         
# 9 HeadingDate mu_sp       realistic cor       0.559   0.0882  0.00702   0.381     0.726    26 *         
# 10 HeadingDate mu_sp       relevant  cor       0.635   0.0837  0.00307   0.460     0.782    26 *         
# 11 HeadingDate variance    realistic cor       0.391   0.191   0.0251    0.0334    0.767    26 *         
# 12 HeadingDate variance    relevant  cor       0.449   0.168   0.0316    0.123     0.776    26 *         
# 13 PlantHeight family_mean realistic cor       0.528   0.128  -0.00212   0.260     0.742    26 *         
# 14 PlantHeight family_mean relevant  cor       0.438   0.132  -0.0109    0.135     0.659    26 *         
# 15 PlantHeight mu_sp       realistic cor       0.621   0.106   0.00100   0.386     0.801    26 *         
# 16 PlantHeight mu_sp       relevant  cor       0.524   0.122  -0.00686   0.252     0.732    26 *         
# 17 PlantHeight variance    realistic cor       0.482   0.137  -0.00319   0.181     0.702    26 *         
# 18 PlantHeight variance    relevant  cor       0.487   0.140  -0.00207   0.176     0.732    26 *

## Create and write a table
pred_table <- pred_acc %>% 
  mutate(annotation = str_c(round(base, 2), " (", round(ci_lower, 2), ", ", round(ci_upper, 2), ")")) %>% 
  select(trait, n_fam, tp_set, parameter, annotation) %>% 
  mutate(parameter = str_replace_all(parameter, param_replace), 
         parameter = factor(parameter, levels = param_replace)) %>% 
  spread(parameter, annotation) %>%
  split(.$tp_set)

# for (df in pred_table) {
#   write_csv(x = select(df, -tp_set), path = file.path(fig_dir, str_c("pred_accuracy_", unique(df$tp_set), ".csv")))
# }




# Plot the same results
g_pred_acc <- popvar_pred_obs %>%
  split(.$tp_set) %>%
  map(~ggplot(., aes(x = prediction, y = estimate)) +
        geom_smooth(method = "lm", se = FALSE) + 
        geom_point() + 
        geom_text(data = mutate(subset(pred_acc, tp_set == unique(.$tp_set)), annotation = str_c("r = ", round(base, 2), annotation)),
                  aes(x = Inf, y = -Inf, label = annotation), size = 3, hjust = 1.2, vjust = -1) + 
        ylab("Observation") +
        xlab("Prediction") + 
        facet_wrap(~ trait + parameter, ncol = 3, scales = "free") + 
        theme_acs() )

# ## Plot for the expectation
# g_pred_acc_exp <- popvar_pred_obs %>%
#   split(.$tp_set) %>%
#   map(~ggplot(., aes(x = prediction, y = expectation)) +
#         geom_smooth(method = "lm", se = FALSE) + 
#         geom_point() + 
#         geom_text(data = mutate(subset(pred_acc_exp, tp_set == unique(.$tp_set)), annotation = str_c("r = ", round(base, 2), annotation)),
#                   aes(x = Inf, y = -Inf, label = annotation), size = 3, hjust = 1.2, vjust = -1) + 
#         ylab("Observation") +
#         xlab("Prediction") + 
#         facet_wrap(~ trait + parameter, ncol = 3, scales = "free") + 
#         theme_acs() )


# Save the plots
for (i in seq_along(g_pred_acc)) {
  filename <- str_c(names(g_pred_acc)[i], "_comb_pred_acc.jpg")
  ggsave(filename = filename, plot = g_pred_acc[[i]], path = fig_dir, height = 6, width = 6, dpi = 1000)
}

## Edit the manuscript-ready plot
g_pred_acc_use <- popvar_pred_obs %>%
  filter(tp_set == "realistic") %>%
  split(.$tp_set) %>%
  map(~{
    df <- .
    df1 <- df %>% 
      left_join(., pred_acc, by = c("trait", "parameter", "tp_set")) %>%
      mutate(parameter = str_replace_all(parameter, param_replace),
             parameter = factor(parameter, levels = param_replace),
             annotation = str_c("r[MP]==", round(base, 2), "^'", annotation, "'"))
    
    # Split by trait and parameter
    plot_list <- df1 %>%
      split(list(.$trait, .$parameter)) %>%
      map(~{
        df2 <- .
        ggplot(df2, aes(x = prediction, y = estimate)) +
          geom_smooth(method = "lm", se = FALSE) + 
          geom_point(size = 1) + 
          geom_text(data = distinct(df2, trait, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
                    parse = TRUE, size = 3, hjust = 1.1, vjust = -0.5) + 
          ylab("Observation") +
          xlab("Prediction") + 
          facet_grid(trait ~ parameter, scales = "free", labeller = labeller(parameter = label_parsed), switch = "y") + 
          scale_y_continuous(breaks = scales::pretty_breaks(), labels = function(x) str_pad(x, width = 2, pad = "0")) + 
          theme_acs() +
          theme(strip.placement = "outside", axis.title = element_blank())
        
      })
    
    ## Re-order
    plot_list <- plot_list[c(1,4,7,2,5,8,3,6,9)]
    
    # # Extract the y axis with strips
    # y_axis_plot <- plot_list[c(1,4,7)] %>% 
    #   map(~ . + theme(axis.title = element_blank())) %>%
    #   plot_grid(plotlist = ., ncol = 1)
    
    
    ## Edit parts of the grid
    # Remove strips
    plot_list[c(2,3,5,6,8,9)] <- plot_list[c(2,3,5,6,8,9)] %>%
      map(~ . + theme(strip.text.y = element_blank(), strip.background.y = element_blank()))
    plot_list[4:9] <- plot_list[4:9] %>%
      map(~ . + theme(strip.text.x = element_blank(), strip.background.x = element_blank()))
    
    # Create the grid by rows
    top <- plot_grid(plotlist = plot_list[1:3], ncol = 3, align = "h", rel_widths = c(1, 0.9, 0.9))
    middle <- plot_grid(plotlist = plot_list[4:6], ncol = 3, align = "h", rel_widths = c(1, 0.9, 0.9))
    bottom <- plot_grid(plotlist = plot_list[7:9], ncol = 3, align = "h", rel_widths = c(1, 0.9, 0.9))
    
    
    # First plot
    plotgrid <- plot_grid(top, middle, bottom, ncol = 1, rel_heights = c(1, 0.9, 0.9))
    
    ## Add axis
    y_axis <- grid::textGrob(label = "Observation", gp = grid::gpar(fontsize = 8), rot = 90)
    x_axis <- grid::textGrob(label = "Prediction", gp = grid::gpar(fontsize = 8))
    
    # Plot again
    grid.arrange(arrangeGrob(plotgrid, left = y_axis, bottom = x_axis))
    
    })

# Save
ggsave(filename = "realistic_comb_pred_acc_use.jpg", plot = grid.arrange(g_pred_acc_use$realistic), path = fig_dir,
       height = 4.5, width = 5, dpi = 1000)



## Filter out some outliers for HD variance
pred_acc_filt <- popvar_pred_obs %>% 
  filter(!(trait == "HeadingDate" & parameter == "variance" & (estimate > 4 | prediction > 0.8))) %>%
  group_by(trait, parameter, tp_set) %>% 
  do(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = boot_reps, alpha = alpha)) %>%
  rowwise() %>%
  mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", "")) %>%
  ungroup()

# Plot the same results
g_pred_acc_filt <- popvar_pred_obs %>%
  filter(!(trait == "HeadingDate" & parameter == "variance" & (estimate > 4 | prediction > 0.8))) %>%
  split(.$tp_set) %>%
  map(~ggplot(., aes(x = prediction, y = estimate)) +
        geom_smooth(method = "lm", se = FALSE) + 
        geom_point() + 
        geom_text(data = mutate(subset(pred_acc_filt, tp_set == unique(.$tp_set)), annotation = str_c("r = ", round(base, 2), annotation)),
                  aes(x = Inf, y = -Inf, label = annotation), size = 3, hjust = 1.2, vjust = -1) + 
        ylab("Observation") +
        xlab("Prediction") + 
        facet_wrap(~ trait + parameter, ncol = 3, scales = "free") + 
        theme_acs() )

# Save the plots
for (i in seq_along(g_pred_acc_filt)) {
  filename <- str_c(names(g_pred_acc_filt)[i], "_comb_pred_acc_filt.jpg")
  ggsave(filename = filename, plot = g_pred_acc_filt[[i]], path = fig_dir, height = 6, width = 6, dpi = 1000)
}


## Are the prediction accuracies using the realistic versus relevant data different?
pred_acc %>% 
  filter(alpha == 0.05) %>% 
  group_by(trait, parameter) %>% 
  mutate(different = !between(base[1], ci_lower[2], ci_upper[2]))

# No







## Analyze bias
# First calculate bias on a per-family basis
popvar_bias <- popvar_pred_obs %>% 
  mutate(bias = (prediction - estimate) / estimate) %>%
  # filter(tp_set == "realistic") %>%
  # Filter out extremely large bias (i.e. very small estimate) %>%
  filter(estimate > 1e-5) %>%
  select(family, trait, tp_set, parameter, prediction:bias)

# Next summarize over traits
popvar_trait_bias <- popvar_bias %>% 
  group_by(tp_set, trait, parameter) %>% 
  summarize(bias = mean(prediction - estimate) / mean(estimate))

# Plot
popvar_bias %>%
  filter(tp_set == "realistic") %>%
  ggplot(aes(x = family, y = bias, color = parameter)) + 
  geom_point() + 
  facet_wrap(~ trait + parameter, scales = "free")

## Just plot bias for variance
g_variance_bias <- popvar_bias %>%
  filter(parameter == "variance", tp_set == "realistic") %>%
  ggplot(aes(x = family, y = bias)) + 
  geom_point(size = 0.5) + 
  facet_wrap(~ trait, scales = "free") +
  ylab("Bias") +
  xlab("Family") + 
  theme_acs() +
  theme(axis.text.x = element_blank())

ggsave(filename = "predicted_variance_bias.jpg", plot = g_variance_bias, path = fig_dir, height = 1.5, width = 3.5, dpi = 1000)


# Plot bias versus heritability
popvar_bias_herit <- popvar_bias %>% 
  filter(tp_set == "realistic") %>%
  left_join(., select(vp_family_varG_method1, trait, family, heritability))
  
popvar_bias_herit %>%
  ggplot(aes(x = heritability, y = bias)) + 
  geom_point() + 
  facet_wrap(~ trait + parameter, scales = "free")

## Correlate and plot for variance
popvar_bias_herit_cor <- popvar_bias_herit %>% 
  filter(parameter == "variance") %>%
  group_by(trait) %>% 
  do(neyhart::bootstrap(x = .$heritability, y = .$bias, fun = "cor", boot.reps = 1000, alpha = alpha))
  


## View the summary over traits
popvar_trait_bias %>% 
  spread(parameter, bias)

# trait       family_mean    mu_sp variance
# 1 FHBSeverity     -0.371  -0.00395   -0.945
# 2 HeadingDate     -0.0453 -0.00751   -0.834
# 3 PlantHeight     -0.102  -0.0356    -0.963

pred_bias_table <- pred_table %>%
  map(~left_join(., subset(popvar_trait_bias, parameter == "variance", c(tp_set, trait, bias)) ))

## Combine bias with prediction accuracy and write tables

for (df in pred_bias_table) {
  write_csv(x = select(df, -tp_set), path = file.path(fig_dir, str_c("pred_accuracy_", unique(df$tp_set), ".csv")))
}



## Ranges of bias on a per-family basis
# How many families were excluded based on low estimated variance?
popvar_pred_obs %>% 
  filter(parameter == "variance") %>%
  group_by(tp_set, trait, parameter) %>%
  summarize(n_removed = sum(estimate <= 1e-5))


popvar_bias %>% 
  filter(parameter == "variance") %>%
  group_by(tp_set, trait) %>% 
  do(tidy(summary(.$bias)))

## FHB Severity
popvar_bias %>% 
  filter(parameter == "variance", tp_set == "realistic") %>% 
  filter(trait == "FHBSeverity") %>%
  summarize(n = n(), n_0.9 = sum(bias > -0.90))

## Heading Date
popvar_bias %>% 
  filter(parameter == "variance", tp_set == "realistic") %>% 
  filter(trait == "HeadingDate") %>%
  summarize(n = n(), n_0.75 = sum(bias > -0.75))

## Plant height
popvar_bias %>% 
  filter(parameter == "variance", tp_set == "realistic") %>% 
  filter(trait == "PlantHeight") %>%
  summarize(n = n(), n_9 = sum(bias > -0.9))















## Prediction accuracy for FHB
load(file.path(result_dir, "prediction_results_FHB.RData"))


# First subset the relevant columns
popvar_pred_FHB <- pred_results_FHB %>%
  left_join(., cross_list, by = c("Par1" = "parent1", "Par2" = "parent2")) %>%
  select(parent1 = Par1, parent2 = Par2, family, trait, location, note, pred_mu = pred.mu, pred_varG = pred.varG,
         musp_high = mu.sp_high, musp_low = mu.sp_low) %>%
  mutate(family = str_c("4", family))

# Filter for the empirical crosses
# Then select only the relevant parameters and tidy it up
popvar_pred_FHB_cross <- popvar_pred_FHB %>%
  select(family, parent1, parent2, trait, location, note, family_mean = pred_mu, variance = pred_varG, mu_sp = musp_low) %>%
  gather(parameter, prediction, family_mean:mu_sp)

## Edit the vp results
## Edit the bootstrapping results
vp_family_results_FHB <- vp_family_boot_FHB %>% 
  mutate(parameter = case_when(parameter == "mu" ~ "family_mean", parameter == "varG" ~ "variance", TRUE ~ "mu_sp"),
         expectation = base + bias,
         trait = "FHBSeverity") %>%
  select(trait, location, family, parameter, estimate = base, expectation)


# Combine the predictions with the estimates - remove NAs
popvar_pred_obs_FHB <- left_join(popvar_pred_FHB_cross, vp_family_results_FHB) %>%
  filter(!is.na(estimate))






### Measure prediction accuracy
# Calculate the correlation between predictions and observations
set.seed(242)
pred_acc <- popvar_pred_obs_FHB %>% 
  group_by(trait, location, parameter) %>% 
  do(cbind(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = boot_reps, alpha = alpha), n_fam = length(.$prediction))) %>%
  rowwise() %>%
  mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", "")) %>%
  ungroup()

# ## Do predictions of the variance improve when using the unbiased expectation?
# pred_acc_exp <- popvar_pred_obs_FHB %>%
#   group_by(trait, location, parameter) %>%
#   do(bootstrap(x = .$prediction, y = .$expectation, fun = "cor", boot.reps = boot_reps, alpha = alpha)) %>%
#   rowwise() %>%
#   mutate(annotation = ifelse(!between(0, ci_lower, ci_upper), "*", "")) %>%
#   ungroup()


# trait       location parameter   statistic    base    se      bias ci_lower ci_upper n_fam annotation
# 1 FHBSeverity CRM      family_mean cor        0.546  0.257 -0.0187    -0.0568    0.916    14 ""        
# 2 FHBSeverity CRM      mu_sp       cor        0.354  0.204 -0.0120    -0.0823    0.711    14 ""        
# 3 FHBSeverity CRM      variance    cor        0.669  0.190 -0.0303     0.203     0.904    14 *         
# 4 FHBSeverity STP      family_mean cor        0.0778 0.177  0.00562   -0.271     0.433    14 ""        
# 5 FHBSeverity STP      mu_sp       cor        0.0912 0.213  0.000154  -0.338     0.491    14 ""        
# 6 FHBSeverity STP      variance    cor       -0.0205 0.267  0.0349    -0.423     0.622    14 "" 

## Create and write a table
pred_table <- pred_acc %>% 
  mutate(annotation = str_c(round(base, 2), " (", round(ci_lower, 2), ", ", round(ci_upper, 2), ")")) %>% 
  select(trait, location, n_fam, parameter, annotation) %>% 
  mutate(parameter = str_replace_all(parameter, param_replace), 
         parameter = factor(parameter, levels = param_replace)) %>% 
  spread(parameter, annotation)

# for (df in pred_table) {
#   write_csv(x = select(df, -tp_set), path = file.path(fig_dir, str_c("pred_accuracy_", unique(df$tp_set), ".csv")))
# }




# Plot the same results
g_pred_acc <- popvar_pred_obs_FHB %>%
  ggplot(., aes(x = prediction, y = estimate)) +
    geom_smooth(method = "lm", se = FALSE) + 
    geom_point(size = 1) + 
    geom_text(data = mutate(pred_acc, annotation = str_c("r = ", round(base, 2), annotation)),
              aes(x = Inf, y = -Inf, label = annotation), size = 3, hjust = 1.2, vjust = -1) + 
    ylab("Observation") +
    xlab("Prediction") + 
    facet_wrap(~ location + parameter, ncol = 3, scales = "free") + 
    theme_acs()

# ## Plot for the expectation (unbiased)
# g_pred_acc_exp <- popvar_pred_obs %>%
#   split(.$tp_set) %>%
#   map(~ggplot(., aes(x = prediction, y = expectation)) +
#         geom_smooth(method = "lm", se = FALSE) + 
#         geom_point() + 
#         geom_text(data = mutate(subset(pred_acc_exp, tp_set == unique(.$tp_set)), annotation = str_c("r = ", round(base, 2), annotation)),
#                   aes(x = Inf, y = -Inf, label = annotation), size = 3, hjust = 1.2, vjust = -1) + 
#         ylab("Observation") +
#         xlab("Prediction") + 
#         facet_wrap(~ trait + parameter, ncol = 3, scales = "free") + 
#         theme_acs() )


# Save the plots
ggsave(filename = "comb_pred_acc_FHB.jpg", plot = g_pred_acc, path = fig_dir, height = 4, width = 5, dpi = 1000)




## Edit the manuscript-ready plot
g_pred_acc_use <- list(popvar_pred_obs_FHB) %>%
  map(~{
    df <- .
    df1 <- df %>% 
      left_join(., pred_acc, by = c("trait", "parameter", "location")) %>%
      mutate(parameter = str_replace_all(parameter, param_replace),
             parameter = factor(parameter, levels = param_replace),
             annotation = str_c("r[MP]==", round(base, 2), "^'", annotation, "'"),
             facet_y = str_c(trait, ": ", location))
    
    # Split by trait and parameter
    plot_list <- df1 %>%
      split(list(.$location, .$parameter)) %>%
      map(~{
        df2 <- .
        ggplot(df2, aes(x = prediction, y = estimate)) +
          geom_smooth(method = "lm", se = FALSE) + 
          geom_point(size = 1) + 
          geom_text(data = distinct(df2, trait, parameter, annotation), aes(x = Inf, y = -Inf, label = annotation), 
                    parse = TRUE, size = 3, hjust = 1.1, vjust = -0.5) + 
          ylab("Observation") +
          xlab("Prediction") + 
          facet_grid(facet_y ~ parameter, scales = "free", labeller = labeller(parameter = label_parsed), switch = "y") + 
          scale_y_continuous(breaks = scales::pretty_breaks(), labels = function(x) str_pad(x, width = 3, pad = "0")) + 
          theme_acs() +
          theme(strip.placement = "outside", axis.title = element_blank())
        
      })
    
    ## Re-order
    plot_list <- plot_list[c(1,3,5,2,4,6)]
    
    # # Extract the y axis with strips
    # y_axis_plot <- plot_list[c(1,4,7)] %>% 
    #   map(~ . + theme(axis.title = element_blank())) %>%
    #   plot_grid(plotlist = ., ncol = 1)
    
    
    ## Edit parts of the grid
    # Remove strips
    plot_list[c(2,3,5,6)] <- plot_list[c(2,3,5,6)] %>%
      map(~ . + theme(strip.text.y = element_blank(), strip.background.y = element_blank()))
    plot_list[4:6] <- plot_list[4:6] %>%
      map(~ . + theme(strip.text.x = element_blank(), strip.background.x = element_blank()))
    
    # Create the grid by rows
    top <- plot_grid(plotlist = plot_list[1:3], ncol = 3, align = "h", rel_widths = c(1, 0.9, 0.9))
    middle <- plot_grid(plotlist = plot_list[4:6], ncol = 3, align = "h", rel_widths = c(1, 0.9, 0.9))

    
    # First plot
    plotgrid <- plot_grid(top, middle, ncol = 1, rel_heights = c(1, 0.9))
    
    ## Add axis
    y_axis <- grid::textGrob(label = "Observation", gp = grid::gpar(fontsize = 8), rot = 90)
    x_axis <- grid::textGrob(label = "Prediction", gp = grid::gpar(fontsize = 8))
    
    # Plot again
    grid.arrange(arrangeGrob(plotgrid, left = y_axis, bottom = x_axis))
    
  })

# Save
ggsave(filename = "comb_pred_acc_FHB_use.jpg", plot = grid.arrange(g_pred_acc_use[[1]]), path = fig_dir,
       height = 3.5, width = 5, dpi = 1000)



## Analyze bias
# First calculate bias on a per-family basis
popvar_bias <- popvar_pred_obs_FHB %>% 
  mutate(bias = (prediction - estimate) / estimate) %>%
  # filter(tp_set == "realistic") %>%
  # Filter out extremely large bias (i.e. very small estimate) %>%
  filter(estimate > 1e-5) %>%
  select(family, trait, location, parameter, prediction:bias)

# Next summarize over traits
popvar_trait_bias <- popvar_bias %>% 
  group_by(trait, location, parameter) %>% 
  summarize(bias = mean(prediction - estimate) / mean(estimate))

# Plot
popvar_bias %>%
  ggplot(aes(x = family, y = bias, color = parameter)) + 
  geom_point() + 
  facet_wrap(~ location + parameter, scales = "free")

# Plot bias versus heritability
popvar_bias_herit <- popvar_bias %>% 
  left_join(., select(vp_family_varG_FHB_method1, family, location, heritability))

popvar_bias_herit %>%
  ggplot(aes(x = heritability, y = bias)) + 
  geom_point() + 
  facet_wrap(~ location + parameter, scales = "free")

## Correlate and plot for variance
popvar_bias_herit_cor <- popvar_bias_herit %>% 
  filter(parameter == "variance") %>%
  group_by(location) %>% 
  do(neyhart::bootstrap(x = .$heritability, y = .$bias, fun = "cor", boot.reps = 1000, alpha = alpha))



## View the summary over traits
popvar_trait_bias %>% 
  spread(parameter, bias)

# trait       family_mean    mu_sp variance
# 1 FHBSeverity     -0.371  -0.00395   -0.945
# 2 HeadingDate     -0.0453 -0.00751   -0.834
# 3 PlantHeight     -0.102  -0.0356    -0.963

pred_bias_table <- pred_table %>%
  left_join(., subset(popvar_trait_bias, parameter == "variance", c(location, trait, bias)))

## Combine bias with prediction accuracy and write tables
write_csv(x = pred_bias_table, path = file.path(fig_dir, "pred_accuracy_FHB.csv"))



## Ranges of bias on a per-family basis
# How many families were excluded based on low estimated variance?
popvar_pred_obs_FHB %>% 
  filter(parameter == "variance") %>%
  group_by(location, trait, parameter) %>%
  summarize(n_removed = sum(estimate <= 1e-5))


popvar_bias %>% 
  filter(parameter == "variance") %>%
  group_by(location, trait) %>% 
  do(broom::tidy(summary(.$bias)))







# ### Analyze discriminatory power
# # Subset the low/high variance families
# popvar_pred_obs_discrim <- popvar_pred_obs %>% 
#   filter(str_detect(note, "VarG")) %>%
#   mutate(discrim_trait = ifelse(str_detect(note, "FHB"), "FHBSeverity", "PlantHeight"),
#          discrim_group = str_extract(note, "low|high")) %>%
#   filter(trait == discrim_trait, parameter != "mu_sp")
# 
# ## Significance test
# popvar_pred_obs_discrim_sig <- popvar_pred_obs_discrim %>% 
#   group_by(trait, parameter, tp_set) %>% 
#   do(test = t.test(estimate ~ discrim_group, data = .)) %>%
#   ungroup() %>%
#   mutate(pvalue = map_dbl(test, "p.value"))
# 
# 
# # Plot
# g_pred_discrim <- popvar_pred_obs_discrim %>%
#   ggplot(aes(x = discrim_group, y = estimate, color = family, shape = trait)) +
#   geom_point(size = 3) +
#   facet_wrap(~ trait + parameter + tp_set, ncol = 2, scales = "free") + 
#   ylab("Observation") +
#   xlab("Predicted variance group") +
#   theme_acs()
# 
# ## Save
# ggsave(filename = "predicted_variance_discrim.jpg", plot = g_pred_discrim, path = fig_dir, height = 5, width = 5, dpi = 1000)
# 
# 
# g_pred_discrim_alt <- popvar_pred_obs_discrim %>% 
#   select(-prediction) %>% 
#   spread(parameter, estimate) %>%
#   ggplot(aes(x = family_mean, y = variance, shape = discrim_group)) +
#   geom_point(size = 3) +
#   facet_wrap(~ trait, ncol = 2, scales = "free") + 
#   scale_shape_discrete(name = "Predicted\nvariance\ngroup") +
#   ylab("Family variance") +
#   xlab("Family mean") +
#   theme_acs()
# 
# ## Save
# ggsave(filename = "predicted_variance_discrim_alt.jpg", plot = g_pred_discrim_alt, path = fig_dir, height = 3, width = 5, dpi = 1000)
# 
# 
# 
