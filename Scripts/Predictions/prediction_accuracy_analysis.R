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

# Combine and tidy the empirical results
vp_family_results <- left_join(x = select(vp_family_musp, -means), y = select(vp_family_varG_method1, -family_mean)) %>%
  gather(parameter, estimate, family_mean:variance)

# Combine the predictions with the estimates - remove NAs
popvar_pred_obs <- left_join(popvar_pred_cross, vp_family_results) %>%
  filter(!is.na(estimate))






### Measure prediction accuracy
# Calculate the correlation between predictions and observations
set.seed(242)
pred_acc <- popvar_pred_obs %>% 
  group_by(trait, parameter, tp_set) %>% 
  do(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = 10000, alpha = c(0.01, 0.05, 0.10))) %>%
  ungroup() %>%
  mutate(alpha = rep(c(0.01, 0.05, 0.10), length.out = nrow(.))) %>%
  rowwise() %>%
  mutate(annotation = case_when(
    !between(0, ci_lower, ci_upper) & alpha == 0.01 ~ "***",
    !between(0, ci_lower, ci_upper) & alpha == 0.05 ~ "**",
    !between(0, ci_lower, ci_upper) & alpha == 0.1 ~ "*",
    TRUE ~ "")) %>%
  ungroup()

# Select the most significant results
pred_acc_sig <- pred_acc %>% 
  group_by(trait, parameter, tp_set) %>% 
  filter(annotation != "" | alpha == max(alpha)) %>%
  slice(1)

# trait       parameter   tp_set    statistic  base     se     bias ci_lower ci_upper alpha annotation
# 1 FHBSeverity family_mean realistic cor       0.504   0.196   0.0126    0.110     0.862  0.05 **        
# 2 FHBSeverity family_mean relevant  cor       0.502   0.200   0.0106    0.0912    0.858  0.05 **        
# 3 FHBSeverity mu_sp       realistic cor       0.693   0.163  -0.0159    0.0513    0.926  0.01 ***       
# 4 FHBSeverity mu_sp       relevant  cor       0.690   0.163  -0.0158    0.0278    0.916  0.01 ***       
# 5 FHBSeverity variance    realistic cor       0.00693 0.230   0.0376   -0.301     0.485  0.1  ""        
# 6 FHBSeverity variance    relevant  cor       0.00454 0.238   0.0389   -0.313     0.495  0.1  ""        
# 7 HeadingDate family_mean realistic cor       0.619   0.0874  0.00293   0.345     0.810  0.01 ***       
# 8 HeadingDate family_mean relevant  cor       0.673   0.0772  0.00386   0.450     0.846  0.01 ***       
# 9 HeadingDate mu_sp       realistic cor       0.559   0.0892  0.00559   0.291     0.771  0.01 ***       
# 10 HeadingDate mu_sp       relevant  cor       0.635   0.0827  0.00421   0.389     0.820  0.01 ***       
# 11 HeadingDate variance    realistic cor       0.391   0.197   0.0293    0.0452    0.778  0.05 **        
# 12 HeadingDate variance    relevant  cor       0.449   0.169   0.0306    0.0267    0.850  0.01 ***       
# 13 PlantHeight family_mean realistic cor       0.524   0.128  -0.00417   0.121     0.793  0.01 ***       
# 14 PlantHeight family_mean relevant  cor       0.433   0.135  -0.00698   0.0202    0.734  0.01 ***       
# 15 PlantHeight mu_sp       realistic cor       0.621   0.109  -0.00151   0.257     0.835  0.01 ***       
# 16 PlantHeight mu_sp       relevant  cor       0.524   0.126  -0.00293   0.146     0.781  0.01 ***       
# 17 PlantHeight variance    realistic cor       0.482   0.141  -0.00496   0.0112    0.762  0.01 ***       
# 18 PlantHeight variance    relevant  cor       0.487   0.140  -0.00316   0.0611    0.770  0.01 ***


# Plot the same results
g_pred_acc <- popvar_pred_obs %>%
  split(.$tp_set) %>%
  map(~ggplot(., aes(x = prediction, y = estimate)) +
        geom_smooth(method = "lm", se = FALSE) + 
        geom_point() + 
        geom_text(data = mutate(subset(pred_acc_sig, tp_set == unique(.$tp_set)), annotation = str_c("r = ", round(base, 2), annotation)),
                  aes(x = Inf, y = -Inf, label = annotation), size = 3, hjust = 1.2, vjust = -1) + 
        ylab("Observation") +
        xlab("Prediction") + 
        facet_wrap(~ trait + parameter, ncol = 3, scales = "free") + 
        theme_acs() )

# Save the plots
for (i in seq_along(g_pred_acc)) {
  filename <- str_c(names(g_pred_acc)[i], "_comb_pred_acc.jpg")
  ggsave(filename = filename, plot = g_pred_acc[[i]], path = fig_dir, height = 6, width = 6, dpi = 1000)
}


## Filter out some outliers for HD variance
pred_acc_filt <- popvar_pred_obs %>% 
  filter(!(trait == "HeadingDate" & parameter == "variance" & (estimate > 4 | prediction > 0.8))) %>%
  group_by(trait, parameter, tp_set) %>% 
  do(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = 10000, alpha = c(0.01, 0.05, 0.10))) %>%
  ungroup() %>%
  mutate(alpha = rep(c(0.01, 0.05, 0.10), length.out = nrow(.))) %>%
  rowwise() %>%
  mutate(annotation = case_when(
    !between(0, ci_lower, ci_upper) & alpha == 0.01 ~ "***",
    !between(0, ci_lower, ci_upper) & alpha == 0.05 ~ "**",
    !between(0, ci_lower, ci_upper) & alpha == 0.1 ~ "*",
    TRUE ~ "")) %>%
  ungroup()

pred_acc_filt_sig <- pred_acc_filt %>% 
  group_by(trait, parameter, tp_set) %>% 
  filter(annotation != "" | alpha == max(alpha)) %>%
  slice(1)

# Plot the same results
g_pred_acc_filt <- popvar_pred_obs %>%
  filter(!(trait == "HeadingDate" & parameter == "variance" & (estimate > 4 | prediction > 0.8))) %>%
  split(.$tp_set) %>%
  map(~ggplot(., aes(x = prediction, y = estimate)) +
        geom_smooth(method = "lm", se = FALSE) + 
        geom_point() + 
        geom_text(data = mutate(subset(pred_acc_filt_sig, tp_set == unique(.$tp_set)), annotation = str_c("r = ", round(base, 2), annotation)),
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





### Analyze discriminatory power
# Subset the low/high variance families
popvar_pred_obs_discrim <- popvar_pred_obs %>% 
  filter(str_detect(note, "VarG")) %>%
  mutate(discrim_trait = ifelse(str_detect(note, "FHB"), "FHBSeverity", "PlantHeight"),
         discrim_group = str_extract(note, "low|high")) %>%
  filter(trait == discrim_trait, parameter != "mu_sp")

## Significance test
popvar_pred_obs_discrim_sig <- popvar_pred_obs_discrim %>% 
  group_by(trait, parameter, tp_set) %>% 
  do(test = t.test(estimate ~ discrim_group, data = .)) %>%
  ungroup() %>%
  mutate(pvalue = map_dbl(test, "p.value"))


# Plot
g_pred_discrim <- popvar_pred_obs_discrim %>%
  ggplot(aes(x = discrim_group, y = estimate, color = family, shape = trait)) +
  geom_point(size = 3) +
  facet_wrap(~ trait + parameter + tp_set, ncol = 2, scales = "free") + 
  ylab("Observation") +
  xlab("Predicted variance group") +
  theme_acs()

## Save
ggsave(filename = "predicted_variance_discrim.jpg", plot = g_pred_discrim, path = fig_dir, height = 5, width = 5, dpi = 1000)


g_pred_discrim_alt <- popvar_pred_obs_discrim %>% 
  select(-prediction) %>% 
  spread(parameter, estimate) %>%
  ggplot(aes(x = family_mean, y = variance, shape = discrim_group)) +
  geom_point(size = 3) +
  facet_wrap(~ trait, ncol = 2, scales = "free") + 
  scale_shape_discrete(name = "Predicted\nvariance\ngroup") +
  ylab("Family variance") +
  xlab("Family mean") +
  theme_acs()

## Save
ggsave(filename = "predicted_variance_discrim_alt.jpg", plot = g_pred_discrim_alt, path = fig_dir, height = 3, width = 5, dpi = 1000)



