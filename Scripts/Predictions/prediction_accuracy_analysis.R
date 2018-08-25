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
load(file.path(result_dir, "PVV_prediction_results.RData"))
# Results
load(file.path(result_dir, "vp_family_analysis.RData"))

# Filter the predictions on the relevant traits
popvar_pred1 <- PVV_all_pred %>%
  filter(trait %in% traits)

# Filter for the empirical crosses
# Then select only the relevant parameters and tidy it up
popvar_pred_cross <- popvar_pred1 %>%
  filter(cross %in% vp_family_musp$family) %>%
  select(family = cross, parent1, parent2, trait, family_mean = pred.mu, variance = pred.varG, mu_sp = mu.sp_low) %>%
  gather(parameter, prediction, family_mean:mu_sp)

# Combine and tidy the empirical results
vp_family_results <- left_join(x = select(vp_family_musp, -means), y = select(vp_family_varG_method1, -family_mean)) %>%
  gather(parameter, estimate, family_mean:variance)

# Combine the predictions with the estimates - remove NAs
popvar_pred_obs <- left_join(popvar_pred_cross, vp_family_results) %>%
  filter(!is.na(estimate)) %>%
  ## Add experimental annotation
  left_join(x = ., y = entry_list %>% filter(Group == "Experimental") %>% distinct(Family, Note) %>% 
              mutate(Family = str_c("4", Family)) %>% rename_all(str_to_lower))






### Measure prediction accuracy
# Calculate the correlation between predictions and observations
set.seed(242)
pred_acc <- popvar_pred_obs %>% 
  # filter(estimate > 1e-10) %>%
  group_by(trait, parameter) %>% 
  do(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = 10000)) %>%
  mutate(annotation = ifelse(between(0, ci_lower, ci_upper), "", "*")) %>%
  ungroup()

# trait       parameter   statistic    base     se      bias ci_lower ci_upper annotation
# 1 FHBSeverity family_mean cor        0.502  0.207   0.00997   0.0702     0.862 *         
# 2 FHBSeverity mu_sp       cor        0.690  0.170  -0.0219    0.229      0.893 *         
# 3 FHBSeverity variance    cor       -0.0257 0.234   0.0334   -0.410      0.536 ""        
# 4 HeadingDate family_mean cor        0.662  0.0828  0.00537   0.487      0.811 *         
# 5 HeadingDate mu_sp       cor        0.626  0.0880  0.00633   0.445      0.787 *         
# 6 HeadingDate variance    cor        0.435  0.172   0.0350    0.135      0.794 *         
# 7 PlantHeight family_mean cor        0.303  0.149  -0.00564  -0.0118     0.568 ""        
# 8 PlantHeight mu_sp       cor        0.347  0.163   0.00113   0.00296    0.639 *         
# 9 PlantHeight variance    cor        0.411  0.137  -0.000953  0.112      0.650 *


# Plot the same results
g_pred_acc <- popvar_pred_obs %>%
  filter(!(trait == "HeadingDate" & parameter == "variance" & (estimate > 4 | prediction > 0.8))) %>%
  ggplot(aes(x = prediction, y = estimate)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  geom_point() + 
  geom_text(data = mutate(pred_acc, annotation = str_c("r = ", round(base, 2), annotation)),
            aes(x = Inf, y = -Inf, label = annotation), size = 3, hjust = 1.2, vjust = -1) + 
  ylab("Observation") +
  xlab("Prediction") + 
  facet_wrap(~ trait + parameter, ncol = 3, scales = "free") + 
  theme_acs()

## Save
ggsave(filename = "comb_pred_acc.jpg", plot = g_pred_acc, path = fig_dir, height = 8, width = 8, dpi = 1000)


## Filter out some outliers for HD variance
pred_acc_filt <- popvar_pred_obs %>% 
  filter(!(trait == "HeadingDate" & parameter == "variance" & (estimate > 4 | prediction > 0.8))) %>%
  group_by(trait, parameter) %>% 
  do(bootstrap(x = .$prediction, y = .$estimate, fun = "cor", boot.reps = 10000)) %>%
  mutate(annotation = ifelse(between(0, ci_lower, ci_upper), "", "*")) %>%
  ungroup()

# Plot the same results
g_pred_acc_filt <- popvar_pred_obs %>%
  filter(!(trait == "HeadingDate" & parameter == "variance" & (estimate > 4 | prediction > 0.8))) %>%
  ggplot(aes(x = prediction, y = estimate)) + 
  geom_smooth(method = "lm", se = FALSE) + 
  geom_point() + 
  geom_text(data = mutate(pred_acc_filt, annotation = str_c("r = ", round(base, 2), annotation)),
            aes(x = Inf, y = -Inf, label = annotation), size = 3, hjust = 1.2, vjust = -1) + 
  ylab("Observation") +
  xlab("Prediction") + 
  facet_wrap(~ trait + parameter, ncol = 3, scales = "free") + 
  theme_acs()

## Save
ggsave(filename = "comb_pred_acc_filt.jpg", plot = g_pred_acc_filt, path = fig_dir, height = 8, width = 8, dpi = 1000)





### Analyze discriminatory power
# Subset the low/high variance families
popvar_pred_obs_discrim <- popvar_pred_obs %>% 
  filter(str_detect(note, "VarG")) %>%
  mutate(discrim_trait = ifelse(str_detect(note, "FHB"), "FHBSeverity", "PlantHeight"),
         discrim_group = str_extract(note, "low|high")) %>%
  filter(trait == discrim_trait, parameter != "mu_sp")

## Significance test
popvar_pred_obs_discrim_sig <- popvar_pred_obs_discrim %>% 
  group_by(trait, parameter) %>% 
  do(test = t.test(estimate ~ discrim_group, data = .)) %>%
  ungroup() %>%
  mutate(pvalue = map_dbl(test, "p.value"))


# Plot
g_pred_discrim <- popvar_pred_obs_discrim %>%
  ggplot(aes(x = discrim_group, y = estimate, color = family, shape = trait)) +
  geom_point(size = 3) +
  facet_wrap(~ trait + parameter, ncol = 2, scales = "free") + 
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



