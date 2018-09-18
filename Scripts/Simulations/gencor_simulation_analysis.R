## PopVarValidation - genetic correlation simulation analysis
## 
## Author: Jeff Neyhart
## Last modified: September 3, 2018
## 

# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other libraries
library(cowplot)

# Significance level
alpha <- 0.05






# Load the simulation results
load(file.path(result_dir, "popvar_gencor_simulation_results.RData"))


## Tidy the results
sim_summary_tidy <- popvar_simulation_out %>% 
  mutate(results = map(results, "summary")) %>%
  unnest(results) %>% 
  mutate_at(vars(trait1_h2, trait2_h2, nQTL, tp_size, gencor), as.factor)

sim_meta_tidy <- popvar_simulation_out %>% 
  mutate(results = map(results, "other")) %>%
  unnest(results) %>% 
  mutate_at(vars(trait1_h2, trait2_h2, nQTL, tp_size, gencor), as.factor)


## What combinations are missing?
popvar_simulation_out %>% 
  mutate(complete = T) %>% 
  complete(trait1_h2, trait2_h2, nQTL, gencor, tp_size, iter, fill = list(complete = F)) %>%
  filter(!complete)

## What are the genetic correlations in the tp for each tp size and intented gencor?

# Plot
sim_meta_tidy %>% 
  filter(variable == "tp_gencor") %>% 
  ggplot(aes(x = value, fill = tp_size)) + 
  geom_density() + 
  facet_grid(~ gencor)



## Summarize
sim_out_summ <- sim_summary_tidy %>%
  gather(variable, value, accuracy, bias) %>% 
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, trait, parameter) %>% 
  summarize_at(vars(value), funs(mean, sd), na.rm = TRUE) %>% 
  ungroup()




### Fit models

## Number of QTL does not seem to have a strong effect
models <- sim_summary_tidy %>%
  group_by(trait, parameter) %>%
  do(fit = lm(accuracy ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor, data = .))

# anovas
models$fit %>% map(anova)
     
# Effect plot
plot(effects::allEffects(models$fit[[7]]))



### Interestingly, the heritability of trait2 affects the accuracy to predict mu_sp of trait1 (and vice versa), 
### especially if the correlation between the two traits is stronger


## Summarize
sim_out_summ1 <- sim_summary_tidy %>%
  gather(statistic, estimate, accuracy, bias) %>%
  # Remove NA
  filter(!is.na(estimate)) %>%
  # Remove bad bias estimates
  filter(!(statistic == "bias" & (estimate > 1 | estimate < -1))) %>%
  group_by_at(vars(trait1_h2, trait2_h2, nQTL, tp_size, gencor, trait, parameter, statistic)) %>%
  summarize_at(vars(estimate), funs(mean = mean(., na.rm = T), n = n(), se = sd(., na.rm = T) / n)) %>%
  mutate(stat = qt(p = 1 - (alpha / 2), df = n - 1) * se,
         lower = mean - stat, upper = mean + stat,
         round_mean = round(mean, 2)) %>%
  ungroup()

### Accuaracy to predict genetic correlation
## Plot, with emphasis on tp size
g_accuracy <- sim_out_summ1 %>% 
  filter(parameter == "corG",
         statistic == "accuracy") %>%
  ggplot(aes(x = trait1_h2, y = trait2_h2, fill = mean, label = round_mean)) + 
  geom_tile() + 
  geom_text() + 
  scale_fill_gradient(low = "white", high = "green") +
  facet_grid(gencor + nQTL ~ tp_size) +
  theme_acs()

## Plot using points
g_accuracy1 <- sim_out_summ1 %>% 
  filter(parameter == "corG", statistic == "accuracy", nQTL == 100) %>%
  mutate(trait2_h2 = factor(trait2_h2, levels = rev(levels(trait2_h2)))) %>%
  ggplot(aes(x = tp_size, y = mean, ymin = lower, ymax = upper, color = gencor, group = gencor)) + 
  geom_point(size = 0.5) +
  geom_line() + 
  geom_errorbar(width = 0.25) +
  facet_grid(trait2_h2 ~ trait1_h2, labeller = label_both) +
  theme_acs() +
  theme(legend.position = c(0.90, 0.75), legend.key.height = unit(0.5, "lines"))

# Save
ggsave(filename = "gen_cor_simulation_accuracy.jpg", plot = g_accuracy1, path = fig_dir,
       height = 5, width = 5, dpi = 1000)

## Bias
g_bias <- sim_out_summ1 %>% 
  filter(parameter == "corG", statistic == "bias", nQTL == 100) %>%
  ggplot(aes(x = trait1_h2, y = trait2_h2, fill = mean, label = round_mean)) + 
  geom_tile() + 
  geom_text(size = 2) + 
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  facet_grid(gencor + nQTL ~ tp_size, labeller = label_both) +
  theme_acs()

# Save
ggsave(filename = "gen_cor_simulation_bias.jpg", plot = g_bias, path = fig_dir,
       height = 3, width = 5, dpi = 1000)


# ### Accuaracy to predict mu_sp
# 
# ## Plot using points
# g_accuracy1 <- sim_out_summ1 %>% 
#   filter(parameter == "mu_sp", statistic == "accuracy") %>%
#   mutate(trait2_h2 = factor(trait2_h2, levels = rev(levels(trait2_h2)))) %>%
#   ggplot(aes(x = tp_size, y = mean, ymin = lower, ymax = upper, color = gencor, group = gencor)) + 
#   geom_point(size = 0.5) +
#   geom_line() + 
#   geom_errorbar(width = 0.25) +
#   facet_grid(trait2_h2 ~ trait1_h2, labeller = label_both) +
#   theme_acs() +
#   theme(legend.position = c(0.90, 0.75), legend.key.height = unit(0.5, "lines"))
# 
# # Save
# ggsave(filename = "gen_cor_simulation_accuracy.jpg", plot = g_accuracy1, path = fig_dir,
#        height = 5, width = 5, dpi = 1000)










### Genetic correlation genetic architecture simulation
# Load the simulation results
load(file.path(result_dir, "popvar_gencor_space_simulation_results.RData"))

# Tidy the results and extract the probability of linkage and the degree of linkage
sim_results_tidy <- popvar_corG_space_simulation_out %>%
  rename(probcor = probor) %>%
  mutate(probcor = map(probcor, ~as.data.frame(.) %>% `names<-`(., c("dLinkage", "pLinkage")) %>% tail(., 1))) %>% # The tail is used to remove the probabilities of pleiotropy
  unnest(probcor) %>%
  ## Extract the results
  mutate(predictions = map(results, "summary"),
         correlations = map(results, "other")) %>%
  select(-results) %>%
  mutate_at(vars(trait1_h2:gencor, dLinkage, pLinkage), as.factor)
  


## Extract the training population genetic correlation
base_cor_summ <- sim_results_tidy %>%
  unnest(correlations) %>%
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, dLinkage, pLinkage, variable) %>%
  summarize_at(vars(value), funs(mean, sd)) %>%
  ## Fill-in missing combinations when dLinkage == 0
  bind_rows(., 
            filter(., dLinkage == 0) %>%
            group_by(variable, add = T) %>% 
            complete(trait1_h2, trait2_h2, nQTL, tp_size, gencor, pLinkage, variable) %>%
            group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, variable) %>% 
            mutate_at(vars(mean, sd), funs(mean), na.rm = T) %>%
            ungroup() %>%
            distinct(trait1_h2, trait2_h2, nQTL, tp_size, gencor, pLinkage, dLinkage, variable, mean, sd) %>%
            filter(pLinkage != 1)
  ) %>% ungroup()





## Plot
g_base_cor <- base_cor_summ %>%
  filter(variable == "tp_gencor") %>%
  ggplot(aes(x = pLinkage, y = dLinkage, fill = mean)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red") +
  facet_wrap(~ gencor, nrow = 1) +
  theme_acs()

# save
ggsave(filename = "gencor_arch_space_base_corG.jpg", plot = g_base_cor, path = fig_dir,
       height = 3, width = 6, dpi = 1000)



  
## Extract the prediction results
# Summarize
pred_results_summ <- sim_results_tidy %>%
  unnest(predictions) %>%
  gather(variable, value, accuracy, bias) %>%
  filter(!(variable == "bias" & abs(value) > 2)) %>%
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, dLinkage, pLinkage, trait, parameter, variable) %>%
  summarize_at(vars(value), funs(mean, sd, n())) %>%
  ungroup()


## Plot results for genetic correlation
g_pred_acc_corG <- pred_results_summ %>%
  filter(dLinkage != 0, parameter == "corG", variable == "accuracy") %>%
  ggplot(aes(x = pLinkage, y = dLinkage, fill = mean)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "green") +
  facet_grid(~ gencor) +
  theme_acs()


# bias
g_pred_bias_corG <- pred_results_summ %>%
  filter(dLinkage != 0, parameter == "corG", variable == "bias") %>%
  ggplot(aes(x = pLinkage, y = dLinkage, fill = mean)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", high = "blue") +
  facet_grid(~ gencor) +
  theme_acs()

# Combine
g_pred_corG <- plot_grid(g_pred_acc_corG, g_pred_bias_corG, ncol = 1, align = "hv")
ggsave(filename = "gencor_arch_space_pred.jpg", plot = g_pred_corG, path = fig_dir, width = 6, height = 5, dpi = 1000)







## Fit a model
fit <- lm(accuracy ~ gencor + dLinkage + pLinkage, data = sim_pred_results, subset = parameter == "corG" & dLinkage != 0)
anova(fit)
plot(effects::allEffects(fit))

## Nothing is significant










### Genetic correlation selection simulation


# Load the simulation results
load(file.path(result_dir, "popvar_gencor_selection_simulation_results.RData"))

## Are there any missing trials?
popvar_gencor_selection_simulation_out %>%
  distinct(trait1_h2, trait2_h2, gencor, arch, iter) %>%
  mutate_all(as.factor) %>% 
  mutate(obs = T) %>% 
  complete(trait1_h2, trait2_h2, gencor, arch, iter, fill = list(obs = F)) %>% 
  filter(!obs) %>%
  nrow()



# Tidy
popvar_selection_sim_tidy <- popvar_gencor_selection_simulation_out %>% 
  mutate(response = map(results, "response"), 
         common_cross = map(results, "common_cross"),
         correlations = map(results, "correlations")) %>% 
  select(-input, -results)

## Unnest the correlations
correlations_tidy <- popvar_selection_sim_tidy %>%
  unnest(correlations) %>%
  rename(base_corG = corG)

response_tidy <- popvar_selection_sim_tidy %>%
  unnest(response) %>%
  left_join(., filter(correlations_tidy, type != "tp_select_cor"), by = c("trait1_h2", "trait2_h2", "gencor", "arch", "iter")) %>%
  mutate(delta_corG = corG - base_corG,
         intensity = round(intensity, 2)) %>%
  # Sum the responses to selection for both traits as an index
  group_by(trait1_h2, trait2_h2, gencor, arch, iter, intensity, selection) %>%
  mutate_at(vars(relative_response, stand_response), funs(index = sum)) %>%
  ungroup()
 

## Plot the genetic correlation in the TP
correlations_tidy %>% 
  filter(type == "tp_base_cor", trait1_h2 == 1, trait2_h2 == 1) %>% 
  ggplot(aes(x = base_corG, fill = as.factor(gencor))) + 
  geom_density() + 
  facet_grid( ~ arch) +
  theme_acs()

## We have a problem with the close linkage and loose linkage groups.




## Summarize the response results
response_tomodel <- response_tidy %>%
  select(trait1_h2:trait, contains("response"), contains("index"), delta_corG) %>%
  gather(variable, value, contains("response"), contains("index"), delta_corG) %>%
  # Change the trait name to "trait" for the indices and the correlation
  mutate(trait = ifelse(variable %in% c("relative_response_index", "stand_response_index", "delta_corG"), "trait", trait)) %>%
  distinct(trait1_h2, trait2_h2, gencor, arch, iter, intensity, selection, trait, variable, value) %>%
  mutate_at(vars(trait1_h2, trait2_h2, gencor), as.factor)
  
  
## Summarize over iterations
response_summary <- response_tomodel %>% 
  group_by(trait1_h2, trait2_h2, gencor, arch, trait, intensity, selection, variable) %>%
  summarize_at(vars(value), funs(mean, sd, n())) %>%
  ungroup() %>%
  mutate(stat = qt(p = 1 - (alpha / 2), df = n - 1) * (sd / sqrt(n) ),
         lower = mean - stat, upper = mean + stat)
    

### Plotting 


## Response to selection of the index - relative to base population
g_stand_resp <- response_summary %>%
  filter(variable == "stand_response_index") %>%
  # filter(gencor == 0.75) %>%
  ggplot(aes(x = intensity, y = mean, shape = selection, fill = selection)) +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  ylab("Standardized response (relative to base population)") +
  xlab("Proportion of selected individuals") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs()

# Save
ggsave(filename = "gencor_sim_stand_response.jpg", plot = g_stand_resp, path = fig_dir,
       height = 6, width = 14, dpi = 1000)

## Take a subset
g_stand_resp_sub <- response_summary %>%
  filter(variable == "stand_response_index") %>%
  filter(trait1_h2 == 1, gencor %in% c(-0.75, 0.75)) %>%
  ggplot(aes(x = intensity, y = mean, shape = selection, fill = selection)) +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  ylab("Standardized response (relative to base population)") +
  xlab("Proportion of selected individuals") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs()

# Save
ggsave(filename = "gencor_sim_stand_response_subset.jpg", plot = g_stand_resp_sub, path = fig_dir,
       height = 5, width = 8, dpi = 1000)


## Model
fit_resp <- response_tomodel %>% 
  filter(variable == "stand_response_index", intensity %in% c(0.01)) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(intensity = as.factor(intensity)) %>%
  lm(value ~ trait1_h2 + trait2_h2 + gencor + arch + selection + gencor:arch, data = .)

anova(fit_resp)
plot(effects::allEffects(fit_resp))

## Notes:
## 1. There is a strong interaction between architecture and intended genetic correlation, where pleiotropic architecture
## is highly influenced by the genetic correlation
## 2. There is little interaction of the selection strategy with architecture or intended genetic correlation



# Response to selection relative to that based on the mean
g_relative_resp <- response_summary %>%
  filter(variable == "relative_response_index") %>%
  ggplot(aes(x = intensity, y = mean, shape = selection, fill = selection)) +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  ylab("Standardized response (relative to mean selection)") +
  xlab("Proportion of selected individuals") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs()

# Save
ggsave(filename = "gencor_sim_relative_response.jpg", plot = g_relative_resp, path = fig_dir,
       height = 6, width = 14, dpi = 1000)

# Response relative to selection on the mean - subset
g_relative_resp_sub <- response_summary %>%
  filter(variable == "relative_response_index") %>%
  filter(trait1_h2 == 1, gencor %in% c(-0.75, 0.75)) %>%
  # filter(selection != "musp") %>%
  ggplot(aes(x = intensity, y = mean, shape = selection, fill = selection)) +
  geom_hline(yintercept = 0, size = 0.25) + 
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  ylab("Standardized response (relative to mean selection)") +
  xlab("Proportion of selected individuals") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs() +
  theme(axis.text.x = element_text(size = 5))

# Save
ggsave(filename = "gencor_sim_relative_response_sub.jpg", plot = g_relative_resp_sub, path = fig_dir,
       height = 5, width = 8, dpi = 1000)


## Model
fit_resp <- response_tomodel %>% 
  filter(variable == "relative_response_index", intensity %in% c(0.01), selection %in% c("musp", "muspC")) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(intensity = as.factor(intensity)) %>%
  lm(value ~ trait1_h2 + trait2_h2 + gencor + arch + selection + gencor:arch, data = .)

anova(fit_resp)
plot(effects::allEffects(fit_resp))

## Notes
## 1. The advantage of selecting on musp or muspC is greater when heritability is lower






## Plot the change in genetic correlation
g_delta_corG <- response_summary %>%
  filter(variable == "delta_corG") %>%
  ggplot(aes(x = intensity, y = mean, shape = selection, fill = selection)) +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  ylab("Change in genetic correlation (relative to base population)") +
  xlab("Proportion of selected individuals") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs()  +
  theme(axis.text.x = element_text(size = 5))

# Save
ggsave(filename = "gencor_sim_delta_corG.jpg", plot = g_delta_corG, path = fig_dir,
       height = 6, width = 14, dpi = 1000)


## Plot the change in genetic correlation - subset
g_delta_corG_sub <- response_summary %>%
  filter(variable == "delta_corG") %>%
  filter(trait1_h2 == 1, gencor %in% c(-0.75, 0.75)) %>%
  ggplot(aes(x = intensity, y = mean, shape = selection, fill = selection)) +
  geom_point(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  ylab("Change in genetic correlation (relative to base population)") +
  xlab("Proportion of selected individuals") +
  facet_grid(trait1_h2 + trait2_h2 ~ arch + gencor, 
             labeller = labeller(arch = label_value, .default = label_both)) +
  theme_acs()  +
  theme(axis.text.x = element_text(size = 5))

# Save
ggsave(filename = "gencor_sim_delta_corG_sub.jpg", plot = g_delta_corG_sub, path = fig_dir,
       height = 5, width = 8, dpi = 1000)



## Model
fit_corG <- response_tomodel %>% 
  filter(variable == "delta_corG", intensity %in% c(0.01)) %>%
  mutate_if(is.character, as.factor) %>%
  mutate(intensity = as.factor(intensity)) %>%
  lm(value ~ trait1_h2 + trait2_h2 + gencor + arch + selection, data = .)

anova(fit_corG)
plot(effects::allEffects(fit_corG))





