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
  unnest(probcor)
  
## Extract the prediction results
sim_pred_results <- sim_results_tidy %>%
  mutate(predictions = map(results, "summary")) %>% 
  unnest(predictions) %>%
  mutate_at(vars(trait1_h2:gencor, dLinkage, pLinkage), as.factor)

# Summarize
sim_pred_summ <- sim_pred_results %>%
  group_by(trait1_h2, trait2_h2, nQTL, tp_size, gencor, dLinkage, pLinkage, trait, parameter) %>%
  mutate(n = n()) %>%
  summarize_at(vars(accuracy, bias, n), funs(mean, sd)) %>%
  select(-n_sd) %>%
  ungroup()


## Plot results for genetic correlation
sim_pred_summ %>%
  filter(dLinkage != 0, parameter == "corG") %>%
  ggplot(aes(x = pLinkage, y = dLinkage, fill = accuracy_mean)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red") +
  facet_wrap(~ gencor, ncol = 2) +
  theme_acs()


## Fit a model
fit <- lm(accuracy ~ gencor + dLinkage + pLinkage, data = sim_pred_results, subset = parameter == "corG" & dLinkage != 0)
anova(fit)
plot(effects::allEffects(fit))

## Nothing is significant





### Genetic correlation selection simulation


# Load the simulation results
load(file.path(result_dir, "popvar_gencor_selection_simulation_results.RData"))

# Tidy
popvar_selection_sim_tidy <- popvar_simulation_out %>% 
  mutate(response = map(results, "response"), 
         common_cross = map(results, "common_cross"),
         correlations = map(results, "correlations")) %>% 
  select(-input, -results)

## Unnest the correlations
correlations_tidy <- popvar_selection_sim_tidy %>%
  unnest(correlations)

response_tidy <- popvar_selection_sim_tidy %>%
  unnest(response) %>%
  left_join(., filter(correlations_tidy, type == "tp_select_cor"), by = c("trait1_h2", "trait2_h2", "gencor", "arch", "iter")) %>%
  mutate(delta_corG = corG.x - corG.y,
         relative_delta_corG = corG.x - corG_mean) %>%
  select(-corG.x, -corG.y) %>%
  mutate_at(vars(trait1_h2, trait2_h2, gencor, arch, trait, selection), as.factor)


n_iter <- n_distinct(response_tidy$iter)


## Summarize the response results
response_summary <- response_tidy %>%
  select(-mean_gv, -corG_mean, -mean_gv_mean, -par_mean_gv, -base_sdG, -type) %>%
  gather(variable, value, stand_response, relative_response, delta_corG, relative_delta_corG) %>%
  group_by(trait1_h2, trait2_h2, gencor, arch, trait, selection, variable) %>%
  summarize_at(vars(value), funs(mean, sd)) %>%
  ungroup() %>%
  mutate(stat = qt(p = 1 - (alpha / 2), df = n_iter - 1) * (sd / sqrt(n_iter) ),
         lower = mean - stat, upper = mean + stat)
    

  
# Plot the response to selection - relative to the base population
response_summary %>%
  filter(variable == "stand_response") %>%
  mutate(gencor = parse_number(gencor)) %>%
  ggplot(aes(x = gencor, y = mean, shape = selection, color = trait)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  facet_grid(trait1_h2 + trait2_h2 ~ arch) +
  theme_acs()

# Response relative to selection on the mean
response_summary %>%
  filter(selection != "mean", variable == "relative_response") %>%
  mutate(gencor = parse_number(gencor)) %>%
  ggplot(aes(x = gencor, y = mean, shape = selection, color = trait)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point() +
  geom_line() +
  facet_grid(trait1_h2 + trait2_h2 ~ arch) +
  theme_acs()


## Plot the change in genetic correlation
response_summary %>% 
  mutate(gencor = parse_number(gencor)) %>%
  filter(trait == "trait1", variable == "delta_corG") %>%
  ggplot(aes(x = gencor, y = mean, shape = selection, color = selection)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point() +
  facet_grid(trait1_h2 + trait2_h2 ~ arch) +
  theme_acs()

# Relative change in genetic correlation
response_summary %>% 
  mutate(gencor = parse_number(gencor)) %>%
  filter(trait == "trait1", selection != "mean") %>%
  ggplot(aes(x = gencor, y = mean, shape = selection, color = selection)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point() +
  facet_grid(trait1_h2 + trait2_h2 ~ arch) +
  theme_acs()



## Create an index of response over traits
response_summary_index <- response_summary %>% 
  group_by(trait1_h2, trait2_h2, gencor, arch, selection) %>% 
  summarize_at(vars(stand_response_mean, relative_response_mean), mean) %>%
  # mutate(gencor = parse_number(gencor)) %>%
  ungroup()
 


# Plot the response to selection - relative to the base population
response_summary_index %>%
  ggplot(aes(x = gencor, y = stand_response_mean, shape = selection, color = selection)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_point() +
  geom_line() +
  facet_grid(trait1_h2 + trait2_h2 ~ arch) +
  theme_acs()

# Response relative to selection on the mean
response_summary_index %>%
  filter(!selection %in% c("mean", "rand")) %>%
  ggplot(aes(x = gencor, y = relative_response_mean, shape = selection, color = selection)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point() +
  # geom_line() +
  facet_grid(trait1_h2 + trait2_h2 ~ arch) +
  theme_acs()


# Plot
response_summary_index %>%
  filter(selection != "mean") %>%
  ggplot(aes(x = gencor, y = index_response_mean, shape = selection, color = selection)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point() +
  facet_grid(trait1_h2 + trait2_h2 ~ arch) +
  theme_acs()

response_tomodel <- response_tidy %>%
  mutate(selection = factor(selection, levels = c("mean", "corG", "musp", "rand")))

## Fit a model - trait2 response
fit <- lm(stand_response ~ trait1_h2 + trait2_h2 + gencor + arch + selection, data = response_tomodel,
          subset = trait == "trait2")

anova(fit)

## Effects plot
plot(effects::allEffects(fit))



# Plot the genetic correlation
response_summary %>%
  filter(trait == "trait1") %>%
  ggplot(aes(x = selection, y = delta_corG_mean, shape = selection)) +
  geom_point() +
  facet_grid(trait1_h2 + trait2_h2 ~ gencor + arch) +
  theme_acs() +
  theme(axis.text.x = element_blank())

## Fit a model - change in correlation relative to selecting on the mean
fit_corG <- lm(delta_corG ~ trait1_h2 + trait2_h2 + gencor + arch + selection, data = response_tomodel,
               subset = trait == "trait1")

anova(fit_corG)

## Effects plot
plot(effects::allEffects(fit_corG))

fit_corG2 <- lm(delta_corG ~ trait1_h2 + trait2_h2 + selection:gencor + arch, data = response_tomodel,
                subset = trait == "trait1")

plot(effects::allEffects(fit_corG2))




