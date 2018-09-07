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
  group_by_at(vars(trait1_h2, trait2_h2, nQTL, tp_size, gencor, trait, parameter)) %>% 
  summarize_at(vars(accuracy, bias), funs(mean, sd), na.rm = TRUE) %>% 
  ungroup()




### Accuracy for genetic correlation


# Fit a model
fit <- lm(accuracy_mean ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor, 
          data = sim_out_summ, subset = parameter == "corG")

## Anova
anova(fit)
# Effect plot
plot(effects::allEffects(fit))


## Number of QTL does not seem to have a strong effect



### Accuracy for genetic variance
fit_trait1 <- lm(accuracy_mean ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor, 
                 data = sim_out_summ, subset = (parameter == "varG" & trait == "trait1"))
anova(fit_trait1)
# Effect plot
plot(effects::allEffects(fit_trait1))

### Accuracy for genetic variance
fit_trait2 <- lm(accuracy_mean ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor, 
                 data = sim_out_summ, subset = (parameter == "varG" & trait == "trait2"))
anova(fit_trait2)
# Effect plot
plot(effects::allEffects(fit_trait2))

## As expected, the heritability of trait1 does not affect the accuracy to predict trait2, and vice versa


### Accuracy of mu_sp
fit_trait1 <- lm(accuracy_mean ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor + trait2_h2:gencor, 
                 data = sim_out_summ, subset = (parameter == "musp" & trait == "trait1"))
anova(fit_trait1)
 # Effect plot
plot(effects::allEffects(fit_trait1))

### Accuracy for genetic variance
fit_trait2 <- lm(accuracy_mean ~ trait1_h2 + trait2_h2 + nQTL + tp_size + gencor + trait1_h2:gencor, 
                 data = sim_out_summ, subset = (parameter == "musp" & trait == "trait2"))
anova(fit_trait2)
# Effect plot
plot(effects::allEffects(fit_trait2))

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









