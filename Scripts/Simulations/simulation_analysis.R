## PopVarValidation - simulation analysis
## 
## Author: Jeff Neyhart
## Last modified: August 31, 2018
## 

# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other libraries
library(cowplot)

# Significance level
alpha <- 0.05

# Load the simulation results
load(file.path(result_dir, "popvar_suitability_simulation_results.RData"))

## Tidy the results
sim_out_tidy <- popvar_simulation_out %>% 
  unnest(results) %>% 
  mutate_at(vars(h2, nQTL, tp_size, map_error), as.factor)

## What combinations are missing?
popvar_simulation_out %>% 
  mutate(complete = T) %>% 
  complete(h2, nQTL, map_error, tp_size, iter, fill = list(complete = F)) %>%
  filter(!complete)
  
## Summarize
sim_out_summ <- sim_out_tidy %>%
  group_by_at(vars(h2:tp_size, parameter)) %>% 
  summarize_at(vars(accuracy, bias), funs(mean, sd)) %>% 
  ungroup()
  


## Plot, with emphasis on map error
sim_out_summ %>% 
  filter(parameter == "varG") %>%
  mutate_at(vars(accuracy_mean:bias_sd), ~round(., 3)) %>%
  ggplot(aes(x = map_error, y = h2, fill = accuracy_mean, label = accuracy_mean)) + 
  geom_tile() + 
  geom_text() + 
  scale_fill_gradient(low = "white", high = "green") +
  facet_grid(nQTL ~ tp_size) +
  theme_acs()

## Fit models
models <- sim_out_tidy %>%
  group_by(parameter) %>%
  do(fit = lm(accuracy ~ h2 + nQTL + map_error + tp_size, data = .))

models$fit %>% map(anova)


## Clearly map error has little effect (probably because most of the 
## genetic variance is accounted by the genetic variance at each locus
sim_out_summ1 <- sim_out_tidy %>%
  filter(map_error == 0) %>%
  gather(statistic, estimate, accuracy, bias) %>%
  group_by_at(vars(h2, nQTL, tp_size, parameter, statistic)) %>%
  summarize_at(vars(estimate), funs(mean = mean(.), n = n(), se = sd(.) / n)) %>%
  mutate(stat = qt(p = 1 - (alpha / 2), df = n - 1) * se,
         lower = mean - stat, upper = mean + stat,
         round_mean = round(mean, 2))

## Plot, with emphasis on tp size
g_accuracy <- sim_out_summ1 %>% 
  filter(statistic == "accuracy") %>%
  ggplot(aes(x = tp_size, y = h2, fill = mean, label = round_mean)) + 
  geom_tile() + 
  geom_text() + 
  scale_fill_gradient(low = "white", high = "green") +
  facet_grid(parameter ~ nQTL) +
  theme_acs()

## Plot using points
g_accuracy1 <- sim_out_summ1 %>% 
  filter(statistic == "accuracy") %>%
  ggplot(aes(x = tp_size, y = mean, ymin = lower, ymax = upper, color = nQTL, group = nQTL)) + 
  geom_point(size = 0.5) +
  geom_line() + 
  geom_errorbar(width = 0.25) +
  facet_grid(parameter ~ h2) +
  theme_acs() +
  theme(legend.position = c(0.90, 0.77), legend.key.height = unit(0.5, "lines"))

## Save
ggsave(filename = "simulation_prediction_accuracy.jpg", plot = g_accuracy1, path = fig_dir,
       height = 4, width = 4, dpi = 1000)




## Plot, with emphasis on tp size
g_bias <- sim_out_summ1 %>% 
  filter(statistic == "bias", parameter == "varG") %>%
  ggplot(aes(x = tp_size, y = h2, fill = mean, label = round_mean)) + 
  geom_tile() + 
  geom_text(size = 3) + 
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  facet_grid(parameter ~ nQTL, labeller = labeller(nQTL = label_both)) +
  theme_acs()

g_bias1 <- sim_out_summ1 %>% 
  filter(statistic == "bias", parameter == "varG") %>%
  ggplot(aes(x = tp_size, y = mean, ymin = lower, ymax = upper, color = nQTL, group = nQTL)) + 
  geom_point() +
  geom_line() + 
  geom_errorbar(width = 0.25) +
  facet_grid(parameter ~ h2) +
  theme_acs()


## Save
ggsave(filename = "simulation_bias.jpg", plot = g_bias, path = fig_dir,
       height = 3, width = 5, dpi = 1000)

