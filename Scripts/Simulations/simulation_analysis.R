## PopVarValidation - simulation analysis
## 
## Author: Jeff Neyhart
## Last modified: August 31, 2018
## 

# Load the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Significance level
alpha <- 0.05

# Load the simulation results
load(file.path(result_dir, "popvar_suitability_simulation_results.RData"))

## Tidy the results
sim_out_tidy <- popvar_simulation_out %>% 
  unnest(results) %>% 
  mutate_at(vars(h2, nQTL, tp_size, map_error), as.factor)

## What combinations are missing?

  
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

## Clearly map error has little effect (probably because most of the 
## genetic variance is accounted by the genetic variance at each locus
sim_out_summ1 <- sim_out_tidy %>%
  filter(map_error == 0) %>%
  group_by_at(vars(h2:tp_size, parameter)) %>% 
  summarize_at(vars(accuracy, bias), funs(mean, sd)) %>% 
  ungroup() %>%
  mutate_at(vars(accuracy_mean:bias_sd), ~round(., 3))

## Plot, with emphasis on tp size
sim_out_summ1 %>% 
  mutate_at(vars(accuracy_mean:bias_sd), ~round(., 3)) %>%
  ggplot(aes(x = tp_size, y = h2, fill = accuracy_mean, label = accuracy_mean)) + 
  geom_tile() + 
  geom_text() + 
  scale_fill_gradient(low = "white", high = "green") +
  facet_grid(parameter ~ nQTL) +
  theme_acs()

## Plot, with emphasis on tp size
sim_out_summ1 %>% 
  filter(parameter == "varG") %>%
  ggplot(aes(x = tp_size, y = h2, fill = bias_mean, label = bias_mean)) + 
  geom_tile() + 
  geom_text() + 
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  facet_grid(parameter ~ nQTL) +
  theme_acs()




