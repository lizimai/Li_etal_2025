rm(list = ls())
# load libraries
library(tidyverse)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(ggfortify)

# load data ----
# ! please load the most recently generated data from step 1 ! 
load("data/exp_processed/processed_data_20240212.RData")

# process data ----
# calculate the total number of encounters for each ant during each trial
encounters_ind <- data_formatted_cut %>%
  filter(behavior == "meet") %>%
  # group by individuals
  group_by(colony, trial, group_size, brood, color) %>%
  count() %>%
  rename(totalEncounters = n) %>%
  ungroup() %>%
  right_join(., complete_ant_list, by = join_by(colony, group_size, brood, color, trial)) %>% 
  mutate(totalEncounters = ifelse(is.na(totalEncounters), 0, totalEncounters))

# calculate complete ant list across all trials
complete_ant_list_across_trials <- complete_ant_list %>% 
  dplyr::select(-trial) %>% 
  arrange() %>% 
  unique()

# calculate the total number of encounters for each ant across all trials
encounters_ind_across_trials <- data_formatted_cut %>%
  filter(behavior == "meet") %>%
  # group by individuals
  group_by(colony, group_size, brood, color) %>%
  count() %>%
  rename(totalEncounters = n) %>%
  ungroup() %>%
  right_join(., complete_ant_list_across_trials, by = join_by(colony, group_size, brood, color)) %>% 
  mutate(totalEncounters = ifelse(is.na(totalEncounters), 0, totalEncounters))

# calculate individual aggression score in each trial
defence_score_ind <- left_join(data_formatted_cut, encounters_ind, by = c("group_size", "brood", "color", "colony", "trial")) %>%
  filter(behavior == "sting") %>%
  group_by(group_size, brood, color, colony, trial, totalEncounters) %>%
  summarise(totalStings_ind = n(),
            totalStingsDur_ind = sum(duration)) %>%
  ungroup() %>%
  left_join(complete_ant_list, ., by = join_by(colony, brood, group_size, color, trial)) %>% 
  mutate(
    totalEncounters = ifelse(is.na(totalEncounters), 0, totalEncounters),
    totalStings_ind = ifelse(is.na(totalStings_ind), 0, totalStings_ind),
    totalStingsDur_ind = ifelse(is.na(totalStingsDur_ind), 0, totalStingsDur_ind),
    totalEncounterNoSting_ind = totalEncounters - totalStings_ind, 
    normalised_stings = totalStings_ind / totalEncounters,
    ID = paste0(colony,color),
    rank_totalEncounters = percent_rank(totalEncounters),
    rank_totalStings_ind = percent_rank(totalStings_ind),
    rank_totalStingsDur_ind = percent_rank(totalStingsDur_ind),
    defence_score = ifelse(totalEncounters == 0, 0, normalised_stings*log1p(totalStingsDur_ind)),
    )

# calculate individual aggression score across all trials
defence_score_ind_across_trials <- left_join(data_formatted_cut, encounters_ind_across_trials, by = c("group_size", "brood", "color", "colony")) %>%
  filter(behavior == "sting") %>%
  group_by(group_size, brood, color, colony, totalEncounters) %>%
  summarise(totalStings_ind = n(),
            totalStingsDur_ind = sum(duration)) %>%
  ungroup() %>%
  left_join(complete_ant_list_across_trials, ., by = join_by(colony, brood, group_size, color)) %>% 
  mutate(
    totalEncounters = ifelse(is.na(totalEncounters), 0, totalEncounters),
    totalStings_ind = ifelse(is.na(totalStings_ind), 0, totalStings_ind),
    totalStingsDur_ind = ifelse(is.na(totalStingsDur_ind), 0, totalStingsDur_ind),
    totalEncounterNoSting_ind = totalEncounters - totalStings_ind, 
    normalised_stings = totalStings_ind / totalEncounters,
    ID = paste0(colony,color),
    rank_totalEncounters = percent_rank(totalEncounters),
    rank_totalStings_ind = percent_rank(totalStings_ind),
    rank_totalStingsDur_ind = percent_rank(totalStingsDur_ind),
    defence_score = ifelse(totalEncounters == 0, 0, normalised_stings*log1p(totalStingsDur_ind)),
  )

# Stats on DOL of defence behaviour ----
# Compare DOL of defence behaviour between group size 4 and 8 using resampling method
defence_score_ind_4 <- defence_score_ind %>% 
  filter(group_size == 4)

defence_score_ind_8 <- defence_score_ind %>% 
  filter(group_size == 8)

# Function to process the data
process_data <- function(data) {
  data %>%
    group_by(colony, trial, brood) %>%
    sample_n(replace = FALSE, size = 4) %>%
    summarise(sd_defence_score = sd(defence_score), .groups = 'drop') %>%
    summarise(mean_sd_defence_score = mean(sd_defence_score))
}

# Store the mean standard deviations
mean_sd_values <- numeric(1000)

# Number of iterations
n_iterations <- 1000

# Iterate 1000 times
set.seed(123) # for reproducibility
for (i in 1:n_iterations) {
  processed_data <- process_data(defence_score_ind_8)
  mean_sd_values[i] <- processed_data$mean_sd_defence_score
}

# Calculate 95% CI
ci <- quantile(mean_sd_values, probs = c(0.025, 0.975))

# Create a tibble to display the result
ci_tibble <- tibble(lower_CI = ci[1], upper_CI = ci[2])

# calculate mean +- SE for group size of 4
defence_score_mean_sd_group_4 <- defence_score_ind_4 %>% 
  summarise(mean_aggression = mean(defence_score),
            sd_aggression = sd(defence_score),
            n = n(),
            se_aggression = sd_aggression / sqrt(n)) %>% 
  dplyr::select(mean_aggression, se_aggression)

# compare differences of DOL between group size 4 and 8 (simulated colonies of size 4 resampled from colonies of group size 8)
defence_score_mean_sd_group_4
ci_tibble

# calculate standard deviation of defence scores as metric for behavioural varaition
# for each trial
sd_defence_score_ind <- defence_score_ind %>%
  group_by(colony, trial, brood, group_size) %>%
  summarise(sd_defence_score = sd(defence_score))

# for all trials
sd_defence_score_ind_across_trials <- defence_score_ind_across_trials %>%
  group_by(colony, brood, group_size) %>%
  summarise(sd_defence_score = sd(defence_score))

# GLMMs
# Start with model with all possible interactions and remove when interactions are not significant
lmer_dol <- lmer(sd_defence_score ~ brood * group_size * trial + (1|colony), data = sd_defence_score_ind)
drop1(lmer_dol, test = "Chisq") # interaction not significant
lmer_dol <- lmer(sd_defence_score ~ brood * group_size + trial + (1|colony), data = sd_defence_score_ind)
drop1(lmer_dol, test = "Chisq") # interaction not significant
lmer_dol <- lmer(sd_defence_score ~ brood + group_size * trial + (1|colony), data = sd_defence_score_ind)
drop1(lmer_dol, test = "Chisq") # interaction not significant
lmer_dol <- lmer(sd_defence_score ~ group_size + brood * trial + (1|colony), data = sd_defence_score_ind)
drop1(lmer_dol, test = "Chisq") # interaction not significant
lmer_dol_final <- lmer(sd_defence_score ~ group_size + brood + trial + (1|colony), data = sd_defence_score_ind)
# Check model assumptions
# simulationOutput <- simulateResiduals(fittedModel = lmer_dol_final, n = 250)
# plot(simulationOutput)
drop1(lmer_dol_final, test = "Chisq")
summary(lmer_dol_final)

# LMs for scores across trials
lmer_dol <- lm(sd_defence_score ~ brood * group_size, data = sd_defence_score_ind_across_trials)
drop1(lmer_dol, test = "F") # interaction not significant
lmer_dol <- lm(sd_defence_score ~ brood + group_size, data = sd_defence_score_ind_across_trials)
drop1(lmer_dol, test = "F") # interaction not significant
summary(lmer_dol)
autoplot(lmer_dol)

# Save data
current_date <- format(Sys.Date(), "%Y%m%d")
file_name <- paste0("data/exp_processed/individual_defence_score", current_date, ".RData")
save(defence_score_ind, sd_defence_score_ind, sd_defence_score_ind_across_trials, file = file_name)
