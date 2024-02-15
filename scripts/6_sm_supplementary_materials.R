rm(list = ls())
# load libraries
library(tidyverse)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(multcomp)
library(ggfortify)

# load data
# ! please load the most recently generated data from step 1 ! 
load("data/exp_processed/processed_data_20240215.RData")

# calculate the encounter rate for each individual (used to normalise defence later)
encounters_ind <- data_formatted_cut %>%
  filter(behavior == "meet") %>%
  # group by individuals
  group_by(colony, trial, group_size, brood, color) %>%
  count() %>%
  rename(totalEncounters = n) %>%
  ungroup() %>%
  right_join(., complete_ant_list, by = join_by(colony, group_size, brood, color, trial)) %>% 
  mutate(totalEncounters = ifelse(is.na(totalEncounters), 0, totalEncounters))

# all trials
complete_ant_list_across_trials <- complete_ant_list %>% 
  dplyr::select(-trial) %>% 
  arrange() %>% 
  unique

encounters_ind_across_trials <- data_formatted_cut %>%
  filter(behavior == "meet") %>%
  # group by individuals
  group_by(colony, group_size, brood, color) %>%
  count() %>%
  rename(totalEncounters = n) %>%
  ungroup() %>%
  right_join(., complete_ant_list_across_trials, by = join_by(colony, group_size, brood, color)) %>% 
  mutate(totalEncounters = ifelse(is.na(totalEncounters), 0, totalEncounters))

# calculate individual defence score
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
    defence_score0 = ifelse(totalEncounters == 0, 0, normalised_stings), # simple normalised stings
    defence_score = ifelse(totalEncounters == 0, 0, normalised_stings*log1p(totalStingsDur_ind)),  # defence score
    )

# correlation between number of sting normalised by number of encounters and defence score
cor.test(defence_score_ind$defence_score0, defence_score_ind$defence_score)

# plotting
defence_score_ind %>% 
  drop_na(defence_score0, defence_score) %>% 
  ggplot(., aes(x = defence_score0, y = defence_score)) +
  geom_point(pch = 21) +
  #geom_jitter(aes(color=brood, alpha = 0.4), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1))+
  #facet_grid(~ trial) +
  #facet_wrap(~ brood, ncol = 3) +
  xlab("Defence score") +
  ylab("Normalised number of stinging attempt") +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_blank(),
        text = element_text(size=15))

# calculate sd of defence score 0 (normalised number of stinging attempts) of each trial
sd_defence_score_ind <- defence_score_ind %>%
  group_by(colony, trial, brood, group_size) %>%
  summarise(sd_defence_score = sd(defence_score0))

# calculate sd over all trials
sd_defence_score_ind_all_trial <- defence_score_ind %>%
  group_by(colony, brood, group_size) %>%
  summarise(sd_defence_score = sd(defence_score0))

defence_score_ind %>%
  group_by(colony, trial, brood, group_size) %>%
  summarise(range = max(defence_score0) - min(defence_score0)) %>% 
  ungroup() %>% 
  summarise(mean = mean(range),
            sd = sd(range),
            n = n(),
            se = sd / sqrt(n))

lmer_dol <- lmer(sd_defence_score ~ brood * group_size * trial + (1|colony), data = sd_defence_score_ind)
drop1(lmer_dol, test = "Chisq") # interaction not significant
lmer_dol <- lmer(sd_defence_score ~ brood * group_size + trial + (1|colony), data = sd_defence_score_ind)
drop1(lmer_dol, test = "Chisq") # interaction not significant
lmer_dol <- lmer(sd_defence_score ~ brood + group_size * trial + (1|colony), data = sd_defence_score_ind)
drop1(lmer_dol, test = "Chisq") # interaction not significant
lmer_dol <- lmer(sd_defence_score ~ group_size + brood * trial + (1|colony), data = sd_defence_score_ind)
drop1(lmer_dol, test = "Chisq") # interaction not significant
lmer_dol <- lmer(sd_defence_score ~ group_size + brood + trial + (1|colony), data = sd_defence_score_ind)
drop1(lmer_dol, test = "Chisq")

# Behavioural consistency using normalised number of stinging attempts to measure defence behaviour
# Rank the data within each colony and trial
ranked_defence <- defence_score_ind %>%
  group_by(trial, colony) %>%
  mutate(rank_defence_score = min_rank(defence_score0)) %>% # I used min rank technique
  ungroup()

# Reshape data for correlation tests
reshaped_defence <- reshape2::dcast(ranked_defence, group_size + brood + ID + colony ~ trial, value.var = "rank_defence_score")

# correlation test
reshaped_defence4 <- reshaped_defence[reshaped_defence$group_size == "4", ]
reshaped_defence8 <- reshaped_defence[reshaped_defence$group_size == "8", ]

# remove colonies where there were no defence in trial 2 and 3
no_agg_col_trial2 <- defence_score_ind %>%
  filter(trial == "Trial 2") %>% 
  group_by(colony, trial) %>% 
  summarise(defence_by_colonies = sum(totalStings_ind)) %>% 
  filter(defence_by_colonies == 0) %>%
  dplyr::select(colony)

no_agg_col_trial3 <- defence_score_ind %>%
  filter(trial == "Trial 3") %>% 
  group_by(colony, trial) %>% 
  summarise(defence_by_colonies = sum(totalStings_ind)) %>% 
  filter(defence_by_colonies == 0) %>% 
  dplyr::select(colony)

# group sizes
reshaped_defence_1 <- anti_join(reshaped_defence, no_agg_col_trial2)
reshaped_defence_2 <- anti_join(reshaped_defence_1, no_agg_col_trial3)

reshaped_defence_1_4 <- reshaped_defence_1[reshaped_defence_1$group_size == "4", ]
reshaped_defence_1_8 <- reshaped_defence_1[reshaped_defence_1$group_size == "8", ]

reshaped_defence_2_4 <- reshaped_defence_2[reshaped_defence_2$group_size == "4", ]
reshaped_defence_2_8 <- reshaped_defence_2[reshaped_defence_2$group_size == "8", ]

# brood conditions
reshaped_defence_1_n <- reshaped_defence_1[reshaped_defence_1$brood == "No brood", ]
reshaped_defence_1_l <- reshaped_defence_1[reshaped_defence_1$brood == "With larvae", ]
reshaped_defence_1_p <- reshaped_defence_1[reshaped_defence_1$brood == "With pupae", ]

reshaped_defence_2_n <- reshaped_defence_2[reshaped_defence_1$brood == "No brood", ]
reshaped_defence_2_l <- reshaped_defence_2[reshaped_defence_1$brood == "With larvae", ]
reshaped_defence_2_p <- reshaped_defence_2[reshaped_defence_1$brood == "With pupae", ]

# Trial 1 v.s. Trial 2
test_result <- cor.test(reshaped_defence_1_n$`Trial 1`, reshaped_defence_1_n$`Trial 2`, method = "spearman",exact=FALSE )
cor.test(reshaped_defence_1_l$`Trial 1`, reshaped_defence_1_l$`Trial 2`)
cor.test(reshaped_defence_1_p$`Trial 1`, reshaped_defence_1_p$`Trial 2`)

# Trial 2 v.s. Trial 3
cor.test(reshaped_defence_2_n$`Trial 2`, reshaped_defence_2_n$`Trial 3`)
cor.test(reshaped_defence_2_l$`Trial 2`, reshaped_defence_2_l$`Trial 3`)
cor.test(reshaped_defence_2_p$`Trial 2`, reshaped_defence_2_p$`Trial 3`)


# Output ----
# Function to perform correlation test and extract results
extract_cor_test_results <- function(data1, data2) {
  test_result <- cor.test(data1, data2, method = "spearman", exact=FALSE)
  return(data.frame(
    cor = test_result$estimate,
    p_value = test_result$p.value
  ))
}

# Perform correlation tests and store results
result_1_4 <- extract_cor_test_results(reshaped_defence_1_4$`Trial 1`, reshaped_defence_1_4$`Trial 2`)
result_1_8 <- extract_cor_test_results(reshaped_defence_1_8$`Trial 1`, reshaped_defence_1_8$`Trial 2`)
result_2_4 <- extract_cor_test_results(reshaped_defence_2_4$`Trial 2`, reshaped_defence_2_4$`Trial 3`)
result_2_8 <- extract_cor_test_results(reshaped_defence_2_8$`Trial 2`, reshaped_defence_2_8$`Trial 3`)

# Create a tibble with all results and include the dataset information
results_tibble <- tibble(
  trial_comparison = c("Trial 1 vs 2", "Trial 2 vs 3", "Trial 1 vs 2", "Trial 2 vs 3"),
  group_size = c(rep("4", 2), rep("8", 2)),  # Add dataset info
  cor = c(result_1_4$cor, result_2_4$cor, result_1_8$cor, result_2_8$cor),
  p_value = c(result_1_4$p_value, result_2_4$p_value, result_1_8$p_value, result_2_8$p_value)
)

# Print the results tibble
print(results_tibble)
clipr::write_clip(results_tibble)

# Perform correlation tests and store results for brood
result_1_n <- extract_cor_test_results(reshaped_defence_1_n$`Trial 1`, reshaped_defence_1_n$`Trial 2`)
result_1_l <- extract_cor_test_results(reshaped_defence_1_l$`Trial 1`, reshaped_defence_1_l$`Trial 2`)
result_1_p <- extract_cor_test_results(reshaped_defence_1_p$`Trial 1`, reshaped_defence_1_p$`Trial 2`)
result_2_n <- extract_cor_test_results(reshaped_defence_2_n$`Trial 2`, reshaped_defence_2_n$`Trial 3`)
result_2_l <- extract_cor_test_results(reshaped_defence_2_l$`Trial 2`, reshaped_defence_2_l$`Trial 3`)
result_2_p <- extract_cor_test_results(reshaped_defence_2_p$`Trial 2`, reshaped_defence_2_p$`Trial 3`)

# Create a tibble with all results
results_tibble <- tibble(
  trial_comparison = c("Trial 1 vs 2", "Trial 2 vs 3", "Trial 1 vs 2", "Trial 2 vs 3",
                       "Trial 1 vs 2", "Trial 2 vs 3"),
  brood = c(rep("n", 2), rep("l", 2), rep("p", 2)),
  cor = c(result_1_n$cor, result_2_n$cor, result_1_l$cor, result_2_l$cor,
          result_1_p$cor, result_2_p$cor),
  p_value = c(result_1_n$p_value, result_2_n$p_value, result_1_l$p_value, result_2_l$p_value,
              result_1_p$p_value, result_2_p$p_value)
)

# Print the results tibble
print(results_tibble)
clipr::write_clip(results_tibble)

# DOL (sd of normalised number of stinging attempts) and efficiency in colony defence
defence_efficiency_col <- defence_score_ind %>%
  mutate(group_size_num = as.integer(as.character(group_size))) %>% 
  group_by(colony, group_size, brood, group_size_num) %>% 
  summarise(total_defence = sum(totalStings_ind)) %>% 
  mutate(defence_efficiency = total_defence) %>% 
  ungroup()

eff_dol_sd <- left_join(defence_efficiency_col, sd_defence_score_ind_all_trial)

lm_eff_dol_int <- lm(defence_efficiency ~ sd_defence_score * brood * group_size, data = eff_dol_sd)
drop1(lm_eff_dol_int, test = "F")
lm_eff_dol_int <- lm(defence_efficiency ~ sd_defence_score * brood + group_size, data = eff_dol_sd)
drop1(lm_eff_dol_int, test = "F")
lm_eff_dol_int <- lm(defence_efficiency ~ sd_defence_score + brood * group_size, data = eff_dol_sd)
drop1(lm_eff_dol_int, test = "F")
lm_eff_dol_int <- lm(defence_efficiency ~ sd_defence_score * group_size + brood , data = eff_dol_sd)
drop1(lm_eff_dol_int, test = "F")
lm_eff_dol_final <- lm(defence_efficiency ~ sd_defence_score + group_size + brood, data = eff_dol_sd)
drop1(lm_eff_dol_final, test = "F")
autoplot(lm_eff_dol_final)
summary(lm_eff_dol_final)

# posthoc analysis
post_hoc_result <- emmeans(lm_eff_dol_int, pairwise ~ brood, adjust = "tukey")
# Viewing the results
print(post_hoc_result$contrasts)

