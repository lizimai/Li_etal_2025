rm(list = ls())
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggExtra)
library(ggh4x)

load(file = "data/exp_processed/individual_defence_score20240215.RData")

# data preparation ----
# Rank the data within each colony and trial
ranked_defence <- defence_score_ind %>%
  group_by(trial, colony) %>%
  mutate(rank_defence_score = min_rank(defence_score)) %>% # I used min rank technique
  ungroup()

# reshape data for correlation tests
reshaped_defence <- reshape2::dcast(ranked_defence, group_size + brood + ID + colony ~ trial, value.var = "rank_defence_score")

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

# update dataset
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

# plotting
reshaped_defence_1 %>% 
  drop_na(`Trial 1`,`Trial 2`) %>% 
  ggplot(., aes(x = `Trial 1`, y = `Trial 2`, fill = group_size)) +
  geom_count(pch = 21) +
  #geom_jitter(aes(color=brood, alpha = 0.4), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1))+
  stat_smooth(method='lm', aes(`Trial 1`, `Trial 2`), color = "gray", size = 0.5, inherit.aes = FALSE, level = 0.95, se = FALSE) +
  #facet_grid(~ trial) +
  #facet_wrap(~ brood, ncol = 3) +
  xlab("Trial 1 defence score rank") +
  ylab("Trial 2 defence score rank") +
  facet_grid2(col = vars(group_size), scales = "free", independent = "y", axes = "all") +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_blank(),
        text = element_text(size=20))
ggsave("plots/grp_t1t2_corr.pdf", width = 8, height = 4)

reshaped_defence_2 %>% 
drop_na(`Trial 2`,`Trial 3`) %>% 
  ggplot(., aes(x = `Trial 2`, y = `Trial 3`, fill = group_size)) +
  geom_count(pch = 21) +
  #geom_jitter(aes(color=brood, alpha = 0.4), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1))+
  stat_smooth(method='lm', aes(`Trial 2`, `Trial 3`), color = "gray", size = 0.5, inherit.aes = FALSE, level = 0.95, se = FALSE) +
  #facet_grid(~ trial) +
  #facet_wrap(~ brood, ncol = 3) +
  xlab("Trial 2 defence score rank") +
  ylab("Trial 3 defence score rank") +
  facet_grid2(col = vars(group_size), scales = "free", independent = "y", axes = "all") +
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_blank(),
        text = element_text(size=20))
ggsave("plots/grp_t2t3_corr.pdf", width = 8, height = 4)

reshaped_defence_1 %>% 
  drop_na(`Trial 1`,`Trial 2`) %>% 
  ggplot(., aes(x = `Trial 1`, y = `Trial 2`, fill = brood)) +
  geom_count(pch = 21) +
  #geom_jitter(aes(color=brood, alpha = 0.4), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1))+
  stat_smooth(method='lm', aes(`Trial 1`, `Trial 2`), color = "gray", size = 0.5, inherit.aes = FALSE, level = 0.95, se = FALSE) +
  #facet_grid(~ trial) +
  #facet_wrap(~ brood, ncol = 3) +
  xlab("Rank of individual defence score during trial 1") +
  ylab("Rank of individual defence score during trial 2") +
  facet_grid2(col = vars(brood), scales = "free", independent = "y", axes = "all") +
  theme_classic() +
  theme(legend.position = "none")
ggsave("plots/bro_t1t2_corr.pdf", width = 10, height = 4)

reshaped_defence_2 %>% 
  drop_na(`Trial 2`,`Trial 3`) %>% 
  ggplot(., aes(x = `Trial 2`, y = `Trial 3`, fill = brood)) +
  geom_count(pch = 21) +
  #geom_jitter(aes(color=brood, alpha = 0.4), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 1))+
  stat_smooth(method='lm', aes(`Trial 2`, `Trial 3`), color = "gray", size = 0.5, inherit.aes = FALSE, level = 0.95, se = FALSE) +
  #facet_grid(~ trial) +
  #facet_wrap(~ brood, ncol = 3) +
  xlab("Rank of individual defence score during trial 2") +
  ylab("Rank of individual defence score during trial 3") +
  facet_grid2(col = vars(brood), scales = "free", independent = "y", axes = "all") +
  theme_classic() +
  theme(legend.position = "none")
ggsave("plots/bro_t2t3_corr.pdf", width = 10, height = 4)

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

# Save data
current_date <- format(Sys.Date(), "%Y%m%d")
file_name <- paste0("data/exp_processed/ranked_defence_score", current_date, ".RData")
save(ranked_defence, file = file_name)

