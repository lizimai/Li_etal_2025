# calculate the individual consistency in aggression
rm(list = ls())
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggExtra)
library(ggh4x)

#load(file = "data/exp_processed/individual_aggression_score20240129.RData")
load(file = "data/exp_processed/individual_aggression_score20240203.RData")

# Rank the data within each colony and trial
ranked_aggression <- aggression_score_ind %>%
  group_by(trial, colony) %>%
  mutate(rank_aggression_score = min_rank(aggression_score0)) %>% # I used min rank technique
  ungroup()

# Reshape data for correlation tests
reshaped_aggression <- reshape2::dcast(ranked_aggression, group_size + brood + ID + colony ~ trial, value.var = "rank_aggression_score")

# correlation test
reshaped_aggression4 <- reshaped_aggression[reshaped_aggression$group_size == "4", ]
reshaped_aggression8 <- reshaped_aggression[reshaped_aggression$group_size == "8", ]

# remove colonies where there were no aggression in trial 2 and 3
no_agg_col_trial2 <- aggression_score_ind %>%
  filter(trial == "Trial 2") %>% 
  group_by(colony, trial) %>% 
  summarise(aggression_by_colonies = sum(totalStings_ind)) %>% 
  filter(aggression_by_colonies == 0) %>%
  dplyr::select(colony)

no_agg_col_trial3 <- aggression_score_ind %>%
  filter(trial == "Trial 3") %>% 
  group_by(colony, trial) %>% 
  summarise(aggression_by_colonies = sum(totalStings_ind)) %>% 
  filter(aggression_by_colonies == 0) %>% 
  dplyr::select(colony)

# group sizes
reshaped_aggression_1 <- anti_join(reshaped_aggression, no_agg_col_trial2)
reshaped_aggression_2 <- anti_join(reshaped_aggression_1, no_agg_col_trial3)

reshaped_aggression_1_4 <- reshaped_aggression_1[reshaped_aggression_1$group_size == "4", ]
reshaped_aggression_1_8 <- reshaped_aggression_1[reshaped_aggression_1$group_size == "8", ]

reshaped_aggression_2_4 <- reshaped_aggression_2[reshaped_aggression_2$group_size == "4", ]
reshaped_aggression_2_8 <- reshaped_aggression_2[reshaped_aggression_2$group_size == "8", ]

# brood conditions
reshaped_aggression_1_n <- reshaped_aggression_1[reshaped_aggression_1$brood == "No brood", ]
reshaped_aggression_1_l <- reshaped_aggression_1[reshaped_aggression_1$brood == "With larvae", ]
reshaped_aggression_1_p <- reshaped_aggression_1[reshaped_aggression_1$brood == "With pupae", ]

reshaped_aggression_2_n <- reshaped_aggression_2[reshaped_aggression_1$brood == "No brood", ]
reshaped_aggression_2_l <- reshaped_aggression_2[reshaped_aggression_1$brood == "With larvae", ]
reshaped_aggression_2_p <- reshaped_aggression_2[reshaped_aggression_1$brood == "With pupae", ]

# Trial 1 v.s. Trial 2
test_result <- cor.test(reshaped_aggression_1_n$`Trial 1`, reshaped_aggression_1_n$`Trial 2`, method = "spearman",exact=FALSE )
cor.test(reshaped_aggression_1_l$`Trial 1`, reshaped_aggression_1_l$`Trial 2`)
cor.test(reshaped_aggression_1_p$`Trial 1`, reshaped_aggression_1_p$`Trial 2`)

# Trial 2 v.s. Trial 3
cor.test(reshaped_aggression_2_n$`Trial 2`, reshaped_aggression_2_n$`Trial 3`)
cor.test(reshaped_aggression_2_l$`Trial 2`, reshaped_aggression_2_l$`Trial 3`)
cor.test(reshaped_aggression_2_p$`Trial 2`, reshaped_aggression_2_p$`Trial 3`)


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
result_1_4 <- extract_cor_test_results(reshaped_aggression_1_4$`Trial 1`, reshaped_aggression_1_4$`Trial 2`)
result_1_8 <- extract_cor_test_results(reshaped_aggression_1_8$`Trial 1`, reshaped_aggression_1_8$`Trial 2`)
result_2_4 <- extract_cor_test_results(reshaped_aggression_2_4$`Trial 2`, reshaped_aggression_2_4$`Trial 3`)
result_2_8 <- extract_cor_test_results(reshaped_aggression_2_8$`Trial 2`, reshaped_aggression_2_8$`Trial 3`)

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
result_1_n <- extract_cor_test_results(reshaped_aggression_1_n$`Trial 1`, reshaped_aggression_1_n$`Trial 2`)
result_1_l <- extract_cor_test_results(reshaped_aggression_1_l$`Trial 1`, reshaped_aggression_1_l$`Trial 2`)
result_1_p <- extract_cor_test_results(reshaped_aggression_1_p$`Trial 1`, reshaped_aggression_1_p$`Trial 2`)
result_2_n <- extract_cor_test_results(reshaped_aggression_2_n$`Trial 2`, reshaped_aggression_2_n$`Trial 3`)
result_2_l <- extract_cor_test_results(reshaped_aggression_2_l$`Trial 2`, reshaped_aggression_2_l$`Trial 3`)
result_2_p <- extract_cor_test_results(reshaped_aggression_2_p$`Trial 2`, reshaped_aggression_2_p$`Trial 3`)

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
file_name <- paste0("data/exp_processed/ranked_aggression_score", current_date, ".RData")
save(ranked_aggression, file = file_name)

