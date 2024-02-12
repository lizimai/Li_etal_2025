rm(list = ls())
# load libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggfortify)
library(ggExtra)
library(ggh4x)

load("data/exp_processed/individual_defence_score20240212.RData")

# data preparation ----
defence_efficiency_col <- defence_score_ind %>%
  mutate(group_size_num = as.integer(as.character(group_size))) %>% 
  group_by(colony, group_size, brood, group_size_num) %>% 
  summarise(total_defence = sum(totalStings_ind)) %>% 
  mutate(defence_efficiency = total_defence) %>% 
  ungroup()

# plotting ----
sd_defence_score_ind_across_trials  %>% 
  ggplot(aes(x = group_size, y = sd_defence_score, fill = brood)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.1)) +
  coord_flip() +
  xlab("Group size") +
  ylab("Variation in defence behaviour") +
  theme_classic() +
  theme(legend.position = "none")
ggsave("plots/sd.pdf", width = 5, height = 2)

defence_efficiency_col  %>% 
  ggplot(aes(x = group_size, y = defence_efficiency, fill = brood)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, position = position_jitterdodge(jitter.width = 0.1)) +
  xlab("Group size") +
  ylab("Defence efficiency") +
  theme_classic() +
  theme(legend.position = "none")
ggsave("plots/de_ef.pdf", width = 2, height = 5)

eff_dol_sd_ac <- left_join(defence_efficiency_col, sd_defence_score_ind_across_trials)
eff_dol_sd_ac %>% 
  ggplot(aes(x = sd_defence_score, y = defence_efficiency)) +
  geom_point(aes(color = brood, shape = group_size), size = 3) +
  xlab("Variation in defence behaviour") +
  ylab("Defence efficiency") +
  theme_classic() +
  theme(legend.position = "none")
ggsave("plots/eff_dol_sum_trial_ac.pdf", width = 5, height = 5)

eff_dol_sd_ac %>% 
  ggplot(aes(x = sd_defence_score, y = total_defence)) +
  geom_point(aes(color = brood, shape = group_size), size = 3) +
  xlab("Variation in defence behaviour") +
  ylab("Defence efficiency") +
  theme_classic() +
  theme(legend.position = "none") +
  facet_wrap(~ group_size + brood)
ggsave("plots/eff_dol_sum_trial_subplot_ac.pdf", width = 5, height = 5)


# Stats ----
lm_eff_dol_int <- lm(defence_efficiency ~ sd_defence_score * brood * group_size, data = eff_dol_sd_ac)
drop1(lm_eff_dol_int, test = "F")
lm_eff_dol_int <- lm(defence_efficiency ~ sd_defence_score * brood + group_size, data = eff_dol_sd_ac)
drop1(lm_eff_dol_int, test = "F")
lm_eff_dol_int <- lm(defence_efficiency ~ sd_defence_score + brood * group_size, data = eff_dol_sd_ac)
drop1(lm_eff_dol_int, test = "F")
lm_eff_dol_int <- lm(defence_efficiency ~ sd_defence_score * group_size + brood , data = eff_dol_sd_ac)
drop1(lm_eff_dol_int, test = "F")
lm_eff_dol_final <- lm(defence_efficiency ~ sd_defence_score + group_size + brood, data = eff_dol_sd_ac)
drop1(lm_eff_dol_final, test = "F")
autoplot(lm_eff_dol_final)
summary(lm_eff_dol_final)

# posthoc analysis
post_hoc_result <- emmeans(lm_eff_dol_final, pairwise ~ brood, adjust = "tukey")
# Viewing the results
print(post_hoc_result$contrasts)

# Save data
current_date <- format(Sys.Date(), "%Y%m%d")
file_name <- paste0("data/exp_processed/eff_dol", current_date, ".RData")
save(eff_dol, file = file_name)

