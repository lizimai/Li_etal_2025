rm(list = ls())
# load libraries
library(tidyverse)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggfortify)
library(ggExtra)

# load data ----
load("data/exp_processed/processed_data_20240212.RData")
load("data/exp_processed/individual_defence_score20240212.RData")
load("data/exp_processed/ranked_defence_score20240212.RData")
track4before <- read.csv("data/exp_processed/updatedAntTable_4before.csv")
track8before <- read.csv("data/exp_processed/updatedAntTable_8Before.csv")
track8before$antColour <- tolower(track8before$antColour)

# data preparation ----
manual_trial1 <- ranked_defence %>% 
  filter(trial == "Trial 1")

track <- rbind(track4before, track8before) %>% 
  dplyr::select(colony, antColour, spatialEntropies) %>% 
  rename(color = antColour)

parts <- str_split(track$colony, pattern = "_", n = 5)

# Create the new columns
track$header <- sapply(parts, function(x) paste(x[2:3], collapse = "_"))
track$group_size <- as.factor(sapply(parts, function(x) str_extract(x[5], "[0-9]+")))
track$brood <- sapply(parts, function(x) str_extract(x[5], "[a-zA-Z]+"))

track_formatted <- track %>% 
  mutate(brood = case_when(
    brood == "n"  ~ "No brood",
    brood == "wl" ~ "With larvae",
    brood == "wp" ~ "With pupae",
    TRUE          ~ brood),  # This leaves the value unchanged if it's not one of the above
  colony = paste0(header, "_",brood ,"_", group_size)) %>% 
  dplyr::select(-header)

str(manual_trial1)
str(track_formatted)
combined_manual_track <- left_join(manual_trial1, track_formatted) %>% 
  group_by(colony) %>% 
  mutate(rank_encounter =  min_rank(totalEncounters),
         rank_spatialEntropies = min_rank(spatialEntropies))

# corr rank between foraging and defence
combined_manual_track4 <- combined_manual_track %>% 
  filter(group_size == 4)
combined_manual_track8 <- combined_manual_track %>% 
  filter(group_size == 8)

# relationship between encounter and spatial entropy
lmer_encounter_entropy <- glmmTMB(totalEncounters ~ spatialEntropies * group_size * brood + (1|colony), data = combined_manual_track)
drop1(lmer_encounter_entropy, test = "Chisq")
lmer_encounter_entropy <- glmmTMB(totalEncounters ~ spatialEntropies + group_size * brood + (1|colony), data = combined_manual_track)
drop1(lmer_encounter_entropy, test = "Chisq")
lmer_encounter_entropy <- glmmTMB(totalEncounters ~ spatialEntropies + group_size + brood + (1|colony), data = combined_manual_track)
summary(lmer_encounter_entropy)
drop1(lmer_encounter_entropy, test = "Chisq") 

# relationship between defence score and spatial entropy
lmer_SpaEntro <- glmmTMB(defence_score ~ spatialEntropies + group_size + brood + (1|colony), ziformula = ~brood, data = combined_manual_track)
summary(lmer_SpaEntro)
drop1(lmer_SpaEntro, test = "Chisq")

# Save data
current_date <- format(Sys.Date(), "%Y%m%d")
file_name <- paste0("data/exp_processed/combined_manual_track", current_date, ".RData")
save(combined_manual_track, file = file_name)
