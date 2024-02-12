rm(list = ls())
# 1. Load libraries ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggsignif)
library(lme4)
library(lmerTest)
library(DHARMa)
library(emmeans)
library(glmmTMB)
library(multcomp)
library(car)

# 2. Load data ----
data_files <- list(
  data29_3_1 = "data/exp_raw/qi_cam29_exp3_a1.csv",
  data34_3_1 = "data/exp_raw/qi_cam34_exp3_a1.csv",
  data42_3_1 = "data/exp_raw/qi_cam42_exp3_a1.csv",
  data47_3_1 = "data/exp_raw/qi_cam47_exp3_a1.csv",
  data54_3_1 = "data/exp_raw/qi_cam54_exp3_a1.csv",
  data61_3_1 = "data/exp_raw/qi_cam61_exp3_a1.csv",
  
  data29_3_2 = "data/exp_raw/qi_cam29_exp3_a2.csv",
  data34_3_2 = "data/exp_raw/qi_cam34_exp3_a2.csv",
  data42_3_2 = "data/exp_raw/qi_cam42_exp3_a2.csv",
  data47_3_2 = "data/exp_raw/qi_cam47_exp3_a2.csv",
  data54_3_2 = "data/exp_raw/qi_cam54_exp3_a2.csv",
  data61_3_2 = "data/exp_raw/qi_cam61_exp3_a2.csv",
  
  data29_3_3 = "data/exp_raw/qi_cam29_exp3_a3.csv",
  data34_3_3 = "data/exp_raw/qi_cam34_exp3_a3.csv",
  data42_3_3 = "data/exp_raw/qi_cam42_exp3_a3.csv",
  data47_3_3 = "data/exp_raw/qi_cam47_exp3_a3.csv",
  data54_3_3 = "data/exp_raw/qi_cam54_exp3_a3.csv",
  data61_3_3 = "data/exp_raw/qi_cam61_exp3_a3.csv"
)

# 3. Process and combine data ----
process_data <- function(files, trial_number) {
  data_list <- lapply(files, read.csv, header = TRUE, sep = ",")
  data_combined <- bind_rows(data_list)
  return(data_combined)
}

data <- process_data(data_files)

# 4. Clear everything except the final combined "data" dataset ----
rm(list=setdiff(ls(), "data"))

# 5. Format data ----
data$number.of.ants <- as.factor(data$number.of.ants)
str(data)
head(data)

# change variable names
data_formatted <- data %>%
  separate(Observation.id, c("qi", "camera", "experiment", "trial", "treatment"), sep = "_") %>%
  mutate(
    obsID = paste(camera, treatment, experiment, trial, sep = "_"),
    antID = paste(camera, treatment, Subject, sep = "_")
  ) %>%
  dplyr::select(
    obsID, antID, camera, experiment, trial, treatment, number.of.ants, 
    brood.or.not, Subject, Behavior, Behavior.type, Start..s., Stop..s., Duration..s.
  ) %>%
  rename(
    group_size = number.of.ants, brood = brood.or.not, color = Subject, 
    behavior = Behavior, behaviorType = Behavior.type, 
    start = Start..s., stop = Stop..s., duration = Duration..s.
  ) %>% 
  mutate(brood = case_when(
    brood == "n"  ~ "No brood",
    brood == "wl" ~ "With larvae",
    brood == "wp" ~ "With pupae",
    TRUE          ~ brood  # This leaves the value unchanged if it's not one of the above
  ))

data_formatted$trial <- paste0("Trial ", sub("a", "", data_formatted$trial))
data_formatted$colony <- with(data_formatted, paste(camera, experiment, brood, group_size, sep="_"))

# Filter based on intruder's stop time:
intruder <- data_formatted %>%
  filter(behavior == "intruder") %>%
  mutate(endTime = start + 1200) %>%
  dplyr::select(obsID, endTime)

data_formatted_cut <- left_join(data_formatted, intruder, by = "obsID") %>%
  filter(stop <= endTime) %>%
  dplyr::select(-endTime)

# 6. Create complete ant list and complete colony list - for downstream analysis (in case some ants does not sting and were excluded) ----
complete_ant_list <- data_formatted %>%
  dplyr::select(colony, brood, group_size, color, trial) %>%
  filter(color != "No focal subject") %>%
  unique()

complete_colony_trial_list <- data_formatted %>%
  dplyr::select(colony, trial, camera, experiment, trial, treatment, group_size, brood) %>%
  arrange() %>% 
  unique()

# 7. Save the data with the appended date in the filename ----
current_date <- format(Sys.Date(), "%Y%m%d")
file_name <- paste0("data/exp_processed/processed_data_", current_date, ".RData")
save(data_formatted, data_formatted_cut, complete_ant_list, complete_colony_trial_list, file = file_name)
str(data_formatted_cut)

