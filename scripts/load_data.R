# Installing (if needed) and loading the packages ----

if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
  
}

p_load(here, dplyr, tidyr)

setwd(here::here())


# Loading data from Garre-Frutos et al. (2024): ----
# Loading raw data of the previous experiment
raw_m <- read.csv("Input/raw_data_multi.csv")

# There are participants with less observations than the max number of observations?
check <- raw_m %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  ungroup() %>%
  complete(ID, Phase) %>% filter(Phase %in% c("Rewarded", "Unrewarded")) %>%
  filter(n < 288 | is.na(n))

length(unique(check$ID)) # N = 23 participants with missing data

# There are some participants with missing data. We will exclude those participants

raw_fm <- raw_m %>% filter(!ID %in% check$ID)

length(unique(raw_fm$ID)) # N = 193 participants

# Participants with less than .7 of accuracy?

acc_ex <- raw_fm %>%
  filter(Phase %in% c("Rewarded", "Unrewarded")) %>%
  group_by(ID) %>%
  dplyr::summarise(mean_ACC = mean(Accuracy)) %>% filter(mean_ACC < .7) %>% pull(ID)

length(unique(acc_ex)) # 11 participants with less than .7 of accuracy that will be excluded from the analysis

# How many participants?
raw_fm %>%
  filter(!ID %in% acc_ex, Phase %in% c("Rewarded", "Unrewarded")) %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  group_by(Phase) %>%
  dplyr::count() # 182 participants with complete observations and ACC > .7

raw_fm$trial_num <- raw_fm$Trial
raw_fm$rt <- raw_fm$RT
raw_fm$correct <- raw_fm$Accuracy

# Experiment 1 ----

raw <- read.csv("Input/exp_1.csv") %>%
  dplyr::rename("ID" = "Code") %>%
  mutate(Phase = ifelse(Phase == "Reward", "Rewarded", Phase))


# Experiment 2 ----

raw2 <- read.csv("Input/exp_2.csv") %>%
  dplyr::rename("ID" = "Code") %>%
  mutate(Phase = ifelse(Phase == "Reward", "Rewarded", Phase))

# Check
# There are participants with less observations than the max number of observations?
check <- raw2 %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  ungroup() %>%
  complete(ID, Phase) %>% filter(Phase %in% c("Rewarded", "Extinction")) %>%
  filter(n < 288 | is.na(n))

length(unique(check$ID)) # N = 0 participants with missing data

length(unique(raw2$ID)) # N = 71 participants

# Participants with less than .7 of accuracy?

acc_ex <- raw2 %>%
  filter(Phase %in% c("Rewarded", "Extinction")) %>%
  group_by(ID) %>%
  dplyr::summarise(mean_ACC = mean(correct)) %>% filter(mean_ACC < .7) %>% pull(ID)

length(unique(acc_ex)) # 0 participants with less than .7 of accuracy that will be excluded from the analysis

# How many participants?
raw %>%
  filter(!ID %in% acc_ex, Phase %in% c("Rewarded", "Extinction")) %>%
  group_by(ID, Phase) %>%
  dplyr::count() %>%
  group_by(Phase) %>%
  dplyr::count() # 111 participants with complete observations and ACC > .7

# Combination ----
comp <- rbind(raw %>% mutate(Exp = .5), raw2 %>% mutate(Exp = -.5)) %>%
  mutate(Experiment = ifelse(Exp == .5, "Punishment", "Omission"))
