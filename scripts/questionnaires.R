#Correct questionnaires (UPPSP, FCQT, MULTICAGE)
#Note - as part of a broader project, participants completed two more questionnaires:
        #FCQT (craving trait) and MULTICAGE-CAD4 (screening for eating disorders)

if (!require(tidyverse)) {
  install.packages("tidyverse")
}

if (!require(tidyverse)) {
  install.packages("psych")
}

library(psych)


#### EXPERIMENT 1 ----
qs_e1 <- read.csv("Input/QsData_e1.csv") #Load data
sociodemographic <- subset(qs_e1, select = c("ID", "Sexo", "Edad", "estudios"))


## UPPS-P 
  #select data 
uppsp = subset(qs_e1, select = c(ID, UPPSPItems.UPPS1.:UPPSPItems.UPPS20.))

  # convert reversed items
uppsp_reversed = uppsp %>% mutate(across(c("UPPSPItems.UPPS2.", "UPPSPItems.UPPS3.", "UPPSPItems.UPPS4.", 
                                           "UPPSPItems.UPPS7.", "UPPSPItems.UPPS9.", "UPPSPItems.UPPS10.", 
                                           "UPPSPItems.UPPS12.", "UPPSPItems.UPPS14.","UPPSPItems.UPPS15.", 
                                           "UPPSPItems.UPPS17.", "UPPSPItems.UPPS18.", "UPPSPItems.UPPS20."), 
                                         ~ 5 - .))
  # calculate factors
uppsp_factors = uppsp_reversed %>%
  group_by(ID) %>%
  summarise(
    urg_neg = sum(c_across(UPPSPItems.UPPS4.|UPPSPItems.UPPS7.|UPPSPItems.UPPS12.|UPPSPItems.UPPS17.)),
    urg_pos = sum(c_across(UPPSPItems.UPPS2.|UPPSPItems.UPPS10.|UPPSPItems.UPPS15.|UPPSPItems.UPPS20.)),
    prem = sum(c_across(UPPSPItems.UPPS1.|UPPSPItems.UPPS6.|UPPSPItems.UPPS13.|UPPSPItems.UPPS19.)),
    preserv = sum(c_across(UPPSPItems.UPPS5.|UPPSPItems.UPPS8.|UPPSPItems.UPPS11.|UPPSPItems.UPPS16.)),
    sen_seek = sum(c_across(UPPSPItems.UPPS3.|UPPSPItems.UPPS9.|UPPSPItems.UPPS14.|UPPSPItems.UPPS18.)),
    UPPSP = rowSums(across(c(urg_neg, urg_pos, prem, preserv, sen_seek)))
  )

  # Cronbach's alpha (reliability)
urg_neg <- select(uppsp_reversed, UPPSPItems.UPPS4.|UPPSPItems.UPPS7.|UPPSPItems.UPPS12.|UPPSPItems.UPPS17.)
a_urg_neg <- alpha(urg_neg)$total$raw_alpha #0.8099037
urg_pos <- select(uppsp_reversed, UPPSPItems.UPPS2.|UPPSPItems.UPPS10.|UPPSPItems.UPPS15.|UPPSPItems.UPPS20.)
a_urg_pos <- alpha(urg_pos)$total$raw_alpha #0.5938761


## FCQ-T 
  #select data 
fcqt = subset(qs_e1, select = c(ID, FCQTItems.FCQT1.:FCQTItems.FCQT18., FCQTItems2.FCQT19.:FCQTItems2.FCQT39.))

  #calculate factors
fcqt_factors = fcqt %>%
  group_by(ID) %>%
  summarise(
    pos_ref = sum(c_across(FCQTItems.FCQT9.|FCQTItems.FCQT10.|FCQTItems.FCQT15.|FCQTItems2.FCQT24.|FCQTItems2.FCQT37.)),
    relief_neg = sum(c_across(FCQTItems.FCQT16.|FCQTItems2.FCQT19.|FCQTItems2.FCQT21.)),
    plans = sum(c_across(FCQTItems.FCQT5.|FCQTItems.FCQT18.|FCQTItems2.FCQT23.)),
    cues = sum(c_across(FCQTItems.FCQT1.|FCQTItems2.FCQT34.|FCQTItems2.FCQT35.|FCQTItems2.FCQT36.)),
    preocup = sum(c_across(FCQTItems.FCQT6.|FCQTItems.FCQT8.|FCQTItems2.FCQT27.|FCQTItems2.FCQT28.|FCQTItems2.FCQT30.|FCQTItems2.FCQT31.|FCQTItems2.FCQT32.)),
    crav_hung = sum(c_across(FCQTItems.FCQT11.|FCQTItems.FCQT12.|FCQTItems.FCQT13.|FCQTItems.FCQT14.)),
    lack_ctrl = sum(c_across(FCQTItems.FCQT2.|FCQTItems.FCQT3.|FCQTItems2.FCQT22.|FCQTItems2.FCQT25.|FCQTItems2.FCQT26.|FCQTItems2.FCQT29.)),
    emotion = sum(c_across(FCQTItems2.FCQT20.|FCQTItems2.FCQT33.|FCQTItems2.FCQT38.|FCQTItems2.FCQT39.)),
    guilt =sum(c_across(FCQTItems.FCQT4.|FCQTItems.FCQT7.|FCQTItems.FCQT17.)),
    FCQT = rowSums(across(c(pos_ref, relief_neg, plans, cues, preocup, crav_hung, lack_ctrl, emotion, guilt)))
  )
  

## MULTICAGE-CAD4
  # select data 
multicage <- subset(qs_e1, select = c(ID, MulticageItems.MC1.:MulticageItems.MC4.))

  # columns to be corrected
columns_to_correct <- c("MulticageItems.MC1.", "MulticageItems.MC2.", 
                        "MulticageItems.MC3.", "MulticageItems.MC4.")

  # define a function to correct the responses
correct_responses <- function(response) {
  ifelse(response == "SI", 1, ifelse(response == "NO", 0, response))
}

  # apply corrections to the specified columns
multicage_corrected <- multicage %>%
  mutate(across(all_of(columns_to_correct), correct_responses)) %>%
  mutate(across(all_of(columns_to_correct), as.numeric)) %>%
  rename(
    induceVomiting = MulticageItems.MC1.,
    lossOfControl = MulticageItems.MC2.,
    bodyPerception = MulticageItems.MC3.,
    foodObsession = MulticageItems.MC4.
  ) %>%
  arrange(ID)

  # calculate the total score 
multicage_factor <- multicage_corrected %>%
  group_by(ID) %>%
  summarise(
    MulticageTotal = sum(c(induceVomiting, lossOfControl, bodyPerception, foodObsession), na.rm = TRUE)
  )

# Join all Qs in one data frame
dfs <- list(sociodemographic, uppsp_factors, fcqt_factors, multicage_factor)
combined_df <- reduce(dfs, full_join, by = "ID") %>%
  arrange(ID)

write.csv(combined_df, "Output/Qs_e1.csv")



#### EXPERIMENT 2 ----
qs_e2 <- read.csv("Input/QsData_e2.csv") #Load data
sociodemographic <- subset(qs_e2, select = c("ID", "Sexo", "Edad", "estudios"))


## UPPS-P 
  #select data 
uppsp = subset(qs_e2, select = c(ID, UPPSPItems.UPPS1.:UPPSPItems.UPPS20.))

  # convert reversed items
uppsp_reversed = uppsp %>% mutate(across(c("UPPSPItems.UPPS2.", "UPPSPItems.UPPS3.", "UPPSPItems.UPPS4.", 
                                           "UPPSPItems.UPPS7.", "UPPSPItems.UPPS9.", "UPPSPItems.UPPS10.", 
                                           "UPPSPItems.UPPS12.", "UPPSPItems.UPPS14.","UPPSPItems.UPPS15.", 
                                           "UPPSPItems.UPPS17.", "UPPSPItems.UPPS18.", "UPPSPItems.UPPS20."), 
                                         ~ 5 - .))
  # calculate factors
uppsp_factors = uppsp_reversed %>%
  group_by(ID) %>%
  summarise(
    urg_neg = sum(c_across(UPPSPItems.UPPS4.|UPPSPItems.UPPS7.|UPPSPItems.UPPS12.|UPPSPItems.UPPS17.)),
    urg_pos = sum(c_across(UPPSPItems.UPPS2.|UPPSPItems.UPPS10.|UPPSPItems.UPPS15.|UPPSPItems.UPPS20.)),
    prem = sum(c_across(UPPSPItems.UPPS1.|UPPSPItems.UPPS6.|UPPSPItems.UPPS13.|UPPSPItems.UPPS19.)),
    preserv = sum(c_across(UPPSPItems.UPPS5.|UPPSPItems.UPPS8.|UPPSPItems.UPPS11.|UPPSPItems.UPPS16.)),
    sen_seek = sum(c_across(UPPSPItems.UPPS3.|UPPSPItems.UPPS9.|UPPSPItems.UPPS14.|UPPSPItems.UPPS18.)),
    UPPSP = rowSums(across(c(urg_neg, urg_pos, prem, preserv, sen_seek)))
  )

  # Cronbach's alpha (reliability)
urg_neg <- select(uppsp_reversed, UPPSPItems.UPPS4.|UPPSPItems.UPPS7.|UPPSPItems.UPPS12.|UPPSPItems.UPPS17.)
a_urg_neg <- alpha(urg_neg)$total$raw_alpha #0.8007531
urg_pos <- select(uppsp_reversed, UPPSPItems.UPPS2.|UPPSPItems.UPPS10.|UPPSPItems.UPPS15.|UPPSPItems.UPPS20.)
a_urg_pos <- alpha(urg_pos)$total$raw_alpha #0.7411139


## FCQ-T 
  #select data 
fcqt = subset(qs_e2, select = c(ID, FCQTItems.FCQT1.:FCQTItems.FCQT18., FCQTItems2.FCQT19.:FCQTItems2.FCQT39.))

  #calculate factors
fcqt_factors = fcqt %>%
  group_by(ID) %>%
  summarise(
    pos_ref = sum(c_across(FCQTItems.FCQT9.|FCQTItems.FCQT10.|FCQTItems.FCQT15.|FCQTItems2.FCQT24.|FCQTItems2.FCQT37.)),
    relief_neg = sum(c_across(FCQTItems.FCQT16.|FCQTItems2.FCQT19.|FCQTItems2.FCQT21.)),
    plans = sum(c_across(FCQTItems.FCQT5.|FCQTItems.FCQT18.|FCQTItems2.FCQT23.)),
    cues = sum(c_across(FCQTItems.FCQT1.|FCQTItems2.FCQT34.|FCQTItems2.FCQT35.|FCQTItems2.FCQT36.)),
    preocup = sum(c_across(FCQTItems.FCQT6.|FCQTItems.FCQT8.|FCQTItems2.FCQT27.|FCQTItems2.FCQT28.|FCQTItems2.FCQT30.|FCQTItems2.FCQT31.|FCQTItems2.FCQT32.)),
    crav_hung = sum(c_across(FCQTItems.FCQT11.|FCQTItems.FCQT12.|FCQTItems.FCQT13.|FCQTItems.FCQT14.)),
    lack_ctrl = sum(c_across(FCQTItems.FCQT2.|FCQTItems.FCQT3.|FCQTItems2.FCQT22.|FCQTItems2.FCQT25.|FCQTItems2.FCQT26.|FCQTItems2.FCQT29.)),
    emotion = sum(c_across(FCQTItems2.FCQT20.|FCQTItems2.FCQT33.|FCQTItems2.FCQT38.|FCQTItems2.FCQT39.)),
    guilt =sum(c_across(FCQTItems.FCQT4.|FCQTItems.FCQT7.|FCQTItems.FCQT17.)),
    FCQT = rowSums(across(c(pos_ref, relief_neg, plans, cues, preocup, crav_hung, lack_ctrl, emotion, guilt)))
  )


## MULTICAGE-CAD4
  # select data 
multicage <- subset(qs_e2, select = c(ID, MulticageItems.MC1.:MulticageItems.MC4.))

  # apply corrections to the specified columns
multicage_corrected <- multicage %>%
  mutate(across(all_of(columns_to_correct), correct_responses)) %>%
  mutate(across(all_of(columns_to_correct), as.numeric)) %>%
  rename(
    induceVomiting = MulticageItems.MC1.,
    lossOfControl = MulticageItems.MC2.,
    bodyPerception = MulticageItems.MC3.,
    foodObsession = MulticageItems.MC4.
  ) %>%
  arrange(ID)

  # calculate the total score 
multicage_factor <- multicage_corrected %>%
  group_by(ID) %>%
  summarise(
    MulticageTotal = sum(c(induceVomiting, lossOfControl, bodyPerception, foodObsession), na.rm = TRUE)
  )

# Join all Qs in one data frame
dfs <- list(sociodemographic, uppsp_factors, fcqt_factors, multicage_factor)
combined_df <- reduce(dfs, full_join, by = "ID") %>%
  arrange(ID)

write.csv(combined_df, "Output/Qs_e2.csv")
