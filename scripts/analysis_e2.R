# Loading packages----

source("scripts/functions.R") # This will load functions and data into the workspace

if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(dplyr,# Data wrangling
       lme4,# LMMs
       lmerTest, # P-values in summary
       marginaleffects,# Model predictions and conditional effects
       hypr,# Contrast matrix
       sjPlot,# Tables
       afex,# ANOVA
       Rmisc,# Function to averaged the data in within-subject designs
       effectsize, # Cohen's d
       BayesFactor,
       ggplot2, # Plotting
       emmeans
)

# RT analysis ----
# Rewarded phase:
d_RT <-
  filter_data(raw2, f.absent = F, fixed = T, phase = "Rewarded")

d_RT$log_RT <- log(d_RT$rt)
d_RT$Singleton <- factor(d_RT$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_RT$Singleton) <- HcRep


fit_r <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d_RT
  )



fit_r1 <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d_RT
  )

summary(fit_r1) # Maximal feasible model

fit_r1_power <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton + scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d_RT
  )

summary(fit_r1_power) # power function model

anova(fit_r1, fit_r1_power) # AIC

## |Difference in AIC|: 340.7 in favor of power function model

# Selected model predictions
predictions(
  fit_r1_power,
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  transform = \(x) exp(x + (sigma(fit_r1_power) ^ 2) / 2),
  by = "Singleton"
)

# VMAC and AC effects:
avg_comparisons(
  fit_r1_power,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(fit_r1_power) ^ 2) / 2) - exp(lo + (sigma(fit_r1_power) ^ 2) / 2)
)

# Extinction phase: 
d_RTe <- filter_data(raw2, f.absent = F, phase = "Extinction", fixed = T)
d_RTe$log_RT <- log(d_RTe$rt)
d_RTe$Singleton <- factor(d_RTe$Singleton, levels = c("High", "Low", "Absent"))
contrasts(d_RTe$Singleton) <- HcRep

fit_e <-
  lmer(
    log_RT ~ Singleton * scale(Block_num) + (Singleton * scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d_RTe
  )

summary(fit_e) # Singular fit

fit_e1 <-
  lmer(
    log_RT ~ Singleton * scale(Block_num)  + (Singleton + scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d_RTe
  )

summary(fit_e1) # Maximal model that converges

fit_e1_power <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num)) + (Singleton + scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d_RTe
  )

summary(fit_e1_power) # power function model

anova(fit_e1, fit_e1_power) # AIC

## |Difference in AIC|: 12 in favor of the linear model

# Selected model predictions
predictions(
  fit_e1,
  newdata = datagrid(
    Singleton = unique,
    Block_num = unique,
    ID = NA
  ),
  re.form = NA,
  transform = \(x) exp(x + (sigma(fit_r1_power) ^ 2) / 2),
  by = "Singleton"
)

# VMAC and AC effects:
avg_comparisons(
  fit_e1,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(Block_num = unique,
                     Singleton = unique,
                     ID = NA),
  re.form = NA,
  comparison = \(hi, lo) exp(hi + (sigma(fit_r1_power) ^ 2) /
                                  2) - exp(lo + (sigma(fit_r1_power) ^ 2) / 2)
)


# Test if the VMAC effect is present in the last epoch (blocks 23 and 24).
 d_RTe[d_RTe$Block_num %in% 23:24 & d_RTe$Singleton != "Absent", ] %>%
  dplyr::summarise(RT = mean(rt), .by = c("ID", "Singleton")) %>%
  t.test(RT ~ Singleton, data = ., paired = T)

d_RTe[d_RTe$Block_num %in% 23:24 & d_RTe$Singleton != "Absent", ] %>%
  dplyr::summarise(RT = mean(rt), .by = c("ID", "Singleton")) %>%
  repeated_measures_d(RT ~ Singleton | ID, data = ., paired = T, method = "z") # cohen's d

# Bayes factor:

d_bayes <- d_RTe[d_RTe$Block_num %in% 23:24 & d_RTe$Singleton != "Absent", ] %>%
  dplyr::summarise(RT = mean(rt), .by = c("ID", "Singleton")) 

1/ttestBF(d_bayes$RT[d_bayes$Singleton == "High"], 
        d_bayes$RT[d_bayes$Singleton == "Low"], paired = T) # Bayes factor for the null

# A bayes factor of 4.18 in favor of the null


# Difference between phases? ----

phase_aovE1 <- filter_data(raw, f.absent = T, fixed = T, phase = c("Rewarded", "Extinction")) %>%
.[.$Block_num %in% 11:14, ] %>%
  dplyr::summarise(RT = mean(rt), .by = c("Phase", "ID", "Singleton")) %>%
  aov_ez(dv = "RT", within = c("Singleton", "Phase"), id = "ID", anova_table = list(es = "pes"))

emmeans(phase_aovE1, ~ Singleton | Phase) %>% pairs()
emmeans(phase_aovE1, ~ Singleton | Phase) %>% pairs() %>% confint()

emm_options(opt.digits = F)
t_to_d(7.019, df_error = 102, paired = T) # cohen's d for Acquisition
t_to_d(5.671, df_error = 102, paired = T) # cohen's d for Acquisition

phase_aovE2 <- filter_data(raw2, f.absent = T, fixed = T, phase = c("Rewarded", "Extinction")) %>%
  .[.$Block_num %in% 11:14, ] %>%
  dplyr::summarise(RT = mean(rt), .by = c("Phase", "ID", "Singleton")) %>%
  aov_ez(dv = "RT", within = c("Singleton", "Phase"), id = "ID", anova_table = list(es = "pes"))

emmeans(phase_aovE2, ~ Singleton | Phase) %>% pairs()
emmeans(phase_aovE2, ~ Singleton | Phase) %>% pairs() %>% confint()

emm_options(opt.digits = F)
t_to_d(4.883, df_error = 110, paired = T) # cohen's d for Acquisition
t_to_d(3.532, df_error = 110, paired = T) # cohen's d for Acquisition


 # ACC analysis----
# Rewarded phase:
d_acc <-
  filter_data(raw2, f.absent = F, acc = F, fixed = T)
d_acc$Singleton <-
  factor(d_acc$Singleton, levels = c("High", "Low", "Absent"))

# Binomial trick for performance
d_acc <- d_acc %>%
  dplyr::summarise(correct = mean(correct),
                   N = n(), .by = c(ID, Singleton, Block_num))

contrasts(d_acc$Singleton) <- contr.hypothesis(HcRep)

fit.acc <-
  glmer(
    correct ~ Singleton * Block_num + (Singleton * Block_num | ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_acc,
    weights = N, # Number of trials
    family = binomial()
  ) # Singular fit

fit.acc2 <-
  glmer(
    correct ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_acc,
    weights = N,
    family = binomial()
  ) # singular fit

fit.acc3 <-
  glmer(
    correct ~ Singleton * scale(Block_num) + (scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_acc,
    weights = N,
    family = binomial()
  ) # maximal model

summary(fit.acc3)

# Overall acc
avg_predictions(
  fit.acc3,
  newdata = datagrid(
    Singleton = unique,
    Block_num = seq(1, 12, 1),
    ID = NA
  ),
  re.form = NA
)# Predicted accuracy averaged over singleton and blocks

# Predictions for each Singleton:
predictions(
  fit.acc3,
  newdata = datagrid(
    Singleton = unique,
    Block_num = seq(1, 12, 1),
    ID = NA
  ),
  re.form = NA,
  by = "Singleton"
)# Predicted accuracy by Singleton averaged over block

avg_comparisons(
  fit.acc3,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(Block_num = seq(1, 12, 1),
                     ID = NA),
  re.form = NA
  ) 

# Extinction phase
d_acce <-
  filter_data(raw2, f.absent = F, acc = F, phase = "Extinction", fixed = T)
d_acce$Singleton <-
  factor(d_acce$Singleton, levels = c("High", "Low", "Absent"))
d_acce <- d_acce %>%
  dplyr::summarise(correct = mean(correct),
                   N = n(), .by = c(ID, Singleton, Block_num))

contrasts(d_acce$Singleton) <- contr.hypothesis(HcRep)


fit.acc_ex <-
  glmer(
    correct ~ Singleton * scale(Block_num) + (Singleton + scale(Block_num) | ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_acce,
    weights = N,
    family = binomial()
  )

fit.acc_ex2 <-
  glmer(
    correct ~ Singleton * scale(Block_num) + (Singleton | ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_acce,
    weights = N,
    family = binomial()
  )

fit.acc_ex3 <-
  glmer(
    correct ~ Singleton * scale(Block_num) + (1 | ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_acce,
    weights = N,
    family = binomial()
  )

summary(fit.acc_ex3)

avg_predictions(
  fit.acc_ex3,
  newdata = datagrid(
    Singleton = unique,
    Block_num = seq(13, 24, 1),
    ID = NA
  ),
  re.form = NA
)# Predicted accuracy averaged over singleton and blocks

# Predictions for each Singleton:
predictions(
  fit.acc_ex3,
  newdata = datagrid(
    Singleton = unique,
    Block_num = seq(13, 24, 1),
    ID = NA
  ),
  re.form = NA,
  by = c("Singleton")
) # Predicted accuracy by Singleton averaged over block


# Contrasts for high-low and low-absent:
avg_comparisons(
  fit.acc_ex3,
  variables = list(Singleton = "revsequential"),
  newdata = datagrid(Block_num = seq(13, 24, 1),
                     ID = NA),
  re.form = NA
)


# Plots ----
# Model predictions:

# Returns a list with raw2 averaged data and model preds
preds <- get_predictions(d_RT, d_RTe, fit_r1_power, fit_e1)

# Plot predictions
labels <- list(
  'Extinction'="Acquisition",
  'Rewarded'="Extinction"
)
phase_labeller <- function(variable,value){
  return(labels[value])
}

Colors <- c("High" = "#330023", "Low" = "#7d0000", "Absent" = "#df5800") 

predictions <- ggplot(data = preds[["mod"]],
                      aes(
                        y = estimate,
                        x = Block,
                        color = Singleton,
                        fill = Singleton,
                        #linetype = Singleton
                      )) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA) +
  facet_wrap( ~ rev(Phase), scales = "free_x", labeller = phase_labeller) +
  geom_point(data = preds[["raw"]],
             aes(x = as.numeric(Block), shape = Singleton),
             size = 1.5,
             position = position_dodge(.5)) +
  geom_errorbar(
    data = preds[["raw"]],
    aes(
      x = as.numeric(Block),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(600, 900, 50)) +
  scale_color_manual(values = Colors, name = "Singleton type") +
  scale_fill_manual(values = Colors, name = "Singleton type") +
  scale_shape_manual(values = c(16, 17, 18), name = "Singleton type") +
  labs(y = "Response time (ms)", x = "Block") +
  theme_Publication(text_size = 10) + 
  theme(legend.spacing.x = unit(.2, 'cm'))


# Plot model comparisons:
# Returns a list with raw averaged VMAC and ACC effects and the conditional effect of
# Singleton as a function of block
comps <- get_comparisons(d_RT, d_RTe, fit_r1_power, fit_e1)

levels(comps$mod$Effect) <- c("VMAC: High - Low", "Attentional Capture: Low - Absent")
levels(comps$raw$Effect) <- c("VMAC: High - Low", "Attentional Capture: Low - Absent")

darkColors <- c("VMAC: High - Low" = "#252373", "Attentional Capture: Low - Absent" = "#d16606") # Lines
pastelColors <- c("VMAC: High - Low" = "#3f3d8b", "Attentional Capture: Low - Absent" = "#fc963b") # Fill

comparisons <- ggplot(data = comps[["mod"]],
                      aes(
                        y = estimate,
                        x = Block,
                        color = Effect,
                        fill = Effect,
                        shape = Effect,
                        linetype = Effect
                      )) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA) +
  facet_wrap( ~ rev(Phase), scales = "free_x", labeller = phase_labeller) +
  geom_point(data = comps[["raw"]],
             aes(x = as.numeric(Block), shape = Effect),
             position = position_dodge(.5)) +
  geom_errorbar(
    data = comps[["raw"]],
    aes(
      x = as.numeric(Block),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01,
    linetype = "solid",
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_manual(values = darkColors, name = "Contrast") +
  scale_fill_manual(values = pastelColors, name = "Contrast") +
  scale_linetype_manual(values = c("solid", "dashed"), name = "Singleton type") + 
  labs(y = "RT Contrast (ms)", x = "Block") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 10) + # text size is adjusted for DPI
  theme(legend.spacing.x = unit(.2, 'cm')) +
  guides(
    fill = guide_legend(title = "Contrast"),  
    color = guide_legend(title = "Contrast"), 
    shape = "none",  
    linetype = "none"  
  )


#Save both plots in one figure
figure2 <- grid.arrange(predictions, comparisons, nrow = 2)

ggsave(
  "plots/Figure2.png",
  plot = figure2,
  height = 20,
  width = 15,
  dpi = 1200,
  units = "cm"
)