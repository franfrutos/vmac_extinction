
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
       grid,
       gridExtra
)

# Comparison: ----
d_RT <-
  filter_data(comp, f.absent = T, fixed = T, phase = c("Rewarded", "Extinction"))

d_RT$log_RT <- log(d_RT$rt)
d_RT$Singleton <- ifelse(d_RT$Singleton == "High", .5, -.5)

fit_comp <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num))*Exp + (Singleton+scale(log(Block_num)) | ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d_RT[d_RT$Phase == "Rewarded",]
  )

summary(fit_comp)

fit_comp2 <-
  lmer(
    log_RT ~ Singleton * scale(log(Block_num))*Exp + (Singleton+scale(Block_num) | ID),
    control = lmerControl(optimizer = 'bobyqa'),
    data = d_RT[d_RT$Phase == "Extinction",]
  )

summary(fit_comp2)

# Differences in last blocks?

d_RT[d_RT$Block_num %in% 23:24, ] %>%
  mutate(Experiment = ifelse(Exp == .5, "Experiment 1", "Experiment 2"),
         Singleton = ifelse(Singleton == .5, "High", "Low")) %>%
  dplyr::summarise(RT = mean(rt), .by = c("Experiment", "ID", "Singleton")) %>%
  aov_ez(dv = "RT", within = "Singleton", between = "Experiment", id = "ID", anova_table = list(es = "pes"))

# Accuracy: ----
d_accComp <-
  filter_data(comp, f.absent = T, fixed = T, phase = c("Rewarded" , "Extinction"), acc = F) %>%
  dplyr::summarise(correct = mean(correct),
                   N = n(), .by = c(Exp, ID, Phase, Singleton, Block_num))

d_accComp$Singleton <- ifelse(d_accComp$Singleton == "High", .5, -.5)

fit_AccComp<-
  glmer(
    correct ~ Singleton * scale(Block_num)*Exp + (Singleton + scale(Block_num)| ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_accComp[d_accComp$Phase == "Rewarded",],
    weights = N,
    family = binomial()
  )

summary(fit_AccComp)

fit_AccComp2<-
  glmer(
    correct ~ Singleton * scale(Block_num)*Exp + (Singleton | ID),
    control = glmerControl(optimizer = 'bobyqa'),
    data = d_accComp[d_accComp$Phase == "Extinction",],
    weights = N,
    family = binomial()
  )

summary(fit_AccComp2)

# Plots ----

comps1 <- get_comparisons(d1 = filter_data(comp, f.absent = F, phase = c("Rewarded", "Extinction")),
                     m1 = fit_comp, m2 = fit_comp2, both = T)

levels(comps1$mod$Effect) <- c("VMAC: High - Low", "Attentional Capture: Low - Absent")
levels(comps1$raw$Effect) <- c("VMAC: High - Low", "Attentional Capture: Low - Absent")

comps1$mod$Experiment <- ifelse(comps1$mod$Exp == .5, "Punishment", "Omission")
comps1$mod$Experiment <- factor(comps1$mod$Experiment, levels=c("Punishment", "Omission"))
comps1$raw$Experiment <- factor(comps1$raw$Experiment, levels=c("Punishment", "Omission"))


comps2 <- get_comparisons(d1 = filter_data(comp, f.absent = F, phase = c("Rewarded", "Extinction"), acc = F),
                     m1 = fit_AccComp, m2 = fit_AccComp2, both = T, acc = T)


levels(comps2$mod$Effect) <- c("VMAC: High - Low", "Attentional Capture: Low - Absent")
levels(comps2$raw$Effect) <- c("VMAC: High - Low", "Attentional Capture: Low - Absent")

comps2$mod$Experiment <- ifelse(comps2$mod$Exp == -.5, "Omission", "Punishment")
comps2$mod$Experiment <- factor(comps2$mod$Experiment, levels=c("Punishment", "Omission"))
comps2$raw$Experiment <- factor(comps2$raw$Experiment, levels=c("Punishment", "Omission"))


darkColors <- c("Punishment" = "#29285a", "Omission" = "#9592da") 

# Plot predictions
labels <- list(
  'Extinction'="Acquisition",
  'Rewarded'="Extinction"
)
phase_labeller <- function(variable,value){
  return(labels[value])
}

p1 <- ggplot(data = comps1[["mod"]],
                      aes(
                        y = estimate,
                        x = Block,
                        color = Experiment,
                        fill = Experiment,
                        shape = Experiment
                      )) +
  geom_line(aes(linetype = Experiment)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA) +
  facet_wrap( ~ rev(Phase), scales = "free_x", labeller = phase_labeller) +
  geom_point(data = comps1[["raw"]],
             aes(x = as.numeric(Block)),
             position = position_dodge(.5)) +
  geom_errorbar(
    data = comps1[["raw"]],
    aes(
      x = as.numeric(Block),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_manual(values = darkColors, name = "Experiment") +
  scale_fill_manual(values = darkColors, name = "Experiment") +
  labs(y = "RT Contrast (ms)", x = "Block") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 10) + # text size is adjusted for DPI
  theme(legend.position = "none")

p2 <- ggplot(data = comps2[["mod"]],
             aes(
               y = estimate,
               x = Block,
               color = Experiment,
               fill = Experiment,
               shape = Experiment
             )) +
  geom_line(aes(linetype = Experiment)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .3,
              color = NA) +
  facet_wrap( ~ rev(Phase), scales = "free_x", labeller = phase_labeller) +
  geom_point(data = comps2[["raw"]],
             aes(x = as.numeric(Block)),
             position = position_dodge(.5)) +
  geom_errorbar(
    data = comps2[["raw"]],
    aes(
      x = as.numeric(Block),
      ymin = estimate - se,
      ymax = estimate + se
    ),
    position = position_dodge(.5),
    width = .01
  ) +
  scale_x_continuous(breaks = seq(1, 24, 1)) +
  #scale_y_continuous(breaks = seq(-20, 60, 20)) +
  scale_color_manual(values = darkColors, name = "Experiment") +
  scale_fill_manual(values = darkColors, name = "Experiment") +
  labs(y = "Accuracy contrast (difference in proportions)", x = "Block") +
  geom_hline(yintercept = 0, linetype = "dashed")+
  theme_Publication(text_size = 10) + # text size is adjusted for DPI
  theme(legend.spacing.x = unit(.2, 'cm'))

figure4 <- grid.arrange(p1, p2, nrow = 2)

ggsave(
  "Output/plots/Figure4.png",
  plot = figure4,
  height = 20,
  width = 15,
  dpi = 1200,
  units = "cm"
)

# SAT analysis ----

d_RT <-
  filter_data(comp, f.absent = F, fixed = T, phase = c("Rewarded", "Extinction"), acc = F)

df_caf <- d_RT %>%
  dplyr::summarise(Accuracy = mean(correct), RT_mean = mean(rt), .by = c(Experiment, Phase, ID, Singleton))

# Calcular la media de cada sujeto dentro de cada condición

df_caf_mean1 <- summarySEwithin2(
  df_caf,
  measurevar = "Accuracy",
  betweenvars = "Experiment",
  withinvars = c("Phase", "Singleton")
) %>% mutate(ci_lowA = Accuracy - se,
             ci_highA = Accuracy + se) 

df_caf_mean2 <- summarySEwithin2(
  df_caf,
  measurevar = "RT_mean",
  betweenvars = "Experiment",
  withinvars = c("Phase", "Singleton")
) %>% mutate(ci_lowRT = RT_mean - se,
             ci_highRT = RT_mean + se) 

df_caf_mean <- cbind(df_caf_mean1[, c("Experiment", "Phase", "Singleton", "Accuracy", "ci_lowA", "ci_highA")],
                     df_caf_mean2[, c("RT_mean", "ci_lowRT", "ci_highRT")])


df_caf_mean$Phase <- ifelse(df_caf_mean$Phase == "Rewarded", "Acquisition", "Extinction")
df_caf_mean$Singleton <- factor(df_caf_mean$Singleton, levels = c("High", "Low", "Absent"))
df_caf_mean$Experiment <- factor(df_caf_mean$Experiment, levels = c("Punishment", "Omission"))

Colors <- c("High" = "#330023", "Low" = "#7d0000", "Absent" = "#df5800") 

ggplot(df_caf_mean, aes(x = RT_mean, y = Accuracy, color = Singleton)) +
  geom_line(data = df_caf_mean, aes(x = RT_mean, y = Accuracy, group = interaction(Experiment)), color = "gray", size = .5) +
  geom_point(data = df_caf_mean, aes(x = RT_mean, y = Accuracy, color = Singleton, shape = Experiment), size = 3) +
  geom_errorbar(aes(ymin = ci_lowA, ymax = ci_highA), width = 0) +
  geom_errorbarh(aes(xmin = ci_lowRT, xmax = ci_highRT)) +
  theme_minimal() +
  facet_wrap(.~Phase, scales = "free_x") +
  scale_color_manual(values = Colors, name = "Singleton type") +
  labs(x = "Response time (ms)", y = "Accuracy (Proportion)") +
  scale_y_continuous(limits = c(.9, 1), breaks = seq(.9, 1, .02)) +
  theme_Publication(text_size = 10) + # text size is adjusted for DPI
  theme(legend.spacing.x = unit(.2, 'cm'))

ggsave(
  "Output/plots/Figure5.png",
  height = 12,
  width = 15,
  dpi = 1200,
  units = "cm"
)

