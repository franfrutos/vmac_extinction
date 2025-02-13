#### CORRELATIONS BETWEEN VMAC EXTINCTION AND INDIVIDUAL DIFFERENCES (QS)

source("scripts/functions.R") # This will load functions and data into the workspace

library(tidyverse)

# EXPERIMENT 1 ----
e1 <- 
  filter_data(raw, f.absent = T, fixed = T, sd_filter = 3, phase = c("Rewarded", "Extinction"))

vmac_e1 <- e1 %>%
  group_by(ID, Phase, Block_num, Singleton) %>%
  summarise(mean_RT = mean(rt, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Singleton, values_from = mean_RT) %>%
  mutate(VMAC = High - Low) %>%
  group_by(ID, Phase) %>%
  summarise(mean_vmac = mean(VMAC), .groups = "drop") %>%
  pivot_wider(names_from = Phase, values_from = mean_vmac)


# Split-half reliability for VMAC ---- 
library(splithalf)

diff_reliability <- splithalf(data = e1,
                          outcome = "RT",
                          score = "difference",
                          conditionlist = c("Rewarded", "Extinction"),
                          halftype = "random",
                          permutations = 5000,
                          var.RT = "rt",
                          var.condition = "Phase",
                          var.participant = "ID",
                          var.compare = "Singleton",
                          compare1 = "High",
                          compare2 = "Low",
                          average = "mean")



# Read UPPS-P subscales data
qs_e1 <- read.csv("Output/Qs_e1.csv") 
corr_e1 <- vmac_e1 %>%
  merge(qs_e1, by = "ID")

# Check if normal distribution (corr_e1$Rewarded, corr_e1$Extinction, corr_e1$urg_neg, corr_e1$urg_pos)
qqnorm(corr_e1$Rewarded) # Quantile-Quantile
qqline(corr_e1$Rewarded, col = "red")

ggplot(corr_e1, aes(x = Rewarded)) + # Visualize distribution
  geom_histogram(aes(y = ..density..), bins = 30, color = "black", fill = "lightblue") +
  geom_density(alpha = 0.2, fill = "blue")


# Correlations
e1_NU_acq <- cor.test(corr_e1$Rewarded, corr_e1$urg_neg, method = "pearson") # r = 0.006055163; p-value = 0.9516
e1_PU_acq <- cor.test(corr_e1$Rewarded, corr_e1$urg_pos, method = "pearson") # r = -0.01328742; p-value = 0.1807019

e1_NU_ext <- cor.test(corr_e1$Extinction, corr_e1$urg_neg, method = "pearson") # r = -0.008228317; p-value = 0.1855919
e1_PU_ext <- cor.test(corr_e1$Extinction, corr_e1$urg_pos, method = "pearson") # r = 0.04032828; p-value = 0.2320420



# EXPERIMENT 2 ----
e2 <- 
  filter_data(raw2, f.absent = T, fixed = T, sd_filter = 3, phase = c("Rewarded", "Extinction"))

vmac_e2 <- e2 %>%
  group_by(ID, Phase, Block_num, Singleton) %>%
  summarise(mean_RT = mean(rt, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Singleton, values_from = mean_RT) %>%
  mutate(VMAC = High - Low) %>%
  group_by(ID, Phase) %>%
  summarise(mean_vmac = mean(VMAC), .groups = "drop") %>%
  pivot_wider(names_from = Phase, values_from = mean_vmac)


# Split-half reliability for VMAC ---- 
library(splithalf)

diff_reliability <- splithalf(data = e1,
                              outcome = "RT",
                              score = "difference",
                              conditionlist = c("Rewarded", "Extinction"),
                              halftype = "random",
                              permutations = 5000,
                              var.RT = "rt",
                              var.condition = "Phase",
                              var.participant = "ID",
                              var.compare = "Singleton",
                              compare1 = "High",
                              compare2 = "Low",
                              average = "mean")


# Read UPPS-P subscales data
qs_e2 <- read.csv("Output/Qs_e2.csv") 
corr_e2 <- vmac_e2 %>%
  merge(qs_e2, by = "ID")

# Correlations
e2_NU_acq <- cor.test(corr_e2$Rewarded, corr_e2$urg_neg, method = "pearson") # r = -0.04967269; p-value = 0.6063
e2_PU_acq <- cor.test(corr_e2$Rewarded, corr_e2$urg_pos, method = "pearson") # r = 0.04495489; p-value = 0.641

e2_NU_ext <- cor.test(corr_e2$Extinction, corr_e2$urg_neg, method = "pearson") # r = 0.09719144; p-value = 0.3102
e2_PU_ext <- cor.test(corr_e2$Extinction, corr_e2$urg_pos, method = "pearson") # r = 0.2383775; p-value = 0.01175



## FISHER R-TO-Z TRANSFORMATION: comparison between experiments ----
compare_correlations <- function(r1, n1, r2, n2) {
  fisher_r_to_z <- function(r) {  # Fisher r-to-z transformation
    0.5 * log((1 + r) / (1 - r))
  }
  
  z1 <- fisher_r_to_z(r1)
  z2 <- fisher_r_to_z(r2)
  
  se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))   # Combined standard error
  z_diff <- (z1 - z2) / se   # z-scores difference
  p_value <- 2 * (1 - pnorm(abs(z_diff)))   # Two-tailed p-value
  
  list(z_score = z_diff, p_value = p_value)
}


correlations <- data.frame(
  r1 = c(e1_NU_acq$estimate, e1_PU_acq$estimate, e1_NU_ext$estimate, e1_PU_ext$estimate),
  n1 = nrow(corr_e1),
  r2 = c(e2_NU_acq$estimate, e2_PU_acq$estimate, e2_NU_ext$estimate, e2_PU_ext$estimate),
  n2 = nrow(corr_e2)
)


results <- apply(correlations, 1, function(row) {
  compare_correlations(row["r1"], row["n1"], row["r2"], row["n2"])
})

results_df <- do.call(rbind, lapply(results, as.data.frame)) 
print(results_df)


