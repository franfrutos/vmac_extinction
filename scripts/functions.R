source("scripts/load_data.R") # Automatically load data

if (!require(pacman)) {
  install.packages(pacman)
  library(pacman)
}

p_load(dplyr, Hmisc, grid, ggthemes)

# Function to filter data:
filter_data <- function(data,
                        f.absent = T,
                        acc = T,
                        sd_filter = NULL,
                        fixed = F,
                        two_t = F,
                        phase = "Rewarded"
) {

  filt_trial <- NA
  if (two_t) {
    filt_trial <- c(seq(1, 577, 24), seq(2, 577, 24))
  }
  
  out <- data %>%
    filter(Phase %in% phase,
      !(trial_num %in% filt_trial),
      !ID %in% acc_ex,
      !is.na(rt))
  
  if (f.absent)
    out <- out[which(out$Singleton != "Absent"),]
  
  if (is.numeric(sd_filter) & ifelse(is.null(sd_filter), F, sd_filter %in% 1:3)) {
    out <- out %>%
      group_by(ID, Phase) %>%
      dplyr::mutate(high_rt = mean(rt, na.rm = T) + sd(rt, na.rm = T)*sd_filter,
             low_rt = mean(rt, na.rm = T) - sd(rt, na.rm = T)*sd_filter) %>%
      ungroup()  %>%
      filter(rt > low_rt, rt < high_rt) %>%
      select(-c(high_rt, low_rt))
  }
  
  if (fixed) {
    out <- out[which(out$rt > 150 & out$rt < 1800),]
    }
  
  if (acc) {
    out <- out[which(out$correct == 1),]
    }
  
  return(out)
}

create_epochs <- function(blocks, epoch = 2) {
  vapply(blocks, function(x, e = epoch) {
    ceiling(x / e)
  }, FUN.VALUE = numeric(1))
  
}

# Theme used in the plots
theme_Publication <-
  function(base_size = 12,
           base_family = "sans",
           text_size = 11) {
    (
      theme_foundation(base_size = base_size, base_family = base_family)
      + theme(
        plot.title = element_text(
          face = "bold",
          size = rel(1.2),
          hjust = 0.5
        ),
        text = element_text(size = text_size),
        panel.background = element_rect(colour = NA),
        plot.background = element_rect(colour = NA),
        panel.border = element_rect(colour = NA),
        axis.title = element_text(face = "bold", size = rel(1)),
        axis.title.y = element_text(angle = 90, vjust = 2),
        axis.title.x = element_text(vjust = -0.2),
        axis.text = element_text(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(),
        panel.grid.major = element_line(colour = "#f0f0f0"),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size = unit(0.2, "cm"),
        #legend.margin = margin(0, "cm"),
        #legend.title = element_text(face="italic"),
        plot.margin = unit(c(10, 5, 5, 5), "mm"),
        strip.background = element_rect(colour = NA, fill = NA),
        strip.text = element_text(face = "bold")
      )
    )
    
  }

# Functions to plot:

# This function gets predictions for both models presented in the main text to plot (Figure 2)
# d indicate raw data fed in a model for each phase. m indicate models for each phase
# Epoch is set to TRUE, which use epochs of 2 for raw data
# returns a list with 1) averaged raw data by Phase, Block (or Epoch) and Singleton
# and 2) model predictions for both models

get_predictions <- function(d1, d2, m1, m2, epoch = T, acc = F) {
  library(marginaleffects)
  library(dplyr)
  
  if (epoch) {
    d1$Block <- create_epochs(d1$Block_num)
    d2$Block <- create_epochs(d2$Block_num)
  }
  if (!acc) {
    raw_data <- rbind(
      d1 %>% group_by(ID, Phase, Block, Singleton) %>% dplyr::summarise(estimate = mean(rt)) %>%
        Rmisc::summarySEwithin(
          data = .,
          measurevar = "estimate",
          withinvars = c("Phase", "Block", "Singleton"),
          idvar = "ID"
        ) %>% mutate(Phase = "Rewarded"),
      d2 %>% group_by(ID, Phase, Block, Singleton) %>% dplyr::summarise(estimate = mean(rt)) %>%
        Rmisc::summarySEwithin(
          data = .,
          measurevar = "estimate",
          withinvars = c("Phase", "Block", "Singleton"),
          idvar = "ID"
        ) %>% mutate(Phase = "Extinction")
    )
    raw_data$Singleton <-
      factor(raw_data$Singleton, levels = c("High", "Low", "Absent"))
    raw_data$Block <- as.numeric(raw_data$Block) * 2 - .5
  } else{
    raw_data <- rbind(
      d1 %>% group_by(ID, Phase, Block, Singleton) %>% dplyr::summarise(estimate = mean(correct)) %>%
        Rmisc::summarySEwithin(
          data = .,
          measurevar = "estimate",
          withinvars = c("Phase", "Block", "Singleton"),
          idvar = "ID"
        ) %>% mutate(Phase = "Rewarded"),
      d2 %>% group_by(ID, Phase, Block, Singleton) %>% dplyr::summarise(estimate = mean(correct)) %>%
        Rmisc::summarySEwithin(
          data = .,
          measurevar = "estimate",
          withinvars = c("Phase", "Block", "Singleton"),
          idvar = "ID"
        ) %>% mutate(Phase = "Extinction")
    )
    raw_data$Singleton <-
      factor(raw_data$Singleton, levels = c("High", "Low", "Absent"))
    raw_data$Block <- as.numeric(raw_data$Block) * 2 - .5
  }
  
  model_preds <- bind_rows(
    predictions(
      m1,
      newdata = datagrid(
        Singleton = unique,
        Block_num = seq(1, 12, .01),
        ID = NA
      ),
      re.form = NA,
      transform = \(x) exp(x + (sigma(m1) ^ 2) / 2),
      #by = c("Singleton", "Block_num")
    ) %>% mutate(Phase = "Rewarded") %>% dplyr::rename("Block" = "Block_num"),
    predictions(
      m2,
      newdata = datagrid(
        Singleton = unique,
        Block_num = seq(13, 24, .01),
        ID = NA
      ),
      re.form = NA,
      transform = \(x) exp(x + (sigma(m2) ^ 2) / 2)
    ) %>% mutate(Phase = "Extinction") %>% dplyr::rename("Block" = "Block_num")
  )
  
  model_preds$Singleton <-
    factor(model_preds$Singleton, levels = c("High", "Low", "Absent"))
  
  return(list(raw = raw_data,
              mod = model_preds))
}

# Intermediate function to get VMAC and AC effects in get_comparisons()
get_raw_effect <- function(d, epoch = TRUE, acc = FALSE, both = FALSE) {
  if (epoch) {
    if (!"Block_num" %in% names(d)) stop("Column 'Block_num' not found in data")
    d$Block <- create_epochs(d$Block_num)
  }
  
  if (both) {
    vars <- c("ID", "Experiment", "Phase", "Block", "Singleton")
  } else {
    vars <- c("ID", "Phase", "Block", "Singleton")
  }
  
  print(vars)
  
  if (!acc) {
    effs <- d %>%
      group_by(across(all_of(vars))) %>%
      dplyr::summarise(RT = mean(rt, na.rm = TRUE), .groups = "drop") %>%
      spread(Singleton, RT) %>%
      mutate(VMAC = High - Low,
             AC = Low - Absent) %>%
      ungroup() %>%
      drop_na() 
    
  } else {
    effs <- d %>%
      group_by(across(all_of(vars))) %>%
      dplyr::summarise(acc = mean(correct, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = Singleton, values_from = acc) %>%
      mutate(VMAC = coalesce(High, 0) - coalesce(Low, 0),
             AC = coalesce(Low, 0) - coalesce(Absent, 0)) %>%
      drop_na()
  }
  
  if (!both) {
    return(
      rbind(
        summarySEwithin2(
          data = effs,
          measurevar = "VMAC",
          withinvars = c("Phase", "Block"),
          idvar = "ID"
        ) %>%
          dplyr::rename("estimate" = "VMAC") %>%
          mutate(Block = as.numeric(as.character(Block)), Effect = "high-low") %>%
          select(-VMAC_norm),
        summarySEwithin2(
          data = effs,
          measurevar = "AC",
          withinvars = c("Phase", "Block"),
          idvar = "ID"
        ) %>%
          dplyr::rename("estimate" = "AC") %>%
          mutate(Block = as.numeric(as.character(Block)), Effect = "low-absent") %>%
          select(-AC_norm)
      )
    )
  }
  return(
    rbind(
      summarySEwithin2(
        data = effs,
        betweenvars = "Experiment",
        measurevar = "VMAC",
        withinvars = c("Phase", "Block"),
        idvar = "ID"
      ) %>%
        dplyr::rename("estimate" = "VMAC") %>%
        mutate(Block = as.numeric(as.character(Block)), Effect = "high-low") %>%
        select(-VMAC_norm),
      summarySEwithin2(
        data = effs,
        betweenvars = "Experiment",
        measurevar = "AC",
        withinvars = c("Phase", "Block"),
        idvar = "ID"
      ) %>%
        dplyr::rename("estimate" = "AC") %>%
        mutate(Block = as.numeric(as.character(Block)), Effect = "low-absent") %>%
        select(-AC_norm)
    )
  )
}


# This function gets conditional effects for each contrasts in both models presented in the main text to plot (Figure 3)
# d indicate raw data fed in a model for each phase. m indicate models for each phase
# Epoch is set to TRUE, which use epochs of 2 for raw data
# returns a list with 1) averaged raw data by Phase, Block (or Epoch) and Singleton
# and 2) model predictions for both models

get_comparisons <- function(d1, d2, m1, m2, epoch = T, acc = F, both = F) {
  
  if (!both) {
    raw_data <- rbind(get_raw_effect(d1, epoch, acc = acc),
                      get_raw_effect(d2, epoch, acc = acc))
  } else {
    raw_data <- get_raw_effect(d1, both = T, acc = acc) %>% filter(Effect == "high-low")
  }
  raw_data$Effect <-
    factor(raw_data$Effect, levels = c("high-low", "low-absent"))
  #raw_data$Block[raw_data$Phase == "Extinction"] <- raw_data$Block[raw_data$Phase == "Extinction"] + 6
  raw_data$Block <- as.numeric(raw_data$Block) * 2 - .5
  
  
  if (!both) {
    if (!acc) {
      model_comps <- bind_rows(
        comparisons(
          m1,
          variables = list(Singleton = "revsequential"),
          newdata = datagrid(Block_num = seq(1, 12, .01),
                             ID = NA),
          re.form = NA,
          comparison = \(hi, lo) exp(hi + (sigma(m1) ^ 2) / 2) - exp(lo +
                                                                       (sigma(m1) ^ 2) / 2),
        ) %>% mutate(Effect = rep(
          c("high-low", "low-absent"), each = nrow(.) / 2
        ), Phase = "Rewarded") %>% dplyr::rename("Block" = "Block_num"),
        comparisons(
          m2,
          variables = list(Singleton = "revsequential"),
          newdata = datagrid(Block_num = seq(13, 24, .01),
                             ID = NA),
          re.form = NA,
          comparison = \(hi, lo) exp(hi + (sigma(m2) ^ 2) / 2) - exp(lo +
                                                                       (sigma(m2) ^ 2) / 2),
        ) %>% mutate(Effect = rep(
          c("high-low", "low-absent"), each = nrow(.) / 2
        ), Phase = "Extinction") %>% dplyr::rename("Block" = "Block_num")
      )
    } else {
      model_comps <- bind_rows(
        comparisons(
          m1,
          variables = list(Singleton = "revsequential"),
          newdata = datagrid(Block_num = seq(1, 12, .01),
                             ID = unique),
          re.form = NA,
          by = c("Singleton", "Block_num"),
        ) %>% mutate(Effect = rep(
          c("high-low", "low-absent"), each = nrow(.) / 2
        ), Phase = "Rewarded") %>% dplyr::rename("Block" = "Block_num"),
        comparisons(
          m2,
          variables = list(Singleton = "revsequential"),
          newdata = datagrid(Block_num = seq(13, 24, .01),
                             ID = unique),
          re.form = NA,
          by = c("Singleton", "Block_num"),
        ) %>% mutate(Effect = rep(
          c("high-low", "low-absent"), each = nrow(.) / 2
        ), Phase = "Extinction") %>% dplyr::rename("Block" = "Block_num")
      )
    }
  } else {
    if (!acc) {
      model_comps <- bind_rows(
        comparisons(
          m1,
          variables = "Singleton",
          newdata = datagrid(Block_num = seq(1, 12, .01),
                             ID = NA,
                             Exp = unique),
          re.form = NA,
          comparison = \(hi, lo) exp(hi + (sigma(m1) ^ 2) / 2) - exp(lo +
                                                                       (sigma(m1) ^ 2) / 2),
          by = c("Exp", "Block_num")
        ) %>% 
           mutate(Effect = rep(
           c("high-low", "low-absent"), each = nrow(.) / 2
         ), Phase = "Rewarded") %>% dplyr::rename("Block" = "Block_num"),
        comparisons(
          m2,
          variables = "Singleton",
          newdata = datagrid(Block_num = seq(13, 24, .01),
                             ID = NA,
                             Exp = unique),
          re.form = NA,
          comparison = \(hi, lo) exp(hi + (sigma(m2) ^ 2) / 2) - exp(lo +
                                                                       (sigma(m2) ^ 2) / 2),
          by = c("Exp", "Block_num")
        ) %>%
           mutate(Effect = rep(
           c("high-low", "low-absent"), each = nrow(.) / 2
         ), Phase = "Extinction") %>% dplyr::rename("Block" = "Block_num")
      )
    } else {
      model_comps <- bind_rows(
        comparisons(
          m1,
          variables = "Singleton",
          newdata = datagrid(Block_num = seq(1, 12, .01),
                             ID = unique,
                             Exp = unique),
          re.form = NA,
          by = c("Exp", "Block_num"),
        ) %>% mutate(Effect = rep(
          c("high-low", "low-absent"), each = nrow(.) / 2
        ), Phase = "Rewarded") %>% dplyr::rename("Block" = "Block_num"),
        comparisons(
          m2,
          variables = "Singleton",
          newdata = datagrid(Block_num = seq(13, 24, .01),
                             ID = unique,
                             Exp = unique),
          re.form = NA,
          by = c("Exp", "Block_num"),
        ) %>% mutate(Effect = rep(
          c("high-low", "low-absent"), each = nrow(.) / 2
        ), Phase = "Extinction") %>% dplyr::rename("Block" = "Block_num")
      )
    }
  }
  
  model_comps$Effect <-
    factor(model_comps$Effect, levels = c("high-low", "low-absent"))
  
  return(list(raw = raw_data,
              mod = model_comps))
}

# Functions to get average and SE for raw data in mixed designs. The source code is in: :
summarySEwithin2 <- function (data = NULL,
                              measurevar,
                              betweenvars = NULL,
                              withinvars = NULL,
                              idvar = NULL,
                              na.rm = FALSE,
                              conf.interval = 0.95,
                              .drop = TRUE)
{
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop = FALSE], FUN = is.factor, FUN.VALUE = logical(1))
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message(
      "Automatically converting the following non-factors to factors: ",
      paste(nonfactorvars, collapse = ", ")
    )
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  datac <- summarySE(
    data,
    measurevar,
    groupvars = c(betweenvars, withinvars),
    na.rm = na.rm,
    conf.interval = conf.interval,
    .drop = .drop
  )
  
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop = .drop)
  measurevar_n <- paste(measurevar, "_norm", sep = "")
  ndatac <- summarySE(
    ndata,
    measurevar_n,
    groupvars = c(betweenvars, withinvars),
    na.rm = na.rm,
    conf.interval = conf.interval,
    .drop = .drop
  )
  nWithinGroups <- prod(vapply(ndatac[, withinvars, drop = FALSE], FUN = nlevels, FUN.VALUE = numeric(1)))
  correctionFactor <- sqrt(nWithinGroups / (nWithinGroups - 1))
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  merge(datac, ndatac)
}
normDataWithin <- function(data = NULL,
                           idvar,
                           measurevar,
                           betweenvars = NULL,
                           na.rm = FALSE,
                           .drop = TRUE) {
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- plyr::ddply(
    data,
    c(idvar, betweenvars),
    .drop = .drop,
    .fun = function(xx, col, na.rm) {
      c(subjMean = mean(xx[, col], na.rm = na.rm))
    },
    measurevar,
    na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep = "")
  data[, measureNormedVar] <- data[, measurevar] - data[, "subjMean"] +
    mean(data[, measurevar], na.rm = na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}
