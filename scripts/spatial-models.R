#_______________________________________________________________________________
#                           Looming Eye Buoy Work
#                            MATT HORTON (2025)  
#                        Influence of the deterrents
#_______________________________________________________________________________


# Aim:
# To assess the influence of the deterrents on the number of birds caught. 
#    1. Look at the column Tom created as meters from closest deterrent. 
#    2. Change to cumulative mortality plots. Then see what the plot lools like. 



#_______________________________________________________________________________
# plot function
{# for viewing
  titlesize = 10
  textsize = 8
  library(ggplot2); theme_update(text = element_text(family = "Calibri", size = textsize), # all text defaults
                                 title = element_text(family = "Calibri", size = titlesize), # all title defaults
                                 plot.subtitle = element_text(face = "bold"), # all subtitle defaults
                                 axis.text = element_text(margin = margin(0,0,0.3,0.3, unit = "lines")), # axis labels including distance from margin 
                                 axis.title = element_text(margin = margin(0.0,0,0.6,0.6, unit = "lines")), # axis labels including distance from margin 
                                 axis.ticks = element_line(linewidth = 0.1),
                                 axis.ticks.length = unit(1, "points"),
                                 panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.25),
                                 panel.grid = element_blank(),
                                 panel.background = element_rect(fill = NA, colour = "black", linewidth = 0.25),
                                 strip.background = element_rect(fill = NA),
                                 strip.text = element_text(hjust = 0, face = "bold", margin = margin(0.1,0,0.1,0, "lines")),
                                 legend.key=element_blank(),
                                 legend.key.height = unit(0.2, 'cm'), #change legend key height
                                 legend.key.width = unit(0.2, 'cm'),
                                 legend.spacing = unit(5, 'points'),
                                 plot.margin = margin(0,0,1,1, unit = "lines"),
                                 legend.position = "right")
  
  update_geom_defaults("point", aes(stroke = 0.2, size = 1))
  update_geom_defaults("line", list(linewidth = 0.2, stroke = 0.2))
  update_geom_defaults("path", list(linewidth = 0.2, stroke = 0.2))
  update_geom_defaults("col", list(linewidth = 0.2, stroke = 0.2))
  update_geom_defaults("vline", list(linewidth = 0.2, stroke = 0.2))
  update_geom_defaults("hline", list(linewidth = 0.1, stroke = 0.1))
  update_geom_defaults("text", list(size = titlesize))
  
}


# Packages
library(tidyr)
library(dplyr)
library(readxl)
library(rstanarm)
library(rstantools)
library(bayesplot)
library(tidybayes)
library(MASS)
library(rstan)
library(SciViews)
library(lme4)
library(glmmTMB)
library(emmeans)
library(ggplot2)


# Load data: 
leb_dat <- read.csv("./data/bycatch.csv")

# LEB ----
alci_leb <- leb_dat[leb_dat$family %in% "Alcidae" &
                      leb_dat$net_description %in% "LEB", ]

# Plot cumulative 
alci_leb_sorted <- alci_leb %>%
                    filter(!is.na(min_m_treat)) %>%  # remove NAs
                    arrange(min_m_treat) %>%
                    mutate(cumulative_mortality = row_number())  # each row is one bird

alci_leb_sorted <- alci_leb %>%
                    filter(!is.na(min_m_treat)) %>%
                    arrange(min_m_treat) %>%
                    mutate(cumulative_mortality = row_number(),
                           cum_prop = cumulative_mortality / max(cumulative_mortality))

library(minpack.lm)  # for nlsLM, more robust than nls

# exponential saturation
exp_mod <- nlsLM(cum_prop ~ a * (1 - exp(-b * min_m_treat)),
                 data = alci_leb_sorted,
                 start = list(a = 1, b = 0.01))

# logistic model
log_mod <- nlsLM(cum_prop ~ L / (1 + exp(-k * (min_m_treat - x0))),
                 data = alci_leb_sorted,
                 start = list(L = 1, k = 0.01, x0 = median(alci_leb_sorted$min_m_treat)))

# weibull model
weib_mod <- nlsLM(cum_prop ~ 1 - exp(- (min_m_treat / lambda)^k),
                  data = alci_leb_sorted,
                  start = list(lambda = 100, k = 2))

AIC(exp_mod, log_mod, weib_mod)
#Weibull is definitely the best for this data

#  deviance explained:
# Response
y <- alci_leb_sorted$cum_prop

# Residual sum of squares function
rss <- function(model) sum(residuals(model)^2)

# RSS for each model
rss_exp <- rss(exp_mod)
rss_log <- rss(log_mod)
rss_weib <- rss(weib_mod)

# Null model RSS (mean only)
rss_null <- sum((y - mean(y))^2)

# Deviance explained
dev_exp <- 1 - rss_exp / rss_null
dev_log <- 1 - rss_log / rss_null
dev_weib <- 1 - rss_weib / rss_null

# Compile results
data.frame(
  model = c("Exponential","Logistic","Weibull"),
  deviance_explained = c(dev_exp, dev_log, dev_weib)
)



#_______________________________________________________________________________
# KITE ----
# look between the different deterrents 
alci_kite <- leb_dat[leb_dat$group %in% "Bird" &
                      leb_dat$net_description %in% "KITE", ]


# Plot cumulative plot:
alci_kite_sorted <- alci_kite %>%
                    filter(!is.na(min_m_treat)) %>%  # remove NAs
                    arrange(min_m_treat) %>%
                    mutate(cumulative_mortality = row_number())  # each row is one bird



alci_kite_sorted <- alci_kite %>%
                    filter(!is.na(min_m_treat)) %>%
                    arrange(min_m_treat) %>%
                    mutate(cumulative_mortality = row_number(),
                           cum_prop = cumulative_mortality / max(cumulative_mortality))

#________________________

# exponential saturation
exp_modK <- nlsLM(cum_prop ~ a * (1 - exp(-b * min_m_treat)),
                 data = alci_kite_sorted,
                 start = list(a = 1, b = 0.01))

# logistic model
log_modK <- nlsLM(cum_prop ~ L / (1 + exp(-k * (min_m_treat - x0))),
                 data = alci_kite_sorted,
                 start = list(L = 1, k = 0.01, x0 = median(alci_kite_sorted$min_m_treat)))

# weibull model
weib_modK <- nlsLM(cum_prop ~ 1 - exp(- (min_m_treat / lambda)^k),
                  data = alci_kite_sorted,
                  start = list(lambda = 100, k = 2))

AIC(exp_modK, log_modK, weib_modK)


#  deviance explained:
# Response
y <- alci_kite_sorted$cum_prop

# Residual sum of squares function
rss <- function(model) sum(residuals(model)^2)

# RSS for each model
rss_exp <- rss(exp_modK)
rss_log <- rss(log_modK)
rss_weib <- rss(weib_modK)

# Null model RSS (mean only)
rss_null <- sum((y - mean(y))^2)

# Deviance explained
dev_exp <- 1 - rss_exp / rss_null
dev_log <- 1 - rss_log / rss_null
dev_weib <- 1 - rss_weib / rss_null

# Compile results
data.frame(
  model = c("Exponential","Logistic","Weibull"),
  deviance_explained = c(dev_exp, dev_log, dev_weib)
)


# logistic model best for this model 

pred_df <- data.frame(min_m_treat = seq(
  min(alci_kite_sorted$min_m_treat),
  max(alci_kite_sorted$min_m_treat),
  length.out = 200)
)

# Add predictions
pred_df$exp <- predict(exp_modK, newdata = pred_df)
pred_df$log <- predict(log_modK, newdata = pred_df)
pred_df$weib <- predict(weib_modK, newdata = pred_df)

#_______________________________________________________________________________
# Model comparison summary table ----
# Combines all six nls fits above (Exponential, Logistic, Weibull x LEB, Kite)
# into a single results table for reporting:
#   Device, Model, Parameters (fitted values), df, AIC, Diff AIC, Deviance explained
# Diff AIC = AIC - (lowest/best AIC within that Device).
# Rows are ordered within each Device from AIC highest (worst fit) to lowest (best fit).

# Helper: format a model's fitted coefficients as "name = value, name = value"
format_params <- function(model) {
  co <- coef(model)
  paste(paste0(names(co), " = ", signif(co, 4)), collapse = ", ")
}

# Helper: deviance explained = 1 - RSS(model) / RSS(null, mean-only) for an nls model
dev_explained <- function(model, y) {
  rss_model <- sum(residuals(model)^2)
  rss_null  <- sum((y - mean(y))^2)
  1 - rss_model / rss_null
}

model_summary <- data.frame(
  Device = factor(c(rep("LEB", 3), rep("Kite", 3)), levels = c("LEB", "Kite")),
  Model = rep(c("Exponential", "Logistic", "Weibull"), 2),
  Parameters = c(
    format_params(exp_mod),  format_params(log_mod),  format_params(weib_mod),
    format_params(exp_modK), format_params(log_modK), format_params(weib_modK)
  ),
  df = c(
    attr(logLik(exp_mod), "df"),  attr(logLik(log_mod), "df"),  attr(logLik(weib_mod), "df"),
    attr(logLik(exp_modK), "df"), attr(logLik(log_modK), "df"), attr(logLik(weib_modK), "df")
  ),
  AIC = c(
    AIC(exp_mod),  AIC(log_mod),  AIC(weib_mod),
    AIC(exp_modK), AIC(log_modK), AIC(weib_modK)
  ),
  Deviance_explained = c(
    dev_explained(exp_mod,  alci_leb_sorted$cum_prop),
    dev_explained(log_mod,  alci_leb_sorted$cum_prop),
    dev_explained(weib_mod, alci_leb_sorted$cum_prop),
    dev_explained(exp_modK,  alci_kite_sorted$cum_prop),
    dev_explained(log_modK,  alci_kite_sorted$cum_prop),
    dev_explained(weib_modK, alci_kite_sorted$cum_prop)
  ),
  stringsAsFactors = FALSE
)

# Diff AIC relative to the largest (worst) AIC within each Device, i.e. how much
# lower each model's AIC is than the top (worst-fitting) row. The top row itself
# has no comparison and is shown as "-".
# Then order rows within Device from highest AIC to lowest.
model_summary <- model_summary %>%
  group_by(Device) %>%
  mutate(Diff_AIC = ifelse(AIC == max(AIC), NA_real_, max(AIC) - AIC)) %>%
  ungroup() %>%
  arrange(Device, desc(AIC)) %>%
  as.data.frame()

# Round for display
model_summary$AIC <- round(model_summary$AIC, 2)
model_summary$Diff_AIC <- round(model_summary$Diff_AIC, 2)
model_summary$Diff_AIC <- ifelse(is.na(model_summary$Diff_AIC), "-", as.character(model_summary$Diff_AIC))
model_summary$Deviance_explained <- round(model_summary$Deviance_explained, 3)

# Final column order
model_summary <- model_summary[, c("Device", "Model", "Parameters", "df", "AIC", "Diff_AIC", "Deviance_explained")]

print(model_summary)
write.csv(model_summary, "./tables/Table3.csv")

# Export if needed:
# write.csv(model_summary, "./outputs/model_comparison_summary.csv", row.names = FALSE)


