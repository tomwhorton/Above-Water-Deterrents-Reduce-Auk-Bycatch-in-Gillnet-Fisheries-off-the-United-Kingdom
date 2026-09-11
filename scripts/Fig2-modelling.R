#_______________________________________________________________________________
#                           Looming Eye Buoy Work
#                            MATT HORTON (2025)  
#                        Zero-inflated hurdle model
#_______________________________________________________________________________


# Packages
library(tidyr)
library(tidyverse)
library(readxl)
library(tidybayes)
library(MASS)
library(rstan)
library(SciViews)
library(lme4)
library(glmmTMB)
library(emmeans)
library(ggplot2)
library(purrr)
library(DHARMa)
library(patchwork)
library(ggimage)
library(extrafont)

{# for viewing
  titlesize = 10
  textsize = 9
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
                                 plot.margin = margin(0.2,0.1,0.2,0.1, unit = "lines"),
                                 legend.position = "right")
  
  update_geom_defaults("point", aes(stroke = 0.2, size = 1))
  update_geom_defaults("line", list(linewidth = 0.2, stroke = 0.2))
  update_geom_defaults("path", list(linewidth = 0.2, stroke = 0.2))
  update_geom_defaults("col", list(linewidth = 0.2, stroke = 0.2))
  update_geom_defaults("vline", list(linewidth = 0.2, stroke = 0.2))
  update_geom_defaults("hline", list(linewidth = 0.1, stroke = 0.1))
  update_geom_defaults("text", list(size = titlesize))
  
}

# clean environment prior to starting
cat("\014") 
rm(list=ls())

# Data
leb_dat <- read.csv("./data/sets-bycatch.csv") %>% 
           mutate(net_description = case_when(net_description == "CONTROL" ~ "Ctrl.",
                                              net_description == "KITE" ~ "Kite",
                                              net_description == "LEB" ~ "LEB"))      

#_______________________________________________________________________________
# Initial Modeling----
# __ALCI ----
# no covariates in the binomial regression
z1 <- glmmTMB(n50_alci ~ net_description + offset(log(eff)) + (1 | haul_trip_id),
                    ziformula = ~1,
                    family = truncated_nbinom2,
                    data = leb_dat)

# the same covariates between the conditional and binomial model
# although it doesn't contain offset, as you modeling binary data not cpue
z2 <- glmmTMB(n50_alci ~ net_description + offset(log(eff)) + (1 | haul_trip_id),
                    ziformula = ~.,
                    family = truncated_nbinom2,
                    data = leb_dat)

# just net_description as the binomial covariate
z3 <- glmmTMB(n50_alci ~ net_description + offset(log(eff)) + (1 | haul_trip_id),
                    ziformula = ~net_description,
                    family = truncated_nbinom2,
                    data = leb_dat)
AIC(z1,z2,z3)

# AIC lowest for z2<z3<z1 - but lets check model assumptions first. 

# List of models
models <- list(z1 = z1, z2 = z2, z3 = z3)

# Simulate residuals for each model
residuals_list <- map(models, simulateResiduals)

# Plot residuals
#walk2(residuals_list, names(residuals_list), ~ {
#  cat("\n--- Residual plots for model:", .y, "---\n")
#  windows()
#  plot(.x, main = .y)
#  })

# Dispersion and zero-inflation tests
#walk2(residuals_list, names(residuals_list), ~ {
#  cat("\n--- Tests for model:", .y, "---\n")
#  print(testDispersion(.x))
#  print(testZeroInflation(.x))
#})


# It looks like model z3 is the only model to validate its assumptions
# No zero inflation 
# Overdispersion correctly modelled 
# Residuals look good. 
# Despite the AIC being the lowest for z2, z3 is the best, as z2 doesn't validate
# model assumptions and would lead to bias. 

#_______________________________________________________________________________

#============================================================
# model outputs
#=============================================================
# 1 emmeans modelled bycatch rates 
emmeans_alci <- as.data.frame(emmeans(z3, ~ net_description, type = "response"))

# 2 - Pairwise comparisons on the catch for the binomial model
pairwise_zi_l <- pairs(emmeans(z3, ~ net_description, type = "link", component = "zi"))
pairwise_zi <- pairs(emmeans(z3, ~ net_description, type = "response", component = "zi"))

# Get confidence intervals
pairs_df_zi <- as.data.frame(confint(pairwise_zi)) 

# 3- Pairwise comparisons on the log(CPUE) scale for the conditional model
pairwise_contrasts <- pairs(emmeans(z3, ~ net_description, type = "response"))

# Get confidence intervals
pairs_df <- as.data.frame(confint(pairwise_contrasts))


#============================================================
# writes estimates to tables
#=============================================================
odds.mod <- as.data.frame(pairs(emmeans(z3, ~ net_description, type = "response", component = "zi")))  %>% 
              dplyr::select(!null) %>% 
              mutate(reduction = (1-odds.ratio) *100,
                     reduction = if_else(p.value <= 0.3, round(reduction,0), NA),
                     odds.ratio = log(odds.ratio)) %>%   
              mutate_if(is.numeric, round, 3) %>% 
              dplyr::select(!df)

write_csv(odds.mod, file = "./tables/Table1.csv")

rate.mod <- as.data.frame(pairs(emmeans(z3, ~ net_description, type = "response")))  %>% 
              dplyr::select(!null) %>% 
              mutate(reduction = (1-1/ratio) *100,
                     reduction = if_else(p.value <= 0.3, round(reduction,0), NA)) %>%   
              mutate_if(is.numeric, round, 3) %>% 
              dplyr::select(!df)

write_csv(rate.mod, file = "./tables/Table2.csv")
#--------------------------------------------------------------


#============================================================
# plotting 
#=============================================================
pal <- c("#FF6A1A","#00528C","black")

idf.sp <- data.frame(ims = c("./images/auk.png","./images/seal.png"),
                     name = c("alci","seal"),
                     x = 2.75, 
                     y = 1.8)

title_b <- expression("Binomial – " * logit(pi[i,j]) == alpha + beta * x[j] + t[i])
title_c <- expression("BPUE – " * log(mu[i,j]) == alpha * "'" + beta * "'" * x[j] + t[i] + log(E[i,j]))


# Plot the differences and confidence intervals
alci_catch <- ggplot(emmeans_alci, aes(x = net_description, y = response)) +
                     geom_point(size = 3, aes(col = net_description)) +
                     geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL,col = net_description), width=0, alpha = 0.5) +
                     scale_colour_manual(values = c("black","#FF6A1A","#00528C"), guide = "none") + 
                     scale_x_discrete(limits = c("Ctrl.","LEB","Kite")) +
                     labs(y = "Auk Bycatch\nper 100 m\u00B2 d\u207B\u00B9",#expression(atop("Auk Bycatch\nper ", 100, " ", m^2, " d"^-1)),
                          x = "Net Type") +
                     ggpp::geom_text_npc(inherit.aes = FALSE, label = "a)", size = 3, npcx = "right", npcy = "top") +
                     geom_image(data = idf.sp %>% filter(name == "alci"), aes(image = ims, x = x, y = y), size = 0.25)

# add significance
stars <- data.frame(contrast = c("Ctrl. / LEB","Ctrl. / Kite"), odds.ratio = 3.8, label = c("*","**"))

# Plot the differences and confidence intervals
alci_prob <- ggplot(pairs_df_zi, aes(x = contrast, y = odds.ratio)) +
                            geom_point(size = 3) +
                            geom_hline(yintercept = 1, linetype = "dashed") +
                            geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width=0, alpha = 0.5) +
                            geom_text(data = stars, aes(label = label), size = 4) +
                            scale_x_discrete(labels = function(x){str_replace_all(x,c("-" = "/"))}) +
                            labs(y = "Odds ratio of \nzero-bycatch probability",
                                 x = "Net Type Comparison" )+
                            ggpp::geom_text_npc(inherit.aes = FALSE, label = "b)", size = 3, npcx = "right", npcy = "top")

alci_prob
# add significance
stars <- data.frame(contrast = "Ctrl. / LEB", ratio = 10, label = "**")

# Plot the differences and confidence intervals
alci_bpue <-ggplot(pairs_df, aes(x = contrast, y = ratio)) +
                    geom_point(size = 3) +
                    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width=0, alpha = 0.5) +
                    geom_text(data = stars, aes(label = label), size = 4) +
                    geom_hline(yintercept = 1, linetype = "dashed") +
                    labs(y = "\nBycatch Rate Ratio",
                         x = "Net Type Comparison") +
                    ggpp::geom_text_npc(inherit.aes = FALSE, label = "c)", size = 3, npcx = "right", npcy = "top")
  


alci_catch /alci_prob / alci_bpue
ggsave(filename = "./plots/Fig2-alci-mods.png", width = 6, height = 15, dpi = 600, unit = "cm")
#----------------------------------------------------------------------------------------------------------




