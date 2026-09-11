#_______________________________________________________________________________
#                           Looming Eye Buoy Work
#                            MATT HORTON (2025)  
#                     Influence of the deterrents--GLMM
#_______________________________________________________________________________

# Model the cumulative mortality by depth and other covariates so that the 
# influence decay can feed into management recommendations for the LEB work. 

#  Plotting 
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
                                 plot.margin = margin(1,1,1,1, unit = "lines"),
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
library(tidyverse)
library(MASS)
library(SciViews)
library(lme4)
library(emmeans)
library(ggplot2)
library(minpack.lm)
library(boot)
library(patchwork)
library(ggpp)
library(ggimage)
library(extrafont)

# clean environment prior to starting
cat("\014") 
rm(list=ls())

# Data
dat <- read_csv("./data/bycatch.csv")

# subset 
alci <- dat[dat$family %in% "Alcidae", ]

# ______________________________________________________________________________
# Data transformation ----
alci_leb <- alci[alci$net_description %in% "LEB", ]
alci_kite <- alci[alci$net_description %in% "KITE", ]

alci_leb_sorted <- alci_leb %>%
                   filter(!is.na(min_m_treat)) %>%
                   arrange(min_m_treat) %>%
                   mutate(cumulative_mortality = row_number(),   # each row is one alci
                          cum_prop = cumulative_mortality / max(cumulative_mortality))  # proportion

alci_kite_sorted <- alci_kite %>%
                    filter(!is.na(min_m_treat)) %>%  # remove NAs
                    arrange(min_m_treat) %>%
                    mutate(cumulative_mortality = row_number(), # each row is one alci
                           cum_prop = cumulative_mortality / max(cumulative_mortality))  # proportion

# ______________________________________________________________________________
# Functions ----
# --- Weibull function ---
boot_weibull <- function(data, indices) {
  d <- data[indices, ]
  fit <- nls(cum_prop ~ 1 - exp(-(min_m_treat / lambda)^k),
             data = d,
             start = list(lambda = 50, k = 2))
  coef(fit)
}

# --- Logistic bootstrap function ---
boot_logistic <- function(data, indices) {
  d <- data[indices, ]
  fit <- try(nlsLM(cum_prop ~ L / (1 + exp(-k * (min_m_treat - x0))),
                   data = d,
                   start = list(L = 1,
                                k = 0.01,
                                x0 = median(d$min_m_treat))),
             silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(c(NA, NA, NA))  # skip failed fits
  } else {
    return(coef(fit))
  }
}



# ______________________________________________________________________________

#  AUKS ----
#  __leb ----
# --- Fit Weibull NLS ---

alci_leb <- nls(cum_prop ~ 1 - exp(- (min_m_treat / lambda)^k),
                data = alci_leb_sorted,
                start = list(lambda = 50, k = 2)) # adjust starting values as needed

# --- Run bootstrap ---
set.seed(123)
alci_leb_boot <- boot(alci_leb_sorted, boot_weibull, R = 1000)
# --- Generate predicted curves from bootstrap ---
dist_seq <- seq(min(alci_leb_sorted$min_m_treat),
                max(alci_leb_sorted$min_m_treat), length.out = 100)

boot_preds <- apply(alci_leb_boot$t, 1, function(p) {
  1 - exp(-(dist_seq / p[1])^p[2])
})
# --- Compute mean and percentile CI ---
ci_lower <- apply(boot_preds, 1, quantile, probs = 0.025)
ci_upper <- apply(boot_preds, 1, quantile, probs = 0.975)
pred_mean <- apply(boot_preds, 1, mean)



# ______________________________________________________________________________
# __kite ----
# --- Fit logistic NLS ---
alci_kite_fit <- nlsLM(cum_prop ~ L / (1 + exp(-k * (min_m_treat - x0))),
                       data = alci_kite_sorted,
                       start = list(L = 1, 
                                    k = 0.01, 
                                    x0 = median(alci_kite_sorted$min_m_treat)))

# --- Run bootstrap ---
set.seed(123)
logistic_boot <- boot(alci_kite_sorted, boot_logistic, R = 1000)
logistic_params <- logistic_boot$t[complete.cases(logistic_boot$t), ]

# --- Generate predicted curves from bootstrap ---
dist_seq_log <- seq(min(alci_kite_sorted$min_m_treat),
                    max(alci_kite_sorted$min_m_treat),
                    length.out = 100)

boot_preds_log <- apply(logistic_params, 1, function(p) {
  L <- p[1]; k <- p[2]; x0 <- p[3]
  L / (1 + exp(-k * (dist_seq_log - x0)))
})

# --- Compute mean and percentile CI ---
ci_lower_log <- apply(boot_preds_log, 1, quantile, probs = 0.025, na.rm = TRUE)
ci_upper_log <- apply(boot_preds_log, 1, quantile, probs = 0.975, na.rm = TRUE)
pred_mean_log <- apply(boot_preds_log, 1, mean, na.rm = TRUE)

#-------------------------------

# join all datasets for faceting 
alci <- bind_rows(alci_leb_sorted %>% dplyr::select(net_description, family, min_m_treat, cumulative_mortality, cum_prop),
                  alci_kite_sorted %>% dplyr::select(net_description, family, min_m_treat, cumulative_mortality, cum_prop)) #

alci_mods <- bind_rows(data.frame(net_description = "LEB", x = dist_seq, mean = pred_mean, ymin = ci_lower, ymax = ci_upper),
                       data.frame(net_description = "KITE", x = dist_seq_log, mean = pred_mean_log, ymin = ci_lower_log, ymax = ci_upper_log))

# scale the profiles
alci_mods <- alci_mods %>% left_join(.,alci %>% group_by(net_description) %>% summarise(cumulative_mortality = max(cumulative_mortality))) %>% 
                            group_by(net_description) %>% mutate(
                              across(
                                c(mean, ymin, ymax), 
                                ~ .x * max(cumulative_mortality), 
                                .names = "scaled_{.col}"
                              )
                            )


# final line plots
idf.sp <- data.frame(ims = c("./images/auk.png","./images/seal.png"),
                     name = c("alci","seal"),
                     x = 20, 
                     y = 24)

lines <- ggplot() +
           geom_point(data = alci, aes(x = min_m_treat, y = cumulative_mortality, color = net_description), alpha = 0) +
           geom_rect(aes(xmin = -Inf, xmax = 50, ymin = -Inf, ymax = Inf), alpha = 0.1) +
           geom_point(data = alci, aes(x = min_m_treat, y = cumulative_mortality, color = net_description)) +
           geom_line(data = alci_mods, aes(x = x, y = scaled_mean, color = net_description), size = 1) +
           scale_colour_manual(values = c("#FF6A1A","#00528C"), labels = c("Kite", "LEB"), name = "Treatment", guide = "none") +
           scale_fill_manual(values = c("#FF6A1A","#00528C"), labels = c("Kite", "LEB"), name = "Treatment", guide = "none") +
           geom_ribbon(data = alci_mods, aes(x = x, ymin = scaled_ymin, ymax = scaled_ymax, fill = net_description),
                       alpha = 0.2) +
           scale_y_continuous(expand = expansion(mult = 0.01), limits = c(0,28), breaks = seq(0,28,3)) +
           scale_x_continuous(expand = expansion(mult = 0.01)) +
           #ggpp::geom_text_npc(aes(label = "a)"), size = 3, npcx = "right", npcy = "top") +
           geom_image(data = idf.sp %>% filter(name == "alci"), aes(image = ims, x = x, y = y), size = 0.2) +
           labs(x = "Distance from AWD (m)", y = "Cumulative Bycatch") +
           theme(plot.margin = margin(0,0,0,0, unit = "cm"))


#lines

# final polar plots
# calculate cumulative frequency 

# distance categories 
sq <- seq(0,120,10)
dist_cats <- dat %>% filter(net_description != "CONTROL", family == "Alcidae") %>% 
             dplyr::select(net_description, min_m_buoys, family) %>% 
             mutate(dist_cat = cut(min_m_buoys,sq, include.lowest = TRUE, labels = seq(0,115,10))) %>%
             group_by(net_description, dist_cat) %>% filter(min_m_buoys <= 50) %>% 
             tally() %>% 
             group_by(net_description) %>% 
             mutate(cum_bycatch = cumsum(n),
                    prop_bycatch = n/sum(n))

# expand the dataframe 
labs <- expand_grid(net_description = c("KITE","LEB"), dist_cat = sq) %>% filter(dist_cat <= 40) %>% 
        mutate(dist_cat_lab = paste0(dist_cat, "-",dist_cat+5),
               dist_cat = factor(dist_cat, levels = sq)) 
  
# add the values 
dist_cats <- left_join(labs, dist_cats) %>% 
             mutate_at(vars(n,prop_bycatch,cum_bycatch), function(x){if_else(is.na(x),0,x)}) %>% 
             group_by(net_description) %>% #arrange(dist_cat) %>% 
             mutate(cum_bycatch = cumsum(n),
                    cum_prop = cumsum(prop_bycatch),
                    dist_cat = as.numeric(as.character(dist_cat)))  


# plot 
# imposted distance labs
m_labs <- data.frame(labs = c("0 m","50 m"),
                     y = c(-5,45),
                     x = 0.5)

# distance at which 10% (or X) bycatch is exceeded
dist_at_X <- dist_cats %>%
             group_by(net_description) %>%
             summarise(dist_at_X = approx(cum_bycatch, dist_cat, xout = c(5,10))$y)

bullseye_leb <- ggplot(dist_cats %>% filter(net_description == "LEB"), aes(y = dist_cat, x = "x", fill = cum_bycatch)) +  
                 geom_tile() +
                 scale_y_continuous(expand = expansion(mult = 0)) +
                 scale_x_discrete(expand = expansion(mult = 0)) +
                 scale_fill_gradientn(colors = c("#EEF6F9","#00528C"), 
                                      breaks = c(0, 5, 10, 15, 20),
                                      #guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                                      guide = "none") +
                 coord_polar() +
                 geom_text(inherit.aes = FALSE, data = m_labs, aes(label = labs, y = y, x = x), size = 2, col = "grey20", hjust = -0.2, vjust = 1) +
                 geom_segment(aes(x = 0.5, y = -5, yend = 45), col = "grey20", linewidth = 0.1) +
                 geom_hline(yintercept = dist_at_X %>% filter(net_description == "LEB") %>% mutate(dist_at_X = dist_at_X+5) %>% pull(dist_at_X), linetype = "dashed", linewidth = 0.2) +
                 theme_minimal() +
                 theme(axis.title = element_blank(),
                       axis.text = element_blank(),
                       panel.grid = element_blank(),
                       #plot.background = element_rect(fill = "yellow"),
                       plot.margin = margin(0,0,0,0))

bullseye_kite <- ggplot(dist_cats %>% filter(net_description == "KITE"), aes(y = dist_cat, x = "x", fill = cum_bycatch)) +  
                 geom_tile() +
                 scale_y_continuous(expand = expansion(mult = 0)) +
                 scale_x_discrete(expand = expansion(mult = 0)) +
                 scale_fill_gradientn(colors = c("#FFEDE1","#FF6A1A"), 
                                      breaks = c(0, 5, 10, 15, 20),   # 5 ind. increments
                                      name = NULL,
                                      guide = "none") +
                 coord_polar() +
                 geom_hline(yintercept = dist_at_X %>% filter(net_description == "KITE") %>% mutate(dist_at_X = dist_at_X+5) %>% pull(dist_at_X), linetype = "dashed", linewidth = 0.2) +
                 geom_text(inherit.aes = FALSE, data = dist_at_X %>% filter(net_description == "KITE") %>% mutate(dist_at_X = dist_at_X+5), 
                            aes(label = c("5 ind.","10 ind."), y = dist_at_X), x = 0.4, angle = 45, size = 1.5, col = "black", vjust = -0.75) +
                 theme_minimal() +
                 theme(axis.title = element_blank(),
                       axis.text = element_blank(),
                       panel.grid = element_blank(),
                       plot.margin = margin(0,0,0,0))


# final line plots
idf.tr <- data.frame(ims = c("./images/kite.png","./images/leb.png"),
                     net_description = c("KITE","LEB"),
                     x = 1, 
                     y = 15)

leb <- png::readPNG("./images/leb.png", native = TRUE)
kite <- png::readPNG("./images/kite.png", native = TRUE)

lines+(bullseye_leb/bullseye_kite) + plot_layout(widths = c(0.65,0.35), tag_level = "new") + 
                                     inset_element(leb, left = 0.9, bottom = 0.9, right = 1, top = 1, align_to = 'full', ignore_tag = TRUE) +
                                     inset_element(kite, left = 0.9, bottom = 0.47, right = 1, top = 0.55, align_to = 'full', ignore_tag = TRUE) +
                                     plot_annotation(tag_levels = list('a'), tag_suffix = ")") & theme(plot.tag = element_text(size = 8))

ggsave("./plots/Fig3-alci-spatial-efficacy.png", width = 13, height = 9, dpi = 600, unit = "cm") 

