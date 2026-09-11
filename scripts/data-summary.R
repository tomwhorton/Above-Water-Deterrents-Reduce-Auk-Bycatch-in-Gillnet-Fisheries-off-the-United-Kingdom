# Packages
library(tidyverse)
library(patchwork)
library(extrafont)

# clean environment prior to starting
cat("\014") 
rm(list=ls())

# functions
fun.list <- list(mean = mean, sd = sd, min = min, max = max, sum = sum)

# Data
dat <- read_csv("./data/sets-bycatch.csv") %>% 
           mutate(net_description = case_when(net_description == "CONTROL" ~ "Ctrl.",
                                             net_description == "KITE" ~ "Kite",
                                             net_description == "LEB" ~ "LEB"))      
# date range
min(dat$set_date)
max(dat$haul_date)

# n vessels 
length(unique(dat$vessel_id))

# n nets by treatment 
nrow(dat)
table(dat$net_description)

# effort
dat %>% summarise_at(vars(net_length_m), fun.list) # length in m
dat %>% summarise_at(vars(soak_time_hrs), fun.list) # soak time hrs
dat %>% summarise_at(vars(net_height_m), fun.list) # average net height

# overnight sets 
dat$overnight = ifelse(dat$haul_date - dat$set_date == 0, FALSE, TRUE)
dat %>% group_by(overnight) %>% summarise(n = n(), mean_soak = mean(soak_time_hrs), sd_soak = sd(soak_time_hrs), min_soak = min(soak_time_hrs), max_soak = max(soak_time_hrs)) %>% ungroup() %>% mutate(prop = n/sum(n)*100)


# bycatch
# data 
byc <- read_csv("./data/bycatch.csv") %>% 
           mutate(net_description = case_when(net_description == "CONTROL" ~ "Ctrl.",
                                              net_description == "KITE" ~ "Kite",
                                              net_description == "LEB" ~ "LEB"),
                  common_name = case_when(common_name == "Shag/Cormorant Species" ~ "Phalacrocoracid spp.",
                                          common_name == "Unidentified Toothed Whales Species" ~ "Odontocete spp.",
                                          common_name == "Unidentified Seal Species" ~ "Pinniped spp.",
                                          TRUE ~ common_name)) 

#dat 
nrow(byc) # total bycatch
byc %>% group_by(group) %>% tally() %>% mutate(prop = (n/sum(n))*100) %>% arrange(desc(n))
byc %>% group_by(common_name) %>% tally() %>% mutate(prop = (n/sum(n))*100) %>% arrange(desc(n))

# now many dead etc
byc %>% group_by(condition) %>% count() %>% ungroup() %>% mutate(prop = n/sum(n)*100)
byc %>% filter(condition == "A") %>% group_by(common_name) %>% count()

