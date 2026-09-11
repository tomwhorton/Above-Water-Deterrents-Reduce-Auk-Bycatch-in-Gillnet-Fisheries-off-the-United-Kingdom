library(tidyverse)
library(patchwork)
library(ggimage)
library(extrafont)

{# for viewing
  titlesize = 10
  textsize = 8
  library(ggplot2); theme_update(text = element_text(family = "Calibri", size = textsize), # all text defaults
                                 title = element_text(family = "Calibri", size = titlesize), # all title defaults
                                 plot.subtitle = element_text(face = "bold"), # all subtitle defaults
                                 axis.text = element_text(margin = margin(0,0,0.4,0.4, unit = "lines")), # axis labels including distance from margin 
                                 axis.title = element_text(margin = margin(0.0,0,0.8,0.8, unit = "lines")), # axis labels including distance from margin 
                                 axis.ticks = element_line(linewidth = 0.1),
                                 axis.ticks.length = unit(1, "points"),
                                 panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.25),
                                 panel.grid = element_blank(),
                                 panel.background = element_rect(fill = NA, colour = "black", linewidth = 0.25),
                                 strip.background = element_rect(fill = NA),
                                 strip.text = element_text(hjust = 0, face = "bold", margin = margin(0.1,0,0.1,0, "lines")),
                                 legend.key=element_blank(),
                                 legend.key.height = unit(0.3, 'cm'), #change legend key height
                                 legend.key.width = unit(0.3, 'cm'),
                                 legend.spacing = unit(5, 'points'),
                                 plot.margin = margin(1,1,1,1, unit = "lines"),
                                 legend.position = "top")
  
  update_geom_defaults("point", aes(stroke = 0.2, size = 1))
  update_geom_defaults("line", list(linewidth = 0.2, stroke = 0.2))
  update_geom_defaults("path", list(linewidth = 0.2, stroke = 0.2))
  update_geom_defaults("col", list(linewidth = 0.2, stroke = 0.2))
  update_geom_defaults("vline", list(linewidth = 0.2, stroke = 0.2))
  update_geom_defaults("hline", list(linewidth = 0.1, stroke = 0.1))
  update_geom_defaults("text", list(size = titlesize))
  
}


cat("\014") 
rm(list=ls())

# se
se <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x))

# target species lookup 
spp <- read_csv("./data/species_lookup.csv") %>% 
       dplyr::select("Species Code", "Common Name") %>% 
       rename("name" = "Common Name", "target_species" = "Species Code") 

# read in the bycatch data 
dat <- read_csv("./data/sets-bycatch.csv") %>% 
       mutate(month = lubridate::month(haul_date, label = TRUE, abbr = TRUE),
              target_species = case_when(target_species == "SLO" ~ "SOL",
                                         TRUE ~ target_species),
              mesh_size = case_when(mesh_size_mm >=60 & mesh_size_mm <71 ~ "60-70",
                                    mesh_size_mm >=100 & mesh_size_mm <110 ~ "100-110",
                                    mesh_size_mm >=110 & mesh_size_mm <120 ~ "110-120",
                                    mesh_size_mm >=120 & mesh_size_mm <130 ~ "120-130",
                                    mesh_size_mm >=130 & mesh_size_mm <151 ~ "130-150",
                                    mesh_size_mm >=300 & mesh_size_mm <331 ~ "300-330"),
              mesh_size = factor(mesh_size, c("60-70", "100-110", "110-120", "120-130", "130-150", "300-330"))) %>% 
              left_join(.,spp)
   

# mesh sizes 
length(unique(dat$mesh_size_mm)) # number unique
range(dat$mesh_size_mm) # range


# PLOT AS BPUE 
#===============================================================================================================================================================================================================================
# summary of n alcidae caught per month per target species 
ts <- dat %>% group_by(month, name) %>% summarise(alci = sum(n_alci),
                                                  alci_bpue = mean(alci.blue.n50,na.rm = TRUE),
                                                  alci_bpue_se = se(alci.blue.n50)) %>% ungroup() %>% mutate(prop = alci/sum(alci)*100) 
ts.ts <- dat %>% group_by(name) %>% summarise(alci = sum(n_alci),
                                              alci_bpue = mean(alci.blue.n50,na.rm = TRUE),
                                              alci_bpue_se = se(alci.blue.n50)) %>% mutate(prop = alci/sum(alci)*100)
ts.mon <- dat %>% group_by(month) %>% summarise(alci = sum(n_alci),
                                                alci_bpue = mean(alci.blue.n50,na.rm = TRUE),
                                                alci_bpue_se = se(alci.blue.n50)) %>% mutate(prop = alci/sum(alci)*100)

# summary of n alcidae caught per month per mesh size 
ms <- dat %>% group_by(month, mesh_size) %>% summarise(alci = sum(n_alci),
                                                       alci_bpue = mean(alci.blue.n50,na.rm = TRUE),
                                                  alci_bpue_se = se(alci.blue.n50)) %>% ungroup() %>% mutate(prop = alci/sum(alci)*100)
ms.ms <- dat %>% group_by(mesh_size) %>% summarise(alci = sum(n_alci),
                                                   alci_bpue = mean(alci.blue.n50,na.rm = TRUE),
                                                   alci_bpue_se = se(alci.blue.n50)) %>% mutate(prop = alci/sum(alci)*100)


# plot 
# images 
idf.sp <- data.frame(ims = c("./images/auk.png","./images/seal.png"),
                     name = c("alci","seal"),
                     x = 4, 
                     y = 7.5)


ts_grid <- ggplot(data = ts, aes(x = month, y = name)) +
            geom_tile(aes(fill = alci_bpue), col = "black") +
            #scale_fill_gradientn(colors = c("#b9e2e9","#EACC2A","#F21A00"), values = c(0,0.1,0.2,1), na.value = "grey80", guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
            scale_fill_gradientn(colors = c("white","black"), values = c(0,0.3,0.4,1),
                                 guide = "none") +
            labs(x = "Month", y = "Capture Species", fill = "Bycatch") +
            geom_text(data = ts %>% filter(alci_bpue < 0.4), aes(label = round(alci_bpue, 2)), col = "black", size = 1.5) +
            geom_text(data = ts %>% filter(alci_bpue >= 0.4), aes(label = round(alci_bpue, 2)), col = "white", size = 1.5, fontface = "bold") +
            geom_image(data = idf.sp %>% filter(name == "alci"), aes(image = ims, x = x, y = y), size = 0.15) +
            ggpp::geom_text_npc(inherit.aes = FALSE, label = "a)", size = 3, npcx = 0.98, npcy = "top") +
            scale_y_discrete(labels = function(x){case_when(x == "Monkfish species" ~ "Monkfish",
                                                            x == "Skates and Rays" ~ "Skates & Rays",
                                                            x == "Anything" ~ "Mixed",
                                                            TRUE ~ x)}) +
            coord_cartesian(expand = 0) +
            theme(plot.margin = margin(0,0,0.5,0,unit = "line"),
                  axis.title.x = element_blank())

ts_bar <- ggplot(data = ts.ts, aes(x = name, y = alci_bpue)) +
            geom_col(aes(fill = alci_bpue), col = "black") +
            scale_fill_gradientn(colors = c("white","black"), values = c(0,0.3,0.4,1),  
                                 na.value = "grey80", guide = "none") +
            labs(x = "Month", y = "Target Species", fill = "Bycatch\n(ind.)") +
            coord_flip(expand = 0) +
            theme_minimal() +
            theme(plot.margin = margin(0,0,0,0,unit = "line"),
                  panel.grid = element_line(colour = NA), 
                  axis.title = element_blank(),
                  axis.text = element_blank())


mon_bar <- ggplot(data = ts.mon, aes(x = month, y = alci_bpue)) +
            geom_col(aes(fill = alci_bpue), col = "black") +
            scale_fill_gradientn(colors = c("white","black"), values = c(0,0.3,0.4,1), na.value = "grey80",
                                 guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"),
                                 labels = c("Low", "High"),
                                 breaks = range(ts.mon$alci_bpue, na.rm = TRUE)) +
            labs(x = "Month", y = "Target Species", fill = "Bycatch\n") +
            coord_cartesian(expand = expansion(mult = c(0,0))) +
            theme_minimal() +
            theme(plot.margin = margin(0,0,0,0,unit = "line"),
                  panel.grid = element_line(colour = NA), 
                  axis.title = element_blank(),
                  axis.text = element_blank(),,
                  legend.position = c(.8, .5),
                  legend.key.height = unit(0.4, "cm"),
                  legend.key.width = unit(0.3, "cm"),
                  legend.title = element_text(size = 6), 
                  legend.text = element_text(size = 5),   
                  legend.direction = "horizontal")



design <- "AAAA#
           BBBBC
           BBBBC
           BBBBC
           BBBBC"
# add
a <- mon_bar + ts_grid + ts_bar + plot_layout(design = design)
                                        
# plot the table # plot the table colour = 
ms_grid <- ggplot(data = ms, aes(x = month, y = mesh_size)) +
            geom_tile(aes(fill = alci_bpue), col = "black") +
            scale_fill_gradientn(colors = c("white","black"), values = c(0,0.3,0.4,1), na.value = "grey80", 
                                 guide = "none") +
            labs(x = "Month", y = "Mesh Size (mm)", fill = "Bycatch") +
            geom_text(data = ms %>% filter(alci_bpue < 0.4), aes(label = round(alci_bpue,2)), col = "black", size = 1.5) +
            geom_text(data = ms %>% filter(alci_bpue >= 0.4), aes(label = round(alci_bpue,2)), col = "white", size = 1.5, fontface = "bold") +
            ggpp::geom_text_npc(inherit.aes = FALSE, label = "b)", size = 3, npcx = 0.98, npcy = "top", nudge_x = 0.1) +
            #geom_image(data = idf.sp %>% filter(name == "alci"), aes(image = ims, x = x, y = y), size = 0.15) +
            coord_cartesian(expand = 0) +
            theme(plot.margin = margin(0.1,0,0.5,0.5,unit = "line"))
#ms_grid

ms_bar <- ggplot(data = ms.ms, aes(x = mesh_size, y = alci_bpue)) +
            geom_col(aes(fill = alci), col = "black") +
            scale_fill_gradientn(colors = c("white","black"), values = c(0,0.3,0.4,1), na.value = "grey80", guide = "none") +
            labs(x = "Month", y = "Target Species", fill = "Bycatch\n(ind.)") +
            coord_flip(expand = 0) +
            theme_minimal() +
            theme(plot.margin = margin(0,0,0,0,unit = "line"),
                  panel.grid = element_line(colour = NA), 
                  axis.title = element_blank(),
                  axis.text = element_blank())

b <- ms_grid + ms_bar + plot_layout(widths = c(0.8,0.2))


a/b + plot_layout(heights = c(0.68,0.32)) 
ggsave("./plots/Fig4-fishery.jpeg", width = 10, height = 13, dpi = 600, unit = "cm")  
#===============================================================================================================================================================================================================================








