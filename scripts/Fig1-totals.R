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

# clean environment prior to starting
cat("\014") 
rm(list=ls())

# data 
dat <- read_csv("./data/bycatch.csv") %>% 
           mutate(net_description = case_when(net_description == "CONTROL" ~ "Ctrl.",
                                              net_description == "KITE" ~ "Kite",
                                              net_description == "LEB" ~ "LEB"),
                  common_name = case_when(common_name == "Shag/Cormorant Species" ~ "Phalacrocoracid spp.",
                                          common_name == "Unidentified Toothed Whales Species" ~ "Odontocete spp.",
                                          common_name == "Unidentified Seal Species" ~ "Pinniped spp.",
                                          TRUE ~ common_name)) 
           

# summarise 
sum <- dat %>% group_by(group, common_name, net_description) %>% tally()
sum %>% filter(group != "Bird") %>% group_by(group, net_description) %>% summarise(inds = sum(n)) %>% mutate(prop = inds/sum(inds)*100)# analysis 


# plotting 
labs <- sum %>% group_by(group, common_name) %>% summarise(x = max(n), tot = sum(n))
order <- labs %>% arrange(desc(tot)) %>% pull(common_name)
order <- c("Common Guillemot","Razorbill","Grey Seal" ,"Common Dolphin","Common Seal","Pinniped spp.",
           "Shag","Cormorant","Phalacrocoracid spp.","Great Northern Diver","Harbour Porpoise","Odontocete spp.","Gannet")
sum <- sum %>% mutate(common_name = factor(common_name, levels = rev(order)),
                      net_description = factor(net_description, levels = rev(c("Ctrl.","LEB","Kite"))))


ggplot(data = sum, aes(y = common_name, x = n)) +
                  geom_col(aes(fill = net_description), position = "dodge") +
                  labs(x = "Bycatch (ind.)", y = "", fill = "Net Type") +
                  geom_text(inherit.aes = FALSE, data = labs, aes(y = common_name, x=x, label = tot), colour = "grey30", size = 2, fontface = "plain", hjust = -0.2) +
                  scale_x_continuous(expand = expansion(mult = c(0.01, 0.1))) +
                  scale_fill_manual(values = c("#FF6A1A","#00528C","black")) + #, labels = c("Illuminated","Non-illuminated")) +
                  guides(fill = guide_legend(reverse=TRUE, title.position = "top", title.hjust = 0.5)) +
                  facet_grid(rows = vars(group), scales = "free", space = "free_y")+#, labeller = labeller(pot_type = facet.labs)) +
                  coord_cartesian(clip = "off") +
                  theme(strip.background = element_rect(fill = NA, colour = NA),
                        strip.text = element_blank(),
                        plot.margin = margin(0,0,0,0),
                        panel.grid.minor.x = element_line(colour = "grey", linewidth = 0.1),
                        panel.grid.major.x = element_line(colour = "grey", linewidth = 0.1))

ggsave("./plots/Fig1-totals.jpeg", width = 6, height = 8, dpi = 600, unit = "cm")  





