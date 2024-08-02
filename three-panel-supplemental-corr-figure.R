#############Script with statistical analyses for synchrony manuscript########

###
#libraries-------
###

library(tidyverse)
library(ggpubr)
# library(car)
# library(MASS)

#Models libraries
library(mgcv)
library(glmmTMB)
# library(gam)
# library(gamm4)
# library(nlme)

#Model check libraries
library(visreg)
library(DHARMa)
library(performance)
library(lsmeans)
library(emmeans)
library(MuMIn)
library(ggeffects)
library(readr)

######################################################################


#####################Data###############################################

###
#Data with Eadj metric, date, season, and hydrological metrics. Metrics are based
#on monthly observations and divided by water year (May-April)
###

### all data
sim = read_csv("data/spat_sim_allthegoods_01_30_2024.csv")
glimpse(sim)

### sample size info
nsize = read_csv("data/snook_sample_size_eadj.csv") |>  
      dplyr::rename(Year.Month = fYear.Month, eadj_n = n)
glimpse(nsize)

sim_all <- left_join(sim, nsize, by = "Year.Month")
glimpse(sim_all)

### freshwater cont info 
d = read_csv("data/snook_ss01302024_UpdatedwYear_2monthFiltered.csv") |> 
      group_by(year, source) |> 
      summarize(mean = mean(mean)) |> 
      pivot_wider(names_from = source, values_from = mean) |> 
      rename(wYear = year)

dt <- left_join(sim_all, d, by = "wYear")

### summarize data
summ_dt <- dt |> 
      group_by(wYear) |> 
      summarize(
            niche = mean(hv_size),
            niche_n = mean(hv_n),
            sus = mean(Eadj),
            sus_n = mean(eadj_n),
            fw_cont = mean(Freshwater)
      ) |> 
      ungroup()

# generate figures and save with same formating ---------------------------

a <- summ_dt |> 
      ggplot(aes(niche_n, log(niche))) +
      geom_point() + 
      geom_smooth(method = "lm", color = "black") +
      labs(x = "Hypervolume Sample Size", y = "log(Niche Volume)") +
      stat_cor(p.accuracy = 0.01, r.accuracy = 0.01, label.x.npc = .1, , label.y.npc = .9)+
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
            axis.title = element_text(size = 15, face = "bold", colour = "black"), 
            plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) 

c <- summ_dt |> 
      ggplot(aes(fw_cont, log(niche))) +
      geom_point() + 
      geom_smooth(method = "lm", color = "black") +
      labs(y = "log(Niche Volume)", x = "Proportion Freshwater Contribution") +
      stat_cor(p.accuracy = 0.01, r.accuracy = 0.01, label.x.npc = .7, , label.y.npc = 0.9)+
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
            axis.title = element_text(size = 15, face = "bold", colour = "black"), 
            plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) 

b <- summ_dt |> 
      ggplot(aes(niche_n, sus)) +
      geom_point() + 
      geom_smooth(method = "lm", color = "black") +
      labs(y = expression(bold('Mean Annual Space Use Specialization (' *E[adj]* ')')), x = "Hypervolume Sample Size") +
      stat_cor(p.accuracy = 0.01, r.accuracy = 0.01, label.x.npc = .7, label.y.npc = .9)+
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
            axis.title = element_text(size = 15, face = "bold", colour = "black"), 
            plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) 

ggarrange(a, b, c,
          labels = c('a)','b)', 'c)'),
          ncol = 3, vjust = 1, align = "v")

#saving for publication
ggsave("./figs/manuscript/supp-three-panel-corr.tiff", units = "in", width = 15,
       height = 6, dpi =  600, compression = "lzw")

sim_all$Season<-factor(sim_all$Season, levels = c("Wet", "Dry"))

sim_all |> 
      ggplot(aes(x=as.factor(wYear), y=Eadj, fill=Season)) +
      geom_boxplot(width =0.8) +
      scale_fill_manual(values = c("cadetblue", "darkgoldenrod")) +
      labs(x = "Hydrologic Year", 
           y = expression(bold('Space Use Specialization (' *E[adj]* ')'))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=14, face="bold", color = "black")) +
      theme(axis.text = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.x = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.y = element_text(size=16,face="bold", color = "black")) +
      theme(axis.title = element_text(size=16,face="bold", color = "black")) +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size=16, face="bold", color = "black")) +
      theme(legend.position = c(0.95, 0.9))

ggsave("./figs/manuscript/supp-eadj-seasonal-boxplot.tiff", units = "in", width = 10,
       height = 6, dpi =  600, compression = "lzw")
