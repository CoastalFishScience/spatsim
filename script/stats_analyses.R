#############Script with statistical analyses for synchrony manuscript########

###
#libraries-------
###

library(tidyverse)
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

######################################################################


#####################Data###############################################

###
#Data with Eadj metric, date, season, and hydrological metrics. Metrics are based
#on monthly observations and divided by water year (May-April)
###

#all data
sim = read_csv("./data/spat_sim_allthegoods_01_15_2023.csv")
glimpse(sim)

#sample size info
nsize = read_csv("./data/snook_sample_size_eadj.csv") |>  
      dplyr::rename(Year.Month = fYear.Month, eadj_n = n)
glimpse(nsize)

sim_all <- left_join(sim, nsize, by = "Year.Month")
glimpse(sim_all)



###############################################################################


#######################Exploratory Plots#######################################

#colored by year
ggplot(sim_all, aes(x = eadj_n, y = Eadj, color = wYear)) + 
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      labs(x = "Eadj Sample Size", 
           y = "Eadj") +
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
      theme(legend.position = c(0.9, 0.9))

#colored by stage
ggplot(sim_all, aes(x = eadj_n, y = Eadj, color = monthly_mean_stage_cm)) + 
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      labs(x = "Eadj Sample Size", 
           y = "Eadj") +
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
      theme(legend.position = c(0.9, 0.9))

#colored by season
ggplot(sim_all, aes(x = eadj_n, y = Eadj, color = Season)) + 
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      labs(x = "Eadj Sample Size", 
           y = "Eadj") +
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
      theme(legend.position = c(0.9, 0.9))

dat_summary <- sim |> 
      group_by(wYear, Season) |> 
      summarise(mean_e = mean(Eadj),
                min_e = min(Eadj),
                max_e = max(Eadj),
                hv_size = mean(hv_size),
                hv_n = mean(hv_n)) |> 
      filter(Season == "Dry")

# run correlation here for hv_size and hv_size sample size

ggplot(dat_summary, aes(x = hv_n, y = hv_size)) + 
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      labs(x = "Hypervolume Sample Size", 
           y = "Niche Volume") +
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
      theme(legend.position = c(0.9, 0.9))

# run correlation here for hv_size and mean_e

ggplot(dat_summary, aes(x = mean_e, y = hv_size)) + 
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      labs(x = "Eadj Dry Season Mean", 
           y = "Niche Volume") +
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
      theme(legend.position = c(0.9, 0.9))

ggplot(filter(dat_summary, hv_size < 200), aes(x = mean_e, y = hv_size)) + 
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      labs(x = "Eadj Dry Season Mean", 
           y = "Niche Volume") +
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
      theme(legend.position = c(0.9, 0.9))
################################################################################


############################Analyzes##########################################

      
#########Q1 - Temporal trends in habitat use spatial similarity################

###
#GAMs to assess the seasonality of Eadj and Eadj annual trend (inter year var)
###

q1.m1 <- gam(Eadj ~ s(wMonth, bs = "cc", k = 10) + s(wYear, bs = "cr", k = 8), 
          data = sim_all,
          family = betar(link = "logit"), 
          method = "REML")
gam.check(q1.m1)

q1.m2 <- gam(Eadj ~ s(wMonth, bs = "cc", k = 10) + wYear, 
             data = sim_all,
             family = betar(link = "logit"), 
             method = "REML")

q1.m3 <- glmmTMB(Eadj ~ wMonth + wYear, 
             data = sim_all,
             family = beta_family(link = "logit"))

compare_performance(q1.m1, q1.m2, q1.m3)

summary(q1.m1)
plot(q1.m1)
visreg(q1.m1, type = "conditional", scale = "response")

