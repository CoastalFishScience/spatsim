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


################################################################################


############################Analyzes##########################################

      
#########Q1 - Temporal trends in habitat use spatial similarity################

###
#GAMs to assess the seasonality of Eadj and Eadj annual trend (inter year var)
###

#Model selection
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
#Best model - q1.m1 based on AICc and r2 metrics


#Model output - seasonal trend
q1.m1.vis = visreg(q1.m1, type = "conditional", scale = "response")

q1.m1.month<-q1.m1.vis[[1]]$fit

seasonal_trend<-ggplot(q1.m1.month, aes(wMonth, visregFit))+ 
      geom_line(size = 2, colour= "red") + theme_bw()+
      geom_line(linetype = 2, colour = "red", aes(y = visregLwr))+
      geom_line(linetype = 2, colour = "red", aes(y = visregUpr)) +
      # geom_rect(aes(xmin=1999, xmax=2003.99, ymin=-Inf, ymax=Inf), alpha = 0.01, fill = "grey80")+ 
      # geom_rect(aes(xmin=2004, xmax=2008.99, ymin=-Inf, ymax=Inf), alpha = 0.01, fill = "grey66")+
      # geom_rect(aes(xmin=2009, xmax=2013.99, ymin=-Inf, ymax=Inf), alpha = 0.01, fill = "grey34")+ 
      # geom_rect(aes(xmin=2014, xmax=2016, ymin=-Inf, ymax=Inf), alpha = 0.01, fill = "grey17")+
      labs(x = "Water Year Month (May-April)", y = "Eadj")+
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
            axis.title = element_text(size = 16, face = "bold", colour = "black"), 
            plot.title = element_text(size = 16, face = "bold", colour = "black"), 
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) 

#Model ouptut - Yearly variability 
q1.m1.year<-q1.m1.vis[[2]]$fit

yearly_vars<-ggplot(q1.m1.year, aes(wYear, visregFit))+ 
      geom_line(size = 2, colour= "red") + theme_bw()+
      geom_line(linetype = 2, colour = "red", aes(y = visregLwr))+
      geom_line(linetype = 2, colour = "red", aes(y = visregUpr)) +
      labs(x = "Water Year (May-April)", y = "Eadj")+
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
            axis.title = element_text(size = 16, face = "bold", colour = "black"), 
            plot.title = element_text(size = 16, face = "bold", colour = "black"), 
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())


#Plot of fitted model q1.m1

ggarrange(seasonal_trend, yearly_vars,
          labels = c('a)','b)'),
          ncol = 1)

#########Q2 - Hydrological effects############################################

###
#GAMs/GLMs to assess how hydrological parameters influence similarity of habitat use
###
library(PerformanceAnalytics)
library(ggcorrplot)

cross_corr = sim_all |> select(monthly_mean_stage_cm, monthly_mean_stage_cm_prev,
                               daysbelow30) |> 
      chart.Correlation(method = "kendall", histogram=TRUE, pch=19)

sim_all = sim_all |> 
      rename(mean_stage = monthly_mean_stage_cm,
             mean_Lstage = monthly_mean_stage_cm_prev)

q2.glm.stage = glmmTMB(Eadj ~ mean_stage, 
                    data = sim_all,
                    family = beta_family(link = "logit"))

q2.glm.Lstage = glmmTMB(Eadj ~ mean_Lstage, 
                       data = sim_all,
                       family = beta_family(link = "logit"))

q2.glm.days = glmmTMB(Eadj ~ daysbelow30, 
                        data = sim_all,
                        family = beta_family(link = "logit"))

q2.gam.stage = gam(Eadj ~ s(mean_stage), 
                   data = sim_all,
                   family = betar(link = "logit"), 
                   method = "REML")

q2.gam.Lstage = gam(Eadj ~ s(mean_Lstage), 
                   data = sim_all,
                   family = betar(link = "logit"), 
                   method = "REML")

q2.gam.days = gam(Eadj ~ s(daysbelow30), 
                  data = sim_all,
                  family = betar(link = "logit"), 
                  method = "REML")

q2.prelim.models = list(q2.glm.stage, q2.gam.stage,
                        q2.glm.Lstage, q2.gam.Lstage,
                        q2.glm.days, q2.gam.days)

map(q2.prelim.models, check_model)

#Checking/comparing prelim models to decide glm vs gam-----

#Best = glm 
compare_performance(q2.glm.stage, q2.gam.stage)
#Best = glm
compare_performance(q2.glm.Lstage, q2.gam.Lstage)
#Best = glm
compare_performance(q2.glm.days, q2.gam.days)

#Checking which stage variable to use since high correlation
#between them
#Best = glm with stage instead of lagged stage
compare_performance(q2.glm.stage, q2.glm.Lstage)

#Full model--------
q2.glm.m1 = glmmTMB(Eadj ~ daysbelow30*mean_stage, 
                    data = sim_all,
                    family = beta_family(link = "logit")) 

check_model(q2.glm.m1)
options(na.action = "na.fail")
dredge(q2.glm.m1, rank = "AICc")

#checking performance of gam full model considering the non-linear trend
#of residuals
q2.gam.m1 = gam(Eadj ~ s(daysbelow30) + s(mean_stage), data = sim_all,
                family = betar(link = "logit"), 
                method = "REML")

#glm still best model
compare_performance(q2.glm.m1, q2.gam.m1)

#Model simplification of glm 
dredge(q2.glm.m1, rank = "AICc")

#Best model - only with daysbelow30
q2.glm.m2 = glmmTMB(Eadj ~ daysbelow30,
             data = sim_all, 
             family = beta_family(link = "logit"))

#Plot - fitted model for Q2----------------
q2.m2.vis = visreg(q2.glm.m2, type = "conditional", scale = "response")

q2.m2.fit<-q2.m2.vis$fit

Enhydro_trend<-ggplot(q2.m2.fit, aes(daysbelow30, visregFit))+ 
      geom_line(size = 2, colour= "red") + theme_bw()+
      geom_line(linetype = 2, colour = "red", aes(y = visregLwr))+
      geom_line(linetype = 2, colour = "red", aes(y = visregUpr)) +
      labs(x = "Days with stage below 30 cm", y = "Eadj")+
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
            axis.title = element_text(size = 16, face = "bold", colour = "black"), 
            plot.title = element_text(size = 16, face = "bold", colour = "black"), 
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) 

###############################################################################


#######Q3 - correlation between Eadj and Hv#####################

#data wraggling to correlate dry montly Eadj and HV
dat_summary <- sim |> 
      group_by(wYear, Season) |> 
      summarise(mean_e = mean(Eadj, na.rm = "TRUE"),
                min_e = min(Eadj, na.rm = "TRUE"),
                max_e = max(Eadj, na.rm = "TRUE"),
                hv_size = mean(hv_size, na.rm = "TRUE"),
                hv_n = mean(hv_n)) |> 
      filter(Season == "Dry")

# run correlation here for hv_size and hv_size sample size---------

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

# run correlation here for hv_size and mean_e-----------------

#HV failed normality test
shapiro.test(dat_summary$hv_size)
shapiro.test(dat_summary$mean_e)

#nonparametric tests are not significant but pearson is when eliminating outlier
#We need to discuss next steps
cor.test(dat_summary$hv_size, dat_summary$mean_e, method = "pearson")

dat_summary2 = filter(dat_summary, hv_size < 200) 
cor.test(dat_summary2$hv_size, dat_summary2$mean_e, method = "pearson")

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
