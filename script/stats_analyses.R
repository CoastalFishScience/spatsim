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
q1.m1 <- gam(Eadj ~ s(wMonth, bs = "cc", k = 10) + s(wYear), 
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
#saving model comparison
compare_performance(q1.m1, q1.m2, q1.m3) |> write_csv("./tables/q1.compare.models.performance.csv")
compare_performance(q1.m1, q1.m2, q1.m3) |> capture.output(file = "./tables/q1.compare.models.performance.txt")

summary(q1.m1) |> capture.output(file = "./tables/q1.summary.models.txt")
plot(q1.m1)
#Best model - q1.m1 based on AICc and r2 metrics


#Model output - seasonal trend
q1.m1.vis = visreg(q1.m1, type = "conditional", scale = "response")
#extracting fitted model from visreg to allow plotting in ggplot
q1.m1.month<-q1.m1.vis[[1]]$fit

seasonal_trend<-ggplot(q1.m1.month, aes(wMonth, visregFit))+ 
      geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), color = "grey60", alpha = .2)+
      geom_line(linewidth = 2, colour= "black") + theme_bw()+
      geom_line(linetype = 2, colour = "black", aes(y = visregLwr))+
      geom_line(linetype = 2, colour = "black", aes(y = visregUpr)) +
      labs(x = "Water Year Month (May-April)", y = "Eadj", title = "Seasonal trend")+
      scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
      ylim(0.3,0.8) +
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
            axis.title = element_text(size = 15, face = "bold", colour = "black"), 
            plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) 

#Model ouptut - Yearly variability 
q1.m1.year<-q1.m1.vis[[2]]$fit

yearly_vars<-ggplot(q1.m1.year, aes(wYear, visregFit))+ 
      geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), color = "grey60", alpha = .2)+
      geom_line(size = 2, colour= "black") + theme_bw()+
      geom_line(linetype = 2, colour = "black", aes(y = visregLwr))+
      geom_line(linetype = 2, colour = "black", aes(y = visregUpr)) +
      labs(x = "Water Year (May-April)", y = "Eadj", title = "Inter-annual trend")+
      scale_x_continuous(breaks = c(2012, 2014, 2016, 2018, 2020, 2022)) +
      ylim(0.3,0.8)+
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
            axis.title = element_text(size = 15, face = "bold", colour = "black"), 
            plot.title = element_text(size = 16, face = "bold", colour = "black", hjust = 0.5), 
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())


#Plot of fitted model q1.m1

ggarrange(seasonal_trend, yearly_vars,
          labels = c('a)','b)'),
          ncol = 1, vjust = 1, align = "v")

#saving for publication
ggsave("./figs/manuscript/q1.fittedmodels.plot.tiff", units = "in", width = 6,
       height = 10, dpi =  600, compression = "lzw")

#Saving model output for publication
summary.final.q1model<-summary(q1.m1)
capture.output(summary.final.q1model, file="tables/summary.final.q1model.txt")

#########Q2 - Hydrological effects############################################

###
#GAMs/GLMs to assess how hydrological parameters influence similarity of habitat use-----
###

#Checking colinearity------- 
library(PerformanceAnalytics)
library(ggcorrplot)

cross_corr = sim |> select(monthly_mean_stage_cm, monthly_mean_stage_cm_prev,
                               daysbelow30) |> 
      chart.Correlation(method = "kendall", histogram=TRUE, pch=19)

tiff("./figs/manuscript/covariate-corr.tiff",width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
cross_corr = sim |> select(monthly_mean_stage_cm, monthly_mean_stage_cm_prev,
                            daysbelow30) |> 
      chart.Correlation(method = "kendall", histogram=TRUE, pch=19)
dev.off()

#Single variable model comparisons------
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

q2.model.compare = compare_performance(q2.glm.stage, q2.gam.stage,
                        q2.glm.Lstage, q2.gam.Lstage,
                        q2.glm.days, q2.gam.days)

q2.model.compare |> write_csv("./tables/q2.compare.models.performance.csv")
q2.model.compare |> capture.output(file = "./tables/q2.compare.models.performance.txt")


# map(q2.prelim.models, check_model)

#Checking/comparing prelim models to decide glm vs gam-----

#Best = gam
compare_performance(q2.glm.stage, q2.gam.stage) #|> capture.output(file = "./tables/q2.stage.performance.txt")
#Best = gam
compare_performance(q2.glm.Lstage, q2.gam.Lstage) #|> capture.output(file = "./tables/q2.Lstage.performance.txt")
#Best = gam
compare_performance(q2.glm.days, q2.gam.days) #|> capture.output(file = "./tables/q2.days.performance.txt")

#################################################################################################
###Code lines to compare full models; however, we decided to present the two best-------- 
###independently due to colinearity issues
#Checking which stage variable to use since high correlation
#between them
#Best = gam with stage instead of lagged stage
# compare_performance(q2.gam.stage, q2.gam.Lstage)
# 
# #Full model--------
# q2.glm.m1 = glmmTMB(Eadj ~ daysbelow30*mean_stage, 
#                     data = sim_all,
#                     family = beta_family(link = "logit")) 
# 
# check_model(q2.glm.m1)
# options(na.action = "na.fail")
# dredge(q2.glm.m1, rank = "AICc")
# 
# #checking performance of gam full model considering the non-linear trend
# #of residuals
# q2.gam.m1 = gam(Eadj ~ s(daysbelow30) + s(mean_stage), data = sim_all,
#                 family = betar(link = "logit"), 
#                 method = "REML")
# 
# q2.gam.m2 = gam(Eadj ~ daysbelow30 + s(mean_stage), data = sim_all,
#                 family = betar(link = "logit"), 
#                 method = "REML")
# 
# #Model performance similar but GAM still best model
# compare_performance(q2.glm.m1, q2.gam.m1, q2.gam.m2)
# check_model(q2.gam.m1)
# 
# 
# #Model simplification of gam
# dredge(q2.gam.m1, rank = "AICc")
# dredge(q2.glm.m1, rank = "AICc") |> write_csv()
# capture.output(summary.final.q1model, file="tables/summary.final.q1model.txt")
# 
# #Best model - only with daysbelow30
# q2.glm.m2 = glmmTMB(Eadj ~ daysbelow30,
#              data = sim_all, 
#              family = beta_family(link = "logit"))
# 
# #Plot - fitted model for Q2----------------
# q2.m2.vis = visreg(q2.glm.m2, type = "conditional", scale = "response")
# 
# q2.m2.fit<-q2.m2.vis$fit
# 
# Enhydro_trend<-ggplot(q2.m2.fit, aes(daysbelow30, visregFit))+ 
#       geom_line(size = 2, colour= "red") + theme_bw()+
#       geom_line(linetype = 2, colour = "red", aes(y = visregLwr))+
#       geom_line(linetype = 2, colour = "red", aes(y = visregUpr)) +
#       labs(x = "Days with stage below 30 cm", y = "Eadj")+
#       theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
#             axis.title = element_text(size = 16, face = "bold", colour = "black"), 
#             plot.title = element_text(size = 16, face = "bold", colour = "black"), 
#             panel.grid.major = element_blank(),
#             axis.line = element_line(colour = "black"),
#             panel.grid.minor = element_blank(),
#             panel.border = element_blank(),
#             panel.background = element_blank()) 
#             



#Plot - fitted model for Q2 based on single models due to colinearity----------------
#####################################################################################
#####################################################################################

###
#Model presentation to account for collinearity----
###
q2.gam.stage = gam(Eadj ~ s(mean_stage), 
                   data = sim_all,
                   family = betar(link = "logit"), 
                   method = "REML")
summary.final.q2model1<-summary(q2.gam.stage)
capture.output(summary.final.q2model1, file="tables/summary.final.q2model1.txt")


q2.gam.days = gam(Eadj ~ s(daysbelow30), 
                  data = sim_all,
                  family = betar(link = "logit"), 
                  method = "REML")
summary.final.q2model2<-summary(q2.gam.days)
capture.output(summary.final.q2model2, file="tables/summary.final.q2model2.txt")

#daysbelow30 effects
q2.days.vis = visreg(q2.gam.days, type = "conditional", scale = "response")
q2.days.fit<-q2.days.vis$fit

days_vs_E<-ggplot(q2.days.fit, aes(daysbelow30, visregFit))+
      geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), color = "grey60", alpha = .2)+
      geom_line(size = 2, colour= "black") + theme_bw()+
      geom_line(linetype = 2, colour = "black", aes(y = visregLwr))+
      geom_line(linetype = 2, colour = "black", aes(y = visregUpr)) +
      labs(x = "Days with stage below 30 cm", y = "Eadj")+
      scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
      ylim(0.3,0.8)+
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
            axis.title = element_text(size = 16, face = "bold", colour = "black"),
            plot.title = element_text(size = 16, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())
      
      
#water stage effects
q2.stage.vis = visreg(q2.gam.stage, type = "conditional", scale = "response")
q2.stage.fit<-q2.stage.vis$fit

stage_vs_E<-ggplot(q2.stage.fit, aes(mean_stage, visregFit))+
      geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr), color = "grey60", alpha = .2)+
      geom_line(size = 2, colour= "black") + theme_bw()+
      geom_line(linetype = 2, colour = "black", aes(y = visregLwr))+
      geom_line(linetype = 2, colour = "black", aes(y = visregUpr)) +
      labs(x = "River Stage (cm)", y = "Eadj")+
      scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30)) +
      ylim(0.3,0.8)+
      theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),
            axis.title = element_text(size = 16, face = "bold", colour = "black"),
            plot.title = element_text(size = 16, face = "bold", colour = "black"),
            panel.grid.major = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank())

#Plot of fitted model q1.m1
ggarrange(stage_vs_E, days_vs_E,
          labels = c('a)','b)'),
          ncol = 1, vjust = 1, align = "v")

#saving for publication
ggsave("./figs/manuscript/q2.fittedmodels.plot.tiff", units = "in", width = 6,
       height = 10, dpi =  600, compression = "lzw")


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
      geom_smooth(method = "lm", se = TRUE) +
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

corr_all = ggplot(dat_summary, aes(x = mean_e, y = hv_size)) + 
      geom_point() +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      labs(x = "Eadj Dry Season Mean", 
           y = "Niche Volume") +
      stat_cor(p.accuracy = 0.01, r.accuracy = 0.01, label.x.npc = .5)+
      xlim(0.4,0.7) +
      # ylim(0, max(dat_summary$hv_size))+
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

corr_NOoutlier = ggplot(filter(dat_summary, hv_size < 200), aes(x = mean_e, y = hv_size)) + 
      geom_point() +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      labs(x = "Eadj Dry Season Mean", 
           y = "Niche Volume") +
      stat_cor(p.accuracy = 0.01, r.accuracy = 0.01, label.x.npc = .5)+
      xlim(0.4,0.7) +
      # ylim()
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

#Plot of fitted model q1.m1
ggarrange(corr_all, corr_NOoutlier,
          labels = c('a)','b)'),
          ncol = 1, vjust = 1, align = "v")

#saving for publication
ggsave("./figs/manuscript/q3.fittedmodels.plot.tiff", units = "in", width = 6,
       height = 10, dpi =  600, compression = "lzw")

