### testing testing testing
###libraries####
library(ggplot2)
# library(lattice)
library(car)
library(reshape2)
library(reshape)
# library(TSA)
library(plyr)
library(dplyr)
library(tidyr)
library(visreg)
# library(MASS)
library(modEvA) 
# library(BiodiversityR)
library(gridExtra)
library(AICcmodavg)
library(nlme)
library(mgcv)
library(lme4)
# library(VTrack)
# library(igraph)
# library(MARSS)
library(splitstackshape) ###package to use cSplit and do text to column like excel
library(chron)
# library(rgdal)#package for geospatial data
library(RInSp)#package for intraspecific niche variation
library(boot)#package for boostrapping operations

library(cowplot)
library(ggpubr)
library(lsmeans)


setwd("E:/MovementMonday/Paper2/synchrony/2018")

tracks<-readRDS("AllSnook_PeriodOfRecord_08272019.rds")

snook<-subset(tracks, Species %in% c("Snook"))%>%
  subset(., select = c("Transmitter", "Station.Name", "Datetime_UTC",
                       "Receiver", "Species", "Station", 
                       "Distance", "Longitude", "Latitude",
                       "n"))

######Processing################################

snook$DateTime2 <- snook$Datetime_UTC
snook<-cSplit(snook, "DateTime2", sep = " ", type.convert = FALSE)
colnames(snook)[11:12]<-c("Date", "Time")

snook$Date2 <- snook$Date
snook<-cSplit(snook, "Date2", sep = "-", type.convert = FALSE)
colnames(snook)[13:15]<-c("Year", "Month", "Day")

snook<-cSplit(snook, "Transmitter", sep = "-", type.convert = FALSE)
colnames(snook)[15:17]<-c("code1", "serial", "ID")

snook<-cSplit(snook, "Receiver", sep = "-", type.convert = FALSE)
colnames(snook)[17:18]<-c("code", "Receiver.ID")

snook$Home<-cut(snook$Distance, breaks = c(-1,15,23,31.72), 
                labels = c("Lower River", "Bay", "Upper River"))

#Classifying columns
snook2<-snook

snook2$fYear.Month <- paste(snook2$Year, snook2$Month, sep = "/")

# snook2$Month<-as.numeric(snook2$Month)
# snook2$Day<-as.numeric(snook2$Day)

snook2$Year<-as.numeric(snook2$Year)
snook2$ID<-as.numeric(snook2$ID)
snook2$Date<-dates(snook2$Date, format = c(dates = "y-m-d"))

snook2$fMonth<-factor(snook2$Month)
snook2$fDay<-factor(snook2$Day)
snook2$fYear<-factor(snook2$Year)
snook2$Home<-factor(snook2$Home)
snook2$f.Distance<-factor(snook2$Distance)
snook2$ID<-factor(snook2$ID)
snook2$Receiver.ID<-factor(snook2$Receiver.ID)

# zones2018 <- select(snook, Station, Distance, Longitude, Latitude) %>%
#   distinct()

###Repeat analysis but with the reclassification of the receivers by zone similar to Matich et al 2017######
zones<-read.csv("zones2019.csv")
zones$f.Distance<-factor(zones$Distance)
zones<-subset(zones, Zone != "NA")
snook2<-merge(snook2, zones, by = "f.Distance", all.x = TRUE)

snook2<-subset(snook2, Year > 2011 & Year < 2020)




saveRDS(snook2, file = "snook_SR_detection_20122019.rds")




###############################################################################
#######Calculate the range of days of detection period and number of detections#######

snook2<-readRDS("snook_SR_detection_20122019.rds")


snook.detdays<-group_by(snook2, ID)%>%
                summarise(det_f = sum(!is.na(fDay)),
                          min.t = min(Datetime_UTC),
                          max.t = max(Datetime_UTC),
                          resid.t = difftime(max.t, min.t, units = "day"),
                          unique.date = length(unique(fYear.Month)))

#Subset to obtain individuals over 90 days and over 100 detections
snook.90day <- subset(snook.detdays, resid.t >= 90 & det_f > 100)

#Filter tag IDs with more that 90 days 
snook.90day <- filter(snook2, ID %in% snook.90day$ID)


#Creating table of unique tags per year for paper/ppt
snook.90day.x <- filter(snook.90day, Year < 2019)
length(unique(snook.90day.x$ID))

tags.summary <- group_by(snook.90day.x, Year)%>%
  summarise(n = length(unique(ID)))

list.tags <- group_by(snook.90day.x, ID, Year)%>%
  summarise(Total.Det = sum(!is.na(Datetime_UTC)),
            Total.Days = length(unique(fYear.Month)))

write.csv(tags.summary, "tags.summary.csv")
write.csv(list.tags, "list.tags.csv")


#Creating the number of detection as a metric of habitat use per ID|Zone|Year|Year.Month
# snook_yrd_zones<-ddply(snook.90day, c("ID", "Zone", "Year","fYear.Month"), summarise,
#                         freq.v = sum(!is.na(f.Distance), na.rm = TRUE))

#Creating the number of days as a metric of habitat use per ID|Zone|Year|Year.Month
snook_yrd_zones<-ddply(snook.90day, c("ID", "Zone", "Year","fYear.Month"), summarise,
                       freq.v = length(unique(fDay)))


snook_yrd_zones<-subset(snook_yrd_zones, fYear.Month != "NA")
snook_yrd_zones2<-subset(snook_yrd_zones, Zone != "NA")

##Create list of tracks vs receivers listed by year.month
dcast.snook<-function(x){
  dcast(x, ID ~ Zone, sum, value.var = "freq.v")}

list.snook.zones<-dlply(snook_yrd_zones2, .(fYear.Month), dcast.snook)

fc <- function(d, i){
  #d1<-d[i,]
  d2<-import.RInSp(d, row.names = 1, col.header = TRUE)
  EAdj<-Emc(d2, popd.type = "average", replicates = 999)
  return(EAdj$Eadj)
}

E_year.month2<-lapply(list.snook.zones, fc)




###UPDATING DATAFRAME###
library(stringr)
E_year.month2_df<-ldply(E_year.month2, data.frame)
colnames(E_year.month2_df)<-c("Year.Month", "Eadj")


########################################Old code that I used when I was using Year.Month in numeric format##################
# E_year.month2_df$Year.Month<-as.numeric(E_year.month2_df$Year.Month)
# E_year.month2_df[,"Year.Month"] = round(E_year.month2_df[,"Year.Month"],2)
# E_year.month2_df$fYear.Month<-factor(E_year.month2_df$Year.Month)
# E_year.month2_df$fYear.Month<-str_replace_all(E_year.month2_df$fYear.Month, "[.]", "-")
# 
# E_year.month2_df$fYear.Month<-str_replace_all(E_year.month2_df$fYear.Month, c("-08" = "-Feb", "-17" = "-Mar", "-25" = "-Apr","-33" = "-May","-42" = "-Jun",
#                                                                               "-5\\b" = "-Jul", "-58" = "-Aug", "-67" = "-Sep", "-75" = "-Oct", "-83" = "-Nov", "-92" = "-Dec"))
# sub<-subset(E_year.month2_df, fYear.Month %in% c("2012", "2013", "2014", "2015", "2016"))
# sub$fYear.Month<-str_replace_all(sub$fYear.Month, c("2012\\b" = "2012-Jan", "2013\\b" = "2013-Jan", "2014\\b" = "2014-Jan", "2015\\b" = "2015-Jan", "2016\\b" = "2016-Jan"))
# E_year.month2_df<-rbind(E_year.month2_df, sub)
# E_year.month2_df<-E_year.month2_df[ which( ! E_year.month2_df$fYear.Month %in% c("2012", "2013", "2014", "2015", "2016")) , ]
# E_year.month2_df$fYear.Month2<-E_year.month2_df$fYear.Month



####################################################################################################################################


library(zoo)
#using function zoo::yearmon to classify in R fYear.Month as date formated as Year|Month
E_year.month2_df$fYear.Month <- as.yearmon(E_year.month2_df$Year.Month, "%Y/%m")

E_year.month2_df$Date2<-E_year.month2_df$Year.Month
E_year.month2_df<-cSplit(E_year.month2_df, "Date2", sep = "/", type.convert = FALSE)
colnames(E_year.month2_df)[4:5]<-c("Year", "Month")

E_year.month2_df$Year <- as.numeric(E_year.month2_df$Year)
E_year.month2_df$fYear<-factor(E_year.month2_df$Year, levels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"))


E_year.month2_df$fMonth<-factor(E_year.month2_df$Month)
levels(E_year.month2_df$fMonth) <- list("Jan" = c("01"), "Feb" = c("02"), "Mar" = c("03"), 
  "Apr" = c("04"), "May" = c("05"), "Jun" = c("06"), 
  "Jul" = c("07"),"Aug" = c("08"), "Sep" = c("09"), 
  "Oct" = c("10"), "Nov" = c("11"), "Dec" = c("12"))

E_year.month2_df$Season<-factor(E_year.month2_df$fMonth)
levels(E_year.month2_df$Season)<-list("Dry" = c("Jan", "Feb", "Mar", "Apr", "May", "Dec"), 
                                      "Wet" = c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))

# E_year.month2_df$nMonth<-(E_year.month2_df$Month)
# levels(E_year.month2_df$nMonth)<-list("1" = "Jan", "2" = "Feb", "3" = "Mar", "4" = "Apr", 
#                                       "5" = "May", "6" = "Jun", "7" = "Jul", "8" =  "Aug", "9" = "Sep", "10" = "Oct", "11"  = "Nov", "12" = "Dec")

E_year.month2_df$nMonth<-as.numeric(E_year.month2_df$Month)


E_year.month2_df<-subset(E_year.month2_df, Year < 2019)
E_year.month2_df$Date <- as.Date(paste(E_year.month2_df$Year, E_year.month2_df$fMonth, "15", sep = "-"),format = "%Y-%b-%d")
E_year.month2_df <- transform(E_year.month2_df, Date.Num = as.numeric(Date)/1000)

# saveRDS(E_year.month2_df,"E_year.month2_df.rds")
#This version is based on E calc using unique number of days
saveRDS(E_year.month2_df,"E_year.month2_df3.rds")
##############################################################################################################################




###########################Stats Analyses#############################################################

###
#Data Prep
###

# E_year.month2_df <- readRDS("E_year.month2_df3.rds")
# 
# 
# E_year.month2_df<-subset(E_year.month2_df, Year < 2019)
# E_year.month2_df$Date <- as.Date(paste(E_year.month2_df$Year, E_year.month2_df$fMonth, "15", sep = "-"),format = "%Y-%b-%d")
# E_year.month2_df <- transform(E_year.month2_df, Date.Num = as.numeric(Date)/1000)
# 
# # saveRDS(E_year.month2_df,"E_year.month2_df.rds")
# #This version is based on E calc using unique number of days
# saveRDS(E_year.month2_df,"E_year.month2_df3.rds")
# 


E_year.month2_df <- readRDS("E_year.month2_df3.rds")




######################Statistical Tests###################################################################################################################



#################Question 1: Temporal Models with GAMs############################################

####Model testing##################
#Based on Zone as resource and looking at E monthly cycle per year

m1_E2_1 <-gam(Eadj ~ s(nMonth, bs = "cc", k = 12, by = fYear), data = E_year.month2_df, method = "REML")
m1_E2_2 <-gam(Eadj ~ s(nMonth, bs = "cc", k = 12, by = fYear), data = E_year.month2_df, family=gaussian(link = "log"), method = "REML")
m1_E2_3 <-gam(Eadj ~ s(nMonth, bs = "cc", k = 12, by = fYear), data = E_year.month2_df, family=gaussian(link = "logit"), method = "REML")
m1_E2_4 <-gam(Eadj ~ s(nMonth, bs = "cc", k = 12, by = fYear), data = E_year.month2_df,family=betar(link = "logit"), method = "REML")

####Testing difference of models with different faimly/link distribution####
AIC(m1_E2_1, m1_E2_2, m1_E2_3, m1_E2_4)
library(MuMIn)
gam_models <- list(m1_E2_1, m1_E2_2, m1_E2_3, m1_E2_4)

model.sel.q1 <- model.sel(gam_models)
capture.output(model.sel.q1, file = "model.sel.q1.txt")

#best model - m1_E2_4

###
#Extract results and fit to make table/plot for paper
###

summary.model.temporal<-summary(m1_E2_4)
# dput(summary.model.temporal, file = "summary_gam_temp.txt", control = "all")
capture.output(summary.model.temporal, file="summary_gam_temp.txt")

v1<-visreg(m1_E2_4, "nMonth", by = "fYear", type = "conditional", scale = c("response"))
model.temporal.fit<-v1$fit
saveRDS(model.temporal.fit, "model.temporal.fit.rds")


##############################################################################################################


##############Question 2: Models with stage##################################################################
library(MuMIn)

###Using avrg/std flows flows###
daily.stage<-read.csv("stage_2011_2018.csv")

daily.stage$Date2 <- daily.stage$Date
daily.stage<-cSplit(daily.stage, "Date2", sep = "/", type.convert = FALSE)
colnames(daily.stage)[6:8]<-c("Month", "Day", "Year")

daily.stage$Date <- as.Date(daily.stage$Date, format = "%m/%d/%Y")

saveRDS(daily.stage, "daily.stage.rds")
daily.stage<-readRDS("daily.stage.rds")


daily.stage.l<-gather(daily.stage, "Site", "Stage_cm", 2:5)
daily.stage.l$Stage_cm<-as.numeric(daily.stage.l$Stage_cm)

month_stage<-ddply(daily.stage.l, c("Site","Year", "Month"), summarise,
                   N = sum(!is.na(Date)),
                   mean.stage = mean(Stage_cm, na.rm = TRUE),
                   sd.stage  = sd(Stage_cm, na.rm = TRUE),
                   se.stage = sd.stage/sqrt(N))


library(DataCombine)
month_stage<-slide(month_stage, Var = "mean.stage", slideBy = -1)
month_stage<-slide(month_stage, Var = "sd.stage", slideBy = -1)
month_stage<-slide(month_stage, Var = "mean.stage", slideBy = -2)
month_stage<-slide(month_stage, Var = "sd.stage", slideBy = -2)

colnames(month_stage)[8:11]<-c("mean.stage.l1", "sd.stage.l1", "mean.stage.l2", "sd.stage.l2")
month_stage$Year.Month <- paste(month_stage$Year, month_stage$Month, sep = "/")
month_stage$fYear.Month <- as.yearmon(month_stage$Year.Month, "%Y/%m")

saveRDS(month_stage, "month_stage.rds")

month_stage <- readRDS("month_stage.rds")

###
#Merge Stage and Eadj values
###

stage_Eadj2<-merge(subset(month_stage, Year > 2011), E_year.month2_df, by =  "fYear.Month", all.x = TRUE)
stage_Eadj2$Site2<-factor(stage_Eadj2$Site)
levels(stage_Eadj2$Site2)<-list("Upstream" = c("BO_Stage"), "Upstream-Bay" = c("CN_Stage"), "Bay" = c("TE_Stage"), "Downstream" = c("GI_Stage"))
stage_Eadj2$Site.Season <- paste(stage_Eadj2$Site2, stage_Eadj2$Season, sep = ".")
stage_Eadj2$Site.Season<-factor(stage_Eadj2$Site.Season)
stage_Eadj2<-subset(stage_Eadj2, Eadj != "NA")

saveRDS(stage_Eadj2, "stage_Eadj2.rds")

####GAM - E as function of stage in Site|Season
lmc<-lmeControl(niterEM = 5000, msMaxIter = 10000, opt = "optim")


###
#Select the best distribution
###
mstage_E2_1<-gam(Eadj ~ s(mean.stage, bs = "cr", k = 5) + s(sd.stage, bs = "cr", k = 5) + 
                 s(mean.stage.l1, bs = "cr", k = 5) + s(sd.stage.l1, bs = "cr", k = 5) +
                 s(mean.stage.l2, bs = "cr", k = 5) + s(sd.stage.l2, bs = "cr", k = 5), 
               gamma = 1.4, data = stage_Eadj2, select = TRUE, family=gaussian(link = "log"))

mstage_E2_2<-gam(Eadj ~ s(mean.stage, bs = "cr", k = 5) + s(sd.stage, bs = "cr", k = 5) + 
                 s(mean.stage.l1, bs = "cr", k = 5) + s(sd.stage.l1, bs = "cr", k = 5) +
                 s(mean.stage.l2, bs = "cr", k = 5) + s(sd.stage.l2, bs = "cr", k = 5), 
               gamma = 1.4, data = stage_Eadj2, select = TRUE, family=Gamma(link = "log"))

mstage_E2_3<-gam(Eadj ~ s(mean.stage, bs = "cr", k = 5) + s(sd.stage, bs = "cr", k = 5) + 
                 s(mean.stage.l1, bs = "cr", k = 5) + s(sd.stage.l1, bs = "cr", k = 5) +
                 s(mean.stage.l2, bs = "cr", k = 5) + s(sd.stage.l2, bs = "cr", k = 5), 
               gamma = 1.4, data = stage_Eadj2, select = TRUE, family=nb(link = "log"))

AIC(mstage_E2_1, mstage_E2_2, mstage_E2_3)
gams_stage1 <- list(mstage_E2_1, mstage_E2_2, mstage_E2_3)
model.sel(gams_stage1)
#Model with gaussian(link = log) - best model

###
#Best model modification/simplification
###
mstage2_E2_1<-gam(Eadj ~ s(mean.stage, bs = "cr", k = 5, by = Site.Season) + s(sd.stage, bs = "cr", k = 5, by = Site.Season) + 
                   s(sd.stage.l1, bs = "cr", k = 5, by = Site.Season) + s(sd.stage.l2, bs = "cr", k = 5, by = Site.Season), 
                 gamma = 1.4, data = stage_Eadj2, select = TRUE, family=gaussian(link = "log"))

mstage2_E2_2<-gam(Eadj ~ s(mean.stage, bs = "cr", k = 5, by = Site.Season) + s(sd.stage, bs = "cr", k = 5, by = Site.Season) + 
                    s(sd.stage.l1, bs = "cr", k = 5, by = Site.Season), 
                  gamma = 1.4, data = stage_Eadj2, select = TRUE, family=gaussian(link = "log"))

mstage2_E2_3<-gam(Eadj ~ s(mean.stage, bs = "cr", k = 5, by = Site.Season) + s(sd.stage, bs = "cr", k = 5, by = Site.Season), 
                  gamma = 1.4, data = stage_Eadj2, select = TRUE, family=gaussian(link = "log"))

mstage2_E2_4<-gam(Eadj ~ s(mean.stage, bs = "cr", k = 5, by = Site.Season), 
                  gamma = 1.4, data = stage_Eadj2, select = TRUE, family=gaussian(link = "log"))

AIC(mstage2_E2_1, mstage2_E2_2, mstage2_E2_3, mstage2_E2_4)
gams_stage2 <- list(mstage2_E2_1, mstage2_E2_2, mstage2_E2_3, mstage2_E2_4)
model.sel(gams_stage2)

#best model is the one with the lags.sd (mstage2_E2_1), however, I think it is over kill adding 32 variables
#Also the stronger/signficant smoothers concentrated in the first two variables
#Therefore, I think is best just to make sense of mstage2_E2_3

mstage3_E2_3<-gam(Eadj ~ s(mean.stage, bs = "cr", k = 5, by = Season) + s(sd.stage, bs = "cr", k = 5, by = Season), 
                  gamma = 1.4, data = stage_Eadj2, select = TRUE, family=gaussian(link = "log"))


visreg_mean.stage<-visreg(mstage2_E2_3, "mean.stage", by = "Site.Season", type = "conditional", scale = c("response"))
visreg_sd.stage<-visreg(mstage2_E2_3, "sd.stage", by = "Site.Season", type = "conditional", scale = c("response"))

model.mstage.fit<-visreg_mean.stage$fit
model.sdstage.fit<-visreg_sd.stage$fit

saveRDS(model.mstage.fit, "model.mstage.fit.rds")
saveRDS(model.sdstage.fit, "model.sdstage.fit.rds")

###
#After inspecting the visreg plots, I decided to cut the mean/sd stage values to manifest 
#the proper range of values for the predictions - i.e., limit extrapolation at the edges

min.max.stage<-group_by(stage_Eadj2, Site.Season)%>%
  summarise(min.mean = min(mean.stage),
            max.mean = max(mean.stage),
            min.sd = min(sd.stage),
            max.sd = max(sd.stage))

model.mstage.fit2<-merge(model.mstage.fit, min.max.stage, by = "Site.Season", all.x = TRUE)%>%
  group_by(., Site.Season)%>%
  filter(mean.stage >= min.mean & mean.stage <= max.mean)

model.sdstage.fit2<-merge(model.sdstage.fit, min.max.stage, by = "Site.Season", all.x = TRUE)%>%
  group_by(., Site.Season)%>%
  filter(sd.stage >= min.sd & sd.stage <= max.sd)


# tmp1<- expand.grid(mean.flow = c(0:38), Season = levels(flows_Eadj2$Season), 
#                    sd.flow=c(2.5))
# 
# pred.E<-as.data.frame(predict(mflow1_e2d, se = TRUE, type = "response", newdata = tmp1))
# pred.E<-transform(pred.E, up = fit + 1.96 * se.fit, low = fit - 1.96 * se.fit)
# pred.E<-cbind(tmp1, pred.E)

#dry season
pred.E_dry<-subset(pred.E, Season %in% c("Dry"))%>%
  subset(mean.flow < 26)

colorx.1<-c("red")
pred.dry<-ggplot(data = pred.E_dry, aes(mean.flow, fit, colour = "")) + geom_line(size = 2) + 
  theme_bw()+ scale_colour_manual(name = "Dry", values = colorx.1) + 
  geom_ribbon(aes(ymax = up, ymin = low), fill = colorx.1, alpha = 0.3) + 
  scale_fill_manual(name = "Dry", values = colorx.1) + scale_y_continuous(limits = c(0.4, 0.9)) +
  labs(x = "Mean Flow", y = "E-adj") + #geom_hline(yintercept = 3.18, colour = "red", linetype = 2) +
  theme(axis.text = element_text(size = 20, face = "bold", colour = "black"), axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"))

#Wet
pred.E_wet<-subset(pred.E, Season %in% c("Wet"))%>%
  subset(mean.flow < 26)

colorx.2<-c("blue")
pred.wet<-ggplot(data = pred.E_wet, aes(mean.flow, fit, colour = "")) + geom_line(size = 2) + 
  theme_bw()+ scale_colour_manual(name = "Wet", values = colorx.2) + 
  geom_ribbon(aes(ymax = up, ymin = low), fill = colorx.2, alpha = 0.3) + 
  scale_fill_manual(name = "Wet", values = colorx.2) + scale_y_continuous(limits = c(0.4, 0.9)) +
  labs(x = "Mean Flow", y = "E-adj") + #geom_hline(yintercept = 3.18, colour = "red", linetype = 2) +
  theme(axis.text = element_text(size = 20, face = "bold", colour = "black"), axis.title = element_text(size = 20, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"))

library(cowplot)
plot_grid(pred.dry, pred.wet, nrow = 1, align = "h")


########################################################################
#After meeting with Jenn 9/6/2019, we decided to use PCA axis for the stage models
#Variable reduction and axis representing different aspect of stage

###
#Modify month_stage to do PCA
###

month_stage2 <- subset(month_stage, Year > 2011)

mean.stage.l <- month_stage2[,c(1:3,5)]%>%
  spread(., Site, mean.stage)
colnames(mean.stage.l)[3:6]<-c("Up.mean", "UpBay.mean", "Down.mean", "Bay.mean")

sd.stage.l <- month_stage2[,c(1:3,6)]%>%
  spread(., Site, sd.stage)
colnames(sd.stage.l)[3:6]<-c("Up.sd", "UpBay.sd", "Down.sd", "Bay.sd")

mean.stage.l2 <- month_stage2[,c(1:3,8)]%>%
  spread(., Site, mean.stage.l1)
colnames(mean.stage.l2)[3:6]<-c("Up.mean.l1", "UpBay.mean.l1", "Down.mean.l1", "Bay.mean.l1")

sd.stage.l2 <- month_stage2[,c(1:3,9)]%>%
  spread(., Site, sd.stage.l1)
colnames(sd.stage.l2)[3:6]<-c("Up.sd.l1", "UpBay.sd.l1", "Down.sd.l1", "Bay.sd.l1")

mean.stage3<-mean.stage.l%>%
  left_join(sd.stage.l, by = c("Year", "Month"))%>%
  left_join(mean.stage.l2, by = c("Year", "Month"))%>%
  left_join(sd.stage.l2, by = c("Year", "Month"))

write.csv(mean.stage3, "mean.stage.pca.csv")

###
#Bringing PCA results from primer using std/normalize values
###

pca.load <-read.csv("stage.pca.loading.csv")
pca.load.l <-pca.load[,1:5]%>%
  gather(., "PCA", "Score", 2:5)

pca.load.l$Variable<-factor(pca.load.l$Variable)
levels(pca.load.l$Variable)<-list("Up.mean" = c("Up.mean"), "UpBay.mean" = c("UpBay.mean"), "Bay.mean" = c("Bay.mean"), "Down.mean" = c("Down.mean"),
                                "Up.sd" = c("Up.sd"), "UpBay.sd" = c("UpBay.sd"), "Bay.sd" = c("Bay.sd"), "Down.sd" = c("Down.sd"),
                                "Up.mean.l1" = c("Up.mean.l1"), "UpBay.mean.l1" = c("UpBay.mean.l1"), "Bay.mean.l1" = c("Bay.mean.l1"), "Down.mean.l1" = c("Down.mean.l1"),
                                "Up.sd.l1" = c("Up.sd.l1"), "UpBay.sd.l1" = c("UpBay.sd.l1"), "Bay.sd.l1" = c("Bay.sd.l1"), "Down.sd.l1" = c("Down.sd.l1"))

#fix this later - something wrong with geom_hline, cutting the bars in half
ggplot(pca.load.l, aes(Variable, Score))+
  geom_bar(stat = "identity", position="dodge")+
  # geom_hline(aes(yintercept = 3), linetype = 2, size = 1, colour = "red")+
  # geom_hline(aes(yintercept = -3), linetype = 2, size = 1, colour = "red")+
  facet_wrap(~PCA)+
  labs(title = "PCA result - River Stage", x = "Stage Variable", 
       y = "Coefficient")+ 
  theme(axis.text = element_text(size = 14, face = "bold", colour = 'black'), 
        axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold")) 
  # coord_cartesian(ylim = c(-4, 4))
 


pca.score <-read.csv("stage.pca.score.csv")


pca_Eadj2<-merge(pca.score, E_year.month2_df, by.x =  c("Year", "Month"), by.y = c("Year", "nMonth"), all.x = TRUE)
pca_Eadj2<-subset(pca_Eadj2, Eadj != "NA")

saveRDS(pca_Eadj2, "pca_Eadj2.rds")

####GAM - E as function of PCA scores
lmc<-lmeControl(niterEM = 5000, msMaxIter = 10000, opt = "optim")


###
#Select the best distribution
###
library(MuMIn)

pca_Eadj2 <- readRDS("pca_Eadj2.rds")

#Old models after last meeting 07/2020------
# mstage_E2_1<-gam(Eadj ~ s(PC1, bs = "cr", k = 5) + s(PC2, bs = "cr", k = 5) + 
#                    s(PC3, bs = "cr", k = 5) + s(PC4, bs = "cr", k = 5),
#                    gamma = 1.4, data = pca_Eadj2, select = TRUE, family=gaussian(link = "log"))
# 
# mstage_E2_2<-gam(Eadj ~ s(PC1, bs = "cr", k = 5) + s(PC2, bs = "cr", k = 5) + 
#                    s(PC3, bs = "cr", k = 5) + s(PC4, bs = "cr", k = 5),
#                  gamma = 1.4, data = pca_Eadj2, select = TRUE, family=Gamma(link = "log"))
# 
# mstage_E2_3<-gam(Eadj ~ s(PC1, bs = "cr", k = 5) + s(PC2, bs = "cr", k = 5) + 
#                    s(PC3, bs = "cr", k = 5) + s(PC4, bs = "cr", k = 5),
#                  gamma = 1.4, data = pca_Eadj2, select = TRUE, family=nb(link = "log"))
# 
# AIC(mstage_E2_1, mstage_E2_2, mstage_E2_3)
# gams_stage1 <- list(mstage_E2_1, mstage_E2_2, mstage_E2_3)
# model.sel.q2 <- model.sel(gams_stage1)
# capture.output(model.sel.q2, file = "model.sel.q2.txt")

#Model with gaussian(link = log) - best model

###
#Best model modification/simplification
###
# mstage2_E2_1<-gam(Eadj ~ s(PC1, bs = "cr", k = 5) + s(PC2, bs = "cr", k = 5) + 
#                     s(PC3, bs = "cr", k = 5) + s(PC4, bs = "cr", k = 5),
#                   gamma = 1.4, data = pca_Eadj2, select = TRUE, family=gaussian(link = "log"))
# 
# mstage2_E2_2<-gam(Eadj ~ s(PC1, bs = "cr", k = 5) + s(PC2, bs = "cr", k = 5) + 
#                     s(PC3, bs = "cr", k = 5),
#                   gamma = 1.4, data = pca_Eadj2, select = TRUE, family=gaussian(link = "log"))
# 
# mstage2_E2_3<-gam(Eadj ~ s(PC1, bs = "cr", k = 5, by = Season) + s(PC3, bs = "cr", k = 5, by = Season),
#                   gamma = 1.4, data = pca_Eadj2, select = TRUE, family=gaussian(link = "log"))


# AIC(mstage2_E2_1, mstage2_E2_2, mstage2_E2_3)
# gams_stage2 <- list(mstage2_E2_1, mstage2_E2_2, mstage2_E2_3)
# model.sel.q2b <- model.sel(gams_stage2)
# capture.output(model.sel.q2b, file = "model.sel.q2b.txt")


#Model update based on Jenn's suggestion - email 08/08/2020-------
mstage2_E2_PC1 <- gam(Eadj ~ s(PC1, bs = "cr", k = 5, by = Season),
                        gamma = 1.4, data = pca_Eadj2, select = TRUE, family=gaussian(link = "log"))

mstage2_E2_PC2 <- gam(Eadj ~ s(PC2, bs = "cr", k = 5, by = Season),
                      gamma = 1.4, data = pca_Eadj2, select = TRUE, family=gaussian(link = "log"))

mstage2_E2_PC3 <- gam(Eadj ~ s(PC3, bs = "cr", k = 5, by = Season),
                      gamma = 1.4, data = pca_Eadj2, select = TRUE, family=gaussian(link = "log"))

mstage2_E2_PC1n2 <- gam(Eadj ~ s(PC1, bs = "cr", k = 5, by = Season) + s(PC2, bs = "cr", k = 5, by = Season),
                        gamma = 1.4, data = pca_Eadj2, select = TRUE, family=gaussian(link = "log"))


mstage2_E2_PC1n3 <- gam(Eadj ~ s(PC1, bs = "cr", k = 5, by = Season) + s(PC3, bs = "cr", k = 5, by = Season),
                  gamma = 1.4, data = pca_Eadj2, select = TRUE, family=gaussian(link = "log"))

mstage2_E2_PC4 <- gam(Eadj ~ s(PC4, bs = "cr", k = 5, by = Season),
                      gamma = 1.4, data = pca_Eadj2, select = TRUE, family=gaussian(link = "log"))


AIC(mstage2_E2_PC1, mstage2_E2_PC2, mstage2_E2_PC3, mstage2_E2_PC1n2, mstage2_E2_PC1n3, mstage2_E2_PC4)
gams_stage2 <- list(mstage2_E2_PC1, mstage2_E2_PC2, mstage2_E2_PC3, mstage2_E2_PC1n2, mstage2_E2_PC1n3, mstage2_E2_PC4)
model.sel.q2b <- model.sel(gams_stage2)
capture.output(model.sel.q2b, file = "model.sel.q2_08082020.txt")

#Best model - PC1 and PC3

# par(mfrow = c(2, 1))

visreg_pca1<-visreg(mstage2_E2_PC1n3, "PC1", by = "Season",type = "conditional", scale = c("response"))
visreg_pca3<-visreg(mstage2_E2_PC1n3, "PC3", by = "Season",type = "conditional", scale = c("response"))

model.pca1.fit<-visreg_pca1$fit
model.pca3.fit<-visreg_pca3$fit

saveRDS(model.pca1.fit, "model.pca1.fit.rds")
saveRDS(model.pca3.fit, "model.pca3.fit.rds")

summary.model.stage<-summary(mstage2_E2_PC1n3)
capture.output(summary.model.stage, file = "summary.model.stage.txt")

# AIC(mstage2_E2_1, mstage2_E2_2, mstage2_E2_3)
# gams_stage2 <- list(mstage2_E2_1, mstage2_E2_2, mstage2_E2_3)
# model.sel.q2b <- model.sel(gams_stage2)
# capture.output(model.sel.q2b, file = "model.sel.q2b.txt")

###
#After inspecting the visreg plots, I decided to cut the mean/sd stage values to manifest 
#the proper range of values for the predictions - i.e., limit extrapolation at the edges

min.max.pca<-group_by(pca_Eadj2, Season)%>%
  summarise(min.PC1 = min(PC1),
            max.PC1 = max(PC1),
            min.PC3 = min(PC3),
            max.PC3 = max(PC3))

model.pca1.fit2<-merge(model.pca1.fit, min.max.pca, by = "Season", all.x = TRUE)%>%
  group_by(., Season)%>%
  filter(PC1 >= min.PC1 & PC1 <= max.PC3)

model.pca3.fit2<-merge(model.pca3.fit, min.max.pca, by = "Season", all.x = TRUE)%>%
  group_by(., Season)%>%
  filter(PC3 >= min.PC3 & PC3 <= max.PC3)






################################################################################


##########################################################################
#############Question 3: Temporal Variability in Spatial Use############

snook_zone_use<-group_by(snook.90day, ID, Month, Year, Zone)%>% 
                summarise(freq.v = length(unique(fDay)))
snook_zone_use<-subset(snook_zone_use, Zone != "NA")

snook_zone_use.w<-dcast(snook_zone_use, ID + Month + Year ~ Zone, sum, value.var = "freq.v")
snook_zone_use.w$Total<-rowSums(snook_zone_use.w[,4:15])

snook_zone_use.l<-gather(snook_zone_use.w, key = Zone, value = Num.Days, c(-ID, -Month, -Year, -Total))%>%
  transform(., Prop.Use = Num.Days/Total)

snook_zone_use.l$fMonth<-factor(snook_zone_use.l$Month)
levels(snook_zone_use.l$fMonth) <- list("Jan" = c("01"), "Feb" = c("02"), "Mar" = c("03"), 
                                        "Apr" = c("04"), "May" = c("05"), "Jun" = c("06"), 
                                        "Jul" = c("07"),"Aug" = c("08"), "Sep" = c("09"), 
                                        "Oct" = c("10"), "Nov" = c("11"), "Dec" = c("12"))

snook_zone_use.l$Season<-factor(snook_zone_use.l$fMonth)
levels(snook_zone_use.l$Season)<-list("Dry" = c("Jan", "Feb", "Mar", "Apr", "May", "Dec"), 
                                      "Wet" = c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))


saveRDS(snook_zone_use.l, "snook_zone_use.l.rds")
snook_zone_use.l <- readRDS("snook_zone_use.l.rds")

snook_zuse_avrg <- group_by(snook_zone_use.l, Year, Season, Zone)%>%
  summarise(N = sum(!is.na(Year)),
            mean.prop = mean(Prop.Use),
            sd.prop  = sd(Prop.Use),
            se.prop = sd.prop/sqrt(N))

snook_zuse_avrg$fZone<-factor(snook_zuse_avrg$Zone)
levels(snook_zuse_avrg$fZone) <- list("1" = "1", "2" = "2", "3" = "3", "4" = "4", "5" = "5",
                                      "6" = "6", "7" = "7", "8" = "8", "9" = "9", "10" = "10",
                                      "11"= "11", "12" = "12")

snook_zuse_avrg<-filter(snook_zuse_avrg, Year < 2019)

library(viridis)
ggplot(snook_zuse_avrg, aes(Year, fZone))+
  geom_tile(alpha = 0.8, aes(fill = mean.prop))+
  theme_bw()+	
  facet_wrap(~Season)+
  scale_y_discrete(limits = rev(levels(snook_zuse_avrg$fZone))) +
  # scale_fill_gradientn(colours = c("red", "green", "blue4"),				
                       # limits = c(0,0.60), breaks = seq(0,0.60,by=0.10), name = "Proportion")+				
  scale_fill_viridis(option = "B", limits = c(0,0.60), breaks = seq(0,0.60,by=0.10), name = "Proportion") +
  labs(title = "Snook - Proportion of Days Detected", x = "Year", y = "Zone")+				
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),				
        axis.title = element_text(size = 16, face = "bold"),				
        plot.title = element_text(size = 16, face = "bold"))				
ggsave("PropDetection_DryWet.png", width = 20, height = 18, units = "cm")


ggplot(snook_zuse_avrg, aes(fZone, mean.prop, colour = Year, group=Year))+
  geom_point(size = 2)+
  geom_line()+
  theme_bw()+	
  facet_wrap(~Season)+
  # scale_y_discrete(limits = rev(levels(snook_zuse_avrg$fZone))) +
  # scale_fill_gradientn(colours = c("red", "green", "blue4"),				
  # limits = c(0,0.60), breaks = seq(0,0.60,by=0.10), name = "Proportion")+				
  scale_color_viridis(option = "B") +
  labs(title = "Snook - Proportion of Days Detected", x = "Zone", y = "Proportion of Days")+				
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),				
        axis.title = element_text(size = 16, face = "bold"),				
        plot.title = element_text(size = 16, face = "bold"))				
ggsave("PropDetection_DryWet_linear.png", width = 20, height = 18, units = "cm")


snook_zuse_avrg.w<-dcast(snook_zuse_avrg, Year + fZone ~ Season, sum, value.var = "mean.prop")
# snook_zuse_avrg.w<-transform(snook_zuse_avrg.w, prop.r.delta = (Wet - Dry)/Dry)
snook_zuse_avrg.w<-transform(snook_zuse_avrg.w, prop.r.delta = (Wet - Dry)/((Wet + Dry)/2), 
                             prop.r.delta2 = log((Wet+1)/(Dry+1)),
                             raw.diff = Wet - Dry)

snook_zuse_avrg.w$prop.r.delta[is.na(snook_zuse_avrg.w$prop.r.delta)] <- 0
# snook_zuse_avrg.w$prop.r.delta2[is.na(snook_zuse_avrg.w$prop.r.delta2)] <- 0

snook_zuse_avrg.w$change.type<-ifelse(snook_zuse_avrg.w$prop.r.delta < 0, "negative", "positive")
snook_zuse_avrg.w$change.type<-ifelse(snook_zuse_avrg.w$prop.r.delta == 0, "no change",snook_zuse_avrg.w$change.type)
snook_zuse_avrg.w<-transform(snook_zuse_avrg.w, prop.r.delta3 = abs(prop.r.delta))


# snook_zuse_avrg2<-gather(snook_zuse_avrg.w, key = Season, value = mean.prop, c(-Year, -fZone, -prop.r.delta))
  

ggplot(snook_zuse_avrg.w, aes(Year, fZone))+
  geom_tile(alpha = 0.8, aes(fill = prop.r.delta))+
  theme_bw()+	
  # facet_wrap(~Season)+
  scale_y_discrete(limits = rev(levels(snook_zuse_avrg$fZone))) +
  # scale_fill_gradientn(colours = c("red", "green", "blue4"),				
  # limits = c(0,0.60), breaks = seq(0,0.60,by=0.10), name = "Proportion")+				
  scale_fill_viridis(option = "D", limits = c(-2,2), breaks = seq(-2,2,by=0.5), name = "d.Proportion") +
  labs(title = "Snook - Relative Change Dry2Wet", x = "Year", y = "Zone")+				
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),				
        axis.title = element_text(size = 16, face = "bold"),				
        plot.title = element_text(size = 16, face = "bold"))				
ggsave("RelativeChange_WetvsDry_1.png", width = 20, height = 20, units = "cm")

ggplot(snook_zuse_avrg.w, aes(Year, fZone))+
  geom_tile(alpha = 0.8, aes(fill = prop.r.delta2))+
  theme_bw()+	
  # facet_wrap(~Season)+
  scale_y_discrete(limits = rev(levels(snook_zuse_avrg$fZone))) +
  # scale_fill_gradientn(colours = c("red", "green", "blue4"),				
  # limits = c(0,0.60), breaks = seq(0,0.60,by=0.10), name = "Proportion")+				
  scale_fill_viridis(option = "D", limits = c(-0.2,0.1), breaks = seq(-0.2,0.1,by=0.1), name = "d.Proportion") +
  labs(title = "Snook - Relative Change Dry2Wet", x = "Year", y = "Zone")+				
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),				
        axis.title = element_text(size = 16, face = "bold"),				
        plot.title = element_text(size = 16, face = "bold"))				
ggsave("RelativeChange_WetvsDry_2.png", width = 20, height = 20, units = "cm")

ggplot(snook_zuse_avrg.w, aes(Year, fZone))+
  geom_tile(alpha = 0.8, aes(fill = raw.diff))+
  theme_bw()+	
  # facet_wrap(~Season)+
  scale_y_discrete(limits = rev(levels(snook_zuse_avrg$fZone))) +
  # scale_fill_gradientn(colours = c("red", "green", "blue4"),				
  # limits = c(0,0.60), breaks = seq(0,0.60,by=0.10), name = "Proportion")+				
  scale_fill_viridis(option = "D", limits = c(-0.2,0.1), breaks = seq(-0.2,0.1,by=0.05), name = "d.Proportion") +
  labs(title = "Wet - Dry", x = "Year", y = "Zone")+				
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),				
        axis.title = element_text(size = 16, face = "bold"),				
        plot.title = element_text(size = 16, face = "bold"))				
ggsave("RelativeChange_WetvsDry_3.png", width = 20, height = 20, units = "cm")


ggplot(snook_zuse_avrg.w, aes(Year, prop.r.delta, colour = fZone))+
  geom_point(size = 2)+
  geom_line()+
  theme_bw()+	
  # facet_wrap(~Season)+
  # scale_y_discrete(limits = rev(levels(snook_zuse_avrg$fZone))) +
  # scale_fill_gradientn(colours = c("red", "green", "blue4"),				
  # limits = c(0,0.60), breaks = seq(0,0.60,by=0.10), name = "Proportion")+				
  scale_color_viridis(option = "B", discrete = TRUE) +
  labs(title = "Snook - Proportion of Days Detected", x = "Year", y = "delta.Prop")+				
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),				
        axis.title = element_text(size = 16, face = "bold"),				
        plot.title = element_text(size = 16, face = "bold"))				


library(ggforce)
ggplot(snook_zuse_avrg.w, aes(Year, prop.r.delta, colour = fZone))+
  geom_point(size = 2)+
  geom_line()+
  theme_bw()+	
  facet_zoom(ylim = c(-1, 10), zoom.data = ifelse(prop.r.delta <= 10, NA, FALSE))+
  scale_color_viridis(option = "B", discrete = TRUE) +
  labs(title = "Snook - Proportion of Days Detected", x = "Year", y = "delta.Prop")+				
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),				
        axis.title = element_text(size = 16, face = "bold"),				
        plot.title = element_text(size = 16, face = "bold"))				

type.change.s <-group_by(snook_zuse_avrg.w, fZone, change.type)%>%
  summarise(count.n = sum(!is.na(prop.r.delta3)),
            mean.prop.delta = mean(prop.r.delta3),
            sd.prop.delta = sd(prop.r.delta3))

ggplot(type.change.s, aes(change.type, fZone, colour = change.type))+
  geom_point(aes(size = count.n))+
  # geom_line()+
  theme_bw()+	
  scale_y_discrete(limits = rev(levels(snook_zuse_avrg$fZone))) +
  scale_size_continuous(range = c(2, 10))+
  scale_color_viridis(option = "D", discrete = TRUE) +
  labs(title = "Snook - Type of Change between Dry-Wet", x = "Type of Change", y = "Zone")+				
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),				
        axis.title = element_text(size = 16, face = "bold"),				
        plot.title = element_text(size = 16, face = "bold"))				
# c("lightskyblue", "blue", "blue4")
ggsave("TypeChange_Count.png", width = 20, height = 20, units = "cm")

ggplot(type.change.s, aes(change.type, fZone, colour = change.type))+
  geom_point(aes(size = mean.prop.delta))+
  # geom_line()+
  theme_bw()+	
  scale_y_discrete(limits = rev(levels(snook_zuse_avrg$fZone))) +
  scale_size_continuous(range = c(2, 10))+
  scale_color_viridis(option = "D", discrete = TRUE) +
  labs(title = "Snook - Type of Change between Dry-Wet", x = "Type of Change", y = "Zone")+				
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"),				
        axis.title = element_text(size = 16, face = "bold"),				
        plot.title = element_text(size = 16, face = "bold"))				
# c("lightskyblue", "blue", "blue4")
ggsave("TypeChange_AvrgAmount.png", width = 20, height = 20, units = "cm")


##########################################################################

##########Cluster Analysis to create heatmap with dendrograms##############################################################
library(vegan)
library(dendextend)
library(heatmaply)

###
#Converting df for cluster analysis
###
seasonal.delta <- snook_zuse_avrg.w[,c(1:2, 5)]%>%
  spread(., Year, prop.r.delta)

#Changing zone as rownames
row.names(seasonal.delta) <-seasonal.delta[,1]
seasonal.delta[,1]<-NULL

# seasonal.delta <- scale(seasonal.delta)

# #Change to as.matrix
# seasonal.delta <- as.matrix(seasonal.delta, rownames = "fZone")

###
#Cluster analysis
##

##For Zones
c.clus.ward <- hclust(dist(seasonal.delta), method = "ward.D")
c.clus.UPGMA <- hclust(dist(seasonal.delta), method = "average")

##For Year
r.clus.ward <- hclust(dist(t(seasonal.delta)), method = "ward.D")
r.clus.UPGMA <- hclust(dist(t(seasonal.delta)), method = "average")

###
#Cophenetic corr
###

#For Zone
cor.c.clus.ward <- cor(dist(seasonal.delta), cophenetic(c.clus.ward))
cor.c.clus.UPGMA <- cor(dist(seasonal.delta), cophenetic(c.clus.UPGMA))

gow.c.clus.ward <- sum((dist(seasonal.delta) - cor.c.clus.ward)^2)
gow.c.clus.UPGMA <- sum((dist(t(seasonal.delta)) - cor.c.clus.UPGMA)^2)

#For Year
cor.r.clus.ward <- cor(dist(t(seasonal.delta)), cophenetic(r.clus.ward))
cor.r.clus.UPGMA <- cor(dist(t(seasonal.delta)), cophenetic(r.clus.UPGMA))

gow.r.clus.ward <- sum((dist(seasonal.delta) - cor.r.clus.ward)^2)
gow.r.clus.UPGMA <- sum((dist(t(seasonal.delta)) - cor.r.clus.UPGMA)^2)

#Decided to use clustering with euclidean distance and average method
#Based on the cor estimates suggested by Numerical Ecology Analysis book

###
#Plot dendrogram
###

#By Zone
plot(c.clus.UPGMA)

#By Year
plot(r.clus.UPGMA)

###
#ID optimal number of groups for cluster analysis
###

#Fusion plot for Zone
plot(c.clus.UPGMA$height, nrow(seasonal.delta) :2, type = "S",
     main = "Fusion Levels - Chord - UPGMA", ylab = "k (number of clusters)",
    xlab = "h (node height)", col = "grey")
text(c.clus.UPGMA$height, nrow(seasonal.delta) :2, col = "red", cex = 0.8)

plot(r.clus.UPGMA$height, nrow(t(seasonal.delta)) :2, type = "S",
     main = "Fusion Levels - Chord - UPGMA", ylab = "k (number of clusters)",
     xlab = "h (node height)", col = "grey")
text(r.clus.UPGMA$height, nrow(t(seasonal.delta)) :2, col = "red", cex = 0.8)


# Optimal number of clusters according to matrix correlation 
# statistic (Pearson)
#source("grpdist) #code is my rnote folder

hc1 <- c.clus.UPGMA
hc2 <- r.clus.UPGMA
spe <- seasonal.delta

kt1 <- data.frame(k = 1:nrow(spe), r = 0)
kt2 <- data.frame(k = 1:nrow(t(spe)), r = 0)

#For Zone
for (i in 2:(nrow(spe) - 1)) 
{
  gr <- cutree(hc1, i)
  distgr <- grpdist(gr)
  mt <- cor(dist(spe), distgr, method = "pearson")
  kt1[i, 2] <- mt
}
k.best <- which.max(kt1$r)
plot(
  kt1$k,
  kt1$r,
  type = "h",
  main = "Zone - Matrix correlation-optimal number of clusters",
  xlab = "k (number of clusters)",
  ylab = "Pearson's correlation"
)
axis(
  1,
  k.best,
  paste("optimum", k.best, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(k.best,
       max(kt1$r),
       pch = 16,
       col = "red",
       cex = 1.5)



###
#For Year
for (i in 2:(nrow(t(spe)) - 1)) 
{
  gr <- cutree(hc2, i)
  distgr <- grpdist(gr)
  mt <- cor(dist(t(spe)), distgr, method = "pearson")
  kt2[i, 2] <- mt
}
k.best <- which.max(kt2$r)
plot(
  kt2$k,
  kt2$r,
  type = "h",
  main = "Year - Matrix correlation-optimal number of clusters",
  xlab = "k (number of clusters)",
  ylab = "Pearson's correlation"
)
axis(
  1,
  k.best,
  paste("optimum", k.best, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(k.best,
       max(kt2$r),
       pch = 16,
       col = "red",
       cex = 1.5)


rcx#Two major clusters for Years and 5 for Zones


###
#Creating dendrogram with optimal groups
###

#For Zone
k1 = 5

re.hc1 <- reorder(hc1, dist(seasonal.delta))

dend1 <- as.dendrogram(re.hc1)

# Plot the dendrogram with colored branches using the dendextend 
# syntax
dev.new(
  title = "Colored dendrogram",
  width = 8,
  height = 6,
  noRStudioGD = TRUE
)
dend1 %>% set("branches_k_color", k = k1) %>% plot

# Use standard colors for clusters
clusters <- cutree(dend1, k1)[order.dendrogram(dend1)]
dend1 %>% set("branches_k_color", 
             k = k1, value = unique(clusters) + 1) %>% plot
# Add a colored bar
colored_bars(clusters + 1,
             y_shift = -0.5,
             rowLabels = paste(k1, "clusters"))

#For Year
k2 = 2

re.hc2 <- reorder(hc2, dist(t(seasonal.delta)))

dend2 <- as.dendrogram(re.hc2)

# Plot the dendrogram with colored branches using the dendextend 
# syntax
dev.new(
  title = "Colored dendrogram",
  width = 8,
  height = 6,
  noRStudioGD = TRUE
)
dend2 %>% set("branches_k_color", k = k2) %>% plot

# Use standard colors for clusters
clusters <- cutree(dend2, k2)[order.dendrogram(dend2)]
dend2 %>% set("branches_k_color", 
             k = k2, value = unique(clusters) + 1) %>% plot
# Add a colored bar
colored_bars(clusters + 1,
             y_shift = -0.5,
             rowLabels = paste(k2, "clusters"))



#####
#Heatmap and dendrogram togehter
###
library(orca)

row_dend <- seasonal.delta %>% dist %>% hclust %>% as.dendrogram %>%
  dendextend::set("branches_k_color", k = 5) %>% dendextend::set("branches_lwd", c(1,5)) %>%
  ladderize

col_dend <- seasonal.delta %>% t %>% dist %>% hclust %>% as.dendrogram %>%
  dendextend::set("branches_k_color", k = 2) %>% dendextend::set("branches_lwd", c(1,5)) %>%
  ladderize

heatmaply(seasonal.delta, Rowv = row_dend, Colv = col_dend, fontsize_col = 14, fontsize_row = 14,
          xlab = "Year", ylab = "Zone")


heatmaply(cor(seasonal.delta), margins = c(40, 40),
          k_col = 2, k_row = 2,
          limits = c(-1,1))

heatmaply(cor(t(seasonal.delta)), margins = c(40, 40),
          k_col = 5, k_row = 5,
          limits = c(-1,1))

heatmaply(cor(t(seasonal.delta)), margins = c(40, 40),
          Rowv = row_dend, Colv = row_dend,
          limits = c(-1,1))

heatmaply(cor(seasonal.delta), margins = c(40, 40),
          Rowv = col_dend, Colv = col_dend,
          limits = c(-1,1))


###
#Add/Plot/Test clusters id above---------
###

###
#Add clusters
###
snook_zuse_avrg.w$Year.Cluster <- factor(snook_zuse_avrg.w$Year)
# levels(snook_zuse_avrg.w$Year.Cluster) <- list("A" = c(2016), "B" = c(2014, 2013), "C" = c(2017, 2015, 2018, 2012))
levels(snook_zuse_avrg.w$Year.Cluster) <- list("A" = c(2016), "B" = c(2014, 2013, 2017, 2015, 2018, 2012))


snook_zuse_avrg.w$Zone.Cluster <- factor(snook_zuse_avrg.w$fZone)
# levels(snook_zuse_avrg.w$Zone.Cluster) <- list("A" = c(1,2,3,8), "B" = c(6,9), "C" = c(5),
#                                                "D" = c(4, 7), "E" = c(10))
levels(snook_zuse_avrg.w$Zone.Cluster) <- list("A" = c(1,2,3,5,9,8,11), "B" = c(6), "C" = c(4,7),
                                               "D" = c(12), "E" = c(10))
###
###
#Plot clusters
###

clusYear_avrg <- group_by(snook_zuse_avrg.w, Year.Cluster)%>%
  summarise(N = sum(!is.na(fZone)),
            mn.delta = mean(prop.r.delta),
            sd.delta = sd(prop.r.delta),
            se.delta = sd.delta/sqrt(N))

clusZone_avrg <- group_by(snook_zuse_avrg.w, Zone.Cluster)%>%
  summarise(N = sum(!is.na(fZone)),
            mn.delta = mean(prop.r.delta),
            sd.delta = sd(prop.r.delta),
            se.delta = sd.delta/sqrt(N))

limits.1<-aes(ymax = mn.delta + se.delta, ymin = mn.delta - se.delta)
colorx.2<-c("#CC476B", "#009681")
ggplot(clusYear_avrg, aes(Year.Cluster, mn.delta, colour = Year.Cluster))+ 
  geom_point(size = 4) +
  geom_errorbar(limits.1, size = 2, width = 0.5)+
  theme_bw()+
  # scale_y_continuous(limits = c(0.5, 0.8), breaks=seq(0.5,0.8,0.1)) + 
  scale_colour_manual(name = "Year Clusters", values = colorx.2) +
  labs(x = "Year Clusters", y = "Dry-Wet Relative Change")+ 
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold"))
ggsave("aov_cluster_Years.png", width = 18, height = 15, units = "cm")

limits.1<-aes(ymax = mn.delta + se.delta, ymin = mn.delta - se.delta)
colorx.2<-c("#A352D1", "#CC4569","#917600", "#009232","#008FB7")
ggplot(clusZone_avrg, aes(Zone.Cluster, mn.delta, colour = Zone.Cluster))+ 
  geom_point(size = 4) +
  geom_errorbar(limits.1, size = 2, width = 0.5)+
  theme_bw()+
  # scale_y_continuous(limits = c(0.5, 0.8), breaks=seq(0.5,0.8,0.1)) + 
  scale_colour_manual(name = "Zone Clusters", values = colorx.2) +
  labs(x = "Zone Clusters", y = "Dry-Wet Relative Change")+ 
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title = element_text(size = 16, face = "bold"), 
        plot.title = element_text(size = 16, face = "bold"))
ggsave("aov_cluster_Zones.png", width = 18, height = 15, units = "cm")

###
#Test diff between clusters
###

aov.year <- aov(prop.r.delta ~ Year.Cluster, data = snook_zuse_avrg.w)
aov.year.out <-summary(aov.year)
capture.output(aov.year.out, file = "aov.year.out.txt")

aov.zone <- aov(prop.r.delta ~ Zone.Cluster, data = snook_zuse_avrg.w)
aov.zone.out <- summary(aov.zone)
capture.output(aov.zone.out, file = "aov.zone.out.txt")


type.change.s$Zone.Cluster <- factor(type.change.s$fZone)
levels(type.change.s$Zone.Cluster) <- list("A" = c(1,2,3,5,9,8,11), "B" = c(6), "C" = c(4,7),
                                           "D" = c(12), "E" = c(10))
###

ggplot(type.change.s, aes(change.type, fZone, colour = change.type))+
  geom_point(aes(size = count.n))+
  # geom_line()+
  theme_bw()+	
  facet_grid(Zone.Cluster~., scales = "free_y")+
  scale_y_discrete(limits = rev(levels(snook_zuse_avrg$fZone))) +
  scale_size_continuous(range = c(2, 10))+
  scale_color_viridis(option = "D", discrete = TRUE) +
  labs(title = "Snook - Type of Change between Dry-Wet", x = "Type of Change", y = "Zone")+				
  theme(axis.text.y = element_text(size = 8, face = "bold", colour = "black"),
        axis.text.x = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 16, face = "bold"),				
        plot.title = element_text(size = 16, face = "bold"))				



####################################################################################

##########################################################################
#############Temporal Variability in LISA############

library(maptools)
library(rgdal)
library(gstat)
library(classInt)
library(spdep)
library(sf)
library(sp)
library(tmap)

#read polygon using sf package and converting to sp object
shark.zones<-st_read("EditedRiverZones_ReprojectedUTM17NCopy.shp")

sz.12d <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2012 & Season %in% c("Dry")),
                by.x = "Region", by.y = "fZone")

sz.12w <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2012 & Season %in% c("Wet")),
                by.x = "Region", by.y = "fZone")

sz.13d <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2013 & Season %in% c("Dry")),
                by.x = "Region", by.y = "fZone")

sz.13w <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2013 & Season %in% c("Wet")),
                by.x = "Region", by.y = "fZone")

sz.14d <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2014 & Season %in% c("Dry")),
                by.x = "Region", by.y = "fZone")

sz.14w <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2014 & Season %in% c("Wet")),
                by.x = "Region", by.y = "fZone")

sz.15d <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2015 & Season %in% c("Dry")),
                by.x = "Region", by.y = "fZone")

sz.15w <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2015 & Season %in% c("Wet")),
                by.x = "Region", by.y = "fZone")

sz.16d <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2016 & Season %in% c("Dry")),
                by.x = "Region", by.y = "fZone")

sz.16w <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2016 & Season %in% c("Wet")),
                by.x = "Region", by.y = "fZone")

sz.17d <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2017 & Season %in% c("Dry")),
                by.x = "Region", by.y = "fZone")

sz.17w <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2017 & Season %in% c("Wet")),
                by.x = "Region", by.y = "fZone")

sz.18d <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2018 & Season %in% c("Dry")),
                by.x = "Region", by.y = "fZone")

sz.18w <- merge(shark.zones, filter(snook_zuse_avrg, Year == 2018 & Season %in% c("Wet")),
                by.x = "Region", by.y = "fZone")

shark.zones2<-merge(shark.zones, snook_zuse_avrg, by.x = "Region", by.y = "fZone")


shark.zones.sp <- as(shark.zones2, "Spatial")



###
#Neighbor connectivity and weight####
###

sz.nb.12d<-poly2nb(as(sz.12d, "Spatial"), queen=T, snap = 100)
sz.nb.12w<-poly2nb(as(sz.12w, "Spatial"), queen=T, snap = 100)
sz.nb.13d<-poly2nb(as(sz.13d, "Spatial"), queen=T, snap = 100)
sz.nb.13w<-poly2nb(as(sz.13w, "Spatial"), queen=T, snap = 100)
sz.nb.14d<-poly2nb(as(sz.14d, "Spatial"), queen=T, snap = 100)
sz.nb.14w<-poly2nb(as(sz.14w, "Spatial"), queen=T, snap = 100)
sz.nb.15d<-poly2nb(as(sz.15d, "Spatial"), queen=T, snap = 100)
sz.nb.15w<-poly2nb(as(sz.15w, "Spatial"), queen=T, snap = 100)
sz.nb.16d<-poly2nb(as(sz.16d, "Spatial"), queen=T, snap = 100)
sz.nb.16w<-poly2nb(as(sz.16w, "Spatial"), queen=T, snap = 100)
sz.nb.17d<-poly2nb(as(sz.17d, "Spatial"), queen=T, snap = 100)
sz.nb.17w<-poly2nb(as(sz.17w, "Spatial"), queen=T, snap = 100)
sz.nb.18d<-poly2nb(as(sz.18d, "Spatial"), queen=T, snap = 100)
sz.nb.18w<-poly2nb(as(sz.18w, "Spatial"), queen=T, snap = 100)

sz.nw.12d<-nb2listw(sz.nb.12d, style="W", zero.policy = TRUE)
sz.nw.12w<-nb2listw(sz.nb.12w, style="W", zero.policy = TRUE)
sz.nw.13d<-nb2listw(sz.nb.13d, style="W", zero.policy = TRUE)
sz.nw.13w<-nb2listw(sz.nb.13w, style="W", zero.policy = TRUE)
sz.nw.14d<-nb2listw(sz.nb.14d, style="W", zero.policy = TRUE)
sz.nw.14w<-nb2listw(sz.nb.14w, style="W", zero.policy = TRUE)
sz.nw.15d<-nb2listw(sz.nb.15d, style="W", zero.policy = TRUE)
sz.nw.15w<-nb2listw(sz.nb.15w, style="W", zero.policy = TRUE)
sz.nw.16d<-nb2listw(sz.nb.16d, style="W", zero.policy = TRUE)
sz.nw.16w<-nb2listw(sz.nb.16w, style="W", zero.policy = TRUE)
sz.nw.17d<-nb2listw(sz.nb.17d, style="W", zero.policy = TRUE)
sz.nw.17w<-nb2listw(sz.nb.17w, style="W", zero.policy = TRUE)
sz.nw.18d<-nb2listw(sz.nb.18d, style="W", zero.policy = TRUE)
sz.nw.18w<-nb2listw(sz.nb.18w, style="W", zero.policy = TRUE)

###
#Moran's I#######
###
sz.mi.12d<-moran.mc(sz.12d$mean.prop, sz.nw.12d, nsim=999) 
sz.mi.12w<-moran.mc(sz.12w$mean.prop, sz.nw.12w, nsim=999) 
sz.mi.13d<-moran.mc(sz.13d$mean.prop, sz.nw.13d, nsim=999) 
sz.mi.13w<-moran.mc(sz.13w$mean.prop, sz.nw.13w, nsim=999) 
sz.mi.14d<-moran.mc(sz.14d$mean.prop, sz.nw.14d, nsim=999) 
sz.mi.14w<-moran.mc(sz.14w$mean.prop, sz.nw.14w, nsim=999) 
sz.mi.15d<-moran.mc(sz.15d$mean.prop, sz.nw.15d, nsim=999) 
sz.mi.15w<-moran.mc(sz.15w$mean.prop, sz.nw.15w, nsim=999) 
sz.mi.16d<-moran.mc(sz.16d$mean.prop, sz.nw.16d, nsim=999) 
sz.mi.16w<-moran.mc(sz.16w$mean.prop, sz.nw.16w, nsim=999) 
sz.mi.17d<-moran.mc(sz.17d$mean.prop, sz.nw.17d, nsim=999) 
sz.mi.17w<-moran.mc(sz.17w$mean.prop, sz.nw.17w, nsim=999) 
sz.mi.18d<-moran.mc(sz.18d$mean.prop, sz.nw.18d, nsim=999) 
sz.mi.18w<-moran.mc(sz.18w$mean.prop, sz.nw.18w, nsim=999) 

###
#LISA#######
###
sz.locm.12d <- localmoran(sz.12d$mean.prop, sz.nw.12d)
sz.locm.12w <- localmoran(sz.12w$mean.prop, sz.nw.12w)
sz.locm.13d <- localmoran(sz.13d$mean.prop, sz.nw.13d)
sz.locm.13w <- localmoran(sz.13w$mean.prop, sz.nw.13w)
sz.locm.14d <- localmoran(sz.14d$mean.prop, sz.nw.14d)
sz.locm.14w <- localmoran(sz.14w$mean.prop, sz.nw.14w)
sz.locm.15d <- localmoran(sz.15d$mean.prop, sz.nw.15d)
sz.locm.15w <- localmoran(sz.15w$mean.prop, sz.nw.15w)
sz.locm.16d <- localmoran(sz.16d$mean.prop, sz.nw.16d)
sz.locm.16w <- localmoran(sz.16w$mean.prop, sz.nw.16w)
sz.locm.17d <- localmoran(sz.17d$mean.prop, sz.nw.17d)
sz.locm.17w <- localmoran(sz.17w$mean.prop, sz.nw.17w)
sz.locm.18d <- localmoran(sz.18d$mean.prop, sz.nw.18d)
sz.locm.18w <- localmoran(sz.18w$mean.prop, sz.nw.18w)

######
#Scaling mean.prop and creating lag variable for each year|season
#To later classify LISA clusters
sz.12d$mean.prop.s <- scale(sz.12d$mean.prop)
sz.12d$lag.prop.s <- lag.listw(sz.nw.12d, sz.12d$mean.prop.s)
sz.12w$mean.prop.s <- scale(sz.12w$mean.prop)
sz.12w$lag.prop.s <- lag.listw(sz.nw.12w, sz.12w$mean.prop.s)

sz.13d$mean.prop.s <- scale(sz.13d$mean.prop)
sz.13d$lag.prop.s <- lag.listw(sz.nw.13d, sz.13d$mean.prop.s)
sz.13w$mean.prop.s <- scale(sz.13w$mean.prop)
sz.13w$lag.prop.s <- lag.listw(sz.nw.13w, sz.13w$mean.prop.s)

sz.14d$mean.prop.s <- scale(sz.14d$mean.prop)
sz.14d$lag.prop.s <- lag.listw(sz.nw.14d, sz.14d$mean.prop.s)
sz.14w$mean.prop.s <- scale(sz.14w$mean.prop)
sz.14w$lag.prop.s <- lag.listw(sz.nw.14w, sz.14w$mean.prop.s)

sz.15d$mean.prop.s <- scale(sz.15d$mean.prop)
sz.15d$lag.prop.s <- lag.listw(sz.nw.15d, sz.15d$mean.prop.s)
sz.15w$mean.prop.s <- scale(sz.15w$mean.prop)
sz.15w$lag.prop.s <- lag.listw(sz.nw.15w, sz.15w$mean.prop.s)

sz.16d$mean.prop.s <- scale(sz.16d$mean.prop)
sz.16d$lag.prop.s <- lag.listw(sz.nw.16d, sz.16d$mean.prop.s)
sz.16w$mean.prop.s <- scale(sz.16w$mean.prop)
sz.16w$lag.prop.s <- lag.listw(sz.nw.16w, sz.16w$mean.prop.s)

sz.17d$mean.prop.s <- scale(sz.17d$mean.prop)
sz.17d$lag.prop.s <- lag.listw(sz.nw.17d, sz.17d$mean.prop.s)
sz.17w$mean.prop.s <- scale(sz.17w$mean.prop)
sz.17w$lag.prop.s <- lag.listw(sz.nw.17w, sz.17w$mean.prop.s)

sz.18d$mean.prop.s <- scale(sz.18d$mean.prop)
sz.18d$lag.prop.s <- lag.listw(sz.nw.18d, sz.18d$mean.prop.s)
sz.18w$mean.prop.s <- scale(sz.18w$mean.prop)
sz.18w$lag.prop.s <- lag.listw(sz.nw.18w, sz.18w$mean.prop.s)


##################################################################################
###Adding LISA clusters using the comparsions of lag and localmoran results

#2012
sz.12d$quad_sig <- NA
sz.12d[(sz.12d$mean.prop.s >= 0 & sz.12d$lag.prop.s >= 0) & (sz.locm.12d[, 5] <= 0.05), "quad_sig"] <- 1
sz.12d[(sz.12d$mean.prop.s <= 0 & sz.12d$lag.prop.s <= 0) & (sz.locm.12d[, 5] <= 0.05), "quad_sig"] <- 2
sz.12d[(sz.12d$mean.prop.s >= 0 & sz.12d$lag.prop.s <= 0) & (sz.locm.12d[, 5] <= 0.05), "quad_sig"] <- 3
sz.12d[(sz.12d$mean.prop.s <= 0 & sz.12d$lag.prop.s >= 0) & (sz.locm.12d[, 5] <= 0.05), "quad_sig"] <- 4
# sz.12d[(sz.12d$mean.prop.s <= 0 & sz.12d$lag.prop.s >= 0) & (sz.locm.12d[, 5] >= 0.05), "quad_sig"] <- 5 
sz.12d[(sz.locm.12d[, 5] >= 0.05), "quad_sig"] <- 5


sz.12w$quad_sig <- NA
sz.12w[(sz.12w$mean.prop.s >= 0 & sz.12w$lag.prop.s >= 0) & (sz.locm.12w[, 5] <= 0.05), "quad_sig"] <- 1
sz.12w[(sz.12w$mean.prop.s <= 0 & sz.12w$lag.prop.s <= 0) & (sz.locm.12w[, 5] <= 0.05), "quad_sig"] <- 2
sz.12w[(sz.12w$mean.prop.s >= 0 & sz.12w$lag.prop.s <= 0) & (sz.locm.12w[, 5] <= 0.05), "quad_sig"] <- 3
sz.12w[(sz.12w$mean.prop.s <= 0 & sz.12w$lag.prop.s >= 0) & (sz.locm.12w[, 5] <= 0.05), "quad_sig"] <- 4
sz.12w[(sz.locm.12w[, 5] >= 0.05), "quad_sig"] <- 5

sz.12d$iI<-sz.locm.12d[,1]
sz.12d$Pr<-sz.locm.12d[,5]
sz.12w$iI<-sz.locm.12w[,1]
sz.12w$Pr<-sz.locm.12w[,5]

#2013
sz.13d$quad_sig <- NA
sz.13d[(sz.13d$mean.prop.s >= 0 & sz.13d$lag.prop.s >= 0) & (sz.locm.13d[, 5] <= 0.05), "quad_sig"] <- 1
sz.13d[(sz.13d$mean.prop.s <= 0 & sz.13d$lag.prop.s <= 0) & (sz.locm.13d[, 5] <= 0.05), "quad_sig"] <- 2
sz.13d[(sz.13d$mean.prop.s >= 0 & sz.13d$lag.prop.s <= 0) & (sz.locm.13d[, 5] <= 0.05), "quad_sig"] <- 3
sz.13d[(sz.13d$mean.prop.s <= 0 & sz.13d$lag.prop.s >= 0) & (sz.locm.13d[, 5] <= 0.05), "quad_sig"] <- 4
sz.13d[(sz.locm.13d[, 5] >= 0.05), "quad_sig"] <- 5

sz.13w$quad_sig <- NA
sz.13w[(sz.13w$mean.prop.s >= 0 & sz.13w$lag.prop.s >= 0) & (sz.locm.13w[, 5] <= 0.05), "quad_sig"] <- 1
sz.13w[(sz.13w$mean.prop.s <= 0 & sz.13w$lag.prop.s <= 0) & (sz.locm.13w[, 5] <= 0.05), "quad_sig"] <- 2
sz.13w[(sz.13w$mean.prop.s >= 0 & sz.13w$lag.prop.s <= 0) & (sz.locm.13w[, 5] <= 0.05), "quad_sig"] <- 3
sz.13w[(sz.13w$mean.prop.s <= 0 & sz.13w$lag.prop.s >= 0) & (sz.locm.13w[, 5] <= 0.05), "quad_sig"] <- 4
sz.13w[(sz.locm.13w[, 5] >= 0.05), "quad_sig"] <- 5

sz.13d$iI<-sz.locm.13d[,1]
sz.13d$Pr<-sz.locm.13d[,5]
sz.13w$iI<-sz.locm.13w[,1]
sz.13w$Pr<-sz.locm.13w[,5]


#2014
sz.14d$quad_sig <- NA
sz.14d[(sz.14d$mean.prop.s >= 0 & sz.14d$lag.prop.s >= 0) & (sz.locm.14d[, 5] <= 0.05), "quad_sig"] <- 1
sz.14d[(sz.14d$mean.prop.s <= 0 & sz.14d$lag.prop.s <= 0) & (sz.locm.14d[, 5] <= 0.05), "quad_sig"] <- 2
sz.14d[(sz.14d$mean.prop.s >= 0 & sz.14d$lag.prop.s <= 0) & (sz.locm.14d[, 5] <= 0.05), "quad_sig"] <- 3
sz.14d[(sz.14d$mean.prop.s <= 0 & sz.14d$lag.prop.s >= 0) & (sz.locm.14d[, 5] <= 0.05), "quad_sig"] <- 4
sz.14d[(sz.locm.14d[, 5] >= 0.05), "quad_sig"] <- 5

sz.14w$quad_sig <- NA
sz.14w[(sz.14w$mean.prop.s >= 0 & sz.14w$lag.prop.s >= 0) & (sz.locm.14w[, 5] <= 0.05), "quad_sig"] <- 1
sz.14w[(sz.14w$mean.prop.s <= 0 & sz.14w$lag.prop.s <= 0) & (sz.locm.14w[, 5] <= 0.05), "quad_sig"] <- 2
sz.14w[(sz.14w$mean.prop.s >= 0 & sz.14w$lag.prop.s <= 0) & (sz.locm.14w[, 5] <= 0.05), "quad_sig"] <- 3
sz.14w[(sz.14w$mean.prop.s <= 0 & sz.14w$lag.prop.s >= 0) & (sz.locm.14w[, 5] <= 0.05), "quad_sig"] <- 4
sz.14w[(sz.locm.14w[, 5] >= 0.05), "quad_sig"] <- 5

sz.14d$iI<-sz.locm.14d[,1]
sz.14d$Pr<-sz.locm.14d[,5]
sz.14w$iI<-sz.locm.14w[,1]
sz.14w$Pr<-sz.locm.14w[,5]

#2015
sz.15d$quad_sig <- NA
sz.15d[(sz.15d$mean.prop.s >= 0 & sz.15d$lag.prop.s >= 0) & (sz.locm.15d[, 5] <= 0.05), "quad_sig"] <- 1
sz.15d[(sz.15d$mean.prop.s <= 0 & sz.15d$lag.prop.s <= 0) & (sz.locm.15d[, 5] <= 0.05), "quad_sig"] <- 2
sz.15d[(sz.15d$mean.prop.s >= 0 & sz.15d$lag.prop.s <= 0) & (sz.locm.15d[, 5] <= 0.05), "quad_sig"] <- 3
sz.15d[(sz.15d$mean.prop.s <= 0 & sz.15d$lag.prop.s >= 0) & (sz.locm.15d[, 5] <= 0.05), "quad_sig"] <- 4
sz.15d[(sz.locm.15d[, 5] >= 0.05), "quad_sig"] <- 5

sz.15w$quad_sig <- NA
sz.15w[(sz.15w$mean.prop.s >= 0 & sz.15w$lag.prop.s >= 0) & (sz.locm.15w[, 5] <= 0.05), "quad_sig"] <- 1
sz.15w[(sz.15w$mean.prop.s <= 0 & sz.15w$lag.prop.s <= 0) & (sz.locm.15w[, 5] <= 0.05), "quad_sig"] <- 2
sz.15w[(sz.15w$mean.prop.s >= 0 & sz.15w$lag.prop.s <= 0) & (sz.locm.15w[, 5] <= 0.05), "quad_sig"] <- 3
sz.15w[(sz.15w$mean.prop.s <= 0 & sz.15w$lag.prop.s >= 0) & (sz.locm.15w[, 5] <= 0.05), "quad_sig"] <- 4
sz.15w[(sz.locm.15w[, 5] >= 0.05), "quad_sig"] <- 5

sz.15d$iI<-sz.locm.15d[,1]
sz.15d$Pr<-sz.locm.15d[,5]
sz.15w$iI<-sz.locm.15w[,1]
sz.15w$Pr<-sz.locm.15w[,5]


#2016
sz.16d$quad_sig <- NA
sz.16d[(sz.16d$mean.prop.s >= 0 & sz.16d$lag.prop.s >= 0) & (sz.locm.16d[, 5] <= 0.05), "quad_sig"] <- 1
sz.16d[(sz.16d$mean.prop.s <= 0 & sz.16d$lag.prop.s <= 0) & (sz.locm.16d[, 5] <= 0.05), "quad_sig"] <- 2
sz.16d[(sz.16d$mean.prop.s >= 0 & sz.16d$lag.prop.s <= 0) & (sz.locm.16d[, 5] <= 0.05), "quad_sig"] <- 3
sz.16d[(sz.16d$mean.prop.s <= 0 & sz.16d$lag.prop.s >= 0) & (sz.locm.16d[, 5] <= 0.05), "quad_sig"] <- 4
sz.16d[(sz.locm.16d[, 5] >= 0.05), "quad_sig"] <- 5

sz.16w$quad_sig <- NA
sz.16w[(sz.16w$mean.prop.s >= 0 & sz.16w$lag.prop.s >= 0) & (sz.locm.16w[, 5] <= 0.05), "quad_sig"] <- 1
sz.16w[(sz.16w$mean.prop.s <= 0 & sz.16w$lag.prop.s <= 0) & (sz.locm.16w[, 5] <= 0.05), "quad_sig"] <- 2
sz.16w[(sz.16w$mean.prop.s >= 0 & sz.16w$lag.prop.s <= 0) & (sz.locm.16w[, 5] <= 0.05), "quad_sig"] <- 3
sz.16w[(sz.16w$mean.prop.s <= 0 & sz.16w$lag.prop.s >= 0) & (sz.locm.16w[, 5] <= 0.05), "quad_sig"] <- 4
sz.16w[(sz.locm.16w[, 5] >= 0.05), "quad_sig"] <- 5

sz.16d$iI<-sz.locm.16d[,1]
sz.16d$Pr<-sz.locm.16d[,5]
sz.16w$iI<-sz.locm.16w[,1]
sz.16w$Pr<-sz.locm.16w[,5]

#2017
sz.17d$quad_sig <- NA
sz.17d[(sz.17d$mean.prop.s >= 0 & sz.17d$lag.prop.s >= 0) & (sz.locm.17d[, 5] <= 0.05), "quad_sig"] <- 1
sz.17d[(sz.17d$mean.prop.s <= 0 & sz.17d$lag.prop.s <= 0) & (sz.locm.17d[, 5] <= 0.05), "quad_sig"] <- 2
sz.17d[(sz.17d$mean.prop.s >= 0 & sz.17d$lag.prop.s <= 0) & (sz.locm.17d[, 5] <= 0.05), "quad_sig"] <- 3
sz.17d[(sz.17d$mean.prop.s <= 0 & sz.17d$lag.prop.s >= 0) & (sz.locm.17d[, 5] <= 0.05), "quad_sig"] <- 4
sz.17d[(sz.locm.17d[, 5] >= 0.05), "quad_sig"] <- 5

sz.17w$quad_sig <- NA
sz.17w[(sz.17w$mean.prop.s >= 0 & sz.17w$lag.prop.s >= 0) & (sz.locm.17w[, 5] <= 0.05), "quad_sig"] <- 1
sz.17w[(sz.17w$mean.prop.s <= 0 & sz.17w$lag.prop.s <= 0) & (sz.locm.17w[, 5] <= 0.05), "quad_sig"] <- 2
sz.17w[(sz.17w$mean.prop.s >= 0 & sz.17w$lag.prop.s <= 0) & (sz.locm.17w[, 5] <= 0.05), "quad_sig"] <- 3
sz.17w[(sz.17w$mean.prop.s <= 0 & sz.17w$lag.prop.s >= 0) & (sz.locm.17w[, 5] <= 0.05), "quad_sig"] <- 4
sz.17w[(sz.locm.17w[, 5] >= 0.05), "quad_sig"] <- 5

sz.17d$iI<-sz.locm.17d[,1]
sz.17d$Pr<-sz.locm.17d[,5]
sz.17w$iI<-sz.locm.17w[,1]
sz.17w$Pr<-sz.locm.17w[,5]

#2018
sz.18d$quad_sig <- NA
sz.18d[(sz.18d$mean.prop.s >= 0 & sz.18d$lag.prop.s >= 0) & (sz.locm.18d[, 5] <= 0.05), "quad_sig"] <- 1
sz.18d[(sz.18d$mean.prop.s <= 0 & sz.18d$lag.prop.s <= 0) & (sz.locm.18d[, 5] <= 0.05), "quad_sig"] <- 2
sz.18d[(sz.18d$mean.prop.s >= 0 & sz.18d$lag.prop.s <= 0) & (sz.locm.18d[, 5] <= 0.05), "quad_sig"] <- 3
sz.18d[(sz.18d$mean.prop.s <= 0 & sz.18d$lag.prop.s >= 0) & (sz.locm.18d[, 5] <= 0.05), "quad_sig"] <- 4
sz.18d[(sz.locm.18d[, 5] >= 0.05), "quad_sig"] <- 5

sz.18w$quad_sig <- NA
sz.18w[(sz.18w$mean.prop.s >= 0 & sz.18w$lag.prop.s >= 0) & (sz.locm.18w[, 5] <= 0.05), "quad_sig"] <- 1
sz.18w[(sz.18w$mean.prop.s <= 0 & sz.18w$lag.prop.s <= 0) & (sz.locm.18w[, 5] <= 0.05), "quad_sig"] <- 2
sz.18w[(sz.18w$mean.prop.s >= 0 & sz.18w$lag.prop.s <= 0) & (sz.locm.18w[, 5] <= 0.05), "quad_sig"] <- 3
sz.18w[(sz.18w$mean.prop.s <= 0 & sz.18w$lag.prop.s >= 0) & (sz.locm.18w[, 5] <= 0.05), "quad_sig"] <- 4
sz.18w[(sz.locm.18w[, 5] >= 0.05), "quad_sig"] <- 5

sz.18d$iI<-sz.locm.18d[,1]
sz.18d$Pr<-sz.locm.18d[,5]
sz.18w$iI<-sz.locm.18w[,1]
sz.18w$Pr<-sz.locm.18w[,5]

sz.all<-rbind(sz.12d, sz.12w, sz.13d, sz.13w, sz.14d, sz.14w,
              sz.15d, sz.15w, sz.16d, sz.16w, sz.17d, sz.17w, sz.18d, sz.18w)

st_write(sz.all, "sz.all.shp")
test<-st_read("sz.all.shp")

###################################################################################
###########Pattern of E across Years and Seasons#####################

###
#Plots
###

E_year.month2_df<-subset(E_year.month2_df, Year < 2019)


#I have to improve sine curve
# colorx.2<-c("orange3", "orange", "orangered", "red", "red4", "lightcyan", "cyan", "deepskyblue", "deepskyblue3", "blue", "blue4", "yellow")
colorx.2 <- c("red", "blue")
# ggplot(E_year.month2_df, aes(fYear.Month, Eadj, colour = fMonth))+ 
#   geom_point(size = 4) +
#   theme_bw()+
#   scale_colour_manual(name = "Months", values = colorx.2) +
#   labs(title = "Individual variation cycle", x = "Year.Month", y = "E index")+ 
#   geom_smooth(method = "lm", formula = y ~ sin(2*pi*x)+cos(2*pi*x), se = FALSE, linetype=2, colour = "black")+
#   theme(axis.text = element_text(size = 18, face = "bold", colour = "black"), 
#         axis.title = element_text(size = 20, face = "bold"),
#         plot.title = element_text(size = 20, face = "bold"))

season.cycle <- ggplot(E_year.month2_df, aes(fYear.Month, Eadj, colour = Season))+ 
  geom_point(size = 4) +
  theme_bw()+
  scale_colour_manual(name = "Seasons", values = colorx.2) +
  labs(x = "Year.Month", y = "E index")+ 
  geom_smooth(method = "lm", formula = y ~ sin(2*pi*x)+cos(2*pi*x), se = FALSE, linetype=2, colour = "black")+
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"))

####Season E means##
meanE_Season<-ddply(E_year.month2_df, c("fYear", "Season"), summarise,
                    N = sum(!is.na(fYear)),
                    mean.E = mean(Eadj),
                    sd.E  = sd(Eadj),
                    se.E = sd.E/sqrt(N))

meanE_Season$Season.Yr <- interaction(meanE_Season$Season, meanE_Season$fYear, drop = TRUE)


limits.1<-aes(ymax = mean.E + se.E, ymin = mean.E - se.E)
colorx.2<-c("red", "blue")
pd <- position_dodge(width = 0.4)
ggplot(meanE_Season, aes(fYear, mean.E, colour = Season))+ 
  geom_point(position = pd, size = 4) +
  geom_errorbar(limits.1, size = 2, width = 0.5 , position = pd)+
  theme_bw()+
  # scale_y_continuous(limits = c(0.5, 0.8), breaks=seq(0.5,0.8,0.1)) + 
  scale_colour_manual(name = "Seasons", values = colorx.2) +
  labs(title = "Inter-annual Variation" , x = "Year", y = "E index")+ 
  theme(axis.text = element_text(size = 18, face = "bold", colour = "black"), 
        axis.title = element_text(size = 20, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold"))


meanE_Season2<-ddply(E_year.month2_df, c("Season"), summarise,
                     N = sum(!is.na(fYear)),
                     mean.E = mean(Eadj),
                     sd.E  = sd(Eadj),
                     se.E = sd.E/sqrt(N))

limits.1<-aes(ymax = mean.E + se.E, ymin = mean.E - se.E)
colorx.2<-c("red", "blue")
ggplot(meanE_Season2, aes(Season, mean.E, colour = Season))+ geom_point(size = 4) +
  geom_errorbar(limits.1, size = 2, width = 0.5)+theme_bw()+
  scale_y_continuous(limits = c(0.5, 0.8), breaks=seq(0.5,0.8,0.1)) + 
  scale_colour_manual(name = "Seasons", values = colorx.2) +
  labs(title = "Seasonal Difference", x = "Year", y = "E index")+ 
  theme(axis.text = element_text(size = 18, face = "bold", colour  = "black"), 
        axis.title = element_text(size = 20, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold"))



################################
#Plots from the models
###############################

###
#gam - monthly cycle of E by Year
###

temporal.fit <- ggplot(model.temporal.fit, aes(nMonth, visregFit))+
  facet_wrap(~fYear, scales = "free_x")+
  theme_bw()+
  geom_line(size = 2, colour= "black", size = 2) + 
  geom_line(linetype = 2, colour = "black", aes(y = visregLwr))+
  geom_line(linetype = 2, colour = "black", aes(y = visregUpr)) +
  labs(x = "Month", y = "E Index")+
  scale_x_continuous(limits = c(2, 12), breaks=seq(2,12,2)) +
  theme(axis.text = element_text(size = 12, face = "bold", colour = "black"), 
        axis.title = element_text(size = 16, face = "bold", colour = "black"), 
        strip.text.x = element_text(size=11, face="bold"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+
  geom_rect(aes(xmin=2, xmax=6, ymin=-Inf, ymax=Inf), alpha = 0.007, fill = "red")+ 
  geom_rect(aes(xmin=6, xmax=11, ymin=-Inf, ymax=Inf), alpha = 0.007, fill = "blue")+
  geom_rect(aes(xmin=11, xmax=12, ymin=-Inf, ymax=Inf), alpha = 0.007, fill = "red")

figure.q1 <- ggarrange(season.cycle, temporal.fit, nrow = 2, labels = c("A", "B"))
ggsave("figure.q1.png", width = 20, height = 25, units = "cm")


###
#GAMs - E as a function of stage at different locations of the river

###
#First approach - using mean and sd of stage across four sites-------
###

model.mstage.fit2$Site.Season2 <- model.mstage.fit2$Site.Season
model.mstage.fit2<-cSplit(model.mstage.fit2, "Site.Season2", sep = ".", type.convert = FALSE)
colnames(model.mstage.fit2)[12:13]<-c("Location", "Season")
model.mstage.fit2$Location<-factor(model.mstage.fit2$Location)
levels(model.mstage.fit2$Location)<-list("Upstream" = c("Upstream"), "Upstream-Bay" = c("Upstream-Bay"), "Bay" = c("Bay"), "Downstream" = c("Downstream"))

model.sdstage.fit2$Site.Season2 <- model.sdstage.fit2$Site.Season
model.sdstage.fit2<-cSplit(model.sdstage.fit2, "Site.Season2", sep = ".", type.convert = FALSE)
colnames(model.sdstage.fit2)[12:13]<-c("Location", "Season")
model.sdstage.fit2$Location<-factor(model.sdstage.fit2$Location)
levels(model.sdstage.fit2$Location)<-list("Upstream" = c("Upstream"), "Upstream-Bay" = c("Upstream-Bay"), "Bay" = c("Bay"), "Downstream" = c("Downstream"))


colorx.3<-c("red", "blue")
plot.mean.stage<-ggplot(model.mstage.fit2, aes(mean.stage, visregFit, colour = Season))+
  facet_grid(scales = "free_x", rows = vars(Location))+
  theme_bw()+
  geom_line(size = 2, size = 2) + 
  geom_line(linetype = 2, aes(y = visregLwr))+
  geom_line(linetype = 2, aes(y = visregUpr)) +
  labs(x = "River Stage (cm)", y = "E Index")+
  # scale_x_continuous(limits = c(2, 12), breaks=seq(2,12,2)) +
  scale_colour_manual(name = "Seasons", values = colorx.3) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title = element_text(size = 16, face = "bold", colour = "black"), 
        strip.text.y = element_text(size=12, face="bold"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank())#,
        #panel.border = element_blank())

colorx.3<-c("red", "blue")
plot.sd.stage<-ggplot(model.sdstage.fit2, aes(sd.stage, visregFit, colour = Season))+
  facet_grid(scales = "free_x", rows = vars(Location))+
  theme_bw()+
  geom_line(size = 2, size = 2) + 
  geom_line(linetype = 2, aes(y = visregLwr))+
  geom_line(linetype = 2, aes(y = visregUpr)) +
  labs(x = "Variation in River Stage (cm)", y = "E Index")+
  # scale_x_continuous(limits = c(2, 12), breaks=seq(2,12,2)) +
  scale_colour_manual(name = "Seasons", values = colorx.3) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title = element_text(size = 16, face = "bold", colour = "black"), 
        strip.text.y = element_text(size=12, face="bold"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank())
        # legend.position = "none")#,
#panel.border = element_blank())

plot_grid(plot.mean.stage, plot.sd.stage,  
          align = "h", ncol = 2, labels = c("a)", "b)"),
          label_size = 16, hjust = -0.3)


###
#Second approach using PCA scores-----
###

colorx.3<-c("red", "blue")
plot.pc1<-ggplot(model.pca1.fit2, aes(PC1, visregFit, colour = Season))+
  # facet_grid(scales = "free_x", rows = vars(Location))+
  theme_bw()+
  geom_line(size = 2, size = 2) + 
  geom_line(linetype = 2, aes(y = visregLwr))+
  geom_line(linetype = 2, aes(y = visregUpr)) +
  labs(x = "Mean Stage Index (PC1)", y = "E Index")+
  # scale_x_continuous(limits = c(2, 12), breaks=seq(2,12,2)) +
  scale_colour_manual(name = "Seasons", values = colorx.3) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title = element_text(size = 16, face = "bold", colour = "black"), 
        strip.text.y = element_text(size=12, face="bold"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank())#,
#panel.border = element_blank())

colorx.3<-c("red", "blue")
plot.pc3<-ggplot(model.pca3.fit2, aes(PC3, visregFit, colour = Season))+
  # facet_grid(scales = "free_x", rows = vars(Location))+
  theme_bw()+
  geom_line(size = 2, size = 2) + 
  geom_line(linetype = 2, aes(y = visregLwr))+
  geom_line(linetype = 2, aes(y = visregUpr)) +
  labs(x = "Stage Variance Index (PC3)", y = "E Index")+
  # scale_x_continuous(limits = c(2, 12), breaks=seq(2,12,2)) +
  scale_colour_manual(name = "Seasons", values = colorx.3) +
  theme(axis.text = element_text(size = 14, face = "bold", colour = "black"), 
        axis.title = element_text(size = 16, face = "bold", colour = "black"), 
        strip.text.y = element_text(size=12, face="bold"),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank())
# legend.position = "none")#,
#panel.border = element_blank())

#cowplot grid
# plot_grid(plot.mean.stage, plot.sd.stage,  
#           align = "h", ncol = 2, labels = c("a)", "b)"),
#           label_size = 16, hjust = -0.3)

figure.q2 <- ggarrange(plot.pc1, plot.pc3, nrow = 2, labels = c("A", "B"))
ggsave("figure.q2.png", width = 20, height = 20, units = "cm")

