# Set up working space (important info, workind directory,  librar --------

### important information
#Author: Mack White
#Date: November 27, 2022
#Project: Spatial Similarity Manuscript
#Description: Trying out generating Eadj values for only peak snook abundance months (i.e., March, April, May)

### set working directory
setwd("/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript/ThreeMonthTest")

### load libraries
library(ggplot2)
library(car)
library(reshape2)
library(reshape)
library(plyr)
library(dplyr)
library(tidyr)
library(visreg)
library(modEvA) 
library(gridExtra)
library(AICcmodavg)
library(nlme)
library(mgcv)
library(lme4)
library(splitstackshape) ###package to use cSplit and do text to column like excel
library(chron)
library(RInSp)#package for intraspecific niche variation
library(boot)#package for boostrapping operations
library(cowplot)
library(ggpubr)
library(lsmeans)
library(MixSIAR)
library(tidyverse)
library(hypervolume)
library(truncnorm)
library(tidyverse)
library(alphahull)
library(dplyr)
library(svMisc)
library(writexl)
library(ggplot2)
library(ggstatsplot)#package for cool correlation plots
library(ggside)#associated package for cool correlation plots

# Calculate new E Values (March,  April,  May - dry months) ---------------

#read in most updated version of rds file from "snook_similarity_HV_MW_V3.R"
snook3 <- readRDS("snook_SR_detection_20122021.rds")

#check out the file
glimpse(snook3)
str(snook3)
head(snook3)

# Calculate POR metrics for each individual -------------------------------

# line 1 = grouping individuals in data set by ID
# line 2 = summary statistics det_f = total number of detections; min.t = first date picked up on array; max.t = last day picked up on array; resid_t = number of days picked up in array; unique.date = number of year/month dates picked up in array

snook.detdays<-group_by(snook3, ID)%>%
      summarise(det_f = sum(!is.na(fDay)),
                min.t = min(Datetime_UTC),
                max.t = max(Datetime_UTC),
                resid.t = difftime(max.t, min.t, units = "day"),
                unique.date = length(unique(fYear.Month)))


# Filter based on >100 detections & >90 days ------------------------------

#Subset to obtain individuals over 90 days and over 100 detections
snook.90day <- subset(snook.detdays, resid.t >= 90 & det_f > 100)

#Filter tag IDs based on cutoff information derived above
snook.90day <- filter(snook3, ID %in% snook.90day$ID)

#Creating table of unique tags per year for paper/ppt
snook.90day.x <- filter(snook.90day, Year < 2022)
length(unique(snook.90day.x$ID))

tags.summary <- group_by(snook.90day.x, Year)%>%
      summarise(n = length(unique(ID)))

list.tags <- group_by(snook.90day.x, ID, Year)%>%
      summarise(Total.Det = sum(!is.na(Datetime_UTC)),
                Total.Yr.Month = length(unique(fYear.Month)))

write.csv(tags.summary, "tags.summary_11_27_2022.csv")
write.csv(list.tags, "list.tags_11_27_2022.csv")

# Creating the number of days as a metric of habitat use per ID|Zo --------

snook_yrd_zones<-ddply(snook.90day, c("ID", "Zone", "Year","fYear.Month"), summarise,
                       freq.v = length(unique(fDay)))

snook_yrd_zones<-subset(snook_yrd_zones, fYear.Month != "NA")
snook_yrd_zones2<-subset(snook_yrd_zones, Zone != "NA")


# Create list of tracks vs receivers listed by year.month -----------------

### custom function one
dcast.snook<-function(x){
      dcast(x, ID ~ Zone, sum, value.var = "freq.v")}

list.snook.zones<-dlply(snook_yrd_zones2, .(fYear.Month), dcast.snook)

### custom function two
fc <- function(d, i){
      #d1<-d[i,]
      d2<-import.RInSp(d, row.names = 1, col.header = TRUE)
      EAdj<-Emc(d2, popd.type = "average", replicates = 999)
      return(EAdj$Eadj)
}

E_year.month3<-lapply(list.snook.zones, fc)


# Update Dataframe with Eadj Values ---------------------------------------

library(stringr)
E_year.month3_df<-ldply(E_year.month3, data.frame)
colnames(E_year.month3_df)<-c("Year.Month", "Eadj")

library(zoo)
#using function zoo::yearmon to classify in R fYear.Month as date formated as Year|Month
E_year.month3_df$fYear.Month <- as.yearmon(E_year.month3_df$Year.Month, "%Y/%m")

E_year.month3_df$Date2<-E_year.month3_df$Year.Month
E_year.month3_df<-cSplit(E_year.month3_df, "Date2", sep = "/", type.convert = FALSE)
colnames(E_year.month3_df)[4:5]<-c("Year", "Month")

E_year.month3_df$Year <- as.numeric(E_year.month3_df$Year)
E_year.month3_df$fYear<-factor(E_year.month3_df$Year, levels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021"))


E_year.month3_df$fMonth<-factor(E_year.month3_df$Month)
levels(E_year.month3_df$fMonth) <- list("Jan" = c("01"), "Feb" = c("02"), "Mar" = c("03"), 
                                        "Apr" = c("04"), "May" = c("05"), "Jun" = c("06"), 
                                        "Jul" = c("07"),"Aug" = c("08"), "Sep" = c("09"), 
                                        "Oct" = c("10"), "Nov" = c("11"), "Dec" = c("12"))

E_year.month3_df$Peak<-factor(E_year.month3_df$fMonth)
levels(E_year.month3_df$Peak)<-list("PeakSnook" = c("Mar", "Apr", "May"), 
                                      "OffPeak" = c("Jan", "Feb", "Jun", "Jul", "Aug", 
                                                    "Sep", "Oct", "Nov", "Dec"))

### save things as is...
write.csv(E_year.month3_df, "Eadj_2012_2021_SRSnook_MW_11_27_2022.csv")
### read into excel and enter water year by hand... kind of confusing and dont know how to do using actual code

### in excel... added water years, included wYear 2011... since it contained e adjusted values corresponded to large # of isotope samples... however, we only have dry 2011 edaj values (four months of dry season values) so wYear 2011 will need to be removed for plotting... further, we only have e adj values for very beginning of wYear 2021...(i.e., wet season) -> completely removed these values because they will plot alone AND gives us no information for trophic questions... once again, remove wYear 2011 for plotting but keep values for future comparisons

### read in revised csv file with water year appropriately categorized
E_year.month3_df <- read.csv("Eadj_wYEAR_2011_2020_SRSnook_MW_11_27_2022.csv")
str(E_year.month3_df)

### remove weird NA rows at the bottom of the dataframe
E_year.month3_df <- E_year.month3_df %>%
      filter(!row_number() %in% c(113,114,115))
view(E_year.month3_df)

### make sure that wYear and season is a factor with levels
E_year.month3_df$wYear<-factor(E_year.month3_df$wYear, levels = c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020"))

E_year.month3_df$Peak<-factor(E_year.month3_df$Peak, levels = c("OffPeak", "PeakSnook"))

Eadj_Peak_Annual = E_year.month3_df %>% 
      group_by(wYear, Peak) %>%
      summarize(Eadj_Seasonal_Mean = mean(Eadj), 
                Eadj_Seasonal_SD = sd(Eadj), 
                Eadj_Seasonal_Min = min(Eadj),
                Eadj_Seasonal_Max = max(Eadj))

write.csv(Eadj_Peak_Annual, "Eadj_Peak_Annual_11_27_2022.csv")

### PLOT SEASONAL EADJ VALUES (Ie Peak VS Off-Peak)

tiff("/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript/ThreeMonthTest/PEAK_withWY2011_Eadj_MW_11_27_2022.tiff", width = 8, height = 8, units = 'in', res = 600, compression = 'lzw')

ggplot(E_year.month3_df, aes(x=Peak, y=Eadj, fill=Peak)) +
      geom_boxplot(width =0.8) +
      scale_fill_manual(values = c("cadetblue", "darkgoldenrod")) +
      labs(x = "Peak vs. Non-Peak Snook Abundance", 
           y = "Individual Specialization (Eadj)") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=14, face="bold", color = "black")) +
      theme(axis.text = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.x = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.y = element_text(size=16,face="bold", color = "black")) +
      theme(axis.title = element_text(size=16,face="bold", color = "black")) +
      theme(legend.position = "none")

dev.off()

### PLOT SEASONAL EADJ VALUES (Ie PEAK VS OFF-PEAK BUT FOR EACH YEAR)

tiff("/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript/ThreeMonthTest/PEAK_ANNUAL_withWY2011_Eadj_MW_11_27_2022.tiff", width = 16, height = 9, units = 'in', res = 600, compression = 'lzw')

ggplot(E_year.month3_df, aes(x=wYear, y=Eadj, fill=Peak)) +
      geom_boxplot(width =0.8) +
      scale_fill_manual(values = c("cadetblue", "darkgoldenrod")) +
      labs(x = "Water Year", 
           y = "Individual Specialization (Eadj)") +
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

dev.off()

### For plotting purposes, get ride of wYear 2011... has good dry season values, but not great wet season... will look funky... remove and plot below using the subset

E_year.month3_df <- E_year.month3_df %>%
      filter(!row_number() %in% c(1,2,3,4))

### looks good! Got ride of wYear 2011 (really just dry 2011)

### PLOT SEASONAL EADJ VALUES (IE Peak VS Offpeak)

tiff("/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript/ThreeMonthTest/PEAK_withoutWY2011_Eadj_MW_11_27_2022.tiff", width = 8, height = 8, units = 'in', res = 600, compression = 'lzw')

ggplot(E_year.month3_df, aes(x=Peak, y=Eadj, fill=Peak)) +
      geom_boxplot(width =0.8) +
      scale_fill_manual(values = c("cadetblue", "darkgoldenrod")) +
      labs(x = "Peak vs. Non-Peak Snook Abundance", 
           y = "Individual Specialization (Eadj)") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=14, face="bold", color = "black")) +
      theme(axis.text = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.x = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.y = element_text(size=16,face="bold", color = "black")) +
      theme(axis.title = element_text(size=16,face="bold", color = "black")) +
      theme(legend.position = "none")

dev.off()

### PLOT SEASONAL EADJ VALUES (IE PEAK VS OFF-PEAK BUT FOR EACH YEAR)

tiff("/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript/ThreeMonthTest/PEAK_ANNUAL_withoutWY2011_Eadj_MW_11_27_2022.tiff", width = 16, height = 9, units = 'in', res = 600, compression = 'lzw')

ggplot(E_year.month3_df, aes(x=wYear, y=Eadj, fill=Peak)) +
      geom_boxplot(width =0.8) +
      scale_fill_manual(values = c("cadetblue", "darkgoldenrod")) +
      labs(x = "Water Year", 
           y = "Individual Specialization (Eadj)") +
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

dev.off()

# Correlation Plots with Updated Eadj Values ------------------------------

# FORMAT AND BROUGHT BACK IN FOR GGPLOT FIGURE #2

long_df <- read.csv("snook_similarity_long_3monthtest_11272022.csv")

long_df$Year<-factor(long_df$Year, levels = c("2011", "2012", "2013", "2014", 
                                              "2016", "2017", "2018", "2019", "2020"))

str(long_df)
#### PRODUCT NUMBER TWO = PLOT OF NICHE SIZE (VIOLIN PLOT?)

################################################################################

tiff("wYearVolByYear_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggplot(long_df, aes(x=Year, y=vol, group=Year)) +
      geom_boxplot(width =0.8) +
      labs(x = "Water Year", 
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

dev.off()

#####

obs_check <- long_df %>%
      group_by(Year) %>%
      summarize(obs = mean(obs), 
                volume = mean(vol), 
                volumesd = sd(vol))

library(ggstatsplot)
library(ggside)

# plot observations alongside volume to look for correlation

tiff("vol_function of obs_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = obs_check,
      x = obs,
      y = volume,
      bf.message = FALSE
)

dev.off()

### plot eadj against volume metrics, but first summarize data

eadj_corr <- long_df %>%
      group_by(Year) %>%
      summarize(volume = mean(vol),
                Eadj_Mean = mean(Eadj_Seasonal_Mean),
                Eadj_SD = mean(Eadj_Seasonal_SD),
                Eadj_Min = mean(Eadj_Seasonal_Min),
                Eadj_Max = mean(Eadj_Seasonal_Max))

### eadj vs mean

tiff("vol_function of eadj_MEAN_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_Mean,
      y = volume,
      bf.message = FALSE
)

dev.off()

### eadj vs sd

tiff("vol_function of eadj_SD_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_SD,
      y = volume,
      bf.message = FALSE
)

dev.off()

### eadj vs min

tiff("vol_function of eadj_MIN_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_Min,
      y = volume,
      bf.message = FALSE
)

dev.off()

### eadj vs max

tiff("vol_function of eadj_MAX_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_Max,
      y = volume,
      bf.message = FALSE
)

dev.off()

# REDO Correlations without Water Year 2012 and 2020 (low sample s --------

long_df <- read.csv("snook_similarity_long_NO OUTLIER_3monthtest_11272022.csv")

long_df$Year<-factor(long_df$Year, levels = c("2011", "2013", "2014", 
                                              "2016", "2017", "2018", "2019"))

str(long_df)
#### PRODUCT NUMBER TWO = PLOT OF NICHE SIZE (VIOLIN PLOT?)

################################################################################

tiff("wYearVolByYear_NO OUTLIERS_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggplot(long_df, aes(x=Year, y=vol, group=Year)) +
      geom_boxplot(width =0.8) +
      labs(x = "Water Year", 
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

dev.off()

#####

obs_check <- long_df %>%
      group_by(Year) %>%
      summarize(obs = mean(obs), 
                volume = mean(vol), 
                volumesd = sd(vol))

library(ggstatsplot)
library(ggside)

# plot observations alongside volume to look for correlation

tiff("vol_function of obs_NO OUTLIERS_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = obs_check,
      x = obs,
      y = volume,
      bf.message = FALSE
)

dev.off()

### plot eadj against volume metrics, but first summarize data

eadj_corr <- long_df %>%
      group_by(Year) %>%
      summarize(volume = mean(vol),
                Eadj_Mean = mean(Eadj_Seasonal_Mean),
                Eadj_SD = mean(Eadj_Seasonal_SD),
                Eadj_Min = mean(Eadj_Seasonal_Min),
                Eadj_Max = mean(Eadj_Seasonal_Max))

### eadj vs mean

tiff("vol_function of eadj_MEAN_NO OUTLIERS_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_Mean,
      y = volume,
      bf.message = FALSE
)

dev.off()

### eadj vs sd

tiff("vol_function of eadj_SD_NO OUTLIERS_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_SD,
      y = volume,
      bf.message = FALSE
)

dev.off()

### eadj vs min

tiff("vol_function of eadj_MIN_NO OUTLIERS_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_Min,
      y = volume,
      bf.message = FALSE
)

dev.off()

### eadj vs max

tiff("vol_function of eadj_MAX_NO OUTLIERS_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_Max,
      y = volume,
      bf.message = FALSE
)

dev.off()
