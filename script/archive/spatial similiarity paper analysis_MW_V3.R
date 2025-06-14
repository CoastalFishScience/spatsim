#Author: Mack White
#Date: July/August 2023
#Project: Spatial Similarity Manuscript - Calculating Eadj Values

# Setting Things Up -------------------------------------------------------
# select softwrap long lines and rainbow parantheses under "Code"

# load in necessary libraries
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

################################################################################
# SKIP TO LINE 122 - DATA SET ALREADY FORMATTED AND PROCESSED CORRECTLY IN LINES 44-119 through September 2021
################################################################################

# Read Data and Subset Snook (though already technically done) ----------------------------------------------
tracks<-readRDS("AllSnook_PeriodOfRecord_09072021.rds")

### subset only variables that you need
snook<-subset(tracks, Species %in% c("Snook"))%>%
      subset(., select = c("Transmitter", "Station.Name", "Datetime_UTC",
                           "Receiver", "Species", "Station", 
                           "Distance", "Longitude", "Latitude",
                           "n"))

# Processing --------------------------------------------------------------

# separate datetime column into date and time columns
snook$DateTime2 <- snook$Datetime_UTC
snook<-cSplit(snook, "DateTime2", sep = " ", type.convert = FALSE)
colnames(snook)[11:12]<-c("Date", "Time")

# separate date column into year, month, and day columns
snook$Date2 <- snook$Date
snook<-cSplit(snook, "Date2", sep = "-", type.convert = FALSE)
colnames(snook)[13:15]<-c("Year", "Month", "Day")

# separate the transmitter column into code1, serial, and ID columns
snook<-cSplit(snook, "Transmitter", sep = "-", type.convert = FALSE)
colnames(snook)[15:17]<-c("code1", "serial", "ID")

# separate receiver column into code and receiver id columns
snook<-cSplit(snook, "Receiver", sep = "-", type.convert = FALSE)
colnames(snook)[17:18]<-c("code", "Receiver.ID")

# generate categorical variable to delineate lower river vs bay vs upper river detections - name "Home"
snook$Home<-cut(snook$Distance, breaks = c(-1,15,23,31.72), 
                labels = c("Lower River", "Bay", "Upper River"))

### SHOULD HAVE ALL COLUMNS OF INTEREST GENERATED FOLLOWING DATA PROCESSING

# Classify each variable appropriately ------------------------------------

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

# Initial Zone Delineation by Rolando Santos ------------------------------

# zones2018 <- select(snook, Station, Distance, Longitude, Latitude) %>%
#   distinct()

# Repeat Zone Delineation According to Matich et al 2017 ------------------

zones<-read.csv("zones2019.csv")
zones$f.Distance<-factor(zones$Distance)
zones<-subset(zones, Zone != "NA")
snook2<-merge(snook2, zones, by = "f.Distance", all.x = TRUE)

# subset for data between 2012 and 2021
snook2<-subset(snook2, Year > 2011 & Year < 2022)

saveRDS(snook2, file = "snook_SR_detection_20122021.rds")

############### START BELOW ######################################################################

snook2 <- readRDS("snook_SR_detection_20122021.rds")

# Calculate POR metrics for each individual -------------------------------

# line 1 = grouping individuals in data set by ID
# line 2 = summary statistics det_f = total number of detections; min.t = first date picked up on array; max.t = last day picked up on array; resid_t = number of days picked up in array; unique.date = number of year/month dates picked up in array

snook.detdays<-group_by(snook2, ID)%>%
      summarise(det_f = sum(!is.na(fDay)),
                min.t = min(Datetime_UTC),
                max.t = max(Datetime_UTC),
                resid.t = difftime(max.t, min.t, units = "day"),
                unique.date = length(unique(fYear.Month)))

#Subset to obtain individuals over 10 detections (based on Massie et al. in review)
### snook.10det <- subset(snook.detdays, det_f > 10)
### this only removes one individual from the data set - 208 inds to 207 indviduals

### decided to revise and use initial cutoff established by Rolo
#Subset to obtain individuals over 90 days and over 100 detections
snook.90day <- subset(snook.detdays, resid.t >= 90 & det_f > 100)

#Filter tag IDs based on cutoff information derived above
snook.90day <- filter(snook2, ID %in% snook.90day$ID)

#Creating table of unique tags per year for paper/ppt
snook.90day.x <- filter(snook.90day, Year < 2022)
length(unique(snook.90day.x$ID))

tags.summary <- group_by(snook.90day.x, Year)%>%
      summarise(n = length(unique(ID)))

list.tags <- group_by(snook.90day.x, ID, Year)%>%
      summarise(Total.Det = sum(!is.na(Datetime_UTC)),
                Total.Yr.Month = length(unique(fYear.Month)))

write.csv(tags.summary, "tags.summary_11_26_2022.csv")
write.csv(list.tags, "list.tags_11_26_2022.csv")


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

library(zoo)
#using function zoo::yearmon to classify in R fYear.Month as date formated as Year|Month
E_year.month2_df$fYear.Month <- as.yearmon(E_year.month2_df$Year.Month, "%Y/%m")

E_year.month2_df$Date2<-E_year.month2_df$Year.Month
E_year.month2_df<-cSplit(E_year.month2_df, "Date2", sep = "/", type.convert = FALSE)
colnames(E_year.month2_df)[4:5]<-c("Year", "Month")

E_year.month2_df$Year <- as.numeric(E_year.month2_df$Year)
E_year.month2_df$fYear<-factor(E_year.month2_df$Year, levels = c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021"))


E_year.month2_df$fMonth<-factor(E_year.month2_df$Month)
levels(E_year.month2_df$fMonth) <- list("Jan" = c("01"), "Feb" = c("02"), "Mar" = c("03"), 
                                        "Apr" = c("04"), "May" = c("05"), "Jun" = c("06"), 
                                        "Jul" = c("07"),"Aug" = c("08"), "Sep" = c("09"), 
                                        "Oct" = c("10"), "Nov" = c("11"), "Dec" = c("12"))

E_year.month2_df$Season<-factor(E_year.month2_df$fMonth)
levels(E_year.month2_df$Season)<-list("Dry" = c("Jan", "Feb", "Mar", "Apr", "May", "Nov", "Dec"), 
                                      "Wet" = c("Jun", "Jul", "Aug", "Sep", "Oct"))

### save things as is...
write.csv(E_year.month2_df, "Eadj_2012_2021_SRSnook_MW_11_26_2022.csv")
### read into excel and enter water year by hand... kind of confusing and dont know how to do using actual code

### in excel... added water years, included wYear 2011... since it contained e adjusted values corresponded to large # of isotope samples... however, we only have dry 2011 edaj values (four months of dry season values) so wYear 2011 will need to be removed for plotting... further, we only have e adj values for very beginning of wYear 2021...(i.e., wet season) -> completely removed these values because they will plot alone AND gives us no information for trophic questions... once again, remove wYear 2011 for plotting but keep values for future comparisons

### read in revised csv file with water year appropriately categorized
E_year.month2_df <- read.csv("Eadj_wYEAR_2011_2020_SRSnook_MW_11_26_2022.csv")
str(E_year.month2_df)
### make sure that wYear and season is a factor with levels
E_year.month2_df$wYear<-factor(E_year.month2_df$wYear, levels = c("2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020"))

E_year.month2_df$Season<-factor(E_year.month2_df$Season, levels = c("Wet", "Dry"))

Eadj_WETvDRY_Annual = E_year.month2_df %>% 
      group_by(wYear, Season) %>%
      summarize(Eadj_Seasonal_Mean = mean(Eadj), 
                Eadj_Seasonal_SD = sd(Eadj), 
                Eadj_Seasonal_Min = min(Eadj),
                Eadj_Seasonal_Max = max(Eadj))

write.csv(Eadj_WETvDRY_Annual, "Eadj_WETvDRY_Annual_11_26_2022.csv")

### PLOT SEASONAL EADJ VALUES (IE DRY VS WET)

tiff("/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript/WETvDRY_withWY2011_Eadj_MW_11_26_2022.tiff", width = 8, height = 8, units = 'in', res = 600, compression = 'lzw')

ggplot(E_year.month2_df, aes(x=Season, y=Eadj, fill=Season)) +
      geom_boxplot(width =0.8) +
      scale_fill_manual(values = c("cadetblue", "darkgoldenrod")) +
      labs(x = "Season", 
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

### PLOT SEASONAL EADJ VALUES (IE DRY VS WET BUT FOR EACH YEAR)

tiff("/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript/WETvDRY_withWY2011_Eadj_Annual_MW_11_26_2022.tiff", width = 16, height = 9, units = 'in', res = 600, compression = 'lzw')


ggplot(E_year.month2_df, aes(x=wYear, y=Eadj, fill=Season)) +
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

E_year.month2_df <- E_year.month2_df %>%
     filter(!row_number() %in% c(1,2,3,4))

### looks good! Got ride of wYear 2011 (really just dry 2011)

### PLOT SEASONAL EADJ VALUES (IE DRY VS WET)

tiff("/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript/WETvDRY_withoutWY2011_Eadj_MW_11_26_2022.tiff", width = 8, height = 8, units = 'in', res = 600, compression = 'lzw')

ggplot(E_year.month2_df, aes(x=Season, y=Eadj, fill=Season)) +
      geom_boxplot(width =0.8) +
      scale_fill_manual(values = c("cadetblue", "darkgoldenrod")) +
      labs(x = "Season", 
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

### PLOT SEASONAL EADJ VALUES (IE DRY VS WET BUT FOR EACH YEAR)

tiff("/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript/WETvDRY_withoutWY2011_Eadj_Annual_MW_11_26_2022.tiff", width = 16, height = 9, units = 'in', res = 600, compression = 'lzw')

ggplot(E_year.month2_df, aes(x=wYear, y=Eadj, fill=Season)) +
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


# NEXT STEPS WOULD BE TO LOOK AT MODELING CHANGES IN EADJ BETWEEN SEASONS AND AS IT RELATES TO STAGE HEIGHT
