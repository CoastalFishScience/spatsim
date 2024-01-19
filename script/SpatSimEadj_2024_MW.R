# background --------------------------------------------------------------

###author: Mack White
###project: SRFCEA - Snook Spatial Similarity Manuscript
###goal of script: Calculate Eadj Values for Common Snook based on space use
###date(s): January 2024
###notes:
# go to line 226 to start Monday, Jan 15

# load in packages --------------------------------------------------------

librarian::shelf(readr, ggplot2, car, reshape, reshape2, plyr, dplyr,
                 tidyr, visreg, modEvA, gridExtra, AICcmodavg, nlme, mgcv,
                 lme4, splitstackshape, chron, RInSp, boot, cowplot, ggpubr,
                 lsmeans, MixSIAR, hypervolume, truncnorm, tidyverse, alphahull,
                 dplyr, svMisc, writexl, stringr, zoo)

# pre-processing and manipulation -----------------------------------------

tracks <- read_rds("~/Library/CloudStorage/Dropbox/R/github/spatsim/data/RehageVUEdatabase_10182023.rds")
glimpse(tracks)
summary(tracks)

tracks$Station.Name <- as.character(tracks$Station.Name)
tracks$Station <- tracks$Station.Name
tracks$VUE_Name <- tracks$Station

distance <- read_csv("data/Station_Distance_Updated07142020.csv")
glimpse(distance)
summary(distance)

tracks_dist <- left_join(tracks, distance,
                         by = "VUE_Name")

snook <- tracks_dist |> 
      filter(Species == "Snook")

# separate datetime column into date and time columns
snook$DateTime2 <- snook$Datetime_UTC
snook<-cSplit(snook, "DateTime2", sep = " ", type.convert = FALSE)
colnames(snook)[16:17]<-c("Date", "Time")

# separate date column into year, month, and day columns
snook$Date2 <- snook$Date
snook<-cSplit(snook, "Date2", sep = "-", type.convert = FALSE)
colnames(snook)[18:20]<-c("Year", "Month", "Day")

# separate the transmitter column into code1, serial, and ID columns
snook<-cSplit(snook, "Transmitter", sep = "-", type.convert = FALSE)
# colnames(snook)[21:23]<-c("code1", "serial", "ID")

snook$ID <- snook$Transmitter_3

# separate receiver column into code and receiver id columns
snook<-cSplit(snook, "Receiver", sep = "-", type.convert = FALSE)
# colnames(snook)[24:25]<-c("code", "Receiver.ID")

snook$Receiver.ID <- snook$Receiver_2

# generate categorical variable to delineate lower river vs bay vs upper river detections - name "Home"
snook$Home<-cut(snook$Distance, breaks = c(-1,15,23,31.72), 
                labels = c("Lower River", "Bay", "Upper River"))
glimpse(snook)

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

zones<- read_csv("data/zones2019.csv")
zones$f.Distance<-factor(zones$Distance)
zones<-subset(zones, Zone != "NA")
snook2<-merge(snook2, zones, by = "f.Distance", all.x = TRUE)

# subset for data between 2012 and 2023
snook2<-subset(snook2, Year > 2011 & Year < 2024)

saveRDS(snook2, file = "spatsim_2012thru2023.rds")

############### START BELOW ######################################################################

snook2 <- readRDS("spatsim_2012thru2023.rds")

# Calculate POR metrics for each individual -------------------------------

# line 1 = grouping individuals in data set by ID
# line 2 = summary statistics det_f = total number of detections; min.t = first date picked up on array; max.t = last day picked up on array; resid_t = number of days picked up in array; unique.date = number of year/month dates picked up in array

snook.detdays<-group_by(snook2, ID) |> 
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
snook.90day.x <- filter(snook.90day, Year < 2024)
length(unique(snook.90day.x$ID))

snook_obs <- snook.90day.x |> 
      group_by(fYear.Month) |> 
      summarise(n = n_distinct(ID))

# write_csv(snook_obs, "snook_sample_size_eadj.csv")

tags.summary <- group_by(snook.90day.x, Year)%>%
      summarise(n = length(unique(ID)))

list.tags <- group_by(snook.90day.x, ID, Year)%>%
      summarise(Total.Det = sum(!is.na(Datetime_UTC)),
                Total.Yr.Month = length(unique(fYear.Month)))

# write.csv(tags.summary, "tags.summary_spatsim_01082024.csv")
# write.csv(list.tags, "list.tags_spatsim_01082024.csv")


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
E_year.month2_df<-ldply(E_year.month2, data.frame)
colnames(E_year.month2_df)<-c("Year.Month", "Eadj")

#using function zoo::yearmon to classify in R fYear.Month as date formated as Year|Month
E_year.month2_df$fYear.Month <- as.yearmon(E_year.month2_df$Year.Month, "%Y/%m")

E_year.month2_df$Date2<-E_year.month2_df$Year.Month
E_year.month2_df<-cSplit(E_year.month2_df, "Date2", sep = "/", type.convert = FALSE)
colnames(E_year.month2_df)[4:5]<-c("Year", "Month")

E_year.month2_df$Year <- as.numeric(E_year.month2_df$Year)
E_year.month2_df$fYear<-factor(E_year.month2_df$Year, levels = c("2012", "2013", "2014", "2015", "2016", 
                                                                 "2017", "2018", "2019", "2020", "2021",
                                                                 "2022", "2023"))


E_year.month2_df$fMonth<-factor(E_year.month2_df$Month)
levels(E_year.month2_df$fMonth) <- list("Jan" = c("01"), "Feb" = c("02"), "Mar" = c("03"), 
                                        "Apr" = c("04"), "May" = c("05"), "Jun" = c("06"), 
                                        "Jul" = c("07"),"Aug" = c("08"), "Sep" = c("09"), 
                                        "Oct" = c("10"), "Nov" = c("11"), "Dec" = c("12"))

E_year.month2_df$Season<-factor(E_year.month2_df$fMonth)
levels(E_year.month2_df$Season)<-list("Dry" = c("Jan", "Feb", "Mar", "Apr", "Nov", "Dec"), 
                                      "Wet" = c("May", "Jun", "Jul", "Aug", "Sep", "Oct"))

### save things as is...
# write.csv(E_year.month2_df, "Eadj_2012_2023_SRSnook_MW_01_08_2024.csv")
### read into excel and enter water year and month by hand... 
### simpler than coding in
### dropped wYear 2011 & 2023 as they only had a handful of observations

### read in revised csv file with water year appropriately categorized
E_year.month2_df <- read.csv("Eadj_wYear_2012_2023_SRSnook_MW_01_08_2024.csv")
str(E_year.month2_df)
### make sure that wYear and season is a factor with levels
E_year.month2_df$wYear<-factor(E_year.month2_df$wYear, levels = c("2012", "2013", "2014", 
                                                                  "2015", "2016", "2017", "2018", 
                                                                  "2019", "2020", "2021", "2022"))

E_year.month2_df$Season<-factor(E_year.month2_df$Season, levels = c("Wet", "Dry"))

Eadj_WETvDRY_Annual = E_year.month2_df %>% 
      group_by(wYear, Season) %>%
      summarize(Eadj_Seasonal_Mean = mean(Eadj), 
                Eadj_Seasonal_SD = sd(Eadj), 
                Eadj_Seasonal_Min = min(Eadj),
                Eadj_Seasonal_Max = max(Eadj))

# write.csv(Eadj_WETvDRY_Annual, "Eadj_WETvDRY_Annual_01_08_2024.csv")

# E_year.month2_df <- read_csv("Eadj_WETvDRY_Annual_01_07_2024.csv")

################################################################################
################################################################################

E_year.month2_df <- read_csv("Eadj_wYear_2012_2023_SRSnook_MW_01_08_2024.csv")

E_year.month2_df$Season<-factor(E_year.month2_df$Season, levels = c("Wet", "Dry"))

E_year.month2_df$wYear<-factor(E_year.month2_df$wYear, levels = c("2012", "2013", "2014", 
                                                                  "2015", "2016", "2017", "2018", 
                                                                  "2019", "2020", "2021", "2022"))
### PLOT SEASONAL EADJ VALUES (IE DRY VS WET)

plot1 <- ggplot(E_year.month2_df, aes(x=Season, y=Eadj, fill=Season)) +
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

# ggsave(filename='figures/WETvDRY_Eadj_MW_01_08_2024.png',
#        plot = last_plot(),
#        scale = 2.5,
#        width = 9,
#        height = 5,
#        units = c("cm"),
#        dpi = 300)


### PLOT SEASONAL EADJ VALUES (IE DRY VS WET BUT FOR EACH YEAR)

plot2 <- ggplot(E_year.month2_df, aes(x=wYear, y=Eadj, fill=Season)) +
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

# ggsave(filename='figures/WETvDRY_SEASONAL_Eadj_Annual_MW_01_08_2024.png',
#        plot = last_plot(),
#        scale = 2.5,
#        width = 9,
#        height = 5,
#        units = c("cm"),
#        dpi = 300)

dat <- E_year.month2_df |> 
      mutate(date = ym(Year.Month))

plot3 <- ggplot(dat, aes(date, Eadj)) +
      geom_point(size = 2) +
      facet_wrap(~wYear, scales = "free_x")+
      geom_smooth()

plot4 <- ggplot(dat, aes(date, Eadj)) +
      geom_point(size = 3, aes(color = Season)) +
      # facet_wrap(~wYear, scales = "free_x")+
      geom_path()

# plot4 <- ggplot(dat, aes(date, Eadj)) +
#       geom_point(size = 3, aes(color = Season)) +
#       # facet_wrap(~wYear, scales = "free_x")+
#       geom_path()

# ggsave(filename='figures/facet_eadj_smoothed_MW_01_08_2024.png',
#        plot = last_plot(),
#        scale = 2.5,
#        width = 12,
#        height = 5,
#        units = c("cm"),
#        dpi = 300)

dat$date_num <- as.numeric(dat$date)
glimpse(dat) #makes date numeric for ts 

# m1 <- gam(Eadj ~ s(date, bs = "cc", k = 12, by = wYear) + s(wYear, bs = "cr") + Season, 
#           data = dat,
#           family = betar(link = "logit"), 
#           method = "REML")
# 
# summary(m1)
# 
# all$wYearFACT <- as.factor(all$wYear)

# all_test = all |> 
#       mutate(wMonth = if_else(
#             Month |> 
#                   between(6,12),
#             Month - 5,
#             Month + 7
#       ))

# move May to Wet season and will need to fix wMonth ifelse statement
# between 5 and 12 and subtract 4, adding 8

### below is correct

# all_test = all |>
#       mutate(wMonth = if_else(
#             Month |>
#                   between(5,12),
#             Month - 4,
#             Month + 8
#       ))

#also going to need to change wYear
# 
# m2 <- gam(Eadj ~ s(date_num, bs = "cr", k = 12, by = wYearFACT) + s(wYear, bs = "cr") + Season, 
#           data = all,
#           family = betar(link = "logit"), 
#           method = "REML")
# 
# summary(m2)
# 
# ggplot(all, aes(date_num, Eadj)) +
#       geom_point(size = 2) +
#       facet_wrap(~wYearFACT, scales = "free_x")+
#       geom_smooth()
# 
# ggplot(all, aes(date, Eadj)) +
#       geom_point(size = 2) +
#       facet_wrap(~wYear, scales = "free_x")+
#       geom_smooth()

# m3 <- gam(Eadj ~ s(wMonth, bs = "cr", k = 12, by = wYear) + s(date_num, bs = "cr") + Season, 
#           data = all,
#           family = betar(link = "logit"), 
#           method = "REML")

### model 3 - we liked last meeting

m3 <- gam(Eadj ~ s(wMonth, bs = "cc", k = 12) + s(date_num), 
          data = dat,
          family = betar(link = "logit"), 
          method = "REML")

summary(m3) #deviance explained increased after revising wYear/wMonth with May in Wet
plot(m3)

### model 4 - we also liked last meeting

m4 <- gam(Eadj ~ s(wMonth, bs = "cc", k = 12) + s(wYear), 
          data = dat,
          family = betar(link = "logit"), 
          method = "REML")

summary(m4)
plot(m4)

all_test$wYear.Month = all_test$wYear + (all_test$wMonth/12) - (1/12)
glimpse(all_test)

dat <- read_csv("spat_sim_allthegoods_01_15_2023.csv")

# meeting next monday
# to-do 
# get ready the dataset of variables for models
# reach out to dr massie -> condition manuscript stuff
# set up with May in Wet Season
# link hydro data with eadj stuff
### average stage at MO215, previous month water stage, # of days <30 cm monthly
# finish HVs

dat <- read_csv("spat_sim_allthegoods_01_15_2023.csv")

###############################################################################
###############################################################################
###############################################################################

# Revisions for Monday (January 15) ---------------------------------------

# - summary stats by dry season for eadj size & hypervolume into new dataset
# - # tack on # of inds that went into each of the hyperpvolume calculations.
# - # look at correlation between eadj and sample size

dat <- read_csv("spat_sim_allthegoods_01_15_2023.csv") 
glimpse(dat)

snook_obs <- read_csv("snook_sample_size_eadj.csv") |> 
      rename(Year.Month = fYear.Month,
             eadj_n = n)
glimpse(snook_obs)

dat_all <- left_join(dat, snook_obs, by = "Year.Month")
glimpse(dat_all)

# run correlation here for eadj and eadj sample size kjhkjg

#colored by year
ggplot(dat_all, aes(x = eadj_n, y = Eadj, color = wYear)) + 
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
ggplot(dat_all, aes(x = eadj_n, y = Eadj, color = monthly_mean_stage_cm)) + 
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
ggplot(dat_all, aes(x = eadj_n, y = Eadj, color = Season)) + 
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

dat_summary <- dat |> 
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
