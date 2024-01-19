setwd("/Users/mack/Dropbox/R/bottlecreek_hydro")

library(ggplot2)
library(tidyverse)
library(lubridate)
library(viridis)
library(plotly)
library(ggmap)
library(xlsx)
library(dplyr)
library(data.table)
library(hrbrthemes)

dat <- read.csv("bottlecreek_waterlevels_02172023.csv")
glimpse(dat)

dat$date <- as.POSIXct(dat$date,"%Y-%m-%d",tz = "UTC")
glimpse(dat)
dat$monthday <- as.factor(dat$monthday)
dat$Year <- as.factor(dat$Year)

ggplot(dat, aes(x = monthday, y = level_cm_corrected, group = Year, colour = Year)) +
      geom_line()

dat_filtered <- dat %>%
      filter(date > '2015-12-31')

tiff("/Users/mack/Desktop/waterlevel_long.tiff", width = 15, height = 5, units = 'in', 
     res = 600, compression = 'lzw')

ggplot(dat, aes(x = date, y = level_cm_corrected, color = Year)) +
      geom_line(size = 1.5) +
      labs(x = "Year", y = "Bottle Creek (ENP) Daily Mean Water Level (cm)") +
      scale_y_continuous(breaks = seq(20, 120, by = 20))+
      theme(panel.background = element_blank()) +
      theme(axis.text = element_text(size=12,face="bold", color = "black")) +
      theme(axis.line = element_line(size = 1)) +
      theme(axis.title = element_text(size=12,face="bold", color = "black")) +
      theme(legend.position = "none") + 
      geom_hline(yintercept = 30.0, size = 1.0, linetype = 2, color = "black")

dev.off()
