
library(tidyverse)
library(lubridate)
library(viridis)
library(plotly)
library(ggmap)
library(writexl)
library(dplyr)
library(data.table)
library(hrbrthemes)

setwd("~/Desktop/RESEARCH/SPATIAL SYNCHRONY")

dat <- read.csv("shark river snook_w dates_June.csv")
View(dat)

dat <- dat[-c(124:150), ]

### determine monthly counts by year

datSNOOK = dat %>% 
  group_by(Month, Year) %>%
  summarize(samples = sum(Count))

datSNOOk <- datSNOOK[-c(124:150), ]

#write_xlsx(datSNOOK, "~/Desktop/RESEARCH/SPATIAL SYNCHRONY/datSNOOKmonthly.xlsx")

### load in water levels

m015_water <- read.csv("m015_2012_2022_waterlevels.csv")

### summarize monthly median levels by year, similiar to snook isotopes

m015_water_median = m015_water %>%
  group_by(Month, Year) %>%
  summarize(monthly_median_water_level = median(median_water_level))

### 

dat = dat[ ,c(2,6:10,12:17)] #wrong


