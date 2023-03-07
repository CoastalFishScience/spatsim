# Author: Mack White
# Date: 11/23/2022
# load libraries:
library(ggplot2)
library(car)
library(reshape2)
library(reshape)
library(plyr)
library(dplyr)
library(tidyr)
# set working directory
setwd("/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript")
# code settings -> rainbow parentheses + soft wrap long lines

# read in most recent raw data from MS Teams
dat <- read.csv("FCE_SI_data_xls_11_30_21_most_recent-PC.csv")
dat <- na.omit(dat)

# check structure
str(dat)

# subset to only the information you need - source file for snook, so only need things they would be likely to eat during the dry season.. talking with ryan think they are feeding at the trophic level immediately below them, but will include small fishes and inverts in source file for bayesian mixing models

dat1 <- subset(dat, functional_grp %in% c("Benthic_invert", "Crab", 
                                          "demersal_fish","shrimp")) %>%
      subset(., select = c("site", "common_name", "species_name", "functional_grp",
                           "d13C", "d15N", "d34S"))

# rename common name to source, what we will need in the end
dat1 <- rename(dat1, Source = common_name)
str(dat1)

# summarize data to get mean and sd for each source at each site
dat_summary <- dat1 %>%
      group_by(site, Source)%>%
      summarise(n = n(),
                Meand13C = mean(d13C),
                SDd13C = sd(d13C),
                Meand15N = mean(d15N),
                SDd15N = sd(d15N),
                Meand34S = mean(d34S),
                SDd34S = sd(d34S))

unique(dat_summary$site)

### group by freshwater (i.e., RB10, SRS3), estuarine (i.e., SRS4, SRS6), and seagrass/bay (i.e., all TS sites)

dat1$Source1 <- dat1$site %>%
      recode(RB10 = "Freshwater", SRS3 = "Freshwater", SRS4 = "Estuarine", 
             SRS6 = "Estuarine", TS10 = "Seagrass", TS11 = "Seagrass", TS3 = "Seagrass", 
             TS7 = "Seagrass", TS9 = "Seagrass")

# summarize data to get mean and sd for each source1 and source (habitat + species)
dat_summary1 <- dat1 %>%
      group_by(Source1, Source)%>%
      summarise(n = n(),
                Meand13C = mean(d13C),
                SDd13C = sd(d13C),
                Meand15N = mean(d15N),
                SDd15N = sd(d15N),
                Meand34S = mean(d34S),
                SDd34S = sd(d34S))

# summarize data to get mean and sd for each source1 (habitat)
dat_summary2 <- dat1 %>%
      group_by(Source1)%>%
      summarise(n = n(),
                Meand13C = mean(d13C),
                SDd13C = sd(d13C),
                Meand15N = mean(d15N),
                SDd15N = sd(d15N),
                Meand34S = mean(d34S),
                SDd34S = sd(d34S))

write.csv(dat_summary2, "ss_snook_source_agg_UPDATED.csv")
