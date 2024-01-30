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

# read in raw data from period of investigation
#dat <- read.csv("FCE_SI_data_xls_11_30_21_most_recent-PC.csv")
dat <- FCE_SI_data_xls_11_30_21_most_recent_PC_snookfood07262023 #new source file w/ only snook food and no missing data - 07262023 (406 rows)
rm(FCE_SI_data_xls_11_30_21_most_recent_PC_snookfood07262023)

# check structure
glimpse(dat)

### NO LONGER NEED CODE BELOW, SORTED BY HAND IN NEW EXCEL FILE
# # subset to only the information you need - source file for snook, so only need things they would be likely to eat during the dry season.. talking with ryan think they are feeding at the trophic level immediately below them, but will include small fishes and inverts in source file for bayesian mixing models
# 
# dat1 <- subset(dat, functional_grp %in% c("Benthic_invert", "Crab", 
#                                           "demersal_fish","shrimp")) %>%
#       subset(., select = c("site", "common_name", "species_name", "functional_grp",
#                            "d13C", "d15N", "d34S"))

# summarize data to get mean and sd for each source at each site
dat_summary <- dat |> 
      group_by(site, Source) |> 
      summarise(n = n(),
                Meand13C = mean(d13C),
                SDd13C = sd(d13C),
                Meand15N = mean(d15N),
                SDd15N = sd(d15N),
                Meand34S = mean(d34S),
                SDd34S = sd(d34S))

unique(dat_summary$site)

# recode habitat types w/ updated info 07262023 ---------------------------

dat$Source1 <- dat$site %>%
      recode(RB10 = "Freshwater", SRS3 = "Freshwater", SRS4 = "Estuarine", 
             SRS6 = "Estuarine", TS10 = "Seagrass", TS11 = "Seagrass", TS3 = "Freshwater", 
             TS7 = "Estuarine", TS9 = "Seagrass")

glimpse(dat)

# summarize data to get mean and sd for each source1 and source (habitat + species)
dat_summary1 <- dat %>%
      group_by(Source1, Source)%>%
      summarise(n = n(),
                Meand13C = mean(d13C),
                SDd13C = sd(d13C),
                Meand15N = mean(d15N),
                SDd15N = sd(d15N),
                Meand34S = mean(d34S),
                SDd34S = sd(d34S))
glimpse(dat_summary1)

# summarize data to get mean and sd for each source1 (habitat)
dat_summary2 <- dat |> 
      drop_na() |> 
      group_by(Source1) |> 
      summarise(n = n(),
                Meand13C = mean(d13C),
                SDd13C = sd(d13C),
                Meand15N = mean(d15N),
                SDd15N = sd(d15N),
                Meand34S = mean(d34S),
                SDd34S = sd(d34S))

glimpse(dat_summary2)

write.csv(dat_summary2, "ss_snook_source_agg_UPDATED_07262023.csv")
