# background --------------------------------------------------------------

#author: Mack White
#project: SRFCEA - Snook Spatial Similarity Manuscript
#goal of script: calculate hydrologic variables of interest
#date(s): January 2024


# load in packages --------------------------------------------------------
librarian::shelf(tidyverse, dataRetrieval, ggplot2, lubridate, 
                 reshape, reshape2, ggeffects)


# #Bottle Creek Gage Height near Rookery Branch, Everglades NPS 
# siteNumber <- "022908295" # found here https://waterdata.usgs.gov/monitoring-location/022908295/#parameterCode=00065&period=P7D
# ChoptankInfo <- readNWISdata(022908295) # don't mind warning, or error, still works
# parameterCd <- "00065" # parameter codes listed here -> http://water.nv.gov/hearings/past/Spring%20Valley%202006/exhibits/SNWA/5__Hydrochemistry_and_Geochemistry/Data/USGS/USGS_NWIS/ParameterCodes.htm
# 
# BC_GH_data <- readNWISdv(
#       siteNumber, parameterCd,
#       "2012-05-01", "2023-04-30"
# )
# 
# head(BC_GH_data)
# 
# botcreekGH <- BC_GH_data |> 
#       mutate(w_lev = X_00065_00003) |> 
#       mutate(s_date = Date) |>
#       mutate(w_lev_m = w_lev*0.3048) |> 
#       select(s_date, w_lev, w_lev_m) 
# glimpse(botcreekGH)
# summary(botcreekGH)
#

# read in csv file which was corrected for elevation on previous project
mo215 <- read_csv("data/mo215_elevation_corrected_water_levels.csv") |> 
      select(-season)
glimpse(mo215)

# specify the date range for analysis
start_date <- as.Date("2012-05-01") #start date
end_date <- as.Date("2023-04-30") #end date

# 
mo215 <- mo215 %>%
      filter(date >= start_date, date <= end_date) |> # filters dates
      mutate(year = year(date),
             month = month(date),
             day = day(date)) #generates year, month, and day column

glimpse(mo215)

mo215$YearMonth <- paste(year(mo215$date),
                             month(mo215$date,
                                   label = FALSE, abbr = FALSE),
                             sep = "-") #create YearMonth column

mo215$YearMonth_Date <- as.Date(paste(mo215$YearMonth, "-01", sep = ""), 
                                 format = "%Y-%m-%d")

glimpse(mo215)

hydro <- mo215 |>
      group_by(YearMonth) |> 
      mutate(monthly_mean_stage_cm = mean(stage_cm)) |> 
      ungroup() |> 
      # arrange(YearMonth) |> 
      # group_by(YearMonth) |> 
      # mutate(monthly_mean_stage_cm_prev = lag(monthly_mean_stage_cm)) |> 
      # ungroup() |> 
      group_by(YearMonth) |> 
      mutate(DaysBelow30 = sum(stage_cm < 30, na.rm = TRUE)) |> 
      ungroup()

hydro_sum <- hydro |> 
      group_by(YearMonth_Date) |> 
      summarise(daysbelow30 = mean(DaysBelow30),
                monthly_mean_stage_cm = mean(monthly_mean_stage_cm))

hydro_sum <- hydro_sum |> 
      arrange(YearMonth_Date) |> 
      group_by(YearMonth_Date) |> 
      mutate(monthly_mean_stage_cm_prev = lag(monthly_mean_stage_cm))
      

write_csv(hydro_sum, "data/mo215_hydro_summary_01_11_2024.csv")
