# background --------------------------------------------------------------

#author: Mack White
#project: SRFCEA - Snook Spatial Similarity Manuscript
#goal of script: calculate hydrologic variables of interest
#date(s): January 2024


# load in packages --------------------------------------------------------
librarian::shelf(tidyverse, dataRetrieval, ggplot2, lubridate, 
                 reshape, reshape2, ggeffects)


#Bottle Creek Gage Height near Rookery Branch, Everglades NPS 
siteNumber <- "022908295" # found here https://waterdata.usgs.gov/monitoring-location/022908295/#parameterCode=00065&period=P7D
ChoptankInfo <- readNWISdata(022908295) # don't mind warning, or error, still works
parameterCd <- "00065" # parameter codes listed here -> http://water.nv.gov/hearings/past/Spring%20Valley%202006/exhibits/SNWA/5__Hydrochemistry_and_Geochemistry/Data/USGS/USGS_NWIS/ParameterCodes.htm

BC_GH_data <- readNWISdv(
      siteNumber, parameterCd,
      "2012-05-01", "2023-04-30"
)

head(BC_GH_data)

botcreekGH <- BC_GH_data |> 
      mutate(w_lev = X_00065_00003) |> 
      mutate(s_date = Date) |>
      mutate(w_lev_m = w_lev*0.3048) |> 
      select(s_date, w_lev, w_lev_m) 
glimpse(botcreekGH)
summary(botcreekGH)

dat <- botcreekGH |> 
      mutate(date = as.Date(s_date), format = "%Y-%m-%d",
             year = year(s_date),
             month = month(s_date),
             day = day(s_date))

dat$YearMonth <- paste(year(dat$date), 
                             month(dat$date, 
                                   label = FALSE, abbr = FALSE), 
                             sep = "-")
glimpse(dat)

# need to impute missing data?
# confirm using bottle creek? maybe use both bottle creek and MO-215
# need to correct from feet to cm

hydro <- dat |>
      group_by(YearMonth) |> 
      mutate(monthly_mean_wl = mean(w_lev)) |> 
      ungroup() |> 
      arrange(YearMonth) |> 
      group_by(YearMonth) |> 
      mutate(monthly_mean_wl_previous = lag(mean(w_lev))) |> 
      ungroup() |> 
      group_by(YearMonth) |> 
      mutate(DaysBelow30 = sum(w_lev < 30, na.rm = TRUE)) |> 
      ungroup()

# BC_gheight_plot <- ggplot(botcreekGH, aes(x=as.Date(s_date), y = as.numeric(w_lev_m), group = 1)) +
#       geom_line(color = "black", linewidth = 0.5) +
#       geom_point(size = 1.0) +
#       geom_smooth() +
#       # facet_wrap(~Site) +
#       labs(x = "Sample Date",
#            y = "Gage Height (m)",
#            title = "Marsh Water Levels") +
#       theme(panel.grid.major = element_blank(),
#             panel.background = element_blank(),
#             axis.line = element_line(colour = "black"),
#             plot.title = element_text(hjust = 0.5, size=14, face="bold", color = "black"),
#             # axis.text = element_text(size=12,face="bold", color = "black"),
#             axis.title = element_text(size=12,face="bold", color = "black"),
#             axis.text.x = element_text(angle = 45, hjust = 1., vjust = 1.1, face = "bold"),axis.text = element_text(color="black"),
#             panel.grid.minor = element_blank(),legend.position = "none") +
#       scale_x_date(date_labels = "%Y",breaks ='1 year')
# BC_gheight_plot

# ggsave(filename='lter cnd wg projects/plots/bottlecreek_GH.png',
#        plot = last_plot(),
#        scale = 2.5,
#        width = 11,
#        height = 5,
#        units = c("cm"),
#        dpi = 300)