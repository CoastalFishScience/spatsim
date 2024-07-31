### goals: generate mean and sd for isotope sizes

librarian::shelf(tidyverse, readxl, readr, writexl)

dat <- read_csv("data/FCE_master_MAP_AT_Isotopes_FORSIZES_June302024.csv") |> 
      ### filter for snookies
      filter(species_name == "Snook") |> 
      ### format dat for filtering of seasonal-lag isotope data
      mutate(Date = as.Date(Date)) |> 
      select(-Year) |> 
      janitor::clean_names() |> 
      mutate(year = year(date),
             month = month(date)) |> 
      ### generate water year information
      # wyear = how we did it/how historically thought about it: align with SpatSimMixing Models_2024_MW.R datasets
      # usace_wyear = how we think about it since SSR
      mutate(wyear = if_else(month <= 4, year-1, year),
             usace_wyear = if_else(month >= 5, year+1, year)) |> 
      ### filtering out months of interest: align with SpatSimMixing Models_2024_MW.R script
      filter(month %in% c(1:6))

size_summary <- dat |> 
      group_by(wyear) |> 
      summarize(tlcm_m = mean(tl_cm, na.rm = TRUE),
                tlcm_sd = sd(tl_cm, na.rm = TRUE))
             