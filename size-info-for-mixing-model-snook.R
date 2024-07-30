### goals: generate mean and sd for isotope sizes

librarian::shelf(tidyverse, readxl, readr, writexl)

dat <- read_csv("data/FCE_master_MAP_AT_Isotopes_FORSIZES_June302024.csv") |> 
      filter(species_name == "Snook")
glimpse(dat)

