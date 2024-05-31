librarian::shelf(readr, ggplot2, car, reshape, reshape2, plyr, dplyr,
                 tidyr, visreg, modEvA, gridExtra, AICcmodavg, nlme, mgcv,
                 lme4, splitstackshape, chron, RInSp, boot, cowplot, ggpubr,
                 lsmeans, MixSIAR, hypervolume, truncnorm, tidyverse, alphahull,
                 dplyr, svMisc, writexl, stringr, zoo)

dat <- read_csv("data/snook_ss01302024_UpdatedwYear_2monthFiltered.csv")

dat1 <- read_csv("data/spat_sim_allthegoods_01_30_2024.csv")

ggplot(dat, aes(as.factor(year), mean, fill = source)) +
      geom_boxplot()

summ <- dat |> 
      group_by(year, source) |> 
      summarize(mean_source = mean(mean),
                sd = sd(mean),
                n = n()) 

wider_df <- summ |> 
      pivot_wider(
            names_from = source,
            values_from = c(mean_source, sd),
            names_sep = "_"
      )

wider_df_clean <- wider_df |> 
      select(year,
             mean_source_Freshwater, sd_Freshwater,
             mean_source_Estuarine, sd_Estuarine, 
             mean_source_Seagrass, sd_Seagrass,
             n) |> 
      mutate(across(starts_with("mean_source"), ~round(., 2))) |> 
      mutate(across(starts_with("sd"), ~round(., 2)))

write_csv(wider_df_clean, "tables/wider_df_clean_niche_table_05302024.csv")
