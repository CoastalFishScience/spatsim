# background --------------------------------------------------------------

#author: Mack White
#project: SRFCEA - Snook Spatial Similarity Manuscript
#goal of script: update hypervolumes through 2023
#date(s): January 2024


# load in packages --------------------------------------------------------
librarian::shelf(tidyverse, hypervolume, MixSIAR, readr, 
                 writexl, ggpubr, ggthemes, ggplot2)

df = read_csv("data/snook_ss01292024_UpdatedwYear_3monthFiltered.csv") |> 
      pivot_wider(names_from = source,
                  values_from = mean) |> 
      mutate(across(Estuarine:Seagrass, scale)) |> 
      select(-ID) |> 
      group_by(year) |> 
      nest() |> 
      mutate(hv = map(data, ~hypervolume_gaussian(.x, name = paste(year), # change to label columns in data
                                                  samples.per.point = 1000,
                                                  kde.bandwidth = estimate_bandwidth(.x), 
                                                  sd.count = 3, 
                                                  quantile.requested = 0.95, 
                                                  quantile.requested.type = "probability", 
                                                  chunk.size = 1000, 
                                                  verbose = F)),
             hv_size = get_volume(hv[[1]]))

head(df)

df$hv_size

df |> 
      select(year, hv_size) |> 
      write_csv("hv_sss_01082024_UpdatedwYear.csv")

dat <- read_csv("hv_sss_01082024.csv") |> 
      rename(wYear = year) |> 
      left_join(mix_summary) |> 
      write_csv("hv_size_with_sample_size_01082024_UpdatedwYear.csv")

dat <- read_csv("hv_size_with_sample_size_01082024_UpdatedwYear.csv")

ggplot(dat, aes(wYear, hv_size))+
      geom_point(size = 2.5)+
      geom_line(linewidth = 1)+
      labs(x = 'Year', y = 'Volume')+
      theme_bw()+
      #scale_y_log10()+
      #scale_color_viridis_d()+
      theme(axis.title = element_text(size = 14), 
            axis.text = element_text(size = 14, colour = "gray0"), 
            plot.title = element_text(size = 14, hjust=0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'right',
            legend.title = element_text(size = 14),
            strip.text.x = element_text(size = 14),
            legend.text = element_text(size = 12))

ggplot(dat |> filter(!wYear%in%c(2011,2020)), aes(wYear, hv_size))+
      geom_point(size = 2.5)+
      geom_smooth(method = "lm")+
      # geom_line(linewidth = 1)+
      labs(x = 'Year', y = 'Volume')+
      theme_bw()+
      #scale_y_log10()+
      #scale_color_viridis_d()+
      theme(axis.title = element_text(size = 14), 
            axis.text = element_text(size = 14, colour = "gray0"), 
            plot.title = element_text(size = 14, hjust=0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'right',
            legend.title = element_text(size = 14),
            strip.text.x = element_text(size = 14),
            legend.text = element_text(size = 12))

#dont take out 2020? Just a big outlier - keep in for analysis of e for dry season and one without
#

ggplot(dat|> filter(!wYear%in%c(2011)), aes(wYear, hv_size))+
      geom_point(size = 2.5)+
      geom_smooth(method = "lm")+
      # geom_line(linewidth = 1)+
      labs(x = 'Year', y = 'Volume')+
      theme_bw()+
      #scale_y_log10()+
      #scale_color_viridis_d()+
      theme(axis.title = element_text(size = 14), 
            axis.text = element_text(size = 14, colour = "gray0"), 
            plot.title = element_text(size = 14, hjust=0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'right',
            legend.title = element_text(size = 14),
            strip.text.x = element_text(size = 14),
            legend.text = element_text(size = 12))

ggplot(dat|> filter(!wYear%in%c(2011,2012,2019,2021,2022)), aes(wYear, hv_size))+
      geom_point(size = 2.5)+
      geom_smooth(method = "lm")+
      # geom_line(linewidth = 1)+
      labs(x = 'Year', y = 'Volume')+
      theme_bw()+
      #scale_y_log10()+
      #scale_color_viridis_d()+
      theme(axis.title = element_text(size = 14), 
            axis.text = element_text(size = 14, colour = "gray0"), 
            plot.title = element_text(size = 14, hjust=0.5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = 'right',
            legend.title = element_text(size = 14),
            strip.text.x = element_text(size = 14),
            legend.text = element_text(size = 12))
