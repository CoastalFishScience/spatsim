df <- read.csv("df_with_year_sss_01082024_correct.csv")
glimpse(df)

df$type <- as.character(df$type)
str(df)
# number or iterations
n = 20

# create tibble to store data
df_hvAll = tibble(rep = rep(seq(1,n,1)),
                  vol_2011 = NA, vol_2012 = NA, vol_2013 = NA, vol_2014 = NA,
                  vol_2016 = NA, vol_2017 = NA, vol_2018 = NA, vol_2019 = NA, 
                  vol_2020 = NA, vol_2021 = NA, vol_2022 = NA)


# run hypervolumes --------------------------------------------------------

for (i in 1:n){
      # generate random from mean and sd that are z-scored
      ss = HVvalues(df = df, ID_rows = c('name','type'),names = c('source'), 
                    mean = 'mean', sd = 'sd', n = 20,
                    end_points = T, low = 'lowend', high = 'highend', 
                    z_score = T)
      #2011 hv
      # subset data
      df1 = ss %>% 
            filter(type == "2011") %>%
            select(-name, -type)
      
      # generate hypervolume
      hv_2011 = hypervolume_gaussian(df1, name = '2011',
                                     #samples.per.point = ceiling((10^(3 + sqrt(ncol(wet_df))))/nrow(wet_df)),
                                     samples.per.point = 5000,
                                     kde.bandwidth = estimate_bandwidth(df1), 
                                     sd.count = 3, 
                                     quantile.requested = 0.95, 
                                     quantile.requested.type = "probability", 
                                     chunk.size = 1000, 
                                     verbose = F)
      
      # get volume
      df_hvAll$vol_2011[i] = get_volume(hv_2011)
      
      #2012 hv
      # subset data
      df2 = ss %>% 
            filter(type == "2012") %>%
            select(-name, -type)
      
      # generate hypervolume
      hv_2012 = hypervolume_gaussian(df2, name = '2012',
                                     #samples.per.point = ceiling((10^(3 + sqrt(ncol(wet_df))))/nrow(wet_df)),
                                     samples.per.point = 5000,
                                     kde.bandwidth = estimate_bandwidth(df2), 
                                     sd.count = 3, 
                                     quantile.requested = 0.95, 
                                     quantile.requested.type = "probability", 
                                     chunk.size = 1000, 
                                     verbose = F)
      
      # get volume
      df_hvAll$vol_2012[i] = get_volume(hv_2012)
      
      #2013 hv
      # subset data
      df3 = ss %>% 
            filter(type == "2013") %>%
            select(-name, -type)
      
      # generate hypervolume
      hv_2013 = hypervolume_gaussian(df3, name = '2013',
                                     #samples.per.point = ceiling((10^(3 + sqrt(ncol(wet_df))))/nrow(wet_df)),
                                     samples.per.point = 5000,
                                     kde.bandwidth = estimate_bandwidth(df3), 
                                     sd.count = 3, 
                                     quantile.requested = 0.95, 
                                     quantile.requested.type = "probability", 
                                     chunk.size = 1000, 
                                     verbose = F)
      
      # get volume
      df_hvAll$vol_2013[i] = get_volume(hv_2013)
      
      #2014 hv
      # subset data
      df4 = ss %>% 
            filter(type == "2014") %>%
            select(-name, -type)
      
      # generate hypervolume
      hv_2014 = hypervolume_gaussian(df4, name = '2014',
                                     #samples.per.point = ceiling((10^(3 + sqrt(ncol(wet_df))))/nrow(wet_df)),
                                     samples.per.point = 5000,
                                     kde.bandwidth = estimate_bandwidth(df4), 
                                     sd.count = 3, 
                                     quantile.requested = 0.95, 
                                     quantile.requested.type = "probability", 
                                     chunk.size = 1000, 
                                     verbose = F)
      
      # get volume
      df_hvAll$vol_2014[i] = get_volume(hv_2014)
      
      #2016 hv
      # subset data
      df5 = ss %>% 
            filter(type == "2016") %>%
            select(-name, -type)
      
      # generate hypervolume
      hv_2016 = hypervolume_gaussian(df5, name = '2016',
                                     #samples.per.point = ceiling((10^(3 + sqrt(ncol(wet_df))))/nrow(wet_df)),
                                     samples.per.point = 5000,
                                     kde.bandwidth = estimate_bandwidth(df5), 
                                     sd.count = 3, 
                                     quantile.requested = 0.95, 
                                     quantile.requested.type = "probability", 
                                     chunk.size = 1000, 
                                     verbose = F)
      
      # get volume
      df_hvAll$vol_2016[i] = get_volume(hv_2016)
      
      #2017 hv
      # subset data
      df6 = ss %>% 
            filter(type == "2017") %>%
            select(-name, -type)
      
      # generate hypervolume
      hv_2017 = hypervolume_gaussian(df6, name = '2017',
                                     #samples.per.point = ceiling((10^(3 + sqrt(ncol(wet_df))))/nrow(wet_df)),
                                     samples.per.point = 5000,
                                     kde.bandwidth = estimate_bandwidth(df6), 
                                     sd.count = 3, 
                                     quantile.requested = 0.95, 
                                     quantile.requested.type = "probability", 
                                     chunk.size = 1000, 
                                     verbose = F)
      
      # get volume
      df_hvAll$vol_2017[i] = get_volume(hv_2017)
      
      #2018 hv
      # subset data
      df7 = ss %>% 
            filter(type == "2018") %>%
            select(-name, -type)
      
      # generate hypervolume
      hv_2018 = hypervolume_gaussian(df7, name = '2018',
                                     #samples.per.point = ceiling((10^(3 + sqrt(ncol(wet_df))))/nrow(wet_df)),
                                     samples.per.point = 5000,
                                     kde.bandwidth = estimate_bandwidth(df7), 
                                     sd.count = 3, 
                                     quantile.requested = 0.95, 
                                     quantile.requested.type = "probability", 
                                     chunk.size = 1000, 
                                     verbose = F)
      
      # get volume
      df_hvAll$vol_2018[i] = get_volume(hv_2018)
      
      #2019 hv
      # subset data
      df8 = ss %>% 
            filter(type == "2019") %>%
            select(-name, -type)
      
      # generate hypervolume
      hv_2019 = hypervolume_gaussian(df8, name = '2019',
                                     #samples.per.point = ceiling((10^(3 + sqrt(ncol(wet_df))))/nrow(wet_df)),
                                     samples.per.point = 5000,
                                     kde.bandwidth = estimate_bandwidth(df8), 
                                     sd.count = 3, 
                                     quantile.requested = 0.95, 
                                     quantile.requested.type = "probability", 
                                     chunk.size = 1000, 
                                     verbose = F)
      
      # get volume
      df_hvAll$vol_2019[i] = get_volume(hv_2019)
      
      #2020 hv
      # subset data
      df9 = ss %>% 
            filter(type == "2020") %>%
            select(-name, -type)
      
      # generate hypervolume
      hv_2020 = hypervolume_gaussian(df9, name = '2020',
                                     #samples.per.point = ceiling((10^(3 + sqrt(ncol(wet_df))))/nrow(wet_df)),
                                     samples.per.point = 5000,
                                     kde.bandwidth = estimate_bandwidth(df9), 
                                     sd.count = 3, 
                                     quantile.requested = 0.95, 
                                     quantile.requested.type = "probability", 
                                     chunk.size = 1000, 
                                     verbose = F)
      
      # get volume
      df_hvAll$vol_2020[i] = get_volume(hv_2020)
      
      #2021 hv
      # subset data
      df10 = ss %>% 
            filter(type == "2021") %>%
            select(-name, -type)
      
      # generate hypervolume
      hv_2021 = hypervolume_gaussian(df10, name = '2021',
                                     #samples.per.point = ceiling((10^(3 + sqrt(ncol(wet_df))))/nrow(wet_df)),
                                     samples.per.point = 5000,
                                     kde.bandwidth = estimate_bandwidth(df10), 
                                     sd.count = 3, 
                                     quantile.requested = 0.95, 
                                     quantile.requested.type = "probability", 
                                     chunk.size = 1000, 
                                     verbose = F)
      
      # get volume
      df_hvAll$vol_2021[i] = get_volume(hv_2021)
      
      #2022 hv
      # subset data
      df11 = ss %>% 
            filter(type == "2022") %>%
            select(-name, -type)
      
      # generate hypervolume
      hv_2022 = hypervolume_gaussian(df3, name = '2022',
                                     #samples.per.point = ceiling((10^(3 + sqrt(ncol(wet_df))))/nrow(wet_df)),
                                     samples.per.point = 5000,
                                     kde.bandwidth = estimate_bandwidth(df11), 
                                     sd.count = 3, 
                                     quantile.requested = 0.95, 
                                     quantile.requested.type = "probability", 
                                     chunk.size = 1000, 
                                     verbose = F)
      
      # get volume
      df_hvAll$vol_2022[i] = get_volume(hv_2022)
      
      # show done
      require(svMisc)
      for(i in 0:101) {
            progress(i, progress.bar = TRUE)
            Sys.sleep(0.01)
            if (i == 101) cat ("Done!\n")
      }
}
# save that sucker
write_csv(df_hvAll, "/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript/snook_similarity_iterations_1126022.csv")

#### PRODUCT NUMBER ONE = BOOTSTRAPPED HYPERVOLUME FIGURE
hvBS= hypervolume_join(hv_2011, hv_2012, hv_2013, hv_2014, hv_2016,
                       hv_2017, hv_2018, hv_2019, hv_2020, hv_2021, hv_2022)

tiff("~hvfigureone_11262022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

plot(hvBS, pairplot = T,
     names= c(expression(italic("Seagrass"),italic("Estuarine"), italic("Freshwater"))),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-6,6), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)
dev.off()

# FORMAT AND BROUGHT BACK IN FOR GGPLOT FIGURE #2

long_df <- read.csv("snook_similarity_long_11262022.csv")

long_df$Year<-factor(long_df$Year, levels = c("2011", "2012", "2013", "2014", 
                                              "2016", "2017", "2018", "2019", "2020"))

str(long_df)
#### PRODUCT NUMBER TWO = PLOT OF NICHE SIZE (VIOLIN PLOT?)

################################################################################

tiff("wYearVolByYear_11262022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggplot(long_df, aes(x=Year, y=vol, group=Year)) +
      geom_boxplot(width =0.8) +
      labs(x = "Water Year", 
           y = "Niche Volume") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=14, face="bold", color = "black")) +
      theme(axis.text = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.x = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.y = element_text(size=16,face="bold", color = "black")) +
      theme(axis.title = element_text(size=16,face="bold", color = "black")) +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size=16, face="bold", color = "black")) +
      theme(legend.position = c(0.9, 0.9))

dev.off()

#####

obs_check <- long_df %>%
      group_by(Year) %>%
      summarize(obs = mean(obs), 
                volume = mean(vol), 
                volumesd = sd(vol))

library(ggstatsplot)
library(ggside)

# plot observations alongside volume to look for correlation

tiff("vol_function of obs_11262022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = obs_check,
      x = obs,
      y = volume,
      bf.message = FALSE
)

dev.off()

# seems to be a pretty strong correlation.. pearons r = -0.72

### plot eadj against volume metrics, but first summarize data

eadj_corr <- long_df %>%
      group_by(Year) %>%
      summarize(volume = mean(vol),
                Eadj_Mean = mean(Eadj_Seasonal_Mean),
                Eadj_SD = mean(Eadj_Seasonal_SD),
                Eadj_Min = mean(Eadj_Seasonal_Min),
                Eadj_Max = mean(Eadj_Seasonal_Max))

### eadj vs mean

tiff("vol_function of eadj_MEAN_11262022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_Mean,
      y = volume,
      bf.message = FALSE
)

dev.off()

### eadj vs sd

tiff("vol_function of eadj_SD_11262022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_SD,
      y = volume,
      bf.message = FALSE
)

dev.off()

### eadj vs min

tiff("vol_function of eadj_MIN_11262022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_Min,
      y = volume,
      bf.message = FALSE
)

dev.off()

### eadj vs max

tiff("vol_function of eadj_MAX_11262022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_Max,
      y = volume,
      bf.message = FALSE
)

dev.off()

# REDO Correlations without Water Year 2012 and 2020 (low sample s --------

long_df <- read.csv("snook_similarity_long_NO OUTLIERS_11262022.csv")

long_df$Year<-factor(long_df$Year, levels = c("2011", "2013", "2014", 
                                              "2016", "2017", "2018", "2019"))

str(long_df)
#### PRODUCT NUMBER TWO = PLOT OF NICHE SIZE (VIOLIN PLOT?)

################################################################################

tiff("wYearVolByYear_NO OUTLIERS_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggplot(long_df, aes(x=Year, y=vol, group=Year)) +
      geom_boxplot(width =0.8) +
      labs(x = "Water Year", 
           y = "Niche Volume") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.title = element_text(size=14, face="bold", color = "black")) +
      theme(axis.text = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.x = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.y = element_text(size=16,face="bold", color = "black")) +
      theme(axis.title = element_text(size=16,face="bold", color = "black")) +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size=16, face="bold", color = "black")) +
      theme(legend.position = c(0.9, 0.9))

dev.off()

#####

obs_check <- long_df %>%
      group_by(Year) %>%
      summarize(obs = mean(obs), 
                volume = mean(vol), 
                volumesd = sd(vol))

library(ggstatsplot)
library(ggside)

# plot observations alongside volume to look for correlation

tiff("vol_function of obs_NO OUTLIERS_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = obs_check,
      x = obs,
      y = volume,
      bf.message = FALSE
)

dev.off()

### plot eadj against volume metrics, but first summarize data

eadj_corr <- long_df %>%
      group_by(Year) %>%
      summarize(volume = mean(vol),
                Eadj_Mean = mean(Eadj_Seasonal_Mean),
                Eadj_SD = mean(Eadj_Seasonal_SD),
                Eadj_Min = mean(Eadj_Seasonal_Min),
                Eadj_Max = mean(Eadj_Seasonal_Max))

### eadj vs mean

tiff("vol_function of eadj_MEAN_NO OUTLIERS_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_Mean,
      y = volume,
      bf.message = FALSE
)

dev.off()

### eadj vs sd

tiff("vol_function of eadj_SD_NO OUTLIERS_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_SD,
      y = volume,
      bf.message = FALSE
)

dev.off()

### eadj vs min

tiff("vol_function of eadj_MIN_NO OUTLIERS_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_Min,
      y = volume,
      bf.message = FALSE
)

dev.off()

### eadj vs max

tiff("vol_function of eadj_MAX_NO OUTLIERS_11272022.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

ggscatterstats(
      data = eadj_corr,
      x = Eadj_Max,
      y = volume,
      bf.message = FALSE
)

dev.off()