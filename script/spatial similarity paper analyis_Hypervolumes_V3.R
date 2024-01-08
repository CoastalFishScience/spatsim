#Author: Mack White
#Date: July/August 2023
#Project: Spatial Similarity Manuscript - Hyper-Volumes

setwd("/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript")
# load packages necessary for generating hypervolumes
library(tidyverse)
library(hypervolume)
library(truncnorm)
library(tidyverse)
library(alphahull)
library(svMisc)
library(writexl)
library(ggplot2)

# read in custom functions developed by Ryan James

#convert to csv----
# updated on 12/8/21****
# function to convert MixSIAR output to a table to generate random data
# to use for hypervolumes
# file = name of .txt of summary statistics of MixSIAR
# type = identifying information of the data set
# ind = if true will output with columns as source values
# nest = if nested mixing model, T will return the sources of nested

mixTable = function(file,type,ind = F,nest = F){
      require(tidyverse)
      cn = c('ID', 'Mean', 'SD', '2.5%', '5%', '25%', '50%', '75%', '95%', '97.5%')
      x = read_fwf(file, skip = 8)
      names(x) = cn
      x$source = NA
      x$name = NA
      x$code = NA
      
      if (nest == F){
            for (i in 1:nrow(x)){
                  temp = strsplit(x$ID, split = '.', fixed = T)
                  x$source[i] = temp[[i]][3]
                  x$name[i] = temp[[i]][2]
                  
                  x$type = type
                  x$ymax = x$`75%` + 1.5*(x$`75%` - x$`25%`)
                  x$ymin = x$`25%` - 1.5*(x$`75%` - x$`25%`)
                  
                  df = data.frame(x$name, x$type, x$source, x$Mean, x$SD, x$`2.5%`, x$`97.5%`,
                                  x$`50%`, x$`25%`, x$`75%`, x$ymax, x$ymin)
                  colnames(df) = c('name', 'type', 'source', 'mean', 'sd', 'lowend', 'highend',
                                   'mid', 'low', 'up', 'ymax', 'ymin')
            }
      }else{
            for (i in 1:nrow(x)){
                  temp = strsplit(x$ID, split = '.', fixed = T)
                  x$source[i] = temp[[i]][4]
                  x$code[i] = temp[[i]][3]
                  x$name[i] = temp[[i]][2]
                  
                  x$type = type
                  x$ymax = x$`75%` + 1.5*(x$`75%` - x$`25%`)
                  x$ymin = x$`25%` - 1.5*(x$`75%` - x$`25%`)
                  
                  df = tibble(x$name, x$type, x$source, x$code, x$Mean, x$SD, x$`2.5%`, x$`97.5%`,
                              x$`50%`, x$`25%`, x$`75%`, x$ymax, x$ymin)
                  colnames(df) = c('name', 'type', 'source', 'code', 'mean', 'sd', 'lowend', 'highend',
                                   'mid', 'low', 'up', 'ymax', 'ymin')
            }
      }
      
      for (i in 1:nrow(df)){
            if (df$ymax[i] > df$highend[i]){
                  df$ymax[i] = df$highend[i]
            }
            if (df$ymin[i] < df$lowend[i]){
                  df$ymin[i] = df$lowend[i]
            }
      }
      df = df %>% drop_na %>%
            filter(name != 'global')
      
      
      if (ind == T){
            if (nest == T){
                  df = df %>% select(name, type, code, source, mean) %>%
                        pivot_wider(names_from = 'source', values_from = 'mean')
            }else{
                  df = df %>% select(name, type, source, mean)%>%
                        pivot_wider(names_from = 'source', values_from = 'mean')
            }
      }
      
      return(df)
}

# Function to make random points of n length----
# data from random sample with mean and sd but 
# can be generated between a high and low value if end_points = T
# *** Note chose either column names or column numbers for ID_rows and names
# either work but must be the same
# df = dataframe or tibble with each row containing 
#        unique entry for making random points 
# ID_rows = vector of column names or numbers with id information
# names = column name or number of name of measure variables
# mean = column name or column number of df with mean 
# sd = column name or column number of df with sd 
# n = number of points to randomly generate
# z_score = T or F, if T z-scores values
# end_points = T or F for if random points need to be generated between
#        a high and low end point (e.g. 5% and 95% interval)
#        low and high required if end_points = T
# low = column name or column number of df with lower bound to sample in
# high = column name or column number of df with upper bound to sample in

HVvalues = function(df, ID_rows, names, mean, sd, n, z_score = F,
                    end_points = F, low = NULL, high = NULL){
      require(tidyverse)
      require(truncnorm)
      
      # check to see if information is needed to restrict where points are
      if (end_points){
            if (is_empty(df[,low]) | is_empty(df[,high])){
                  return(cat('Warning: low and/or high columns not specified \n
                  Specific and run again or end_points = F \n'))
            }
      }
      
      # check to see if there are more 
      if(T %in% duplicated(df[,c(ID_rows,names)])){
            return(cat('Warning: some of the rows contain duplicated information \n
                make sure data is correct \n'))
      }
      
      # rename variables to make code work
      if (is.numeric(mean)){
            names(df)[mean] = 'Mean'
      }else {
            df = df %>% rename(Mean = mean)
      }
      
      if (is.numeric(sd)){
            names(df)[sd] = 'SD'
      }else {
            df = df %>% rename(SD = sd)
      }
      
      if (end_points){
            if (is.numeric(low)){
                  names(df)[low] = 'lower'
            }else {
                  df = df %>% rename(lower = low)
            }
            
            if (is.numeric(high)){
                  names(df)[high] = 'upper'
            }else {
                  df = df %>% rename(upper = high)
            }
      }
      
      # make sure the names is not numeric 
      if (is.numeric(names)){
            names = names(df)[names]
      }
      
      # generate random points within bounds
      if(end_points){
            
            df_tot = df %>% slice(rep(1:n(), each=n))%>% 
                  mutate(point = 
                               truncnorm::rtruncnorm(1, a = lower, b = upper,
                                                     mean = Mean, sd = SD),
                         num = rep(1:n, times=nrow(df))) %>%
                  select(-Mean, -SD, -lower, -upper)%>%
                  pivot_wider(names_from = names, values_from = point)%>% 
                  select(-num)
      }else {
            # generate random points outside of bounds
            df_tot = df %>% slice(rep(1:n(), each=n))%>%
                  mutate(point = 
                               truncnorm::rtruncnorm(1, mean = Mean, sd = SD),
                         num = rep(1:n, times=nrow(df))) %>%
                  select(-Mean, -SD)%>%
                  pivot_wider(names_from = names, values_from = point)%>% 
                  select(-num)
      }
      if (z_score){
            df_tot = df_tot %>% 
                  mutate_if(is.numeric, scale)
      }
      
      return(df_tot)
      
}

#df <- mixTable("MixingModels/mm_results/snook_ss11262022.txt",'year',ind = F,nest = F)
#df <- df [-c(8:12)]
#write_csv(df, "/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript/df_with_year_ss.csv")
# exported... year should be type, turn this right in excel
df <- read.csv("df_with_year_ss.csv") # year should be the type - like the actual year itself
# type is an integer, but needs to be a character
df$type <- as.character(df$type)
str(df)
# number or iterations
n = 20

# create tibble to store data
df_hvAll = tibble(rep = rep(seq(1,n,1)),
                  vol_2011 = NA, vol_2012 = NA, vol_2013 = NA, vol_2014 = NA,
                  vol_2016 = NA, vol_2017 = NA, vol_2018 = NA, vol_2019 = NA, vol_2020 = NA)

##### zscore and generate hypervolumes

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
                       hv_2017, hv_2018, hv_2019, hv_2020)

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