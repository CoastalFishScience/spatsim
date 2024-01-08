# background --------------------------------------------------------------

#author: Mack White
#project: SRFCEA - Snook Spatial Similarity Manuscript
#goal of script: update hypervolumes through 2023
#date(s): January 2024


# load in packages --------------------------------------------------------

library(tidyverse)
library(hypervolume)
library(ggpubr)
library(ggthemes)

#Function to make random points of n length
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
            df = df |>  rename(Mean = mean)
      }
      
      if (is.numeric(sd)){
            names(df)[sd] = 'SD'
      }else {
            df = df |>  rename(SD = sd)
      }
      
      if (end_points){
            if (is.numeric(low)){
                  names(df)[low] = 'lower'
            }else {
                  df = df |>  rename(lower = low)
            }
            
            if (is.numeric(high)){
                  names(df)[high] = 'upper'
            }else {
                  df = df |>  rename(upper = high)
            }
      }
      
      # make sure the names is not numeric 
      if (is.numeric(names)){
            names = names(df)[names]
      }
      
      labs = unique(as.vector(df[,names])[[1]])
      # generate random points within bounds
      if(end_points){
            
            df_tot = df |> slice(rep(1:n(), each=n))|> 
                  mutate(point = 
                               truncnorm::rtruncnorm(1, a = lower, b = upper,
                                                     mean = Mean, sd = SD)) |> 
                  ungroup() |> 
                  mutate(num = rep(1:n, times=nrow(df))) |>
                  select(-Mean, -SD, -lower, -upper)|>
                  pivot_wider(names_from = names, values_from = point)|> 
                  select(-num)
      }else {
            # generate random points outside of bounds
            df_tot = df |> slice(rep(1:n(), each=n))|>
                  mutate(point = 
                               truncnorm::rtruncnorm(1, mean = Mean, sd = SD)) |> 
                  ungroup() |> 
                  mutate(num = rep(1:n, times=nrow(df))) |>
                  select(-Mean, -SD)|>
                  pivot_wider(names_from = names, values_from = point)|> 
                  select(-num)
      }
      if (z_score){
            df_tot = df_tot  |>  
                  mutate(across(all_of(labs), scale))
      }
      
      return(df_tot)
      
}


# load code
g = read_csv('data/combined.csv') |> 
      pivot_longer(cols = c(green, brown), 
                   names_to = 'source',
                   values_to = 'value') |> 
      group_by(site, season, fill, source) |> 
      summarize(m = mean(value),
                sd = sd(value))


# run hypervolumes
gb = HVvalues(df = g, ID_rows = c('site','season', 'fill'),
              names = c('source'), 
              mean = 'm', sd = 'sd', n = 20,
              end_points = F, 
              z_score = T) |> 
      group_by(site,season,fill) |> 
      nest() |> 
      mutate(hv = map(data, ~hypervolume_gaussian(.x, name = paste(site,season,sep = '_'), # change to label columns in data
                                                  samples.per.point = 1000,
                                                  kde.bandwidth = estimate_bandwidth(.x), 
                                                  sd.count = 3, 
                                                  quantile.requested = 0.95, 
                                                  quantile.requested.type = "probability", 
                                                  chunk.size = 1000, 
                                                  verbose = F)),
             hv_size = get_volume(hv[[1]]))
