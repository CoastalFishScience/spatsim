# background --------------------------------------------------------------

#author: Mack White
#project: SRFCEA - Snook Spatial Similarity Manuscript
#goal of script: update hypervolumes through 2023
#date(s): January 2024


# load in packages --------------------------------------------------------
librarian::shelf(tidyverse, hypervolume, MixSIAR, readr, 
                 writexl, ggpubr, ggthemes, ggplot2)


# convert MixSIAR output to csv for HVs -----------------------------------

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

df <- mixTable("snook_ss01082024_NORMAL_CRRCT.txt",'year',ind = F,nest = F)
#df <- df [-c(8:12)]
#write_csv(df, "/Users/mack/Desktop/RESEARCH/Manuscripts_R Scripts/movement/spatial similarity manuscript/df_with_year_ss.csv")
# exported... year should be type, turn this right in excel
df <- read.csv("df_with_year_ss.csv") # year should be the type - like the actual year itself
# type is an integer, but needs to be a character
df$type <- as.character(df$type)
str(df)
# number or iterations
n = 20
# Hypervolumes Function ---------------------------------------------------

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



# load csv generated above ------------------------------------------------

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
