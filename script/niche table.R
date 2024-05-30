dat <- read_csv("data/snook_ss01302024_UpdatedwYear_2monthFiltered.csv")

dat1 <- read_csv("data/spat_sim_allthegoods_01_30_2024.csv")

ggplot(dat, aes(as.factor(year), mean, fill = source)) +
      geom_boxplot()

summ <- dat |> 
      group_by(year, source) |> 
      summarize(mean_source = mean(mean),
                sd = sd(mean),
                n = n()) 
