### project: snook spatial similarity - sci reports special issue
### author: mack white
### goal: generate supplemental figure to address Sci Reports reviewer comment

librarian::shelf(tidyverse, readr, zoo)

dat <- read_csv("data/spat_sim_allthegoods_01_30_2024.csv")

dat$date <- as.Date(as.yearmon(paste(dat$Year, dat$Month, sep="-"), "%Y-%m"))

ggplot(dat, aes(x = date)) + 
      geom_line(aes(y = monthly_mean_stage_cm), colour = "black", linetype = "solid", size = 0.8) +
      geom_line(aes(y = daysbelow30), colour = "black", linetype = "dashed", size = 0.8) +
      scale_y_continuous(
            name = "Monthly Mean Stage (cm)",
            sec.axis = sec_axis(~., name="Days Below 30")
      ) +
      labs(x = "Date", y = 'Monthly Mean Stage (cm)') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      # theme(plot.title = element_text(hjust = 0.5)) +
      # theme(plot.title = element_text(size=14, face="bold", color = "black")) +
      theme(axis.text = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.x = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.y = element_text(size=16,face="bold", color = "black")) +
      theme(axis.title = element_text(size=16,face="bold", color = "black")) +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size=16, face="bold", color = "black"))

#saving for publication
# ggsave("figs/manuscript/new_supp_marshstage_07092024.tiff", units = "in", width = 10,
#        height = 6, dpi =  600, compression = "lzw")

### simplified version of supplemental figure for JR talk at IFCT
ggplot(dat, aes(x = date)) + 
      geom_line(aes(y = monthly_mean_stage_cm), colour = "darkblue", linetype = "solid", size = 1.6) +
      # geom_line(aes(y = daysbelow30), colour = "black", linetype = "dashed", size = 0.8) +
      # scale_y_continuous(
      #       name = "Monthly Mean Stage (cm)",
      #       sec.axis = sec_axis(~., name="Days Below 30")
      # ) +
      labs(x = "Date", y = 'Monthly Mean Stage (cm)') +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      # theme(plot.title = element_text(hjust = 0.5)) +
      # theme(plot.title = element_text(size=14, face="bold", color = "black")) +
      theme(axis.text = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.x = element_text(size=16,face="bold", color = "black")) +
      theme(axis.text.y = element_text(size=16,face="bold", color = "black")) +
      theme(axis.title = element_text(size=16,face="bold", color = "black")) +
      theme(legend.title = element_blank()) +
      theme(legend.text = element_text(size=16, face="bold", color = "black"))

#saving for publication
# ggsave("figs/manuscript/new_supp_marshstage_07092024.tiff", units = "in", width = 10,
#        height = 6, dpi =  600, compression = "lzw")

ggsave("../../../../../../Downloads/revstagefig-v1.tiff", units = "in", width = 10,
       height = 6, dpi =  600, compression = "lzw")
