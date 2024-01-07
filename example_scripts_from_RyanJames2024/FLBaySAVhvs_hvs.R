#' """ Create hypervolumes comparing seagrass basins
#'     @author: Ryan James
#'     Date : 8/28/23
#'     edit : 9/15/23

library(tidyverse)
library(hypervolume)
library(vegan)
library(ggpubr)
library(ggthemes)


# prepare data ----
# TSG, TT, HW, SF, Seagrass Richness (TT, HW, SF, RM, HE, HD), TMA 
# z-score for all basins of interest and 
# add tiny amount so not all values same and can make hv
# only 2007-2022 because has all variable (minus 2020)

basins = c("BLK", "CAL", "CRN", "EAG", "JON",
           "MAD", "RAN", "RKB", "TWN", "WHP")

set.seed(14)
# prepare data 
df = read_csv('data/FHAP_BB_QAQCd.csv') |> 
  filter(BASIN %in% basins) |> 
  mutate(sg_rich = specnumber(across(TT:HD))) |> 
  mutate(across(TOT:RHZ, ~ case_when(
    .x == 0 ~ 0,
    .x == 0.1 ~ 0.03,
    .x == 0.5 ~ 0.03,
    .x == 1 ~ 0.03,
    .x == 2 ~ 0.15,
    .x == 3 ~ 0.375,
    .x == 4 ~ 0.625,
    .x == 5 ~ 0.875))) |>
  select(BASIN, YEAR, STATION, TT, HW, SF, TMA, TDR, sg_rich) |> 
  drop_na() |>
  group_by(BASIN, YEAR, STATION) |> 
  summarize(across(where(is.numeric), mean)) |>
  ungroup() |> 
  mutate(across(c(TT:sg_rich), scale), # z score
         across(c(TT:sg_rich), # add tiny amount so when all values the same can make hv
                ~map_dbl(., ~. + rnorm(1, mean = 0, sd = 0.0001)))) |> 
  select(-STATION) |> 
  group_by(BASIN, YEAR) |>
  nest() 

# make hypervolumes ----
# assigned basins
basin = c('JON')

# generate hypervolumes
df = df |> 
  #filter(BASIN %in% basin) |> 
  mutate(hv = map(data, ~hypervolume_gaussian(.x, name = paste(BASIN,YEAR,sep = '_'),
                                              samples.per.point = 1000,
                                              kde.bandwidth = estimate_bandwidth(.x), 
                                              sd.count = 3, 
                                              quantile.requested = 0.95, 
                                              quantile.requested.type = "probability", 
                                              chunk.size = 1000, 
                                              verbose = F)),
         hv_size = map_dbl(hv, \(hv) get_volume(hv)),
         centroid = map(hv, \(hv) get_centroid(hv)))


# save files 
#saveRDS(df, paste0('data/hv',paste(basin, collapse = '_'),'.rds'))
saveRDS(df, 'data/hvAll.rds')

# df = readRDS('data/hvAll.rds')

dp = df |> unnest_wider(centroid) |> select(-hv)
write_csv(dp, 'data/hvAll.csv')

# comparison of across each year
df_y= tibble(y1 = unique(df$YEAR),
            y2 = unique(df$YEAR)) |> 
   expand(y1,y2)

df_y = df_y[!duplicated(t(apply(df_y,1,sort))),] %>% 
  filter(!(y1 == y2))

# years to make 
df1 = df |> 
  select(BASIN, y1 = YEAR, hv1 = hv, hv1_size = hv_size, cent1 = centroid)

df2 = df |> 
  select(BASIN, y2 = YEAR, hv2 = hv, hv2_size = hv_size, cent2 = centroid)


# create data frame of all data and make yearly comparisons
df_ov = tibble(BASIN = rep(unique(df$BASIN),
                         each = nrow(df_y)),
               y1 = rep(df_y$y1, times = length(unique(df$BASIN))),
               y2 = rep(df_y$y2, times = length(unique(df$BASIN)))) |> 
  inner_join(df1, by = c('BASIN', 'y1')) |> 
  inner_join(df2, by = c('BASIN', 'y2')) |> 
  mutate(ychange = y2-y1,
         size_rat = hv2_size/hv1_size,
         set = map2(hv1,hv2, \(hv1, hv2) hypervolume_set(hv1, hv2, check.memory = F, verbose = F)),
         ov = map(set, \(set) hypervolume_overlap_statistics(set)),
         dist_cent = map2_dbl(hv1, hv2, \(hv1,hv2) hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F)),
         dif = map2(cent1, cent2, \(cent1,cent2) cent2 - cent1)) |> 
  unnest_wider(ov) |> 
  unnest_wider(dif) |> 
  select(BASIN, y1, y2, ychange, hv1_size, hv2_size, size_rat, 
         jaccard, sorensen,uniq_y1 = frac_unique_1, uniq_y2 = frac_unique_2, 
         dist_cent, TT, HW, SF, sg_rich, TMA, TDR)

# save outputs
write_csv(df_ov, 'data/hv_ovAll.csv')
# df_ov = read_csv('data/hv_ovAll.csv')

# compare years
# for(i in 1:nrow(df_ov)){
#   # set hvs
#   hv1 = df_ov$hv1[i][[1]]
#   hv2 = df_ov$hv2[i][[1]]
#   
#   #join hvs
#   set = hypervolume_set(hv1, hv2, check.memory = F, verbose = F)
#   
#   # calculate sorensen overlap
#   df_ov$sorensen[i] = hypervolume_overlap_statistics(set)[2]
#   df_ov$uniq_y1[i] = hypervolume_overlap_statistics(set)[3]
#   df_ov$uniq_y2[i] = hypervolume_overlap_statistics(set)[4]
#   
#   # calculate distance
#   df_ov$dist_cent[i] = hypervolume_distance(hv1, hv2, type = 'centroid', check.memory=F)
#   
#   # calculate difference of each time step
#   df_ov$dif[i] = list(df_ov$cent2[i][[1]] - df_ov$cent1[i][[1]])
#   
#   cat(i, '/',nrow(df_ov),'\n')
# }

# modify and save outputs
# dd = df_ov |> 
#   mutate(ychange = y2-y1) |> 
#   select(BASIN, y1, y2, ychange, hv1_size, hv2_size, size_rat, sorensen,
#          uniq_y1, uniq_y2, dist_cent, TT, HW, SF, sg_rich, TMA, TDR)
# 
# write_csv(dd, paste0('data/hvOV_',paste(basin, collapse = '_'),'.csv'))

# plot year to year changes----

# size 
ggplot(df, aes(YEAR, hv_size, color = BASIN))+
  geom_point(size = 2.5)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'Volume')+
  theme_bw()+
  facet_wrap(~BASIN, scales = 'free_y', nrow = 2)+
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

ggsave('figs/hvSizeYearly.tiff', 
       units="in", width=13, height=9, dpi=600,compression = 'lzw')

dfb = df_ov |> 
  filter(ychange == 1 )

ggplot(dfb, aes(y2, sorensen, color = BASIN))+
  geom_point(size = 2.5)+
  geom_line(linewidth = 1)+
  facet_wrap(~BASIN, nrow = 2)+
  labs(x = 'Year', y = 'Overlap')+
  theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/hvOvYearly.tiff', 
       units="in", width=13, height=9, dpi=600,compression = 'lzw')

ggplot(dfb, aes(y2, dist_cent, color = BASIN))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_line(linewidth = 1)+
  facet_wrap(~BASIN, nrow = 2)+
  labs(x = 'Year', y = 'Centroid distance')+theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/hvDistYearly.tiff', 
       units="in", width=13, height=9, dpi=600,compression = 'lzw')
  
# d = ggplot(dfb, aes(y2,uniq_y1, color = BASIN))+
#   geom_point(size = 2.5)+
#   geom_line(linewidth = 1)+
#   labs(x = 'Year', y = 'Year 1 proportion unique')+theme_bw()+
#   theme(axis.title = element_text(size = 14), 
#         axis.text = element_text(size = 14, colour = "gray0"), 
#         plot.title = element_text(size = 14, hjust=0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = 'right',
#         legend.title = element_text(size = 14),
#         strip.text.x = element_text(size = 14),
#         legend.text = element_text(size = 12))
# 
# e = ggplot(dfb, aes(y2,uniq_y2, color = BASIN))+
#   geom_point(size = 2.5)+
#   geom_line(linewidth = 1)+
#   labs(x = 'Year', y = 'Year 2 proportion unique')+theme_bw()+
#   theme(axis.title = element_text(size = 14), 
#         axis.text = element_text(size = 14, colour = "gray0"), 
#         plot.title = element_text(size = 14, hjust=0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = 'right',
#         legend.title = element_text(size = 14),
#         strip.text.x = element_text(size = 14),
#         legend.text = element_text(size = 12))

# ggarrange(a,b,c,d,e,
#           nrow = 2, ncol = 3, common.legend = T)

# ggsave(paste0('figs/yearly_',paste(basin, collapse = '_'),'.tiff'), 
#        units="in", width=12, height=8, dpi=600,compression = 'lzw')

# ggsave('figs/yearly_All.tiff', 
#        units="in", width=12, height=8, dpi=600,compression = 'lzw')

# all years comparisons ----
ggplot(df_ov, aes(ychange, sorensen, color = BASIN))+
  geom_point(size = 2.5)+
  geom_smooth(linewidth = 1)+
  labs(x = 'Years between comparison', y = 'Overlap')+
  theme_bw()+
  facet_wrap(~BASIN, nrow = 2)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/hvOvAll.tiff', 
       units="in", width=13, height=9, dpi=600,compression = 'lzw')

ggplot(df_ov, aes(ychange, dist_cent, color = BASIN))+
  geom_hline(aes(yintercept = 1), linetype = 'dashed')+
  geom_point(size = 2.5)+
  geom_smooth(linewidth = 1)+
  facet_wrap(~BASIN,  nrow = 2)+
  labs(x = 'Years between comparison', y = 'Centroid distance')+theme_bw()+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/hvDistAll.tiff', 
       units="in", width=13, height=9, dpi=600,compression = 'lzw')

# ggplot(df_ov, aes(ychange,size_rat, color = BASIN))+
#   geom_point(size = 2.5)+
#   geom_smooth(linewidth = 1)+
#   facet_wrap(~BASIN,  scales = 'free_y', nrow = 2)+
#   labs(x = 'Years between comparison', y = 'Y2 size/ Y1 size')+theme_bw()+
#   theme(axis.title = element_text(size = 14), 
#         axis.text = element_text(size = 14, colour = "gray0"), 
#         plot.title = element_text(size = 14, hjust=0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = 'right',
#         legend.title = element_text(size = 14),
#         strip.text.x = element_text(size = 14),
#         legend.text = element_text(size = 12))
# 
# d = ggplot(df_ov, aes(ychange,uniq_y2, color = BASIN))+
#   #geom_point(size = 2.5)+
#   geom_smooth(linewidth = 1)+
#   labs(x = 'Years between comparison', y = 'Year 2 proportion unique')+theme_bw()+
#   theme(axis.title = element_text(size = 14), 
#         axis.text = element_text(size = 14, colour = "gray0"), 
#         plot.title = element_text(size = 14, hjust=0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = 'right',
#         legend.title = element_text(size = 14),
#         strip.text.x = element_text(size = 14),
#         legend.text = element_text(size = 12))
# 
# e = ggplot(df_ov, aes(ychange,uniq_y1, color = BASIN))+
#   #geom_point(size = 2.5)+
#   geom_smooth(linewidth = 1)+
#   labs(x = 'Years between comparison', y = 'Year 1 proportion unique')+theme_bw()+
#   theme(axis.title = element_text(size = 14), 
#         axis.text = element_text(size = 14, colour = "gray0"), 
#         plot.title = element_text(size = 14, hjust=0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = 'right',
#         legend.title = element_text(size = 14),
#         strip.text.x = element_text(size = 14),
#         legend.text = element_text(size = 12))
# 
# ggarrange(a,b,c,d,e,
#           nrow = 2, ncol = 3, common.legend = T)
# 
# # ggsave(paste0('figs/ychange_',paste(basin, collapse = '_'),'.tiff'), 
# #        units="in", width=12, height=8, dpi=600,compression = 'lzw')
# 
# ggsave('figs/ychange_All.tiff', 
#        units="in", width=12, height=8, dpi=600,compression = 'lzw')

# centroids for each year
dc = df |> 
  unnest_longer(centroid)

ggplot(dc, aes(YEAR, centroid, color = BASIN))+
  geom_point(size = 2.5)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'z-score')+
  theme_bw()+
  # scale_y_continuous(labels = scales::percent)+
  facet_wrap(~centroid_id)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/zSAVbasin.tiff', 
       units="in", width=12, height=8, dpi=600,compression = 'lzw')

ggplot(dc, aes(YEAR, centroid, color = centroid_id))+
  geom_point(size = 2.5)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'z-score', color = NULL)+
  theme_bw()+
  scale_color_colorblind()+
  facet_wrap(~BASIN, nrow = 2)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave('figs/zBasinSAV.tiff', 
       units="in", width=13, height=9, dpi=600,compression = 'lzw')

# change for each year
dd = dfb |> 
  pivot_longer(TT:TDR, names_to = 'centroid_id', values_to = 'centroid')

ggplot(dd, aes(y2, centroid, color = BASIN))+
  geom_point(size = 2.5)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'z-score')+
  theme_bw()+
  # scale_y_continuous(labels = scales::percent)+
  facet_wrap(~centroid_id)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggsave('figs/zSAVbasinDif.tiff', 
       units="in", width=12, height=8, dpi=600,compression = 'lzw')

ggplot(dd, aes(y2, centroid, color = centroid_id))+
  geom_point(size = 2.5)+
  geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'z-score', color = NULL)+
  theme_bw()+
  scale_color_colorblind()+
  facet_wrap(~BASIN, nrow = 2)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave('figs/zBasinSAVDif.tiff', 
       units="in", width=13, height=9, dpi=600,compression = 'lzw')


# plot hypervolumes
eag = df |> filter(BASIN == 'EAG')
hv= hypervolume_join(eag$hv[[1]], eag$hv[[2]])


tiff("figs/EAG07_08.tiff", width = 9, height = 9, units = 'in',
     res = 600, compression = 'lzw')
plot(hv, pairplot = T, colors=c('darkolivegreen3','deepskyblue3'),
     # names= c("Freshwater\n Benthic",
     #          "Marine\n Benthic\n Detritivore",
     #          "Marine\n Benthic\n Zoobenthivore",
     #          "Marine\n Pelagic",
     #          "Trophic\n Position"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-4,4), 
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


# summary of differences in years ----- 
s = dfb |> 
  group_by(BASIN) |> 
  summarize(across(size_rat:dist_cent, var))

v = df |> 
  group_by(BASIN) |> 
  summarize(vol = var(hv_size))

s = dfb |> 
  group_by(BASIN) |> 
  summarize(across(size_rat:dist_cent, ~sd(.x)/mean(.x)))

m = df |> 
  pivot_longer(TT:sg_rich, names_to = 'name', values_to = 'value')

ggplot(m, aes(YEAR, value, color = BASIN))+
  geom_pointrange(stat = "summary",
                  fun.min = 'min',
                  fun.max = 'max',
                  fun = 'mean')+
  stat_summary(aes(y = value), fun = mean, geom = 'line')+
  # geom_point(size = 2.5)+
  # geom_line(linewidth = 1)+
  labs(x = 'Year', y = 'z-score')+
  theme_bw()+
  # scale_y_continuous(labels = scales::percent)+
  facet_wrap(~BASIN)+
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 14, colour = "gray0"), 
        plot.title = element_text(size = 14, hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'right',
        legend.title = element_text(size = 14),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12))

ggplot(df_ov, aes(log(size_rat), dist_cent))+
  geom_point()+
  geom_smooth()

n = df_ov |> 
  pivot_longer(c(hv1_size, hv2_size), names_to = 'name', values_to = 'size')

# overlap only useful when size > 1 when we have 6 axes
ggplot(n |>  filter(size >1),
       aes(log(size), sorensen))+
  geom_point()+
  geom_smooth()
