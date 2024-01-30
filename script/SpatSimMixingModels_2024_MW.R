# background --------------------------------------------------------------

#author: Mack White
#project: SRFCEA - Snook Spatial Similarity Manuscript
#goal of script: update mixing models through 2023
#date(s): January 2024


# load in packages --------------------------------------------------------

librarian::shelf(tidyverse, MixSIAR, readr, writexl)
options(max.print = 500000) #changed from 6000000 to match expressions on 01/08/2024
options(expressions=500000) #saw note online that this could help with output issues?

# check out mixing file ---------------------------------------------------

mixcheck <- read_csv("data/archive/mix_formatted_01082024_ALL.csv")
glimpse(mixcheck)

mix_summary <- mixcheck |> 
      group_by(wYear) |> 
      summarise(n = length(unique(ID)))

glimpse(mix_summary)

###january 29 -> checking to see if filtered properly... it was not!
### filtering based on 3-month below:
mixcheck_filtered <- mixcheck |> 
      filter(Month %in% c("Mar", "Apr", "May", "Jun", "Jul"))

mix_summary_filtered <- mixcheck_filtered |> 
      group_by(wYear) |> 
      summarise(n = length(unique(ID)))

# write_csv(mixcheck_filtered, "data/archive/mix_formatted_01292024_FILTERED.csv")
### lost 24 samples - break down below
### initial sample size without filtering

#2011 - 12 to 9
#2012 - 7 to 6
#2013 - 11 to 10
#2014 - 22 to 22
#2015 - 0 to 0
#2016 - 8 to 5
#2017 - 14 to 11
#2018 - 14 to 11
#2019 - 4 to 2
#2020 - 10 to 8
#2021 - 4 to 4
#2022 - 7 to 1

# check out source file  --------------------------------------------------

sourcecheck <- read_csv("data/archive/ss_snook_source_agg_UPDATED_07262023.csv")

glimpse(sourcecheck)

# check out tdf file ------------------------------------------------------

TDFcheck <- read_csv("data/archive/snook_agg_nona_UPDATED.csv")

glimpse(TDFcheck)

# set up mixing file ------------------------------------------------------

mix = load_mix_data(file("data/archive/mix_formatted_01292024_FILTERED.csv"),
                    iso_names=c("d13C","d15N","d34S"),
                    factors= c("wYear", "ID"),
                    fac_random=c(F,T), #water year is nonrandom, id is random
                    fac_nested=c(F,T), #water year is nonrandom, id is random
                    cont_effects=NULL)

# set up source file ------------------------------------------------------

source = load_source_data(file("data/archive/ss_snook_source_agg_UPDATED_07262023.csv"),
                          source_factors=NULL,
                          conc_dep=FALSE,
                          data_type="means",
                          mix)
### will likely receive warning message here about incomplete final line - has been an issue for ~ 1 year but models still go

# set up tdf file ---------------------------------------------------------

discr = load_discr_data(file("data/archive/snook_agg_nona_UPDATED.csv"), mix)
### will likely receive warning message here about incomplete final line - has been an issue for ~ 1 year but models still go

# Generate Isospace Plot --------------------------------------------------

# plot_data(filename="isospace_plot_01_08_2024", plot_save_pdf=TRUE, 
#           plot_save_png=TRUE, mix, source, discr)

# Write and Run Jags Model ------------------------------------------------
### set up model
model_filename = "MixSIAR_model.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

###run the model
jags.sss = run_model(run= "long", mix, source, discr, model_filename,
                     alpha.prior = 1, resid_err, process_err)

save(jags.sss, file = "data/SpatSimRJags01292024_UpdatedWwYear.RData") #-> #following filtering the 3 months
load("data/SpatSimRJags01292024_UpdatedWwYear.RData")
# save(jags.sss, file = "SpatSimRJags01082024_UpdatedWwYear.RData") -> prior to filtering based on 3 months rule
#load("SpatSimRJags01082024_UpdatedWwYear.RData") #model run on normal saved since it takes time
# Process JAGS output

### everything besides summary saved
# output_sss_NORMAL = list(summary_save = TRUE,
#                   summary_name = "snook_ss01082024_NORMAL",
#                   sup_post = FALSE,
#                   plot_post_save_pdf = FALSE,
#                   plot_post_name = "lower_posterior_density",
#                   sup_pairs = FALSE,
#                   plot_pairs_save_pdf = FALSE,
#                   plot_pairs_name = "lower_pairs_plot",
#                   sup_xy = FALSE,
#                   plot_xy_save_pdf = FALSE,
#                   plot_xy_name = "lower_xy_plot",
#                   gelman = FALSE, 
#                   heidel = FALSE,
#                   geweke = FALSE, 
#                   diag_save = FALSE,
#                   diag_name = "snook_diag_01082024_NORMAL",
#                   indiv_effect = FALSE,
#                   plot_post_save_png = F,
#                   plot_pairs_save_png = FALSE,
#                   plot_xy_save_png = FALSE)

# output_JAGS(jags.sss, mix, source, output_sss_NORMAL)

###normal output list
output_sss = list(summary_save = TRUE,
                         summary_name = "snook_ss01292024_NORMAL_Updated",
                         sup_post = FALSE,
                         plot_post_save_pdf = FALSE,
                         plot_post_name = "lower_posterior_density",
                         sup_pairs = FALSE,
                         plot_pairs_save_pdf = FALSE,
                         plot_pairs_name = "lower_pairs_plot",
                         sup_xy = FALSE,
                         plot_xy_save_pdf = FALSE,
                         plot_xy_name = "lower_xy_plot",
                         gelman = TRUE,
                         heidel = FALSE,
                         geweke = TRUE,
                         diag_save = TRUE,
                         diag_name = "snook_diag_01292024_NORMAL_Updated",
                         indiv_effect = FALSE,
                         plot_post_save_png = F,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)

output_JAGS(jags.sss, mix, source, output_sss)
#if you get an error, just keep running the same line of code over and over again - it works, trust me
