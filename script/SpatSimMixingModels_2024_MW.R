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

mixcheck <- read_csv("data/mix_formatted_01072024_ALL.csv")
glimpse(mixcheck)

mix_summary <- mixcheck |> 
      group_by(wYear) |> 
      summarise(n = length(unique(ID)))

glimpse(mix_summary)

# check out source file  --------------------------------------------------

sourcecheck <- read_csv("data/ss_snook_source_agg_UPDATED_07262023.csv")

glimpse(sourcecheck)

# check out tdf file ------------------------------------------------------

TDFcheck <- read_csv("data/snook_agg_nona_UPDATED.csv")

glimpse(TDFcheck)

# set up mixing file ------------------------------------------------------

mix = load_mix_data(file("data/mix_formatted_01072024_ALL.csv"),
                    iso_names=c("d13C","d15N","d34S"),
                    factors= c("wYear", "ID"),
                    fac_random=c(F,T), #water year is nonrandom, id is random
                    fac_nested=c(F,T), #water year is nonrandom, id is random
                    cont_effects=NULL)

# set up source file ------------------------------------------------------

source = load_source_data(file("data/ss_snook_source_agg_UPDATED_07262023.csv"),
                          source_factors=NULL,
                          conc_dep=FALSE,
                          data_type="means",
                          mix)
### will likely receive warning message here about incomplete final line - has been an issue for ~ 1 year but models still go

# set up tdf file ---------------------------------------------------------

discr = load_discr_data(file("data/snook_agg_nona_UPDATED.csv"), mix)
### will likely receive warning message here about incomplete final line - has been an issue for ~ 1 year but models still go

# Generate Isospace Plot --------------------------------------------------

plot_data(filename="isospace_plot_01_08_2024", plot_save_pdf=TRUE, 
          plot_save_png=TRUE, mix, source, discr)

# Write and Run Jags Model ------------------------------------------------
### set up model
model_filename = "MixSIAR_model.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

###run the model
jags.sss = run_model(run= "normal", mix, source, discr, model_filename,
                     alpha.prior = 1, resid_err, process_err)

# save(jags.sss, file = "SpatSimRJags01082024.RData")
#load("SpatSimRJags01082024.RData") #model run on normal saved since it takes time
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
output_sss_try = list(summary_save = TRUE,
                         summary_name = "snook_ss01082024_NORMAL_CRRCT",
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
                         diag_name = "snook_diag_01082024_NORMAL_CRRCT",
                         indiv_effect = FALSE,
                         plot_post_save_png = F,
                         plot_pairs_save_png = FALSE,
                         plot_xy_save_png = FALSE)

output_JAGS(jags.sss, mix, source, output_sss_try)

# output_JAGS(jags.sss, mix, source, output_sss_NORMAL)

# output_JAGS(jags.sss, mix, source, output_options = list(summary_save = TRUE))


