#Author: Mack White
#Date: July/August 2023
#Project: Spatial Similarity Manuscript - Mixing Models

# NICHE VOLUME TIME -------------------------------------------------------

# load libraries 
library(MixSIAR)
library(tidyverse)
#library(rjags) - loads with mixsiar
options(max.print = 6000000)

# check out files used for mixing models (mix, source, TDF) ---------------

### mixing file (i.e., common snook values)

mix_formatted_07272023 <- read_csv("data/mix_formatted_07272023.csv") #or should it be simply, mix_formatted.csv
mixcheck <- mix_formatted_07272023 
glimpse(mixcheck)

mix_summary <- mixcheck |> 
      group_by(wYear) |> 
      summarise(n = length(unique(ID)))
glimpse(mix_summary)

### source file (i.e., common snook resources)
sourcecheck <- read_csv("data/ss_snook_source_agg_UPDATED_07262023.csv")
glimpse(sourcecheck)

### trophic discrimination factor file (i.e., amt expct values to change for each trophic step)
TDFcheck <- read_csv("data/snook_agg_nona_UPDATED.csv")
glimpse(TDFcheck)


mix = load_mix_data(file("~/Library/CloudStorage/Dropbox/R/github/spatsim/data/mix_formatted.csv"),
                    iso_names=c("d13C","d15N","d34S"),
                    factors= c("wYear", "ID"),
                    fac_random=c(F,T), #water year is nonrandom, id is random
                    fac_nested=c(F,T),
                    cont_effects=NULL)


# mix = load_mix_data(file("~/Library/CloudStorage/Dropbox/R/github/spatsim/data/mix_formatted.csv"),
#                     iso_names=c("d13C","d15N","d34S"),
#                     factors= c("wYear"),
#                     fac_random=c(F), #water year is nonrandom, id is random
#                     fac_nested=c(F),
#                     cont_effects=NULL)

# load source data
source = load_source_data(file("data/ss_snook_source_agg_UPDATED_07262023.csv"),
                          source_factors=NULL,
                          conc_dep=FALSE,
                          data_type="means",
                          mix)
### received warning upon reading in source file on 07262023 - see below:
# Warning message:
# In read.table(file = file, header = header, sep = sep, quote = quote,  :
# incomplete final line found by readTableHeader on '~/Library/CloudStorage/Dropbox/R/github/spatsim/data/ss_snook_source_agg_UPDATED_07262023.csv'

# load TDF data
discr = load_discr_data(file("data/snook_agg_nona_UPDATED.csv"), mix)

### same warning as noted above 

# Make an isospace plot
plot_data(filename="isospace_plot", plot_save_pdf=TRUE, 
          plot_save_png=FALSE, mix, source, discr)

# Write the JAGS model file
model_filename = "MixSIAR_model.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# run jags
jags.sss = run_model(run="test", mix, source, discr, model_filename,
                     alpha.prior = 1, resid_err, process_err)

# Process JAGS output
output_sss = list(summary_save = TRUE,
                  summary_name = "snook_ss07272023_1",
                  sup_post = FALSE,
                  plot_post_save_pdf = FALSE,
                  plot_post_name = "lower_posterior_density",
                  sup_pairs = FALSE,
                  plot_pairs_save_pdf = FALSE,
                  plot_pairs_name = "lower_pairs_plot",
                  sup_xy = TRUE,
                  plot_xy_save_pdf = FALSE,
                  plot_xy_name = "lower_xy_plot",
                  gelman = TRUE,
                  heidel = FALSE,
                  geweke = TRUE,
                  diag_save = TRUE,
                  diag_name = "snook_diag07272023_1",
                  indiv_effect = FALSE,
                  plot_post_save_png = F,
                  plot_pairs_save_png = FALSE,
                  plot_xy_save_png = FALSE)

output_JAGS(jags.sss, mix, source, output_sss)
# output_JAGS_07_26_2023 <- 
# Redo Analysis with peak months (i.e.,  March, April, May) ---------------
### conducted in separated R file... "ThreeMonthTest_SpatialSimilarityV1"