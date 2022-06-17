library(here)
library(rstudioapi)

# pull pollen

# create prediction grid


# Model speed comparison
# jobRunScript(here::here("R", "speed-testing.R"))


# jobRunScript(here::here("R", "generate_maps_data"))

# Build prediction grid
# this file needs to be resolved

# Matern model
jobRunScript(here::here("R", "estimate_pars_ST.R")) # currently takes nearly 4 days (96 hours) on my server
jobRunScript(here::here("R", "predict_ST.R"))
jobRunScript(here::here("R", "generate-maps_ST.R"))

# Overdispersed Matern model
jobRunScript(here::here("R", "estimate_pars_ST_overdispersed.R"))
jobRunScript(here::here("R", "predict_ST_overdispersed.R"))
jobRunScript(here::here("R", "generate-maps_ST_overdispersed.R"))

# Latent Overdispersed Matern model
jobRunScript(here::here("R", "estimate_pars_ST_latent_overdispersed.R"))
jobRunScript(here::here("R", "predict_ST_latent_overdispersed.R"))
jobRunScript(here::here("R", "generate-maps_ST_latent_overdispersed.R"))

# MRA model -- currently fits slow and has prediction (indexing?) errors
# jobRunScript(here::here("R", "estimate_pars_ST_MRA.R"))
# jobRunScript(here::here("R", "predict_ST_MRA.R"))
# jobRunScript(here::here("R", "generate-maps_ST_MRA.R"))

# Evaluate the models
jobRunScript(here::here("R", "evaluate_fit_ST.R"))