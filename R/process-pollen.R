library(tidyverse)
library(patchwork)
require(rasterVis)
require(fields)
require(rgdal)
require(raster)
require(enmSdm)
require(rgeos)

version <- '3.1'

# load the species names
species_names <- readRDS(here::here('data', 'species-names.RDS'))

# load the preds data
preds <- readRDS(here::here("output", paste0('polya-gamma-predictions_', version, '_overdispersed.RDS')))

# load the preds grid
locs_grid <- readRDS(here::here('data', paste0('grid_', version, '.RDS')))

pi_mean <- apply(preds$pi, c(2, 3, 4), mean)
pi_sd <- apply(preds$pi, c(2, 3, 4), sd)

# modify this based on Alissa's input
dimnames(pi_mean) <- list(
    location = 1:dim(pi_mean)[1],
    species = species_names,
    time = 1:dim(pi_mean)[3]
)


dat_pi_mean <- as.data.frame.table(pi_mean, responseName = "pi_mean") %>%
    mutate(location = as.numeric(location), time = as.numeric(time)) %>%
    left_join(locs_grid %>%
                  mutate(location = 1:dim(locs_grid)[1]))

dimnames(pi_sd) <- list(
    location = 1:dim(pi_sd)[1],
    species = species_names,
    time = 1:dim(pi_sd)[3]
)

dat_pi_sd <- as.data.frame.table(pi_sd, responseName = "pi_sd") %>%
    mutate(location = as.numeric(location), time = as.numeric(time)) %>%
    left_join(locs_grid %>%
                  mutate(location = 1:dim(locs_grid)[1]))


time_bins <- seq(-285, 21000, by=990)[-1]
time_vec <- paste(time_bins[1:21], "ybp")
names(time_vec) <- paste(1:21)
time_vec
