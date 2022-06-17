if(!dir.exists(here::here("output"))) {
  dir.create(here::here("output"))
}
library(raster)
library(geosphere)
library(tidyverse)
library(sp)
library(rgdal)
library(rdist)
library(rgeos)
library(ggplot2)
library(invgamma)
library(mvnfast)
library(splines)
library(pgdraw)
library(fields)
library(geoR)
library(dplyr)
library(data.table)
library(BayesMRA)
# install.packages("/home/adawson/Documents/projects/pgR", repos=NULL, type="source")
library(pgR)

version='5.0'

#### DATA PREP #### -- renamed version 5.0
# y <- readRDS(here::here('data', paste0('paleo_pollen_dat_', version, '.RDS')))
# taxa.keep <- readRDS(here::here('data', paste0('pollen_taxa_', version, '.RDS')))
# locs <- readRDS(here::here('data', paste0('paleo_pollen_locs_', version, '.RDS')))
y <- readRDS(here::here('data', paste0('pollen_dat_', version, '.RDS')))
taxa.keep <- readRDS(here::here('data', paste0('taxa_', version, '.RDS')))
locs <- readRDS(here::here('data', paste0('pollen_locs_', version, '.RDS')))
rescale <- 1e3

#### RUNNING THE MODEL & SAVING OUTPUT####
# scale locations so distance values aren't too large (from 1m to 100km units)
locs_scaled <- locs/rescale
N_locs_dat = nrow(locs_scaled)
X <- matrix(rep(1, N_locs_dat),N_locs_dat, 1)

params <- default_params()
params$n_adapt <- 2000
params$n_mcmc <- 5000
params$n_message <- 10
params$n_thin <- 5
priors <- default_priors_pg_stlm(y, X, corr_fun = "matern")

J <- ncol(y)


# XXX: need this to work around undefined variable in pgSPLM
d <- ncol(y)
Y = y
X = as.matrix(X)
locs = as.matrix(locs_scaled)
n_cores = 1L
n_chain = 1

for (n in 1:dim(Y)[1]){
  for (t in 1:dim(Y)[3]){
    if (!any(is.na(Y[n,,t]))){
      if (sum(Y[n,,t])==0){
        #print(Y[n,,t])
        Y[n,,t] = rep(NA, dim(Y)[2])
      }
    }
  }
}


# check the MRA grid
M <- 3
# M <- 2
n_coarse_grid <- 25
# n_coarse_grid <- 10
MRA <- mra_wendland_2d(locs, M = M, n_coarse_grid = n_coarse_grid)
# number of basis functions
sum(MRA$n_dims)

# plot the MRA grid
# plot_MRA_grid(MRA)

# code to run matern model
if (!file.exists( here::here('output', paste0('polya-gamma-posts_', version, '_MRA.RDS')))) {
  
  out <- pg_stlm_mra(Y = Y,
                     X = X,
                     locs = locs,
                     params,
                     priors,
                     n_cores = n_cores,
                     n_chain = n_chain,
                     M = M,
                     n_coarse_grid = n_coarse_grid)
  saveRDS(out, here::here('output', paste0('polya-gamma-posts_', version, '_MRA.RDS')),
          compress = FALSE)
  
  pushoverr::pushover(message = "Finished fitting MRA model")
}

dat <- list(y = y,
            X = X, 
            locs = locs_scaled,
            rescale = rescale,
            taxa.keep = taxa.keep)

if (!file.exists(here::here('output', paste0('polya-gamma-dat_', version,'.RDS')))) {
  saveRDS(dat, here::here('output', paste0('polya-gamma-dat_', version,'.RDS')))
}

