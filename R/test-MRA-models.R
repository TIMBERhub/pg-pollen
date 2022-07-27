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
params$n_adapt <- 10#0
params$n_mcmc <- 10#0
params$n_message <- 1
params$n_thin <- 1
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
n_coarse_grid <- 25
MRA <- mra_wendland_2d(locs, M = M, n_coarse_grid = n_coarse_grid)
# number of basis functions
sum(MRA$n_dims)

# plot the MRA grid
# plot_MRA_grid(MRA)

# code to run mra model
start1 <- Sys.time()
profvis::profvis(out <- pg_stlm_mra(Y = Y,
                                    X = X,
                                    locs = locs,
                                    params,
                                    priors,
                                    n_cores = n_cores,
                                    n_chain = n_chain,
                                    M = M,
                                    n_coarse_grid = n_coarse_grid, 
                                    verbose = TRUE,
                                    store_R = TRUE),
                 prof_output = here::here("output", "MRA1.Rprof"))
stop1 <- Sys.time()
time1 <- stop1 - start1
## 16 minutes
## 11 minutes after revising sampler for rho



# code to run mra model
start2 <- Sys.time()
profvis::profvis(out <- pg_stlm_mra(Y = Y,
                                    X = X,
                                    locs = locs,
                                    params,
                                    priors,
                                    n_cores = n_cores,
                                    n_chain = n_chain,
                                    M = M,
                                    n_coarse_grid = n_coarse_grid,
                                    verbose = TRUE,
                                    store_R = FALSE),
                 prof_output = here::here("output", "MRA2.Rprof"))
stop2 <- Sys.time()
time2 <- stop2 - start2
## 38 minutes
## 33 minutes after revising sampler for rho


library(tidyverse)
# profvis::profvis(prof_input = here::here("output", "MRA1.Rprof"))
# profvis::profvis(prof_input = here::here("output", "MRA2.Rprof"))
dat <- summaryRprof(here::here("output", "MRA1.Rprof"))
dat2 <- summaryRprof(here::here("output", "MRA2.Rprof"))
# rbind(mutate(dat$by.total, fit="a"), mutate(dat2$by.total, fit = "b")) %>%
#     arrange(desc(total.pct)) %>%
#     head(n  = 10)
dat$by.total %>% arrange(self.pct)
dat$by.total %>% arrange(total.pct)
dat2$by.total %>% arrange(self.pct)
dat2$by.total %>% arrange(total.pct)

pushoverr::pushover(message = "Finished fitting MRA test models")
