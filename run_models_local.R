#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Packages needed
library(rstan)
library(matrixStats)
library(doParallel)


load(paste0('Rout/model_run_setup_',args[1],'.RData'))

writeLines(sprintf('Doing analysis for %s......', args[1]))

for(i in 1:nrow(model_settings)){
  
  writeLines('Doing the following job:')
  print(model_settings[i, ])
  
  options(mc.cores = model_settings$Nchain[i])
  stopifnot(model_settings$Nchain[i]>getDoParWorkers()) # check worker number assigned
  
  mod = stan_model(file = as.character(model_settings$mod[i])) # compile 
  
  stan_input_job = stan_inputs
  
  analysis_data_stan = stan_input_job$analysis_data_stan
  analysis_data_stan$trt_mat = stan_input_job$Trt_matrix
  analysis_data_stan$K_trt = ncol(analysis_data_stan$trt_mat)
  
  x_intercept = stan_input_job$cov_matrices$X_int[[model_settings$cov_matrices[i]]]
  if(ncol(x_intercept)==0) x_intercept = array(0, dim=c(nrow(x_intercept),1))
  analysis_data_stan$x_intercept = x_intercept
  analysis_data_stan$K_cov_intercept= ncol(x_intercept)
  
  
  x_slope = stan_input_job$cov_matrices$X_slope[[model_settings$cov_matrices[i]]]
  if(ncol(x_slope)==0) x_slope = array(0, dim=c(nrow(x_slope),1))
  analysis_data_stan$x_slope = x_slope
  analysis_data_stan$K_cov_slope=ncol(x_slope)
  
  
  # sample posterior
  out = sampling(mod, 
                 data=c(analysis_data_stan,
                        all_priors[[model_settings$prior[i]]]),
                 iter=model_settings$Niter[i],
                 chain=model_settings$Nchain[i],
                 thin=model_settings$Nthin[i],
                 warmup=model_settings$Nwarmup[i],
                 save_warmup = FALSE,
                 seed=i,
                 pars=c('L_Omega'), # we don't save this as it takes up lots of memory!
                 include=FALSE)
  
  
  save(out, file = paste0('Rout/model_fits_',i,'_',args[1],'.RData'))# save output
  
  writeLines('Finished job')
}
