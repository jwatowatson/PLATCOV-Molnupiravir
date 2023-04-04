//   Copyright 2022 James Watson, Mahidol Oxford Tropical Medicine Research Unit
//
//   Licensed under a CC-BY License
/**
Model 3 with the following characteristics:
- analysis on the copies per ml scale (batch effect adjustement done without uncertainty propagation)
- Adjustment for RNaseP
- no covariate adjustment
- non-linear function

**/



data {
  // Patient data
  int<lower=0> Ntot;                           // Number of PCR data points
  int<lower=0,upper=Ntot> N_obs;               // Number of PCR data points with viral load less than LOD
  int<lower=0> n_id;                           // Number of individuals
  int<lower=1,upper=n_id> id[Ntot];            // Patient identifier for each PCR sample
  int<lower=1,upper=Ntot> ind_start[n_id];     // Starting index for each patient: used in generated quantities block
  real<lower=0> Time_max;                      // Max follow-up time
  real<lower=0,upper=Time_max> obs_day[Ntot];  // Time since randomisation for sample
  real log_10_vl[Ntot];                        // log base 10 viral load in copies per mL
  real log10_cens_vl[Ntot];                    // censoring value for censored observation
  real RNaseP[Ntot];                           // Scaled RNaseP CT values (mean 0)
  int<lower=1> K_trt;                          // Number of treatment arms
  matrix[Ntot,K_trt] trt_mat;                  // Trt matrix
  int<lower=0> K_cov_intercept;                // number of columns in covariate design matrix for intercept
  matrix[Ntot,K_cov_intercept] x_intercept;    // covariate design matrix for regression onto intercept
  int<lower=0> K_cov_slope;                    // number of columns in covariate design matrix for slope
  matrix[Ntot,K_cov_slope] x_slope;            // covariate design matrix for regression onto slope
  
  // priors
  real alpha_0_prior_mean; // prior mean intercept
  real alpha_0_prior_sd;   // prior sd intercept
  
  real beta_0_prior_mean;  // prior mean slope
  real beta_0_prior_sd;    // prior sd slope
  
  real trt_effect_sd;   // prior sd on treatment effect
  
  real sigma_logvl_mean;   // prior mean of sd in error model
  real sigma_logvl_sd;     // prior sd of sd in error model
  
  real slope_coefs_sd;     // prior sd for prior on covariate effect on slope
  real intercept_coefs_sd; // prior sd for prior on covariate effect on intercept
  
  real t_max_pop_mean;     // prior on mean time since peak viral load
  real t_max_pop_sd;       // priors on sd for mean time since peak viral load
}

transformed data {
  vector[3] zeros3;
  for(i in 1:3) zeros3[i] = 0;
}

parameters {
  // hyperparameters
  cholesky_factor_corr[3] L_Omega;     // correlation matrix
  vector<lower=0>[3] sigmasq_u;        // variance of random effects
  
  // Measurement error
  real<lower=0> sigma_logvl;
  
  // Population parameters
  real alpha_0;                             // population intercept at tmax
  real tmax_pop;                            // population intercept at tmax
  vector<lower=0>[2] beta_0;                // growth and decline rates
  real gamma_rnasep;                        // Adjustment for RNaseP
  vector[K_trt-1] trt_effect;               // Estimates of the treatment effect
  vector[K_cov_slope] slope_coefs;          // Covariate effects on slope
  vector[K_cov_intercept] intercept_coefs;  // Covariate effects on intercept
  
  // Random effects
  vector[3] theta_rand_id[n_id];            // individual random effects vector
  
  // Degrees of freedom for the t-distribution error model
  real<lower=0> t_dof;
}

transformed parameters {
  real pred_log10_vl[Ntot];
  vector[Ntot] trt_slope;
  vector[Ntot] beta_cov;
  vector[Ntot] alpha_cov;
  
  {
    vector[K_trt] trt_effect_prime;
    
    trt_effect_prime = append_row(0, trt_effect);
    trt_slope = trt_mat * trt_effect_prime;
    
    // make individual covariate transform
    beta_cov = x_slope*slope_coefs;
    alpha_cov = x_intercept*intercept_coefs;
    
    for(i in 1:Ntot){
      real intercept = alpha_0 + theta_rand_id[id[i]][1] + alpha_cov[i];
      real a = beta_0[1] * exp(trt_slope[i]); // no random effect for the first term as does not converge
      real b = beta_0[2] * exp(theta_rand_id[id[i]][2] + trt_slope[i] + beta_cov[i]);
      real tmax = tmax_pop + theta_rand_id[id[i]][3];
      pred_log10_vl[i] = gamma_rnasep*RNaseP[i]+intercept+log(a+b)-
      log_sum_exp(log(b)-a*(obs_day[i]-tmax),log(a)+b*(obs_day[i]-tmax));
    }
  }
}

model {
  //***** Prior *****
  // error model - degrees of freedom
  t_dof ~ exponential(1);
  // error model variance
  sigma_logvl ~ normal(sigma_logvl_mean, sigma_logvl_sd) T[0,];
  
  // random effects
  sigmasq_u[1] ~ exponential(1);
  sigmasq_u[2] ~ exponential(1);
  sigmasq_u[3] ~ normal(1,.5);
  L_Omega ~ lkj_corr_cholesky(3);  // covariance matrix - random effects for individs
  // individual random effects
  for(i in 1:n_id) theta_rand_id[i] ~ multi_normal_cholesky(zeros3, diag_pre_multiply(sigmasq_u, L_Omega));
  
  // Population parameters
  alpha_0 ~ normal(alpha_0_prior_mean,alpha_0_prior_sd);
  beta_0[1] ~ normal(1.5,.5);
  beta_0[2] ~ normal(beta_0_prior_mean,beta_0_prior_sd);
  gamma_rnasep ~ normal(0,1);
  slope_coefs ~ normal(0,slope_coefs_sd);
  intercept_coefs ~ normal(0,intercept_coefs_sd);
  tmax_pop ~ normal(t_max_pop_mean, t_max_pop_sd);
  trt_effect ~ normal(0,trt_effect_sd);   // Treatment effect
  
  //***** Likelihood *****
  // Non censored observations
  log_10_vl[1:N_obs] ~ student_t(t_dof, pred_log10_vl[1:N_obs], sigma_logvl);
  
  // Censored observations
  for(i in (N_obs+1):Ntot){
    target += student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }
}

generated quantities {
  real preds[Ntot]; // For plotting
  vector[Ntot] log_lik;
  vector[n_id] slope;
  
  for(i in 1:N_obs){
    preds[i] = pred_log10_vl[i]-gamma_rnasep*RNaseP[i];
    log_lik[i] = student_t_lpdf(log_10_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }
  for(i in (N_obs+1):Ntot){
    preds[i] = pred_log10_vl[i]-gamma_rnasep*RNaseP[i];
    log_lik[i] = student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }
  for(i in 1:n_id){
    int j = ind_start[i];
    slope[i] = beta_0[2] * exp(theta_rand_id[id[j]][2] + trt_slope[j] + beta_cov[j]);
  }
}
