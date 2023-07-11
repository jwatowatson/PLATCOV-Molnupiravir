# list of functions for plotting data
plot_effect_estimates = function(effect_ests, #list of stan outputs
                                 plot_models, # indices of models to plot in list
                                 my_pch=1,
                                 mod_cols = NULL,
                                 study_threshold){
  
  if (length(my_pch)==1) my_pch = (1:length(plot_models))+15
  if (length(my_pch)!=length(plot_models)) stop('length of my_pch needs to be the same as the number of input models')
  if (is.null(mod_cols)){
    mod_cols = brewer.pal(n = length(plot_models), name = 'Set1')
    if(length(plot_models)>9) writeLines('too many models - supply user defined colors')
  }
  if(length(mod_cols)!=length(plot_models)) stop('number of colors needs to be equal to number of models to be plotted')
  
  K_treatments = nrow(effect_ests[[plot_models[1]]])
  
  xlims = (exp(range(c(0, range(sapply(effect_ests[plot_models], rbind)))) )-1)*100
  x_points = pretty((xlims),6)
  plot(NA, NA, xlim = range(x_points),
       ylim = c(0.75,K_treatments+.25),
       panel.first=grid(), ylab='', yaxt='n', type='n',
       xlab = 'Change in rate of clearance (%)',
       xaxt = 'n')
  
  index_p = rev(seq(-.2,.2, length.out = length(plot_models)))
  abline(v=0,lwd=2)
  polygon(c(-1000, 100*(study_threshold-1), 100*(study_threshold-1), -1000),
          c(-100, -100, 100, 100), border = NA,
          col = adjustcolor('grey',.4))
  
  sort_ind = order(effect_ests[[plot_models[1]]][,'50%'])
  for(i in 1:length(plot_models)){
    points((exp(effect_ests[[plot_models[i]]][sort_ind,'50%'])-1)*100,
           1:K_treatments+index_p[i],pch=my_pch[i],
           col=mod_cols[i],cex=1.5)
    for(j in 1:length(sort_ind)){
      kk = sort_ind[j]
      lines((exp(c(effect_ests[[plot_models[i]]][kk,'2.5%'],
                   effect_ests[[plot_models[i]]][kk,'97.5%']))-1)*100,
            rep(j+index_p[i],2),col=mod_cols[i],lwd=1)
      lines((exp(c(effect_ests[[plot_models[i]]][kk,'10%'],
                   effect_ests[[plot_models[i]]][kk,'90%']))-1)*100,
            rep(j+index_p[i],2),col=mod_cols[i],lwd=3)
    }
  }
  axis(1, at = x_points)
  axis(2, at = 1:K_treatments, labels = rownames(effect_ests[[plot_models[1]]])[sort_ind])
  
}

plot_baseline_data = function(input_data){
  
  baseline_ind = input_data$Timepoint_ID==0
  bvl = aggregate(log10_viral_load ~ ID, input_data[baseline_ind, ], median)
  
  hist(bvl$log10_viral_load,
       breaks = seq(1,9,by=.5),
       xlab='Baseline viral load (SARS CoV2 genomes/mL)',
       ylab ='Number of patients',xlim=c(1,9),
       main='', xaxt ='n')
  axis(1, at = c(2,4,6,8), labels = c(expression(10^2),
                                      expression(10^4),
                                      expression(10^6),
                                      expression(10^8)))
  grid(); 
}


plot_serial_data = function(xx, xlims=c(0,7), plot_points=T, alpha.f=.3){
  
  daily_VL_data = xx %>% group_by(ID, Time) %>%
    mutate(daily_VL = mean(log10_viral_load,na.rm = T))
  
  summary_VL_data = daily_VL_data %>% group_by(Timepoint_ID, Trt) %>%
    summarise(daily_VL = mean(daily_VL,na.rm = T),
              trt_color=unique(trt_color))
  
  summary_dat = daily_VL_data %>% ungroup() %>% distinct(ID, .keep_all = T) %>%
    group_by(Trt) %>%
    summarise(n = n(),
              trt_color = unique(trt_color)) %>% ungroup() %>%
    mutate(legend = paste0(Trt, ' (n=',n,')'))
  
  par(las=1)
  plot(summary_VL_data$Timepoint_ID, summary_VL_data$daily_VL,
       ylab = 'SARS-CoV-2 genomes/mL', panel.first=grid(),
       xlab = 'Time since randomisation (days)',
       xlim = xlims, yaxt='n',type='n',
       ylim = c(0.7, 8))
  axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
                                      expression(10^4),
                                      expression(10^6),
                                      expression(10^8)))
  
  if(plot_points) {
    points(daily_VL_data$Time, daily_VL_data$daily_VL,
           col= adjustcolor(daily_VL_data$trt_color,alpha.f = alpha.f))
  }
  summary_VL_data$Trt_pch = as.numeric(as.factor(summary_VL_data$Trt))
  for(tt in unique(summary_VL_data$Trt)){
    ind = summary_VL_data$Trt==tt
    lines(summary_VL_data$Timepoint_ID[ind],
          summary_VL_data$daily_VL[ind], pch=summary_VL_data$Trt_pch[ind],
          type='b',col = summary_VL_data$trt_color[ind],lwd=3)
  }
  
  legend('topright', col = summary_dat$trt_color, 
         legend = summary_dat$legend,
         lwd=2, pch=summary_VL_data$Trt_pch,
         cex=1, inset = 0.03)
}


bayes_R2 = function(mod_preds, mod_residuals) {
  var_pred = apply(mod_preds, 1, var)
  var_res = apply(mod_residuals, 1, var)
  var_pred / (var_pred + var_res)
}


make_stan_inputs = function(input_data_fit, 
                            int_covs_base,
                            int_covs_full,
                            slope_covs_base,
                            slope_covs_full,
                            trt_frmla,
                            epoch=NA,
                            Dmax
){
  
  ## check censored values come last
  if(!all(diff(input_data_fit$log10_viral_load == input_data_fit$log10_cens_vl)>=0)) stop()
  ind_dup = !duplicated(input_data_fit$ID)
  
  for(ll in unique(input_data_fit$Lab)){
    ind = input_data_fit$Lab==ll
    med_val = median(input_data_fit$CT_RNaseP[ind], na.rm = T)
    input_data_fit$CT_RNaseP[ind & is.na(input_data_fit$CT_RNaseP)]=med_val
    input_data_fit$CT_RNaseP[ind] = input_data_fit$CT_RNaseP[ind]-med_val
  }
  input_data_fit$RnaseP_scaled = input_data_fit$CT_RNaseP
  
  input_data_fit$Age_scaled = (input_data_fit$Age-mean(input_data_fit$Age[ind_dup]))/sd(input_data_fit$Age[ind_dup])
  
  # make the covariate matrix
  # check no missing data
  if(!all(!apply(input_data_fit[, union(int_covs_full,slope_covs_full)], 2, function(x) any(is.na(x))))){
    stop('Missing data in covariate matrix!')
  }
  
  ind_contr = which(apply(input_data_fit[, int_covs_base,drop=F], 2, 
                          function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_intcpt_1 = model.matrix( ~ ., 
                               data = input_data_fit[, int_covs_base[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_intcpt_1 = array(dim = c(nrow(input_data_fit),0))
  }
  
  ind_contr = which(apply(input_data_fit[, int_covs_full,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_intcpt_2 = model.matrix( ~ ., 
                               data = input_data_fit[, int_covs_full[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_intcpt_2 = array(dim = c(nrow(input_data_fit),0))
  }
  
  ind_contr = which(apply(input_data_fit[, slope_covs_base,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_slope_1 = model.matrix( ~ ., 
                              data = input_data_fit[, slope_covs_base[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_slope_1 = array(dim = c(nrow(input_data_fit),0))
  }
  
  ind_contr = which(apply(input_data_fit[, slope_covs_full,drop=F], 2, function(x) length(unique(x))>1))
  if(length(ind_contr)>0){
    X_slope_2 = model.matrix( ~ ., 
                              data = input_data_fit[, slope_covs_full[ind_contr],drop=F])[, -1, drop=F]
  } else {
    X_slope_2 = array(dim = c(nrow(input_data_fit),0))
  }
  
  if(!nrow(X_intcpt_1) == nrow(input_data_fit)) stop()
  if(!nrow(X_slope_1) == nrow(input_data_fit)) stop()
  
  cov_matrices = list(X_int=list(X_intcpt_1, X_intcpt_2),
                      X_slope=list(X_slope_1, X_slope_2))
  
  ID_map = data.frame(ID_key = input_data_fit$ID,
                      ID_stan = as.numeric(as.factor(input_data_fit$ID)))
  writeLines(sprintf('There are a total of %s patients in the database with a total of %s PCRs analysable',
                     max(ID_map$ID_stan),
                     nrow(input_data_fit)))
  
  ind_cens = !input_data_fit$log10_viral_load>
    input_data_fit$log10_cens_vl
  
  writeLines(sprintf('%s%% of samples are below LOD',
                     round(100*mean(ind_cens),digits = 2)))
  
  analysis_data_stan = list(Ntot = nrow(input_data_fit),
                            N_obs = sum(!ind_cens),
                            n_id = max(ID_map$ID_stan),
                            id = ID_map$ID_stan,
                            ind_start = which(!duplicated(ID_map$ID_stan)),
                            obs_day = input_data_fit$Time,
                            log_10_vl = input_data_fit$log10_viral_load,
                            log10_cens_vl = input_data_fit$log10_cens_vl,
                            RNaseP = input_data_fit$RnaseP_scaled,
                            Time_max = Dmax)
  if(!is.na(epoch)){
    analysis_data_stan$epoch = as.numeric(as.factor(input_data_fit$Epoch))
    analysis_data_stan$K_epoch = max(analysis_data_stan$epoch)
  }
  ID_map = ID_map[!duplicated(ID_map$ID_key), ]
  
  writeLines('check stan data formatting:')
  all(analysis_data_stan$log_10_vl[1:analysis_data_stan$N_obs]>
        analysis_data_stan$log10_cens_vl[1:analysis_data_stan$N_obs]) &
    all(analysis_data_stan$log_10_vl[(1+analysis_data_stan$N_obs):analysis_data_stan$Ntot] ==
          analysis_data_stan$log10_cens_vl[(1+analysis_data_stan$N_obs):analysis_data_stan$Ntot])
  
  Trt_matrix = model.matrix(trt_frmla, data = input_data_fit)
  Trt_matrix[,1]=0 # first column is dummy
  
  analysis_inputs = list(cov_matrices=cov_matrices,
                         analysis_data_stan=analysis_data_stan,
                         Trt_matrix=Trt_matrix,
                         ID_map=ID_map)
  return(analysis_inputs)
}


plot_variants = function(platcov_dat){
  platcov_dat = platcov_dat %>%
    mutate(year_month = paste(year(Rand_date), month(Rand_date), sep = '_')) %>%
    distinct(ID, .keep_all = T)
  xx=aggregate(Variant ~ year_month, data = platcov_dat, table)
  
}


get_rates_linear_mod = function(mod_out, # single model fit - not a list
                                analysis_data_stan){
  
  # get the indices of first datapoint for each individual
  ind_id = which(!duplicated(analysis_data_stan$id))
  thetas = 
    rstan::extract(mod_out, 
                   pars = c('beta_0','beta_cov','theta_rand_id','trt_effect'))
  beta_cov = x_slope*slope_coefs;
  
  beta_0*exp(trt_slope[i]+theta_rand_id[id[i]][2]+beta_cov[i])
}

# plot_data_model_fits = function(mod_out, # model fits
#                              models_plot, # which models to plot
#                              K_plots,
#                              mod_cols,
#                              ID_map,
#                              analysis_data_stan
# ){
#   
#   # extract posterior parameters and outputs
#   thetas = list()
#   for(mm in 1:length(mod_out)){
#     thetas[[mm]] = rstan::extract(mod_out[[mm]])
#   }
#   counter = 1
#   
#   ID_map$Trt = gsub(pattern = '\n',
#                     replacement = '',
#                     x = ID_map$Trt,fixed = T)
#   while(counter <= nrow(ID_map)){
#     
#     # draw individual model fit with data
#     ind = analysis_data_stan$id==ID_map$ID_stan[counter]
#     plot(analysis_data_stan$obs_day[ind],
#          analysis_data_stan$log_10_vl[ind],
#          xlab='', ylab='', 
#          xaxt='n', yaxt='n',
#          panel.first=grid(), xlim=c(0,7),
#          ylim = range(analysis_data_stan$log_10_vl))
#     # if(counter %% sqrt(K_plots) == 1){
#     #   mtext(text = 'RNA copies per mL',side = 2,
#     #         line = 3,las = 3)
#     # }
#     axis(1, at = c(0,3,7))
#     axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
#                                         expression(10^4),
#                                         expression(10^6),
#                                         expression(10^8)))
#     # if((counter%%K_plots) >= K_plots - sqrt(K_plots)){
#     #   mtext(text = 'Days',side = 1,line = 2)
#     # }
#     for(mm in models_plot){
#       ix = order(analysis_data_stan$obs_day[ind])
#       my_xs = analysis_data_stan$obs_day[ind][ix]
#       polygon(x = c(my_xs, rev(my_xs)),
#               y = c(apply(thetas[[mm]]$preds[,ind],2,
#                           quantile,probs=0.025)[ix],
#                     rev(apply(thetas[[mm]]$preds[,ind],2,
#                               quantile,probs=0.975)[ix])),
#               border = NA, 
#               col = adjustcolor(mod_cols[mm],alpha.f = .3))
#       lines(my_xs,
#             colMeans(thetas[[mm]]$preds[,ind])[ix],
#             col = mod_cols[mm],lwd=2)
#     }
#     points(analysis_data_stan$obs_day[ind],
#            analysis_data_stan$log_10_vl[ind],pch=16)
#     
#     mtext(text = paste0(ID_map$ID[counter],
#                         '\n',
#                         ID_map$Trt[counter]),
#           side = 3, line = -0.5, cex=0.8)
#     counter=counter+1
#   }
#   
# }

plot_data_model_fits = 
  function(model_list, # list of model fits
           models_to_plot, # which models to plot
           K_plots,
           mod_cols,
           ID_map,
           analysis_data_stan
  ){
    
    # extract posterior parameters and outputs
    thetas = list()
    for(mm in 1:length(model_list)){
      thetas[[mm]] = rstan::extract(model_list[[mm]], pars='preds')
    }
    counter = 1
    
    ID_map$Trt = gsub(pattern = '\n',
                      replacement = '',
                      x = ID_map$Trt,fixed = T)
    while(counter <= nrow(ID_map)){
      
      # draw individual model fit with data
      ind = analysis_data_stan$id==ID_map$ID_stan[counter]
      plot(analysis_data_stan$obs_day[ind],
           analysis_data_stan$log_10_vl[ind],
           xlab='', ylab='', 
           xaxt='n', yaxt='n',
           panel.first=grid(), xlim=c(0,7),
           ylim = range(analysis_data_stan$log_10_vl))
      
      axis(1, at = c(0,3,7))
      axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
                                          expression(10^4),
                                          expression(10^6),
                                          expression(10^8)))
      
      for(mm in models_to_plot){
        ix = order(analysis_data_stan$obs_day[ind])
        my_xs = analysis_data_stan$obs_day[ind][ix]
        polygon(x = c(my_xs, rev(my_xs)),
                y = c(apply(thetas[[mm]]$preds[,ind],2,
                            quantile,probs=0.025)[ix],
                      rev(apply(thetas[[mm]]$preds[,ind],2,
                                quantile,probs=0.975)[ix])),
                border = NA, 
                col = adjustcolor(mod_cols[mm],alpha.f = .3))
        lines(my_xs,
              colMeans(thetas[[mm]]$preds[,ind])[ix],
              col = mod_cols[mm],lwd=2)
      }
      points(analysis_data_stan$obs_day[ind],
             analysis_data_stan$log_10_vl[ind],pch=16)
      
      mtext(text = paste0(ID_map$ID_key[counter],
                          ' ',
                          ID_map$Trt[counter]),
            side = 3, line = 0.5, cex=0.8)
      counter=counter+1
    }
    
  }



plot_individ_curves = function(platcov_dat, IDs, xlims){
  
  if(!all(IDs %in% platcov_dat$ID)) stop('missing IDs!!')
  ylims = range(platcov_dat$log10_viral_load)
  platcov_dat$Day = floor(platcov_dat$Time)
  
  for(id in unique(IDs)){
    
    ind = platcov_dat$ID==id
    plot(platcov_dat$Time[ind],
         platcov_dat$log10_viral_load[ind],
         pch = as.numeric(platcov_dat$log10_cens_vl[ind]==platcov_dat$log10_viral_load[ind])+1,
         xlab='', ylab='', 
         xaxt='n', yaxt='n',
         panel.first=grid(), xlim=xlims,
         ylim = ylims)
    # title(paste(id, platcov_dat$Trt[ind][1]))
    title(id)
    lines(platcov_dat$Time[ind], platcov_dat$daily_VL[ind], lwd=2, lty=2,col='grey')
    points(platcov_dat$Time[ind],
           platcov_dat$log10_viral_load[ind],
           pch = as.numeric(platcov_dat$log10_cens_vl[ind]==platcov_dat$log10_viral_load[ind])+16)
    axis(1, at = c(0,7,14,21))
    axis(2, at = c(2,4,6,8), labels = c(expression(10^2),
                                        expression(10^4),
                                        expression(10^6),
                                        expression(10^8)))
    
  }
  
}


plot_coef_effects = function(stan_out, cov_mat, stan_inputs){
  
  thetas = rstan::extract(stan_out)
  alpha_coefs = apply(thetas$intercept_coefs,2,
                      quantile,probs=c(0.025,.1,.5,.9,0.975))
  xlims=range(alpha_coefs)
  
  cov_names_intercept =
    plyr::mapvalues(x = colnames(stan_inputs$cov_matrices$X_int[[cov_mat]]),
                    from = c('Age_scaled','Antibody_test',
                             'Symptom_onset','N_dose'),
                    to = c('Age','Serology RDT',
                           'Days since\nsymptom onset',
                           'Number of\nvaccine doses'))
  
  cov_names_slope =
    plyr::mapvalues(x = colnames(stan_inputs$cov_matrices$X_slope[[cov_mat]]),
                    from = c('Age_scaled','Antibody_test',
                             'Symptom_onset','N_dose'),
                    to = c('Age','Serology RDT',
                           'Days since\nsymptom onset',
                           'Number of\nvaccine doses'))
  
  plot(alpha_coefs['50%', ], 1:ncol(alpha_coefs),
       xlim=xlims,yaxt='n',ylab='',bty='n',xaxt='n',
       panel.first=grid(), xlab='Intercept (fold change)')
  abline(v=0,lty=2,lwd=2)
  for(i in 1:ncol(alpha_coefs)){
    lines(c(alpha_coefs['10%',i], alpha_coefs['90%',i]),
          c(i,i), lwd=3)
    lines(c(alpha_coefs['2.5%',i], alpha_coefs['97.5%',i]),
          c(i,i), lwd=1)
  }
  axis(2, at =1:ncol(alpha_coefs), labels = cov_names_intercept,tick = F)
  
  x_points = signif(10^seq(xlims[1], xlims[2],length.out = 5),2)
  axis(1, at = log10(x_points), labels = x_points)
  
  beta_coefs = apply(thetas$slope_coefs,2,quantile,
                     probs=c(0.025,.1,.5,.9,0.975))
  xlims=range(beta_coefs)
  plot(beta_coefs['50%', ], 1:ncol(beta_coefs),
       xlim=xlims,yaxt='n',ylab='',bty='n',xaxt='n',
       panel.first=grid(), xlab='Slope (multiplicative effect)')
  abline(v=0,lty=2,lwd=2)
  for(i in 1:ncol(beta_coefs)){
    lines(c(beta_coefs['10%',i], beta_coefs['90%',i]),
          c(i,i), lwd=3)
    lines(c(beta_coefs['2.5%',i], beta_coefs['97.5%',i]),
          c(i,i), lwd=1)
  }
  axis(2, at =1:ncol(beta_coefs), labels = cov_names_slope, tick = F)
  x_points = signif(10^seq(xlims[1], xlims[2],length.out = 5),2)
  axis(1, at = log10(x_points), labels = x_points)
}

# Checks a function for use of global variables
# Returns TRUE if ok, FALSE if globals were found.
checkStrict <- function(f, silent=FALSE) {
  vars <- codetools::findGlobals(f)
  found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
  if (!silent && any(found)) {
    warning("global variables used: ", paste(names(found)[found], collapse=', '))
    return(invisible(FALSE))
  }
  
  !any(found)
}

calculate_fever_clearance = function(temp_dat,
                                     window_clear = 24/24, # look ahead window to define "fever clearance"
                                     threshold=37){
  
  if(!'ax_temperature' %in% colnames(temp_dat)) stop('needs to contain a ax_temperature column')
  
  temp_dat$clearance_time = NA
  # For interval censored data, the status indicator is 0=right censored, 1=event at time, 2=left censored, 3=interval censored. 
  temp_dat$clearance_time_cens = 1
  
  temp_dat$fever_binary = temp_dat$ax_temperature>threshold
  temp_dat = dplyr::arrange(temp_dat, ID, Time) 
  temp_dat = temp_dat[!is.na(temp_dat$ax_temperature), ]
  
  for(id in unique(temp_dat$ID)){
    ind = temp_dat$ID==id
    if(all(!temp_dat$fever_binary[ind])){ # never fever
      temp_dat$clearance_time[ind]=0
    } else if(all(temp_dat$fever_binary[ind])){ # always fever
      writeLines(sprintf('all fever for %s with %s FUP points',id,sum(ind)))
      temp_dat$clearance_time[ind] = max(temp_dat$Time[ind])
      temp_dat$clearance_time_cens[ind] = 0 #censored obs
    } else { # fever cleared
      j_cleared = which(ind & !temp_dat$fever_binary)
      check_ahead=F
      for(j in j_cleared){
        if(!check_ahead){
          ind_check = 
            which(ind & 
                    temp_dat$Time>temp_dat$Time[j] &
                    temp_dat$Time<temp_dat$Time[j] + window_clear)
          if(length(ind_check)>0 & all(!temp_dat$fever_binary[ind_check])){
            temp_dat$clearance_time[ind]=temp_dat$Time[j]
            check_ahead=T
          }
        }
      }
      if(!check_ahead){
        temp_dat$clearance_time[ind]=tail(temp_dat$Time[ind],1)
        temp_dat$clearance_time_cens[ind]=0
      }
    }
  }
  
  return(temp_dat[!duplicated(temp_dat$ID), ])
}



make_slopes_plot = function(stan_out, 
                            analysis_data_stan,
                            ID_map, 
                            data_summary,
                            my_lims = c(5, 72), # hours
                            my_vals = c(7,24,48,72)){
  
  slopes = -abs(rstan::extract(stan_out, pars='slope')$slope)
  
  t12_output = data.frame(t_12_med = 24*log10(.5)/(apply(slopes,2,mean)),
                          t_12_up = 24*log10(.5)/(apply(slopes,2,quantile,.9)),
                          t_12_low = 24*log10(.5)/(apply(slopes,2,quantile,.1)),
                          slope = apply(slopes,2,mean),
                          ID_stan = analysis_data_stan$id[analysis_data_stan$ind_start])
  t12_output = merge(t12_output, ID_map, by = 'ID_stan')
  data_summary = merge(data_summary, t12_output, by.x = 'ID', by.y = 'ID_key')
  
  data_summary = dplyr::arrange(data_summary, Trt, t_12_med)
  
  par(mar=c(5,2,2,2))
  plot(data_summary$t_12_med, 1:nrow(data_summary),
       yaxt='n', xaxt='n',
       xlim=my_lims, panel.first=grid(), xlab='', 
       pch=15, ylab='', col=data_summary$trt_color)
  mtext(text = expression('t'[1/2] ~ ' (hours)'),side = 1,line=3, cex=1.3)
  axis(1, at = my_vals,labels = my_vals)
  
  for(i in 1:nrow(data_summary)){
    lines(c(data_summary$t_12_low[i],
            data_summary$t_12_up[i]),
          rep(i,2),
          col=adjustcolor(data_summary$trt_color[i],alpha.f = .5))
  }
  
  for(kk in which(!duplicated(data_summary$Trt))){
    ind = data_summary$Trt==data_summary$Trt[kk]
    writeLines(sprintf('In %s the median clearance half life was %s (IQR %s to %s)',
                       data_summary$Trt[kk],
                       round(median(data_summary$t_12_med[ind]),1),
                       round(quantile(data_summary$t_12_med[ind],probs=0.25),1),
                       round(quantile(data_summary$t_12_med[ind],probs=0.75),1)))
    abline(v = median(data_summary$t_12_med[ind]), col=data_summary$trt_color[kk],lty=2,lwd=2)
  }
  
  legend('bottomright', legend = unique(data_summary$Trt),
         col = unique(data_summary$trt_color),
         pch=15,lwd=2,inset=0.03)
  return(data_summary)
}

get_itt_population = function(){
  rand.TH58 <- read.csv("~/Dropbox/PLATCOV/rand-TH58.csv")[0:9, ]
  rand.TH57 <- read.csv("~/Dropbox/PLATCOV/rand-TH57.csv")[0:10, ]
  
  data.TH1 <- read.csv("~/Dropbox/PLATCOV/data-TH1.csv")
  data.TH1$Date = as.POSIXct(data.TH1$Date,format='%a %b %d %H:%M:%S %Y')
  
  data.BR3 <- read.csv("~/Dropbox/PLATCOV/data-BR3.csv")
  data.BR3$Date = as.POSIXct(data.BR3$Date,format='%a %b %d %H:%M:%S %Y')
  
  data.LA08 <- read.csv("~/Dropbox/PLATCOV/data-LA08.csv")
  data.LA08$Date = as.POSIXct(data.LA08$Date,format='%a %b %d %H:%M:%S %Y')
  
  data.TH1$ID = paste('PLT-TH1-',data.TH1$randomizationID,sep='')
  rand.TH58$ID = paste('PLT-TH58-',rand.TH58$RandomisationID,sep='')
  rand.TH57$ID = paste('PLT-TH57-',rand.TH57$RandomisationID,sep='')
  data.BR3$ID = paste('PLT-BR3-',data.BR3$randomizationID,sep='')
  data.LA08$ID = paste('PLT-LA08-',data.LA08$randomizationID,sep='')
  
  xx = rbind(data.TH1[, c('ID', 'Treatment')],
             rand.TH58[, c('ID', 'Treatment')],
             rand.TH57[, c('ID', 'Treatment')],
             data.BR3[, c('ID', 'Treatment')])
  
  library(stringr)
  for(i in 1:nrow(xx)){
    id = unlist(strsplit(xx$ID[i],split = '-'))
    id[3] = str_pad(id[3], 3, pad = "0")
    id = paste(id, collapse = '-')
    xx$ID[i]=id
  }
  
  
  return(xx)
}

get_trt_colors = function(plot_cols=F){
  trt_cols = array(dim = 11)
  names(trt_cols) = 
    c("Ivermectin",
      "Regeneron",
      'No study drug',
      "Remdesivir",
      "Favipiravir",
      "Nitazoxanide",           
      "Fluoxetine",
      "Molnupiravir",
      "Nirmatrelvir + Ritonavir",
      "Evusheld",
      'Ensitrelvir')
  trt_cols['No study drug'] = viridis::viridis(n = 10)[8]
  trt_cols['Fluoxetine'] = viridis::viridis(n = 10)[5]
  trt_cols['Nitazoxanide'] = viridis::magma(n = 10)[8]
  trt_cols['Evusheld'] = viridis::magma(n = 10)[1]
  trt_cols['Favipiravir'] = viridis::plasma(n = 100)[92]
  trt_cols['Ivermectin'] = viridis::plasma(n = 10)[4]
  trt_cols['Nirmatrelvir + Ritonavir'] = viridis::plasma(n = 10)[1]
  trt_cols['Regeneron'] = viridis::inferno(n = 10)[5]
  trt_cols['Molnupiravir'] = viridis::inferno(n = 10)[7]
  trt_cols['Remdesivir'] = RColorBrewer::brewer.pal('Dark2',n=8)[8]
  trt_cols['Ensitrelvir'] = RColorBrewer::brewer.pal('Set1',n=8)[1]
  
  if(plot_cols){
    my_labels = gsub(pattern = ' + Ritonavir',replacement = '',fixed = T,x = names(trt_cols))
    plot(1:length(trt_cols), col=trt_cols, pch=16, cex=5, xlim = c(1,15),
         xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
    text(x = 1:length(trt_cols)+2.5, y= 1:length(trt_cols), labels = my_labels)
  }
  return(trt_cols)
}


assess_rebound = function(patient_dat,
                          lower_bound=2,  # lower level such that VL is defined as non-detectable
                          upper_bound=4,  # upper level such that VL is defined as "high"
                          t_window=1.5      # time window during which it has to be undetectable
){
  xx = patient_dat %>% arrange(Time) %>% distinct(Timepoint_ID, .keep_all = T)
  rebound = virus_cleared = F
  if(nrow(xx)>3){
    for(i in 2:nrow(xx)){
      ind = which(xx$Time <= xx$Time[i] & (xx$Time >= (xx$Time[i]-t_window)))
      if(all(xx$daily_VL[ind] <= lower_bound)){
        virus_cleared=T
      }
      if(virus_cleared & xx$daily_VL[i] >= upper_bound){
        rebound = T
        writeLines(sprintf('patient %s treated with %s had a rebound identified on day %s', 
                           xx$ID[1], xx$Trt[1],xx$Timepoint_ID[i]))
      }
      # print(rebound)
    }
  }
  return(rebound)
}

find_rebounds = function(platcov_dat, 
                         lower_bound=2,  # lower level such that VL is defined as non-detectable
                         upper_bound=4,  # upper level such that VL is defined as "high"
                         t_window=2      # time window during which it has to be undetectable
){
  platcov_dat$rebound=NA
  for(id in unique(platcov_dat$ID)){
    ind=platcov_dat$ID==id
    platcov_dat$rebound[ind] = assess_rebound(platcov_dat[ind,],
                                              lower_bound = lower_bound,
                                              upper_bound = upper_bound,
                                              t_window = t_window)
  }
  return(platcov_dat)
}


checkStrict(make_stan_inputs)
checkStrict(plot_serial_data)
checkStrict(plot_effect_estimates)
checkStrict(plot_data_model_fits)
checkStrict(plot_coef_effects)
checkStrict(calculate_fever_clearance)
