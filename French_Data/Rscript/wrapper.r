############################################################
############################################################


wrapper <- function(dat, regions, plot_incidence = TRUE,
                    variants , t_window , 
                    SI , mean_prior ,
                    std_prior , n_sample_R ){
  
  # the function below will divide incidence by 7 as it weekly rolling count
  I <- format_I(dat = dat, variants = variants)
  col=c('black','blue3','red3')
  
  if(plot_incidence==TRUE){
    layout(matrix(1:4,2,2))
    
    for(i in 1:length(regions)){
      ylim = c(0,0)
      for(j in 1:length(I)){
        ylim[2] <- max(c(ylim[2],I[[j]][,i+1]))
      }
      
      plot(I[[1]][,1],I[[1]][,i+1],
           type = 'l',col='black',
           xlab = '',ylab='I',
           main = names(I[[1]])[i+1],
           ylim =ylim )
      
      for (j in 2:length(variants)){
        lines(I[[j]][,1],I[[j]][,i+1],
              type = 'l',col=col[j])
      }
      
      legend('topleft',legend = variants,
             lwd=2,col=col[1:length(variants)],
             bty='n')
    }
  }
  
  EpiEstim_Rt <- get_epiestim_rt_multi_Locations(I = I, 
                                                 variants = variants, 
                                                 regions = regions,
                                                 t_window = t_window, 
                                                 SI = SI,
                                                 mean_prior = mean_prior,
                                                 std_prior = std_prior,
                                                 n_sample_R = n_sample_R)
  
  
  return(list(I = I, EpiEstim_Rt = EpiEstim_Rt))
}

##############################################################################################
##############################################################################################

# wrapper for joint estimate

wrapper_joint_Rt <- function(incid, t_start, t_end, n_intervals, 
                             plot_Rt = TRUE, plot_epsilon_trace = TRUE, 
                             plot_epsilon = TRUE,initial_res ){
  R <- list()
  
  t_int <- round(seq(t_start, t_end, length.out = n_intervals+1))
  
  R[[1]] <- EpiEstim::estimate_joint(incid = incid, si_distr = si_distr, 
                                     priors = priors,mcmc_control = mcmc_control,
                                     t_min = t_start,
                                     t_max = t_end)
  
  for(k in 1:n_intervals){
    R[[k+1]] <- EpiEstim::estimate_joint(incid = incid, si_distr = si_distr, 
                                         priors = priors,mcmc_control = mcmc_control,
                                         t_min = as.integer(t_int[k]),
                                         t_max = as.integer(t_int[k+1] ) )
  }
  
  
  layout(matrix(1:4,2,2))
  cols <- RColorBrewer::brewer.pal(n_intervals+1,name = 'Dark2')
  cols2 <- yarrr::transparent(cols, trans.val = .9)
  Sum_R <- list()
  for(k in 1:(n_intervals+1)){
    Sum_R[[k]] <- list(median = apply(R[[k]]$R,c(1,2),median),
                       low = apply(R[[k]]$R,c(1,2),quantile,0.025,na.rm=TRUE),
                       high = apply(R[[k]]$R,c(1,2),quantile,0.975,na.rm=TRUE))
  } 
  
  for(i in 1:length(regions)){
    for(k in 1:(n_intervals+1)){
      x <- initial_res$I$alpha$date #1:nrow(R[[k]]$R)
      y <- cbind(Sum_R[[k]]$median[,i],Sum_R[[k]]$low[,i],Sum_R[[k]]$high[,i])
      f <- which(!is.na(Sum_R[[k]]$median[,i]))
      x<-x[f]
      y<-y[f,]
      
      if(k==1){
        plot(x,y[,1],
             xlab = '',ylab = 'Rt',
             type = 'l',col=cols[k],
             main = regions[i],
             ylim = c(0, 2) )
      }else{
        lines(x,y[,1],col=cols[k])
      }
      
      abline(h = 1,lty=2,col='red3')
      
      polygon(c(x,rev(x)), c(y[,2],rev(y[,3])),
              col = cols2[k],
              border = NA)
    }
  }
  
  
  Sum_epsilon <- list()
  for(i in 1:(n_intervals+1)){
    Sum_epsilon[[i]]<- apply(R[[i]]$epsilon,1,quantile,c(.5,.025,.975))
  }
  Sum_epsilon <- matrix(unlist(Sum_epsilon),nrow = 3, ncol = (n_intervals+1)*(dim(incid)[3]-1))
  
  if(plot_epsilon_trace){
    layout(matrix(1:4,2,2))
    for(i in 1:(dim(incid)[3]-1)){
      for(k in 1:(n_intervals+1)){
        if (k==1){
          plot(R[[k]]$epsilon[i,],col=cols[k],ylim=c(0,2))
        }else{
          lines(R[[k]]$epsilon[i,],col=cols[k],type = 'p')
        }
      }
    }
  }
  
  if(plot_epsilon){
    Hmisc::errbar(x = 1:((n_intervals+1)*(dim(incid)[3]-1)), y = Sum_epsilon[1,],
                  ylab = 'epsilon',
                  yplus = Sum_epsilon[2,], yminus = Sum_epsilon[3,],
                  xlim = c(0,((n_intervals+1)*(dim(incid)[3]-1))+1),ylim = c(0,2))
    abline(h = 1,lty=2,col='red3')
  }
  
  
  return(list(R = R, Sum_epsilon = Sum_epsilon))
}