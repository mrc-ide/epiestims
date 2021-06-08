###########################################################################
###########################################################################
# get list of matrix, 1 for each variant with 1 column for dates and 'n' columns for each location

format_I <- function(d,age_group,variants){
  f <- which(d$cl_age90 == age_group)
  d <- d[f,]
  # names(d)
  
  I <- list()
  for(i in 1:length(variants)){
    temp <- d[,c(which(names(d) %in% 'region'),
                 which(names(d) %in% variants[i]),
                 which(names(d) %in% 'week_end'))] # 
    names(temp) <- c('region','I','date')
    temp$I <- temp$I/7
    I[[variants[i]]] <- tidyr::spread(data = temp, region,  I)
  }
  return(I)
}


###########################################################################
###########################################################################

# get Rts for each locations and variants
# extract R matrix from epiestim and samples of posterior

get_epiestim_rt_multi_Locations <- function(I, variants, regions, t_window, 
                                            SI, mean_prior, std_prior,
                                            n_sample_R){
  
  t_start <- seq(2, nrow(I[[1]])-(t_window-1)) 
  t_end <- t_start + (t_window-1)
  
  EpiEstim_Rt <- list()
  
  for (k in 1:length(variants)){
    
    EpiEstim_Rt[[variants[k]]] <- list()
    
    for(i in 1:length(regions)){
      EpiEstim_Rt[[variants[k]]][[regions[i]]] <- list()
      
      temp <- EpiEstim::estimate_R(incid = I[[k]][,i+1], 
                                   method = "non_parametric_si", 
                                   config = EpiEstim::make_config(list(si_distr = SI,
                                                                       t_start = t_start ,
                                                                       t_end = t_end,
                                                                       mean_prior = mean_prior,
                                                                       std_prior = std_prior)))
      
      EpiEstim_Rt[[variants[k]]][[regions[i]]]$R <- temp$R
      
      temp2 <- matrix(NA, nrow = nrow(temp$R),
                      ncol = n_sample_R)
      for(j in 1: nrow(temp$R)){
        param <-epitrix::gamma_mucv2shapescale(mu = temp$R$`Mean(R)`[j],
                                               cv = temp$R$`Std(R)`[j]/temp$R$`Mean(R)`[j])
        temp2[j,] <- rgamma(n = n_sample_R, shape = param$shape, scale = param$scale)
      }
      
      EpiEstim_Rt[[variants[k]]][[regions[i]]]$R_samples <- temp2
    }
  }
  
  return(EpiEstim_Rt)
}


###########################################################################
###########################################################################

# plot 2 d hist of variants Rts

plot2d_hist <- function(x,k){
  sub('-vs-.*', '', names(x)[k])
  sub('.*-vs-', '', names(x)[k])
  ggplot(x[[k]], aes(x,y)) +  
    xlab(paste0('Rt (',sub('-vs-.*', '', names(x)[k]),')')) + 
    ylab(paste0('Rt (',sub('.*-vs-', '', names(x)[k]),')')) +
    stat_bin2d(bins=25) + 
    xlim(0,max(x[[k]])) + ylim(0,max(x[[k]])) +
    scale_fill_gradientn(colours=r) + 
    geom_abline(intercept = 0,slope = 1,col='red3') 
  
}


###########################################################################
###########################################################################


# select right Rt based on threshold for 95CrI range


select_Rt_get_median_samples <- function(th, EpiEstim_Rt, 
                                         regions, variants,
                                         SI ,
                                         trim = 0){
  
  trim_time <- which(cumsum(SI)>=trim)[1]
  select_Rt <- list()
  temp <- data.frame(matrix(0,nrow = nrow(EpiEstim_Rt[[1]][[regions[1]]]$R),
                            ncol = length(regions)))
  
  # check combination of 2 variants and see whether 95CrI range bekow threshold at each time steps
  for(i in 1:(length(variants)-1)){
    for(j in (i+1):length(variants)){
      select_Rt[[paste0(variants[i],'-vs-',variants[j])]] <- temp
      for(k in 1:length(regions)){
        range_95_i <- EpiEstim_Rt[[i]][[regions[k]]]$R$`Quantile.0.975(R)` - EpiEstim_Rt[[i]][[regions[k]]]$R$`Quantile.0.025(R)` 
        range_95_j <- EpiEstim_Rt[[j]][[regions[k]]]$R$`Quantile.0.975(R)` - EpiEstim_Rt[[j]][[regions[k]]]$R$`Quantile.0.025(R)` 
        
        f <- which((range_95_i < th) & (range_95_j < th))
        f <- f[which(f>=trim_time)]
        select_Rt[[paste0(variants[i],'-vs-',variants[j])]][f,k] <- 1
      }
      
    }
  }
  
  # number of time step for which conditions is met
  summary_select <- data.frame(matrix(NA,nrow = length(regions), ncol = length(select_Rt)))
  summary_select[,1] <- regions
  c<-1
  for(i in 1:(length(variants)-1)){
    for(j in (i+1):length(variants)){
      c<-c+1
      for(k in 1:length(regions)){
        
        summary_select[k,c] <- sum(select_Rt[[paste0(variants[i],'-vs-',variants[j])]][,k])
      }
      
    }
  }
  names(summary_select) <- c('region',names(select_Rt))
  
  # median of Rts pairs
  median_Rts <- list()
  
  for(i in 1:(length(variants)-1)){
    for(j in (i+1):length(variants)){
      x <- c()
      y <- c()
      for(k in 1:length(regions)){
        f <- which(select_Rt[[paste0(variants[i],'-vs-',variants[j])]][,k]==1) 
        x <- c(x,EpiEstim_Rt[[variants[i]]][[regions[k]]]$R$`Median(R)`[f])
        y <- c(y,EpiEstim_Rt[[variants[j]]][[regions[k]]]$R$`Median(R)`[f])
        # if(length(x)!=length(y)) break
      }
      median_Rts[[paste0(variants[i],'-vs-',variants[j])]] <- data.frame(x = x, y = y)
      
    }
  }
  
  
  # samples 
  samples_Rts <- list()
  
  for(i in 1:(length(variants)-1)){
    for(j in (i+1):length(variants)){
      x <- c()
      y <- c()
      for(k in 1:length(regions)){
        f <- which(select_Rt[[paste0(variants[i],'-vs-',variants[j])]][,k]==1) 
        x <- c(x,c(EpiEstim_Rt[[variants[i]]][[regions[k]]]$R_samples[f,]) )
        y <- c(y,c(EpiEstim_Rt[[variants[j]]][[regions[k]]]$R_samples[f,]) )
        # if(length(x)!=length(y)) break
      }
      samples_Rts[[paste0(variants[i],'-vs-',variants[j])]] <- data.frame(x = x, y = y)
      
    }
  }
  
  return(list(select_Rt = select_Rt,
              summary_select = summary_select,
              median_Rts = median_Rts,
              samples_Rts = samples_Rts))
}




###########################################################################
###########################################################################
# \ plot 2d hist + violin

plot_hist_dist <- function(x, x_sum, keep = FALSE){
  

  p <- q <-list()
  epsi <- list()
  for(k in 1:3){
    # print(k)
    if(sum(x_sum[k+1])>0){
      
      p[[k]] <- plot2d_hist(x = x, k)
      
      epsilon <- data.frame(v = names(x)[k],
                            epsilon = x[[k]]$y/x[[k]]$x)
      
      q[[k]] <- ggplot(epsilon, aes(factor(v),epsilon)) +  
        ylim(0,2.51) +
        geom_violin() + 
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
        xlab('') +
        ylab('epsilon') + 
        ggtitle(paste0(sub('.*-vs-', '', names(x)[k]),'-vs-',sub('-vs-.*', '', names(x)[k]))) +
        theme(plot.title = element_text( size=10),
              axis.ticks.x = element_blank()) +
        scale_x_discrete(labels=c('')) + 
        geom_abline(intercept = 1,slope = 0,col='red3') 
      
      
      epsi[[paste0(sub('.*-vs-', '', names(x)[k]),'-vs-',sub('-vs-.*', '', names(x)[k]))]] <- epsilon$epsilon
    }else{
      
      epsilon <- data.frame(v = names(x)[k],
                            epsilon = NA)
      
      p[[k]] <- ggplot(epsilon, aes(factor(v),epsilon))
      q[[k]] <- ggplot(epsilon, aes(factor(v),epsilon))
      epsi[[paste0(sub('.*-vs-', '', names(x)[k]),'-vs-',sub('-vs-.*', '', names(x)[k]))]] <- NA
    }
  }
  
  
  gridExtra::grid.arrange(p[[1]],q[[1]],p[[2]],q[[2]],p[[3]],q[[3]],
                          nrow = 3,
                          layout_matrix = matrix(c(1,1,1,2,3,3,3,4,5,5,5,6),
                                                 nrow = 3, ncol = 4,
                                                 byrow = TRUE))
  if(keep){
    return(epsi)
  }
}
