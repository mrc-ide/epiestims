############################################################
############################################################


wrapper <- function(age_group,regions, plot_incidence = TRUE,
                    variants , t_window , 
                    SI , mean_prior ,
                    std_prior , n_sample_R, plot_Rt = TRUE ){
  
  # the function below will divide incidence by 7 as it weekly rolling count
  I <- format_I(d = d, age_group = age_group, variants = variants)
  
  if(plot_incidence==TRUE){
    layout(matrix(1:4,2,2))
    
    for(i in 1:length(regions)){
      plot(I[[1]][,1],I[[1]][,i+1],
           type = 'l',col='black',
           xlab = '',ylab='I',
           main = names(I[[1]])[i+1],
           ylim =c(0, max(c(I[[1]][,i+1],I[[2]][,i+1],I[[3]][,i+1]))) )
      lines(I[[2]][,1],I[[2]][,i+1],
            type = 'l',col='blue3')
      lines(I[[3]][,1],I[[3]][,i+1],
            type = 'l',col='red3')
      
      legend('topleft',legend = variants,
             lwd=2,col=c('black','blue3','red3'),
             bty='n')
    }
  }
  
  EpiEstim_Rt <- get_epiestim_rt_multi_Locations(I = I, 
                                                 variants = variants, 
                                                 regions = regions,
                                                 t_window = t_window, 
                                                 SI = SI_assumed,
                                                 mean_prior = mean_prior,
                                                 std_prior = std_prior,
                                                 n_sample_R = n_sample_R)
  
  if (plot_Rt == TRUE){
    cols <- c('black','blue3','red3')
    cols2 <- c(rgb(0,0,0,.1),rgb(0,0,1,.1),rgb(1,0,0,.1))
    
    layout(matrix(1:4,2,2))
    for(i in 1:length(regions)){
      for(j in 1:length(variants)){
        
        x <- I[[1]]$date[EpiEstim_Rt[[j]][[regions[i]]]$R$t_end]
        y <- cbind(EpiEstim_Rt[[variants[j]]][[regions[i]]]$R$`Median(R)`,
                   EpiEstim_Rt[[variants[j]]][[regions[i]]]$R$`Quantile.0.025(R)`,
                   EpiEstim_Rt[[variants[j]]][[regions[i]]]$R$`Quantile.0.95(R)`)
        if (j==1){
          plot(x,y[,1],
               xlab = '',ylab = 'Rt',
               type = 'l',col=cols[j],
               main = regions[i],
               ylim = c(0, 8) )
        }else{
          lines(x,y[,1],type = 'l',col=cols[j])
        }
        polygon(c(x,rev(x)), c(y[,2],rev(y[,3])),
                col = cols2[j],
                border = NA)
      }
      
    }
  }
  
  return(list(I = I, EpiEstim_Rt = EpiEstim_Rt))
}
