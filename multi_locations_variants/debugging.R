library(rstan)
library(dplyr)
library(purrr)
library(ggplot2)

set.seed(1)
data_incidence <- readRDS('../simulations/incidence.rds')
data_incidence <- data_incidence[4:99,]
SI_sim <- readRDS('../simulations/si.rds')
n_location <- 2L
# n_variant <- as.integer(2)
I <- data_incidence[, 1:2]

t_start <- seq(2, nrow(I)-6)
t_end <- t_start + 6

# preprocessed

O_I <- apply(I,2,EpiEstim::overall_infectivity, si_distr = SI_sim)

mI <- array(NA, dim = c( n_location,
                         length(t_start),
                         ncol = t_end[1] - t_start[1]+1) )
m_O_I <- mI
for(k in 1:(n_location)){
  for(i in 1:dim(mI)[2]){
    mI[k,i,] <- I[t_start[i]:t_end[i],k]
    m_O_I[k,i,] <- O_I[t_start[i]:t_end[i],k]
  }
}
useful <- (rep(1,dim(mI)[3]))


## prior R
mean_prior <- c(2)
std_prior <- c(1)

param_gamma <- epitrix::gamma_mucv2shapescale(mu = mean_prior,
                                              cv = std_prior/mean_prior)


## 1.  Turned off likelihood in both stan versions to see what prior
## is being sampled
##

standata <- list(nt = as.integer(dim(mI)[2]),
                 tw = as.integer(dim(mI)[3]),
                 n_location = n_location,
                 I = mI,
                 O_I = m_O_I,
                 U = useful,
                 prior_shape = param_gamma$shape,
                 prior_rate = 1/param_gamma$scale)
fit <- stan(
  file = 'ML_NoVariant_epi_estim_v0.stan', data = standata, seed = 42
)




individual_fits <- map(
  seq_len(n_location), function(i) {
    standata <- list(
      nR = dim(mI)[2],
      tw = dim(mI)[3],
      I = mI[i,,],
      O_I = m_O_I[i,,],
      U = useful,
      prior_alpha = param_gamma$shape,
      prior_beta = 1/param_gamma$scale
    )
    fit <- stan(file = 'epi_estim.stan', data = standata, seed = 42)
    fit
  }
)

individual_fit_samples <- map(individual_fits, rstan::extract)

v2_samples <- rstan::extract(fit)

## dim(v2_samples[["Rt"]])
## [1] 4000   89    2
## dim(individual_fit_samples[[1]][["Rt"]])
## [1] 4000   89

p <- ggplot() +
  geom_density(
    aes(v2_samples[["Rt"]][, 2, 1], fill = "blue"), alpha = 0.3,
    col = NA
  ) +
  geom_density(
    aes(
      individual_fit_samples[[1]][["Rt"]][, 2], fill = "red"
    ), alpha = 0.3, col = NA
  ) +
  scale_fill_identity(
    breaks = c("blue", "red"),
    labels = c("Loop in stan", "Loop in R"),
    guide = "legend"
  ) +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())


cowplot::save_plot("figs/loop_in_stan_vs_in_R.png", p)
