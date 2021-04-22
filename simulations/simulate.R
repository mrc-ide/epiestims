##################################################################################
### Generating time varying reproduction numbers for 2 locations and 2 strains ###
##################################################################################

source("./projections_code.R")

type <- "SEASONAL"
#type <- "CHANGEPOINT"

calc_seasonality <- function(date, seasonality_date_peak, seasonality, 
                             period = 365) {
  delta <- ((date - seasonality_date_peak) %% period) / period
  beta_mult <- 1 + cos(2 * pi * delta) * seasonality
  beta_mult
}

ndays <- 100
date <- seq(1, ndays, 1)
# set up changes over time in location 1 due mostly to seasonality
seasonality_date_peak_loc1 <- 90
seasonality_loc1 <- 0.2 # do + / - 10% of the average at the peak
seas_loc1 <- calc_seasonality(date, seasonality_date_peak_loc1, seasonality_loc1,
                              period = 20)
# set up changes over time in location 2 due mostly to rapid cahnges in policy
seasonality_date_peak_loc2 <- 30
seasonality_loc2 <- 0.3 # do + / - 30% of the average at the peak
seas_loc2 <- calc_seasonality(date, seasonality_date_peak_loc2, seasonality_loc2,
                              period = 100)

# set up changes over time in location 3 
seasonality_date_peak_loc3 <- 15
seasonality_loc3 <- 0.5 # do + / - 50% of the average at the peak
seas_loc3 <- calc_seasonality(date, seasonality_date_peak_loc3, seasonality_loc3,
                              period = 150)

if(type == "SEASONAL") {
  # generate a time varying reproduction number for location 1
  R1_loc1 <- 2.5 * seas_loc1
  R1_loc2 <- 2.5 * seas_loc2
  R1_loc3 <- 2.5 * seas_loc3
  
  # bring min to 1
  R1_loc1 <- R1_loc1 - min(R1_loc1) +1
  R1_loc2 <- 2.5 * seas_loc2 - min(R1_loc2) + 1
  R1_loc3 <- 2.5 * seas_loc3 - min(R1_loc3) + 1
} else if(type == "CHANGEPOINT"){
  # generate a time varying reproduction number for location 1
  R1_loc1 <- c(rep(2.5, round(ndays / 2)), rep(1.5, round(ndays / 2)))
  R1_loc2 <- c(rep(2.5, round(ndays / 4)), rep(1.5, round(3 * ndays / 4)))
  R1_loc3 <- c(rep(2.5, round(ndays / 10)), rep(1.5, round(9 * ndays / 10)))
}

transmission_advantage <- 1.5

R2_loc1 <- transmission_advantage * R1_loc1
R2_loc2 <- transmission_advantage * R1_loc2
R2_loc3 <- transmission_advantage * R1_loc3

n_loc <- 3
n_v <- 2

R <- array(NA, dim= c(ndays, n_loc, n_v))
R[, 1, 1] <- R1_loc1
R[, 2, 1] <- R1_loc2
R[, 3, 1] <- R1_loc3
R[, 1, 2] <- R2_loc1
R[, 2, 2] <- R2_loc2
R[, 3, 2] <- R2_loc3


par(mfrow = c(2, 2))

plot(date, R1_loc1, type = "l", ylim = c(0, 8), ylab = "R")
lines(date, R1_loc2, col = "blue")
lines(date, R1_loc3, col = "red")
lines(date, R2_loc1, lty = 2)
lines(date, R2_loc2, col = "blue", lty = 2)
lines(date, R2_loc3, col = "red", lty = 2)
legend("topright", c("Strain 1, location 1", "Strain 1, location 2", "Strain 1, location 3",
                       "Strain 2, location 1", "Strain 2, location 2", "Strain 2, location 3"),
       col = c("black", "blue", "red", "black", "blue", "red"), 
       lty = c(1, 1, 1, 2, 2, 2), cex = .5)

##################################################################################
### Simulate epidemics with those reproduction numbers ###
##################################################################################

set.seed(1)

## covid serial interval vaguely resembling this
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7448781/
si <- EpiEstim::discr_si(0:30, mu = 5.40, sigma = 2) # I made up the sd
si_no_zero <- si[-1]
plot(seq(0, length(si) - 1, 1), si, 
     type = "h", xlab = "Day", ylab = "PMF of the SI")

## start with 5 infected individual
initial_incidence <- incidence::incidence(rep(1, 5))

incid <- array(NA, dim= c(ndays, n_loc, n_v))

for (loc in 1:n_loc) {
  for (v in 1:n_v) {
    incid[, loc, v] <- rbind(initial_incidence$counts,
                             as.matrix(#projections::project(initial_incidence,
                               debug_project(initial_incidence,
                                             R = R[-1, loc, v], # R in the future so removing time of seeding
                                             si = si_no_zero,
                                             n_sim = 1,
                                             n_days = ndays - 1,
                                             time_change = seq_len(length(R[, loc, v]) - 2) - 1,
                                             instantaneous_R = TRUE)))
  }
}

plot(log(1 + incid[, 1, 1]), type= "l", xlab = "date", ylab = "log(1 + Incidence)")
lines(log(1 + incid[, 2, 1]), col = "blue")
lines(log(1 + incid[, 3, 1]), col = "red")
lines(log(1 + incid[, 1, 2]), lty = 2)
lines(log(1 + incid[, 2, 2]), col = "blue", lty = 2)
lines(log(1 + incid[, 3, 2]), col = "red", lty = 2)
legend("bottomright", c("Strain 1, location 1", "Strain 1, location 2", "Strain 1, location 3",
                        "Strain 2, location 1", "Strain 2, location 2", "Strain 2, location 3"),
       col = c("black", "blue", "red", "black", "blue", "red"), 
       lty = c(1, 1, 2, 2), cex = 0.7)


##################################################################################
### Generate summary objects ###
##################################################################################

saveRDS(incid, "incidence.rds")
saveRDS(si, "si.rds")
saveRDS(R, 'target_Rt.rds')
