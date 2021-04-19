##################################################################################
### Generating time varying reproduction numbers for 2 locations and 2 strains ###
##################################################################################

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
seasonality_loc2 <- 0.3 # do + / - 10% of the average at the peak
seas_loc2 <- calc_seasonality(date, seasonality_date_peak_loc2, seasonality_loc2,
                              period = 100)

# generate a time varying reproduction number for location 1
R1_loc1 <- 2.5 * seas_loc1
R1_loc2 <- 2.5 * seas_loc2

# bring min to 1
R1_loc1 <- R1_loc1 - min(R1_loc1) +1
R1_loc2 <- 2.5 * seas_loc2 - min(R1_loc2) +1

transmission_advantage <- 1.5

R2_loc1 <- transmission_advantage * R1_loc1
R2_loc2 <- transmission_advantage * R1_loc2

par(mfrow = c(2, 2))

plot(date, R1_loc1, type = "l", ylim = c(0, 5), ylab = "R")
lines(date, R1_loc2, col = "blue")
lines(date, R2_loc1, lty = 2)
lines(date, R2_loc2, col = "blue", lty = 2)
legend("bottomleft", c("Strain 1, location 1", "Strain 1, location 2",
                       "Strain 2, location 1", "Strain 2, location 2"),
       col = c("black", "blue", "black", "blue"), lty = c(1, 1, 2, 2), cex = .7)

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

## location 1 strain 1
inc1_loc1 <- rbind(initial_incidence$counts,
                   as.matrix(projections::project(initial_incidence, 
                     R = R1_loc1, 
                     si = si_no_zero,
                     n_sim = 1,
                     n_days = ndays - 1,
                     time_change = seq_len(length(R1_loc1) - 1))))

## location 2 strain 1
inc1_loc2 <- rbind(initial_incidence$counts,
                   as.matrix(projections::project(initial_incidence, 
                                            R = R1_loc2, 
                                            si = si_no_zero,
                                            n_sim = 1,
                                            n_days = ndays - 1,
                                            time_change = seq_len(length(R1_loc2) - 1))))

## location 1 strain 2
inc2_loc1 <- rbind(initial_incidence$counts,
                   as.matrix(projections::project(initial_incidence, 
                                            R = R2_loc1, 
                                            si = si_no_zero,
                                            n_sim = 1,
                                            n_days = ndays - 1,
                                            time_change = seq_len(length(R2_loc1) - 1))))

## location 2 strain 2
inc2_loc2 <- rbind(initial_incidence$counts,
                   as.matrix(projections::project(initial_incidence, 
                                            R = R2_loc2, 
                                            si = si_no_zero,
                                            n_sim = 1,
                                            n_days = ndays - 1,
                                            time_change = seq_len(length(R2_loc2) - 1))))

plot(log(1 + inc1_loc1[, 1]), type= "l", xlab = "date", ylab = "log(1 + Incidence)")
lines(log(1 + inc1_loc2[, 1]), col = "blue")
lines(log(1 + inc2_loc1[, 1]), lty = 2)
lines(log(1 + inc2_loc2[, 1]), col = "blue", lty = 2)
legend("bottomright", c("Strain 1, location 1", "Strain 1, location 2",
                       "Strain 2, location 1", "Strain 2, location 2"),
       col = c("black", "blue", "black", "blue"), lty = c(1, 1, 2, 2), cex = 0.7)


##################################################################################
### Generate summary objects ###
##################################################################################

incidence <- as.data.frame(cbind(inc1_loc1, inc1_loc2, inc2_loc1, inc2_loc2))
names(incidence) <- c("strain1_loc1", "strain1_loc2", "strain2_loc1", "strain2_loc2")
saveRDS(incidence, "incidence.rds")
saveRDS(si, "si.rds")

# target Rt
target_Rt = list(R1_loc1 = R1_loc1,
                 R2_loc1 = R2_loc1,
                 R1_loc2 = R1_loc2,
                 R2_loc2 = R2_loc2)
saveRDS(target_Rt, 'target_Rt.rds')