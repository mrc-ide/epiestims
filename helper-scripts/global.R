library(EpiEstim)
library(projections)
library(incidence)
library(purrr)
library(dplyr)
library(ggplot2)
library(here)
source("simulations/simulation_functions.R")

if (! dir.exists("figures")) dir.create("figures")
seed <- 42
set.seed(seed)
