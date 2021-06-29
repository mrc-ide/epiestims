library(dplyr)
library(EpiEstim)
library(ggplot2)
library(glue)
library(here)
library(incidence)
library(projections)
library(purrr)

source("simulations/simulation_functions.R")
source("simulations/cluster_functions.R")
if (! dir.exists("figures")) dir.create("figures")
seed <- 42
set.seed(seed)
