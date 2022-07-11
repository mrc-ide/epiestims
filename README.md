# A generic method and software to estimate the transmission advantage of pathogen variants in real-time: SARS-CoV-2 as a case-study

This repository contains code to reproduce the figures and analysis for the paper found here: https://doi.org/10.1101/2021.11.26.21266899.

## Reproducing the analysis

The code to perform the analysis for this work is set up as an [`orderly`](https://github.com/vimc/orderly) project.  The directories are:

* `src`: create new reports here
* `archive`: versioned results of running your report
* `data`: copies of data used in the reports

Each step of the analysis is an orderly task and corresponds to a directory within the src directory. To reproduce the analysis, run the tasks included in the script [`workflow.R`](https://github.com/mrc-ide/epiestims/blob/main/workflow.R) in order. 

The incidence data from England and France are available in the
directory `clean_french_england_data`.
A complete description of the French data set is
available
[here](https://www.data.gouv.fr/fr/datasets/donnees-de-laboratoires-pour-le-depistage-indicateurs-sur-les-variants/).


For the data for England, we used the daily number of positive tests
from England's community SARS-CoV-2 testing system stratified by NHS
region. The relevant columns are:

- specimen_date: Date of the specimen
- nhser_name: NHS England region
- s_positive_adj1: Number of cases with S-gene target positive.
- s_negative_adj1: Number of cases with S-gene target failure.
- s_na_adj1: Number of cases where the S-gene target failure was not
  available.

Further details are available in the Supplementary Material
accompanying the manuscript.
