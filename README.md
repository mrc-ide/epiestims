# A generic method and software to estimate the transmission advantage of pathogen variants in real-time: SARS-CoV-2 as a case-study

This repository contains code to reproduce the figures and analysis for the paper found here: https://doi.org/10.1101/2021.11.26.21266899.

## Reproducing the analysis

The code to perform the analysis for this work is set up as an [`orderly`](https://github.com/vimc/orderly) project.  The directories are:

* `src`: create new reports here
* `archive`: versioned results of running your report
* `data`: copies of data used in the reports

Each step of the analysis is an orderly task and corresponds to a directory within the src directory. To reproduce the analysis, run the tasks included in the script [`workflow.R`](https://github.com/mrc-ide/epiestims/blob/main/workflow.R) in order. 
