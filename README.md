# ZIQ-SIR
We provide a tool for doing hypothesis testing for zero-inflated microbiome data given covariates of interest.
## Description
This R document provides a novel zero-inflated semiparametric single-index quantile regression algorithm to conduct hypothesis testing for the relationship between zero-inflated responses and covariates of interest.
## Usage
You first need to load "ZIQSI-Rank.R" and install R packages: "MASS", "aod", "lme4", "PearsonDS". We provide two processes for hypothesis testing, one for small sample size and one for large sample size. We describe the data generation and how to conduct a pseudo-simulation extrapolation method for the given datasets with small and large sample sizes in "main.R". After running the "main.R", you can obtain the p-values for hypothesis testing.

## License
This R file is free and open source software, licensed under GPL (>=2)

## Author
Zirui Wang (wzr23@mails.tsinghua.edu.cn) and Tianying Wang (Tianying.Wang@colostate.edu)
