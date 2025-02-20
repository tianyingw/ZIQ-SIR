# Powerful and robust non-parametric association testing for zero-inflated data via a zero-inflated quantile approach (ZIQ-SIR)
## Introduction
The **ZIQ-SIR** R package is designed to provide advanced statistical tools for association testing for zero-inflated data. This package includes one main function:

- `ZIQSIR`: A function for testing association for zero-inflated data and covariates of interest.

This document provides an overview of the package, how to install it, and examples of how to use these functions.

## Usage
## Installation
To install the **ZIQ-SIR** package from GitHub, you can use the following code:

```{r}
# If devtools is not already installed, install it first
# install.packages("devtools")
# Install the ZIQSI package from GitHub
devtools::install_github("https://github.com/tianyingw/ZIQ-SIR")
```

## Functions in ZIQ-SIR Package
### ZIQSIR
The `ZIQSIR` function provides a novel zero-inflated semiparametric single-index quantile regression algorithm to obtain the p-value for testing associations for zero-inflated response. This function provides two options for hypothesis testing: `large` for large sample cases; `small` for small sample cases.

## Examples and Usage
```{r}
library(ZIQ-SIR)
 #### demo 1:
 # small sample size
 # using Pearson Type III method
 # alternative distribution
 # sample size
 set.seed(10001)
 n = 500
 # the probability of Y>0 given covariate x
 p = function(x1,x2,x3,x4,x5,gam0=-0.4,gam1=-0.480,gam2=-0.022,gam3=0.021,gam4 = 0.015,gam5 = -0.009){
   lc = gam0 + gam1*x1 + gam2*x2 + gam3*x3 + gam4*x4 + gam5*x5
   exp(lc)/(1+exp(lc))
 }
 # beta_tau
 bet1 = function(x){(0.3*sqrt(x)-x)*2}
 bet2 = function(x){x^2*2.2}
 bet3 = function(x){
   (x^2-0.5*x+0.6)*2/3
 }
 bet5 = function(x){-(0.3*x^2-x)*2}
 bet4 = function(x){-sin(x*2*pi)*0.1}
 bet0 = function(x){-147.7*x-50*x^2-20}
 bet = function(x,u){x^4*u*10^(-5)/6+x^2*u*0.2/3}

 # G_tau function
 func <- function(x, tau)
 {
   return(bet(x %*% rbind(bet1(tau),bet2(tau),bet3(tau),bet4(tau),bet5(tau))+bet0(tau),tau))
 }

 # given samples
 x1 = rbinom(n,1,0.5) # Medicament use
 x2 = rnorm(n,28,2) # bmi
 x3 = rnorm(n,92.5,13) # waist
 x4 = rnorm(n,80,12) # diastolic_bp
 x5 = rnorm(n,124,18.5) # systolic_bp
 x0 = rep(1,n)
 X = cbind(x0,x1,x2,x3,x4,x5)
 u = runif(n)
 b = rbinom(n,1,p(x1,x2,x3,x4,x5))
 w = bet(bet1(u)*x1+bet2(u)*x2+bet3(u)*x3+bet4(u)*x4+bet5(u)*x5+bet0(u),u)
 y = b*w

 Combination(y,X,m = 3,test_num = 4,method = "Pearson")

 ### demo 2
 # simulation results under large sample size
 # using Chi-square method
 # under null hypothesis
 # sample size
 n = 2000

 # the probability of Y>0 given covariate x
 p = function(x1,x2,x3,x4,x5,x6,
              gam0=2.32,gam1=-0.06,gam2=-0.03*0,gam3=-0.010*0,gam4 = -0.005,gam5 = 0.0005,gam6 = -0.030){
   lc = gam0 + gam1*x1 + gam2*x2 + gam3*x3 + gam4*x4 + gam5*x5 + gam6*x6
   exp(lc)/(1+exp(lc))
 }

 # beta_tau
 bet1 = function(x){5*x^2+1}
 bet2 = function(x){(0.1*sin(x*2*pi)+0.05)*0}
 bet3 = function(x){
   ( 0.05*(x-0.5)^2+0.04)*0
 }
 bet4 = function(x){(0.1*sqrt(x)+2.9*x)*0.05}
 bet5 = function(x){(0.4*(x-1)^2)*0.1}
 bet6 = function(x){((x-0.6)*(x-1.1))*1.1}
 bet0 = function(x){62.9*x^2+33.4*x}
 bet = function(x,u){x^4*u*10^(-4)*0.4+x^3*u*10^{-3}*0.1}

 # G_tau function
 func <- function(x, tau)
 {
   return(bet(x %*% rbind(bet1(tau),bet2(tau),bet3(tau),bet4(tau),bet5(tau),bet6(u))+bet0(tau),tau))
 }

 # given samples
 x1 = rbinom(n,1,0.5) # Gender
 x2 = rnorm(n,28,2) # bmi
 x3 = 2*x2+rnorm(n,36.5,9) # waist
 x4 = rnorm(n,80,12) # diastolic_bp
 x5 = x4*1.3+rnorm(n,20,7.75) # systolic_bp
 x6 = sample(1:4, n, replace = TRUE, prob = c(0.25,0.025,0.07,0.625))
 x0 = rep(1,n)
 u = runif(n)
 b = rbinom(n,1,p(x1,x2,x3,x4,x5,x6))
 w = bet(bet1(u)*x1+bet2(u)*x2+bet3(u)*x3+bet4(u)*x4+bet5(u)*x5+bet6(u)*x6+bet0(u),u)
 y = b*w
 X = cbind(1,x1,x2,x3,x4,x5,x6)

 Combination(y,X,m = 3,test_num = c(2,3),method = "Chi")

```

## License
This R file is free and open source software, licensed under GPL (>=2)


