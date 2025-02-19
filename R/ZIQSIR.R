#' Combine the marginal p-values
#'
#' @import MASS
#' @import lme4
#' @import PearsonDS
#' @import splines2
#' @import stats
#' @import quantreg
#' @param y n*1 vector, the observed outcome for n samples
#' @param X n*p matrix, the observed p covariates for n samples
#' @param taus k*1 vector, a grid of quantile levels; e.g., 0.5 for the median, 0.75 for the 3rd quartile; default is c(0.1, 0.25, 0.5, 0.75, 0.9)
#' @param m numeric variable, the order of B-spline function; default is 3
#' @param test_num a vector, representing the test corresponds to which covariate(s) in X.
#' @param method different method for calculating p-value: 'Chi' for large sample cases; 'Pearson' for small sample cases
#' @details
#' \itemize{
#'   \item Please choose 'Chi' or 'Pearson' for \code{method}, no other options.
#'   \item \code{taus} must be a subset or equal to the grid used to produce \code{input}.
#' }
#' @return quantiles of a m*k matrix, each row is the estimated quantiles for each new case
#' @export
#' @examples
#'
#' # example code
#' #### demo 1:
#' # small sample size
#' # using Pearson Type III method
#' # alternative distribution
#' # sample size
#' set.seed(10001)
#' n = 500
#' # the probability of Y>0 given covariate x
#' p = function(x1,x2,x3,x4,x5,gam0=-0.4,gam1=-0.480,gam2=-0.022,gam3=0.021,gam4 = 0.015,gam5 = -0.009){
#'   lc = gam0 + gam1*x1 + gam2*x2 + gam3*x3 + gam4*x4 + gam5*x5
#'   exp(lc)/(1+exp(lc))
#' }
#' # beta_tau
#' bet1 = function(x){(0.3*sqrt(x)-x)*2}
#' bet2 = function(x){x^2*2.2}
#' bet3 = function(x){
#'   (x^2-0.5*x+0.6)*2/3
#' }
#' bet5 = function(x){-(0.3*x^2-x)*2}
#' bet4 = function(x){-sin(x*2*pi)*0.1}
#' bet0 = function(x){-147.7*x-50*x^2-20}
#' bet = function(x,u){x^4*u*10^(-5)/6+x^2*u*0.2/3}

#' # G_tau function
#' func <- function(x, tau)
#' {
#'   return(bet(x %*% rbind(bet1(tau),bet2(tau),bet3(tau),bet4(tau),bet5(tau))+bet0(tau),tau))
#' }
#' # the given covariate x
#' X1 = c(0,1) # sex
#' X2 = qnorm(0.5,28,2) # bmi
#' X3 = qnorm(0.5,92.5,13) # waist
#' X4 = qnorm(0.5,80,12) # diastolic_bp
#' X5 = qnorm(0.5,124,18.5) # systolic_bp
#' X0 =cbind(c(rep(X1[1], 1), rep(X1[2], 1)),  rep(X2, 2), rep(X3, 2), rep(X4, 2), rep(X5,2))

#' # given samples
#' x1 = rbinom(n,1,0.5) # Medicament use
#' x2 = rnorm(n,28,2) # bmi
#' x3 = rnorm(n,92.5,13) # waist
#' x4 = rnorm(n,80,12) # diastolic_bp
#' x5 = rnorm(n,124,18.5) # systolic_bp
#' x0 = rep(1,n)
#' X = cbind(x0,x1,x2,x3,x4,x5)
#' u = runif(n)
#' b = rbinom(n,1,p(x1,x2,x3,x4,x5))
#' w = bet(bet1(u)*x1+bet2(u)*x2+bet3(u)*x3+bet4(u)*x4+bet5(u)*x5+bet0(u),u)
#' y = b*w
#'
#' ZIQSIR(y,X,m = 4,test_num = 4,method = "Pearson")
#'
#' ### demo 2
#' # simulation results under large sample size
#' # using Chi-square method
#' # under null hypothesis
#' # sample size
#' n = 2000
#'
#' # the probability of Y>0 given covariate x
#' p = function(x1,x2,x3,x4,x5,x6,
#'              gam0=2.32,gam1=-0.06,gam2=-0.03*0,gam3=-0.010*0,gam4 = -0.005,gam5 = 0.0005,gam6 = -0.030){
#'   lc = gam0 + gam1*x1 + gam2*x2 + gam3*x3 + gam4*x4 + gam5*x5 + gam6*x6
#'   exp(lc)/(1+exp(lc))
#' }
#'
#' # beta_tau
#' bet1 = function(x){5*x^2+1}
#' bet2 = function(x){(0.1*sin(x*2*pi)+0.05)*0}
#' bet3 = function(x){
#'   ( 0.05*(x-0.5)^2+0.04)*0
#' }
#' bet4 = function(x){(0.1*sqrt(x)+2.9*x)*0.05}
#' bet5 = function(x){(0.4*(x-1)^2)*0.1}
#' bet6 = function(x){((x-0.6)*(x-1.1))*1.1}
#' bet0 = function(x){62.9*x^2+33.4*x}
#' bet = function(x,u){x^4*u*10^(-4)*0.4+x^3*u*10^{-3}*0.1}
#'
#' # G_tau function
#' func <- function(x, tau)
#' {
#'   return(bet(x %*% rbind(bet1(tau),bet2(tau),bet3(tau),bet4(tau),bet5(tau),bet6(u))+bet0(tau),tau))
#' }
#'
#' # given samples
#' x1 = rbinom(n,1,0.5) # Gender
#' x2 = rnorm(n,28,2) # bmi
#' x3 = 2*x2+rnorm(n,36.5,9) # waist
#' x4 = rnorm(n,80,12) # diastolic_bp
#' x5 = x4*1.3+rnorm(n,20,7.75) # systolic_bp
#' x6 = sample(1:4, n, replace = TRUE, prob = c(0.25,0.025,0.07,0.625))
#' x0 = rep(1,n)
#' u = runif(n)
#' b = rbinom(n,1,p(x1,x2,x3,x4,x5,x6))
#' w = bet(bet1(u)*x1+bet2(u)*x2+bet3(u)*x3+bet4(u)*x4+bet5(u)*x5+bet6(u)*x6+bet0(u),u)
#' y = b*w
#' X = cbind(x1,x2,x3,x4,x5,x6)
#'
#' ZIQSIR(y,X,m = 3,test_num = c(2,3),method = "Chi")

# Hypothesis testing for the "test_num"th covariate(s)
# choose method according to the sample size: large sample size using "Chi", small sample size using "Pearson"
ZIQSIR = function(y, X, taus = c(0.1,0.25,0.5,0.75,0.9), m=3, test_num,method = "Chi"){
  # getting the p-value for logistic regression
  b = 1*(y>0)
  mod.logistic = glm(b ~ X, family=binomial(link = 'logit'))
  mod.logistic.null = glm(b ~ X[,-c(test_num)], family=binomial(link = 'logit'))
  pvalue.logistic = anova(mod.logistic.null, mod.logistic, test="Rao")$`Pr(>Chi)`[2]

  X = cbind(1,X)
  test_num = test_num+1
  MM = rep(0,length(taus))
  p_value = rep(0,length(taus))
  if (method == "Chi"){
    for (j in 1:length(taus)){
      MM[j] = test_stats(y,X,taus = taus[j],m,test_num)
      p_value[j] = 1 - pchisq(unlist(MM[j]),length(test_num))
    }
  }
  else if (method == "Pearson"){
    for (j in 1:length(taus)){
      MM[j] = fKMQR.test(y,X,tau = taus[j],m,test_num)
      p_value = MM
    }
  }
  else{
    return("Please enter a method that can be used: 'Chi', 'Pearson'")
  }
  # getting the p-value for single-index quantile regression

  # getting the combining proportions for each p-values
  omega = lapply(taus,function(x){
    a = taus
    s = sum(a*(a-0.5<0) + (1-a)*(1 - (a-0.5<0)))
    (x*(x-0.5<0) + (1-x)*(1 - (x-0.5<0)))/s
  })
  rn = length(y[y==0])/length(y)
  omega = (1-rn)*unlist(omega)

  # getting the p-value for hypothesis testing
  Test_stat = rn*tan((0.5-pvalue.logistic)*pi) + sum(omega*tan((0.5-p_value)*pi))
  oo1 = 1-pcauchy(Test_stat)
  return(oo1)
}

