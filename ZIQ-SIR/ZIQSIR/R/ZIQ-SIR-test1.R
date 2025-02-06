#' Test statistic computation for hypothesis testing
#'
#' @param Y n*1 vector, the observed outcome for inference
#' @param X n*p matrix, the observed covariates for inference
#' @param taus k*1 vector, a grid of quantile levels; e.g., 0.5 for the median, 0.75 for the 3rd quartile; default is c(0.1, 0.25, 0.5, 0.75, 0.9)
#' @param m numeric variable, the order of B-spline function; default is 3
#' @param test_num a vector, representing the test corresponds to which covariate(s) in X.
#'
#' @return test statistics for method 'Pearson'
#' @export

# constructing the score for hypothesis testing and approximated variance-covariance for the score
test_stat_separate = function(Y,X,taus,m,test_num){
  rho = function(x,tau){tau*x - x*(x<0)} # score function

  profile_likelihood_beta = function(u,X1,Y1,tau1,Nn,m)
    # to get the score function regarding to single-index quantile regression given beta(x)
  {
    Z1 = vector(length = length(Y1));
    y = vector(length = ncol(X1));
    # normalization of beta
    y = u/rep(sqrt(crossprod(u)),length(u))
    # the single index given normalized beta
    Z1 = X1 %*% y;

    # the structure of single-index B-spline regression
    S = quantreg::rq(Y1~splines::bs(Z1,df = Nn+m+1, degree = m, intercept = T)-1,tau = tau1);
    # the difference between real value and predicted value
    Res = Y1 - predict(S);
    # calculate the score
    Resq = rho(Y1-predict(S),tau1)

    sumit = sum(Resq) # sum the score
    return(sumit)
  }

  # calculating BIC to get the coefficient theta
  BIC = function(Z1,Y1,tau1,Nn,m){
    n = length(Y1)
    S = quantreg::rq(Y1~splines::bs(Z1,df = Nn+m+1, degree = m, intercept = T)-1,tau = tau1);
    Resq = (Y1-predict(S))*tau1-(Y1-predict(S)<0)*(Y1-predict(S))
    sumit = log(sum(Resq)/n) +  log(n)/2/n*(Nn+m)# sum the score
    return(sumit)
  }
  y1 = Y[which(Y>0)]
  y = y1
  #y = quantreg::dither(log(y1))
  x = X[which(Y>0),-c(test_num)]

  n = length(y)
  l = ncol(x)
  u = rep(1/sqrt(l),l)
  beta = matrix(0,nrow = length(taus),ncol = ncol(x))
  hat_Nn = floor(n^{1/(2*m+1)})+1

  # estimation for the coefficient beta_tau
  for (i in 1:length(taus))
  {
    profile_score_beta = function(v)
    {
      return(profile_likelihood_beta(v,x,y,taus[i],hat_Nn,m))
    }
    beta1 = optim(u,profile_score_beta)$par
    beta[i, ] = beta1/sqrt(sum(beta1*beta1))
  }
  # estimation for the G_tau function
  theta = lapply(1:length(taus),function(ii){
    Z = x %*% beta[ii, ]

    # choosing the knots num Nn using BIC defined above
    h1 = ceiling(hat_Nn/2)
    h2 = ceiling(hat_Nn*3/2)
    h = h1:h2
    Res = lapply(1:length(h), function(j){BIC(Z, y, taus[ii], h[j], m)})
    Res = unlist(Res)
    Nn = h[Res==min(Res)]
    basis = splines::bs(Z,df = Nn+m+1, degree = m, intercept = T)
    # using the Nn with the locally minimum BIC to estimate the G function
    gamma = quantreg::rq(y ~ basis-1, tau = taus[ii])
    u = list(gamma = gamma, Nn = Nn, basis = basis)
    return(u)
  })

  result = lapply(1:length(taus),function(ii){
    Z = x %*% beta[ii, ]
    L = theta[[ii]]$gamma
    Nn = as.numeric(theta[[ii]]$Nn)

    D = theta[[ii]]$basis
    expect = D %*% solve(t(D) %*% D) %*% t(D) %*% X[which(Y>0),]
    # calculate the estimated value of quantile curve G_tau
    Gtau =  predict(L,data.frame(Z))

    # calculate the estimated derivative of quantile curve G_tau
    dev = splines2::dbs(Z, df = Nn+m+1, degree = m, intercept = T) %*% L$coefficients
    dev2 = taus[ii] - (y-Gtau<0)
    s1 = 0

    X.star = X[which(Y>0),] - expect
    for (j in 1:n){
      s1 = s1 + (X.star[j,]*dev[j]) %*% t(X.star[j,]*dev[j])
    }
    Omega = s1[test_num,test_num] - s1[test_num,-test_num] %*% MASS::ginv(s1[-test_num,-test_num]) %*% s1[-test_num,test_num]
    dev = matrix(dev,nrow = nrow(X.star),ncol = length(test_num))
    hat_X = dev*as.matrix(X.star[,test_num])
    hat_Omega = hat_X %*% solve(Omega) %*% t(hat_X)
    hat_Omega = (1-taus[ii])^{-1}*taus[ii]^{-1}*hat_Omega
    LLLL = list(score = dev2, W = hat_Omega)
    return(LLLL)
  })
  return(result)
}
