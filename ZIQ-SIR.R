library(quantreg)
library(aod)

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

test_stats = function(Y,X,taus,m,test_num){
  y = Y[which(Y>0)]
  x = X[which(Y>0),-test_num]
  
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
  
  score = lapply(1:length(taus),function(ii){
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
    s = 0
    s1 = 0
    
    X.star = X[which(Y>0),] - expect
    for (j in 1:n){
      s1 = s1 + (X.star[j,]*dev[j]) %*% t(X.star[j,]*dev[j])
      s = s+(X.star[j,test_num])*dev[j]*dev2[j]
    }
    Omega = s1[test_num,test_num] - s1[test_num,-test_num] %*% MASS::ginv(s1[-test_num,-test_num]) %*% s1[-test_num,test_num]
    ss = s %*% MASS::ginv(Omega) %*% s *taus[ii]^(-1) * (1-taus[ii])^(-1)
    return(ss)
  })
  
  return(score)  
}

# constructing the score for hypothesis testing and approximated variance-covariance for the score
test_stat_separate = function(Y,X,taus,m,test_num){
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

# using Pearson Type III distribution to get the p-value for hypothesis testing for single-index quantile regression 
fKMQR.test = function(Y, X, tau, m, test_num, score = NULL, K = NULL) {
  ## fast Kernel machine quantile regression
  ### X nonparamteric var (nxp)
  ## Z parametric var (nxq)
  ## Y response var (nx1)
  ## tau quantile
  
  ## Define some auxiliary functions: tr and KRV function
  KRV=function(K,L){
    n=nrow(K)
    A=scale(K,scale=F) ## that is A=HK
    W=scale(L,scale=F) ## that is W=HL
    Fstar=tr(A%*%W)
    mean.krv=tr(A)*tr(W)/(n-1)	## mean of KRV 
    T=tr(A);T2=tr(A%*%A);S2=sum(diag(A)^2)
    Ts=tr(W);T2s=tr(W%*%W);S2s=sum(diag(W)^2)
    temp1=2*((n-1)*T2-T^2)*((n-1)*T2s-Ts^2)/(n-1)^2/(n+1)/(n-2)
    temp21=n*(n+1)*S2- (n-1)*(T^2+2*T2)
    temp22=n*(n+1)*S2s- (n-1)*(Ts^2+2*T2s)
    temp23=(n+1)*n*(n-1)*(n-2)*(n-3)
    temp2=temp21*temp22/temp23
    variance.krv=temp1+temp2		## variance of KRV
    T3=tr(A%*%A%*%A);S3=sum(diag(A)^3);U=sum(A^3);R=t(diag(A))%*%diag(A%*%A);B=t(diag(A))%*%A%*%diag(A)
    T3s=tr(W%*%W%*%W);S3s=sum(diag(W)^3);Us=sum(W^3);Rs=t(diag(W))%*%diag(W%*%W);Bs=t(diag(W))%*%W%*%diag(W)
    t1=n^2*(n+1)*(n^2+15*n-4)*S3*S3s
    t2=4*(n^4-8*n^3+19*n^2-4*n-16)*U*Us
    t3=24*(n^2-n-4)*(U*Bs+B*Us)
    t4=6*(n^4-8*n^3+21*n^2-6*n-24)*B*Bs
    t5=12*(n^4-n^3-8*n^2+36*n-48)*R*Rs
    t6=12*(n^3-2*n^2+9*n-12)*(T*S2*Rs+R*Ts*S2s)
    t7=3*(n^4-4*n^3-2*n^2+9*n-12)*T*Ts*S2*S2s
    t81=(n^3-3*n^2-2*n+8)*(R*Us+U*Rs);t82=(n^3-2*n^2-3*n+12)*(R*Bs+B*Rs)
    t8=24*(t81+t82)
    t9=12*(n^2-n+4)*(T*S2*Us+U*Ts*S2s)
    t10=6*(2*n^3-7*n^2-3*n+12)*(T*S2*Bs+B*Ts*S2s)
    t11=-2*n*(n-1)*(n^2-n+4)*((2*U+3*B)*S3s+(2*Us+3*Bs)*S3)
    t12=-3*n*(n-1)^2*(n+4)*((T*S2+4*R)*S3s+(Ts*S2s+4*Rs)*S3)
    t13=2*n*(n-1)*(n-2)*((T^3+6*T*T2+8*T3)*S3s+(Ts^3+6*Ts*T2s+8*T3s)*S3)
    t14=T^3*((n^3-9*n^2+23*n-14)*Ts^3+6*(n-4)*Ts*T2s+8*T3s)
    t15=6*T*T2*((n-4)*Ts^3+(n^3-9*n^2+24*n-14)*Ts*T2s+4*(n-3)*T3s)
    t16=8*T3*(Ts^3+3*(n-3)*Ts*T2s+(n^3-9*n^2+26*n-22)*T3s)
    t17=-16*(T^3*Us+U*Ts^3)-6*(T*T2*Us+U*Ts*T2s)*(2*n^2-10*n+16)
    t18=-8*(T3*Us+U*T3s)*(3*n^2-15*n+16)-(T^3*Bs+B*Ts^3)*(6*n^2-30*n+24)
    t19=-6*(T*T2*Bs+B*Ts*T2s)*(4*n^2-20*n+24)-8*(T3*Bs+B*T3s)*(3*n^2-15*n+24)
    t201=24*(T^3*Rs+R*Ts^3)+6*(T*T2*Rs+R*Ts*T2s)*(2*n^2-10*n+24)
    t202=8*(T3*Rs+R*T3s)*(3*n^2-15*n+24)+(3*n^2-15*n+6)*(T^3*Ts*S2s+T*S2*Ts^3)
    t203=6*(T*T2*Ts*S2s+Ts*T2s*T*S2)*(n^2-5*n+6)+48*(T3*Ts*S2s+T3s*T*S2)
    t20=-(n-2)*(t201+t202+t203)
    temp31=t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20
    temp32=n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)
    mom3=temp31/temp32
    skewness.krv= (mom3-3*mean.krv*variance.krv-mean.krv^3)/variance.krv^1.5 ## skewness of KRV
    m1=mean.krv
    m2=variance.krv
    m3=skewness.krv
    shape=4/m3^2
    scale=sqrt(m2)*m3/2
    location=m1-2*sqrt(m2)/m3
    PIIIpars=list(shape,location,scale)
    pv=1-PearsonDS::ppearsonIII(Fstar, params=PIIIpars) 
    return(pv)
  }
  
  tr=function(x){return(sum(diag(x))) }
  #if(!is.null(K)){
  #  K = IBS.kernel(X) ##kernel matrix
  #}
  RRR = test_stat_separate(Y,X,taus = tau,m = m,test_num = test_num)
  w = RRR[[1]]$score
  K1 = RRR[[1]]$W
  
  Kw=w %*% t(w)
  pv=KRV(Kw,K1)
  return(c(pv))
}

# Combining the p-values for logistic regression and single-index quantile regression
# choose method according to the sample size: large sample size using "Chi", small sample size using "Pearson"
Combination = function(y, X, taus = c(0.1,0.25,0.5,0.75,0.9), m, test_num,method = "Chi"){
  # getting the p-value for logistic regression
  b = 1*(y>0)
  mod.logistic = glm(b ~ X[,-c(1)], family=binomial(link = 'logit'))
  mod.logistic.null = glm(b ~ X[,-c(1,test_num)], family=binomial(link = 'logit'))
  pvalue.logistic = anova(mod.logistic.null, mod.logistic, test="Rao")$`Pr(>Chi)`[2]
  
  MM = rep(0,length(taus))
  p_value = rep(0,length(taus))
  if (method == "Chi"){
    for (j in 1:length(taus)){
      MM[j] = test_stats(y,X,tau = taus[j],m,test_num)
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

