#### package to be installed ####
install.packages(CEoptim)
install.packages(corrplot)
install.package(DescTools)
install.packages(glasso)
install.packages(glmnet)
install.packages(meboot)
install.packages(grplasso)
install.packages(pls)
install.packages(EnvStats)
#################################

setwd("C://Users/jhjoo/Desktop/var/final.codes/")
source("codes/var_utils.R")

library(glmnet)
library(DescTools)
library(CEoptim)
library(glasso)
library(meboot)
library(pls)
detach("package:corrplot", unload = TRUE)
library(corrplot)
library(EnvStats)
library(MASS)

#### Parameters ####
Times=53
N=Times-1
#n=N-1 # P=1
K=30

################################
#### Generate Random Thetas ####
################################
# Simulated sales Growth series only (N observation) 
################################
coef.sim=function(q, sparsity, b.lim, m.diag, metric="None"){
  ################################
  # q: dimension of coef matrix
  # sparsity: #valued within the full matrix.
  # b.lim: Given the positive/negative value in coefficient matrix.
  # m.diag: Given the value of diganal(own) effects
  ################################
  if(sparsity < (1/q)){stop("Sparsity error. Should greater than 1/q")}
  beta <- matrix(0, q, q)
  diag(beta) <- m.diag
  valued = sample(which(upper.tri(beta)|lower.tri(beta)), max(1, q*q*sparsity-q))
  beta[valued] <- c(b.lim[2], sample(b.lim, length(valued)-1, T))
  print(unitroot.test(beta)) # generate the reasonable coefficient
  colnames(beta) <- rep(metric, ncol(beta))
  return(beta)
}
sigma.sim=function(q, sd.vec){
  ################################
  # q: dimesnion of the error covariance matrix
  # sd.vec: Standard error of each variables on left-handed side
  ################################
  cor<-list()
  cov.raw<-matrix(sd.vec, q, 1) %*% matrix(sd.vec, 1, q)
  
  # Diagnal: Identity correlations
  diagonal<-diag(q)
  cov.diagonal<-diagonal*cov.raw
  cor$diagonal<-cov2cor(cov.diagonal)
  
  # Sparse: Randomly picked q cross-correlation = 0.3
  sparse<-diagonal
  sparse[sample(which(lower.tri(diag(q))))[1:q]] <- sample(c(-0.3, 0.3), q, T)
  cov.sparse = sparse %*% t(sparse) * cov.raw
  cor$sparse = cov2cor(cov.sparse)
  
  s=round(0.8*q*(q-1)/2)
  # Dense: Randomly picked s cross-correlation = 0.3
  dense<-diagonal
  dense[sample(which(lower.tri(diag(q))))[1:s]] <- sample(c(-0.3, 0.3), s, T)
  cov.dense = dense %*% t(dense) * cov.raw
  cor$dense = cov2cor(cov.dense)
  
  out<-list(diagonal=cov.diagonal, 
            sparse=cov.sparse,
            dense=cov.dense, 
            cor=cor,
            sigma=sd.vec)
}

set.seed(2020)
beta.star = coef.sim(K, 60/900, c(-0.25, 0.25), rep(-0.5, 30), "SALES")
set.seed(2019)
alpha.star = coef.sim(K, 60/900, c(-0.15, 0.15), rep(-0.3, 30), "PZ.REDUCT")
set.seed(2018)
gamma.star = coef.sim(K, 40/900, c(-0.15, 0.15), rep(0.3, 30), "AD.RATE")
set.seed(2017)
delta.star = coef.sim(K, 40/900, c(-0.15, 0.15), rep(0.3, 30), "DI.RATE")
theta.star = cbind(beta.star, alpha.star, gamma.star, delta.star)
set.seed(2016)
sigma.star.dense = sigma.sim(K, rep(0.2, 30))$dense
sigma.star.diagonal = sigma.sim(K, rep(0.2, 30))$diagonal
sigma.star.sparse = sigma.sim(K, rep(0.2, 30))$sparse

save(theta.star, file="data/truth/theta.star.RData")
save(sigma.star.sparse, file="data/truth/sigma.star.sparse.RData")
save(sigma.star.diagonal, file="data/truth/sigma.star.diag.RData")
save(sigma.star.dense, file="data/truth/sigma.star.diag.RData")

#### DGP & Variables for training data ####
pr.sim=function(q, time, min, max, mode){
  ################################
  # Triangle distribution with min/max/mode 
  ################################
  if(length(min) == 1){min = rep(min, q)}
  if(length(max) == 1){max = rep(max, q)}
  if(length(mode) == 1){mode = rep(mode, q)}
  
  pr.data = matrix(0, time, q)
  for(i in 1:q){pr.data[, i] = rtri(time, min[i], max[i], mode[i])}
  colnames(pr.data) <- paste("PZ.REDUCT", 1:ncol(pr.data), sep='.')
  return(pr.data)
}
addi.sim=function(q, time, type=rep(1, K), metric){
  ################################
  # Beta distribution with alpha/beta parameter
  # type==1 for the normal
  # type==2 for high-proprotional 0
  ################################
  if(length(type)<q){type=sample(c(1,k), q, replace=T)}
  fd.data = matrix(0, time, q)
  a = c(2.5, 0.5); b=c(34, 26)
  for(i in 1:q){
    fd.data[, i] = rbeta(time, a[type[i]], b[type[i]])
  }
  colnames(fd.data) <- paste(metric, 1:ncol(fd.data), sep='.')
  return(fd.data)
}
dgp=function(K, time, S, beta, alpha, gamma, delta, Sigma){
  # K: dimensional.
  # Beta, alpha, gamma, delta: Coefficients of Ylag, PR, AD, DI
  # Noise: multivariate normal with Sigma covariance.
  # Return: list of Y and list of X and Residual
  data=list()
  for(s in 1:S){
    pr=pr.sim(K, time, -0.1, 0.01, 0)
    ad=addi.sim(K, time, rep(1, K), "AD.RATE")
    di=addi.sim(K, time, rep(2, K), "DI.RATE")
    pr.std=stdize(pr)
    ad.std=stdize(ad)
    di.std=stdize(di)
    residual=mvrnorm(time, rep(0, K),  Sigma)
    
    y=matrix(0, ncol=K, nrow=time)
    y[1, 1:K] = residual[1,]
    
    for(t in 2:time){
      y[t, 1:K] <- 
        beta %*% matrix(y[t-1, ]) + 
        alpha %*% matrix(pr.std[t, ]) +
        gamma %*% matrix(ad.std[t, ]) + 
        delta %*% matrix(di.std[t, ]) +
        residual[t, ]
    }
    
    colnames(y) = paste("SALESG", 1:K, sep=".")
    colnames(residual) = paste("Residual", 1:K, sep=".")
    
    data[[s]] = list(Y=y[-1, ], 
                     X=cbind(y[(1:time-1), ],
                             pr.std[-1, ], 
                             ad.std[-1, ],
                             di.std[-1, ]),
                     X.ex=cbind(pr, ad, di),
                     R=residual)
  }
  return(data)
}

################################
# S simulation data of train/test
################################
set.seed(2020)
SIM.data=dgp(K, N+N+1, 50,
             beta.star, alpha.star, gamma.star, delta.star,
             sigma.star.sparse)
save(SIM.data, file='data/SIM.data.lag1.RData')


################################
# Simulation follow wilm's setting
################################
N=50
K=5
p=2
q=K*p

b1=diag(0.4, q/p)
b1[, 1]=0.4
B1.star=kronecker(diag(p), b1)
b2=diag(0.2, q/p)
b2[, 1]=0.2
B2.star=kronecker(diag(p), b2)
theta.star=cbind(B1.star, B2.star)
sigma.star=diag(rep(0.1, q))
save(theta.star, file='data/truth/theta.star.wilms.RData')
save(sigma.star, file='data/truth/sigma.star.wilms.RData')

dgp.wilms=function(q, time, S, B1, B2, Sigma){
  data=list()
  p=2
  for(s in 1:S){
    residual=matrix(0, nrow=time, ncol=q)
    y=matrix(0, ncol=q, nrow=time)
    y[1:p, ] = rnorm(q*p)

    for(i in 3:time){
      residual[i, ]=mvrnorm(n=1, mu=rep(0, q), Sigma)
      y[i,]=B1%*%y[i-1,] + B2%*%y[i-2,] + matrix(residual[i, ])
    }
    
    colnames(y) = paste("SALESG", 1:q, sep=".")
    colnames(residual) = paste("Residual", 1:q, sep=".")
    
    data[[s]] = list(Y=y[-c(1:p), ], 
                     X=y[-c(1,time), ],
                     R=residual)
  }
  return(data)
}
WILMS.data=dgp.wilms(q, N+N+2, 50, B1.star, B2.star, sigma.star)
save(WILMS.data, file='data/WILMS.data.RData')
