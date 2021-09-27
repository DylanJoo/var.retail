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

# DGP
# step(1): Retrieve true theta, true error covariance matrix.
# step(2): Generate the exogenous variables if necessary.
# Step(3): Decide the simulation setting of the data, say lag order, simluation times.

#### step(1) ####
# - sim-en-q: Wilms settings (q=5, 10, 30)
q=30
theta.star=diag(0.4, q)
theta.star[, 1]=0.4

sigma.star=diag(rep(0.1, q))
save(theta.star, file=sprintf('data/simulation/sim-en-%d/theta.star.RData', q))
save(sigma.star, file=sprintf('data/simulation/sim-en-%d/sigma.star.RData', q))

# - sim-ex-q: Wilms settings with 1 exogenous (q=5, 10, 30)
q=30
theta.star=cbind(diag(0.4, q), diag(0.2, q))
theta.star[, 1]=0.4
theta.star[, (1+q)]=0.2

sigma.star=diag(rep(0.1, q))
save(theta.star, file=sprintf('data/simulation/sim-ex-%d/theta.star.RData', q))
save(sigma.star, file=sprintf('data/simulation/sim-ex-%d/sigma.star.RData', q))

# -sim-retail-q: IRI setting with 3 exogenous
q=30
coef.sim=function(q, sparsity, b.lim, m.diag, metric="None", seed=8888){
  ################################
  # q: dimension of coef matrix
  # sparsity: #valued within the full matrix.
  # b.lim: Given the positive/negative value in coefficient matrix.
  # m.diag: Given the value of diganal(own) effects
  ################################
  set.seed(seed)
  if(sparsity < (1/q)){stop("Sparsity error. Should greater than 1/q")}
  beta <- matrix(0, q, q)
  diag(beta) <- m.diag
  valued = sample(which(upper.tri(beta)|lower.tri(beta)), max(1, q*q*sparsity-q))
  beta[valued] <- c(b.lim[2], sample(b.lim, length(valued)-1, T))
  print(unitroot.test(beta)) # generate the reasonable coefficient
  colnames(beta) <- rep(metric, ncol(beta))
  return(beta)
}
sigma.sim=function(q, sd.vec, seed=8888){
  set.seed(seed)
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
beta.star = coef.sim(30, 60/900, c(-0.25, 0.25), rep(-0.5, 30), "SALESG", seed=1)
alpha.star = coef.sim(30, 60/900, c(-0.15, 0.15), rep(-0.3, 30), "PZ.CHANGE", seed=2)
gamma.star = coef.sim(30, 40/900, c(-0.15, 0.15), rep(0.3, 30), "AD.RATE", seed=3)
delta.star = coef.sim(30, 40/900, c(-0.15, 0.15), rep(0.3, 30), "DI.RATE", seed=4)
theta.star = cbind(beta.star, alpha.star, gamma.star, delta.star)
sigma.star = sigma.sim(30, rep(0.2, 30), seed=5)$sparse
# sigma.star.dense = sigma.sim(30, rep(0.2, 30))$dense
# sigma.star.diagonal = sigma.sim(30, rep(0.2, 30))$diagonal

save(theta.star, file=sprintf('data/simulation/sim-iri-%d/theta.star.RData', q))
save(sigma.star, file=sprintf('data/simulation/sim-iri-%d/sigma.star.RData', q))

#### step(2) ####
#### step(3) ####
dgp.simple=function(k, time, S, B, A=NA, Sigma){
  data=list()
  p=ncol(B)/k
  #q=length(A)/q
  
  k.left=nrow(B)
  k.right=ncol(B)+ ifelse(is.na(A), 0, ncol(A))
  
  timefull=(time)*2 + p
  
  # simulation process
  for(s in 1:S){
    residual=mvrnorm(n=timefull, mu=rep(0, k.left), Sigma=Sigma)
    Y=matrix(0, ncol=k.left, nrow=timefull)
    Y[1:p, 1:k.left]=residual[1:p, ]
    if(is.matrix(A)){
      X=matrix(rbinom(timefull*ncol(A), 1, 0.5), nrow=timefull, ncol=ncol(A))
      colnames(X)=paste("PZ.CHANGE", 1:ncol(A), sep='.')
    }
    
    # auto-regressive process
    for(t in (p+1):timefull){
      y=B %*% c(t(Y[(t-1):(t-p), ])) + residual[t, ]
      if(is.matrix(A)){y = y + A %*% c(t(X[t:(t-p+1), ]))}else{X=NULL}
      
      Y[t, 1:k.left]=y
    }
    
    colnames(Y)=paste("SALESG", 1:k, sep=".")
    colnames(residual)=paste("Residual", 1:k, sep=".")
    
    data[[s]]=list(Y=Y, X=X, R=residual, p=p)
  }
  return(data)
}
pr.sim=function(k, time, min, max, mode){
  ################################
  # Triangle distribution with min/max/mode 
  ################################
  if(length(min) == 1){min = rep(min, k)}
  if(length(max) == 1){max = rep(max, k)}
  if(length(mode) == 1){mode = rep(mode, k)}
  
  pr.data = matrix(0, time, k)
  for(i in 1:q){pr.data[, i] = rtri(time, min[i], max[i], mode[i])}
  colnames(pr.data) <- paste("PZ.CHANGE", 1:ncol(pr.data), sep='.')
  return(pr.data)
}
addi.sim=function(k, time, type, metric){
  ################################
  # Beta distribution with alpha/beta parameter
  # type==1 for the normal
  # type==2 for high-proprotional 0
  ################################
  if(length(type)<q){type=sample(c(1,k), k, replace=T)}
  fd.data = matrix(0, time, q)
  a = c(2.5, 0.5); b=c(34, 26)
  for(i in 1:k){
    fd.data[, i] = rbeta(time, a[type[i]], b[type[i]])
  }
  colnames(fd.data) <- paste(metric, 1:ncol(fd.data), sep='.')
  return(fd.data)
}
dgp=function(k, time, S, B, A, Sigma){

  data=list()
  p=ncol(B)/k
  #q=length(A)/q
  
  k.left=nrow(B)
  k.right=ncol(B)+ ifelse(is.na(A), 0, ncol(A))
  
  timefull=time*2 + p
  
  for(s in 1:S){
    residual=mvrnorm(n=timefull, mu=rep(0, k.left), Sigma=Sigma)
    Y=matrix(0, ncol=k.left, nrow=timefull)
    Y[1:p, 1:k.left]=residual[1:p, ]
    
    pr=pr.sim(k.left, timefull, -0.1, 0.01, 0)
    ad=addi.sim(k.left, timefull, rep(1, k.left), "AD.RATE")
    di=addi.sim(k.left, timefull, rep(2, k.left), "DI.RATE")
    X=cbind(stdize(pr), stdize(ad), stdize(di))
    # X=cbind(pr, ad, di)

    for(t in (p+1):timefull){
      y=B %*% c(t(Y[(t-1):(t-p), ]))+
        A %*% c(t(X[t:(t-p+1), ]))+
        residual[t, ] 
      Y[t, 1:k.left]=y
    }
    
    colnames(Y) = paste("SALESG", 1:k.left, sep=".")
    colnames(residual) = paste("Residual", 1:k.left, sep=".")
    
    data[[s]] = list(Y=Y, 
                     X=X,
                     R=residual,
                     p=p,
                     X.raw=cbind(pr, ad, di))
  }
  return(data)
}

N=52
# - sim-en-q: Wilms settings (q=5, 10, 30)
q=30
load(file=sprintf("data/simulation/sim-en-%d/theta.star.RData", q))
load(file=sprintf("data/simulation/sim-en-%d/sigma.star.RData", q))
SIM.data=dgp.simple(k=q, time=N, S=50, B=theta.star, A=NA, Sigma=sigma.star)
save(SIM.data, file=sprintf('data/simulation/sim-en-%d/data.RData', q))

# - sim-ex-q: Wilms settings (q=5, 10, 30)
q=30
load(file=sprintf("data/simulation/sim-ex-%d/theta.star.RData", q))
load(file=sprintf("data/simulation/sim-ex-%d/sigma.star.RData", q))
SIM.data=dgp.simple(k=q, time=N, S=50, B=theta.star[, 1:q], A=theta.star[, -(1:q)], Sigma=sigma.star)
save(SIM.data, file=sprintf('data/simulation/sim-ex-%d/data.RData', q))

# -sim-retail-q: IRI setting with 3 exogenous
q=30
load(file=sprintf("data/simulation/sim-iri-%d/theta.star.RData", q))
load(file=sprintf("data/simulation/sim-iri-%d/sigma.star.RData", q))
SIM.data=dgp(k=q, time=N, S=50, B=theta.star[, 1:q], A=theta.star[, -(1:q)], Sigma=sigma.star)
save(SIM.data, file=sprintf('data/simulation/sim-iri-%d/data.RData', q))


# dgp.wilms=function(q, time, S, B1, B2, Sigma){
#   data=list()
#   p=2
#   for(s in 1:S){
#     residual=matrix(0, nrow=time, ncol=q)
#     y=matrix(0, ncol=q, nrow=time)
#     y[1:p, ] = rnorm(q*p)
# 
#     for(i in 3:time){
#       residual[i, ]=mvrnorm(n=1, mu=rep(0, q), Sigma)
#       y[i,]=B1%*%y[i-1,] + B2%*%y[i-2,] + matrix(residual[i, ])
#     }
#     
#     colnames(y) = paste("SALESG", 1:q, sep=".")
#     colnames(residual) = paste("Residual", 1:q, sep=".")
#     
#     data[[s]] = list(Y=y[-c(1:p), ], 
#                      X=y[-c(1,time), ],
#                      R=residual)
#   }
#   return(data)
# }
