#### package to be installed ####
install.packages(CEoptim)
install.packages(corrplot)
install.package(DescTools)
install.packages(glasso)
install.packages(glmnet)
install.packages(meboot)
install.packages(grplasso)
install.packages(pls)
#################################

setwd("C://Users/jhjoo/Desktop/var/final.codes/")
source("codes/var_model.R")
source("codes/var_utils.R")

library(glmnet)
library(DescTools)
library(CEoptim)
library(glasso)
library(meboot)
library(pls)
detach("package:corrplot", unload = TRUE)
library(corrplot)

simulation=function(s, forecast='naive', DATA=SIM.data, estimate='igls'){
  
  # matrix prepration
  N=nrow(DATA[[s]]$Y)/2 #52
  K=ncol(DATA[[s]]$Y)   #30
  M=ncol(DATA[[s]]$X)
  nameleft=colnames(DATA[[1]]$Y)
  nameright=colnames(DATA[[1]]$X)
  
  if(forecast=='naive'){
    data=matrixPrep(DATA[[s]]$Y, DATA[[s]]$X, DATA[[s]]$R, 1, N)
    Xtrain=kronecker(diag(K), data$Xtrain)
    Ytrain=c(data$Ytrain)
    
    # Estimation
    if(estimate=='igls'){
      fit=var.iter.est(Xtrain, Ytrain, K, 
                       1/seq(1e2, K, length.out = 100), 
                       1/seq(1e1, K, length.out = 100)**2, 
                       T, 'BIC2', c(T, T), 'eigen',
                       10^-2, 10^-1, 10, NA, T, 
                       list(nameleft, nameright), 'standard', 1)
      
      theta=theta.fn(nameleft, nameright, fit$beta.est, K)
      sigma=theta.fn(nameleft, nameleft, fit$sigma.est, K)
    
    } else if (estimate=='ceoptim'){
      
      if(length(thetas) >= s){
        pick = which(thetas[[s]]!=0)
        fit=CEoptim(f=ols.l0,
                    f.arg=list(Y=data$Ytrain, X=data$Xtrain, k=0, lambda=0,
                               warm=T, Theta.igls=rep(1, length(pick))),
                    continuous = list(mean=c(t(thetas[[s]]))[pick],
                                      sd=rep(1, length(pick))),
                    discrete= list(probs=rep(list(c(0.5, 0.5)), length(pick))),
                    maximize = F, verbose = T, N = 10000,
                    iterThr = 50, rho = 0.01, noImproveThr = 5)
        theta=c(t(thetas[[s]]))
        theta[pick]=fit$optimizer$continuous * fit$optimizer$discrete
        theta=theta.fn(nameleft, nameright, theta, K)
        
      } else {
        fit=CEoptim(f=ols,
                    f.arg=list(Y=data$Ytrain, X=data$Xtrain),
                    continuous = list(mean=rep(0, K*M), sd=rep(0, K*M)),
                    maximize = F, verbose = T, N = 10000,
                    iterThr = 50, rho = 0.01, noImproveThr = 5)
        theta=theta.fn(nameleft, nameright, fit$optimizer$continuous, K)
      }
      
      
      # fit=CEoptim(f=ols.l0,
      #             f.arg=list(Y=data$Ytrain, X=data$Xtrain, k=50, lambda=1),
      #             continuous = list(mean=rep(0, K*M), sd=rep(1, K*M)),
      #             discrete = list(probs=rep(list(c(0.5, 0.5)), K*M)),
      #             maximize = F, verbose = T, N = 10000,
      #             iterThr = 50, rho = 0.01, noImproveThr = 5)
      # theta.ce=theta.fn(nameleft, nameright, fit$optimizer$continuous * fit$optimizer$discrete, K)
      # 
      # fit=CEoptim(f=bic.l0,
      #             f.arg=list(Y=data$Ytrain, X=data$Xtrain, sep=K*M, k.dim=K),
      #             continuous = list(mean=rep(0, K*M), sd=rep(1, K*M)),
      #             discrete = list(probs=rep(list(c(0.5, 0.5)), K*M)),
      #             maximize = F, verbose = T, N = 10000,
      #             iterThr = 50, rho = 0.01, noImproveThr = 5)
      # num=fit$optimizer$continuous * fit$optimizer$discrete 
      # theta=theta.fn(nameleft, nameright, num, K)
      sigma=theta.fn(nameleft, nameleft, cov(data$Ytrain-data$Xtrain%*%t(theta)), K)
    }
    
    
    # selection/prediction consistent
    perf=performance(theta.star[1:nrow(theta), 1:ncol(theta)], theta, K)
    print(unlist(perf))
    tpr=c(perf$TPR.en, perf$TPR.ex)
    tnr=c(perf$TNR.en, perf$TNR.ex)
    maee=perf$MAEE
    
    eval.train=evaluation(theta, data$Ytrain, data$Xtrain, NULL, NULL)
    eval.test=evaluation(theta, data$Ytest, data$Xtest, NULL, NULL)
    mafe=c(eval.train$mae, eval.test$mae)
    N=1
    i=1
  } else if(forecast=='rolling'){
    
    for(i in 1:N){
      data=matrixPrep(DATA[[s]]$Y, DATA[[s]]$X, DATA[[s]]$R, i, (i+N-1))
      Xtrain=kronecker(diag(K), data$Xtrain)
      Ytrain=c(data$Ytrain)

      theta=matrix(0, K, M)
      sigma=matrix(0, K, K)
      tpr=rep(0, 2)
      tnr=rep(0, 2)
      maee=rep(0, 2)

      
      # Estimation
      fit=var.iter.est(Xtrain, Ytrain, K, 
                       seq(1e-3, 1/K, length.out = 100), 
                       seq(1e-3, 1/K, length.out = 100), 
                       T, 'BIC2', T, 'eigen',
                       10^-2, 10^-1, 10, NA, F, 
                       list(nameleft, nameright), 'standard', 1)
      
      theta.temp=theta.fn(nameleft, nameright, fit$beta.est, K)
      sigma.temp=theta.fn(nameleft, nameleft, fit$sigma.est, K)
      theta=theta+theta.temp
      sigma=sigma+sigma.temp
    
      # selection/prediction consistent
      perf=performance(theta.star, theta.temp, K)
      tpr=tpr+c(perf$TPR.en, perf$TPR.ex)
      tnr=tnr+c(perf$TNR.en, perf$TNR.ex)
      maee=maee+perf$MAEE
      
      eval.train=evaluation(theta.temp, data$Ytrain, data$Xtrain, NULL, NULL)
      eval.test=evaluation(theta.temp, data$Ytest, data$Xtest, NULL, NULL)
      mafe=c(eval.train$mae, eval.test$mae)
    }
  }
  
  return(list(theta=theta/N, 
              sigma=sigma/N,
              tpr=tpr/N,
              tnr=tnr/N,
              maee=maee/N,
              mafe=mafe/N))
}

##########################
# Simulation from RERAIL
##########################
results=c()
results.ce=c()
thetas=list()
thetas.ce=list()
sigmas=list()
sigmas.ce=list()
S=10
load(file="data/SIM.data.lag1.RData")
load(file="data/truth/theta.star.RData")
load(file="data/truth/sigma.star.sparse.RData")

for(s in 1:S){
  cat('##################\nIteration', s, '\n##################\n')
  cat('===IGLS===\n')
  # sim=simulation(s, forecast='naive', SIM.data)
  # thetas[[s]]=sim$theta
  # sigmas[[s]]=sim$sigma
  # results=cbind(results, c(sim$tpr, sim$tnr, sim$maee, sim$mafe))

  cat('===CEOPTIM===\n')
  sim.ce=simulation(s, forecast='naive', SIM.data, estimate='ceoptim')
  thetas.ce[[s]]=sim.ce$theta
  sigmas.ce[[s]]=sim.ce$sigma
  results.ce=cbind(results.ce, c(sim.ce$tpr, sim.ce$tnr, sim.ce$maee, sim.ce$mafe))
}

##########################
# Simulation from WILM's
##########################
results=c()
results.ce=c()
thetas=list()
thetas.ce=list()
sigmas=list()
sigmas.ce=list()
S=50
load(file="data/WILMS.data.RData")
load(file="data/truth/theta.star.wilms.RData")

ols=function(ThetaC, Y, X, warm=F, Theta.igls=NA){
  if(any(warm==T)){
    Theta=Theta.igls
    Theta[which(Theta!=0)]=ThetaC
  }
  Theta = matrix(Theta, nrow=ncol(Y), byrow=T)
  mse=(1/2)*mean( (Y-X%*%t(Theta))**2 )
  return(mse)
}
ols.l0=function(ThetaC, ThetaZ, Y, X, k, lambda=1, warm=F, Theta.igls=NA){
  if(any(warm==T)){
    Theta=Theta.igls
    Theta[which(Theta!=0)]=ThetaC*ThetaZ
  }
  Theta = matrix(Theta, nrow=ncol(Y), byrow=T)
  mse=mean( (Y-X %*% t(Theta)) ** 2)
  return((1/2)*mse+lambda*pmax(0, sum(Theta!=0)-k))
}
bic=function(C, Y, X, sep, k.dim){
  Theta = C[1:sep]
  Theta = matrix(Theta, nrow=ncol(Y), byrow=T)
  
  Sigma = diag(C[sep:(sep+k.dim-1)], k.dim)
  Sigma[upper.tri(Sigma)] = Sigma[lower.tri(Sigma)] = C[-(1:(sep+k.dim-1))]
  return( nrow(Y) * log(det(Sigma)) + (Y-X %*% t(Theta)) ** 2 + sum(Theta!=0) * log(nrow(Y)))
}
bic.l0=function(C, Z, Y, X, sep, k.dim){
  Theta = C[1:sep]*Z[1:sep]
  Theta = matrix(Theta, nrow=ncol(Y), byrow=T)
  resid = Y-X %*% t(Theta)
  #Sigma = diag(apply(resid, 2, sd))
  #Sigma[lower.tri(Sigma)] = C[-(1:sep)] * Z[-(1:sep)]
  #Sigma = Sigma %*% t(Sigma)
  #return( nrow(Y) * log(det(Sigma)) + sum(resid**2) + sum(Theta!=0) * log(nrow(Y)))
  #return( sum(resid**2) + sum(Theta!=0) * log(nrow(Y)))
  return( sum(resid**2) + sum(Theta!=0)*log(nrow(Y)))
}


for(s in 1:S){
  cat('##################\nIteration', s, '\n##################\n')
  cat('===IGLS===\n')
  # sim=simulation(s, forecast='naive', WILMS.data)
  # thetas[[s]]=sim$theta
  # sigmas[[s]]=sim$sigma
  # results=cbind(results, c(sim$tpr, sim$tnr, sim$maee, sim$mafe))
  
  cat('===CEOPTIM===\n')
  sim.ce=simulation(s, forecast='naive', WILMS.data, estimate='ceoptim')
  thetas.ce[[s]]=sim.ce$theta
  sigmas.ce[[s]]=sim.ce$sigma
  results.ce=cbind(results.ce, c(sim.ce$tpr, sim.ce$tnr, sim.ce$maee, sim.ce$mafe))
}


#### (1) Standard estimation on simulation result  #### 
data=matrixPrep(SIM.data[[1]]$Y, SIM.data[[1]]$X, SIM.data[[1]]$R, 1, N)
Xtrain=kronecker(diag(K), data$Xtrain)
Ytrain=c(data$Ytrain)
nameleft=colnames(SIM.data[[1]]$Y)
nameright=colnames(SIM.data[[1]]$X)

fit=var.iter.est(Xtrain, Ytrain, K, 
                 1/seq(1e2, K, length.out = 100), 
                 1/seq(1e1, K, length.out = 100)**2, 
                 T, 'BIC2', 'eigen',
                 10^-2, 10^-1, 10, NA, F, 
                 list(nameleft, nameright), 'standard', 1)
theta=theta.fn(nameleft, nameright, fit$beta.est, K)
unlist(performance(theta.star, theta, K))

#### Visualize the theta effct of simulated result
#########################
# Figure of coefficient matrix.
#########################
pdf(file="figure/SIM-theta-BIC2.pdf")
corrplot(theta[, grep("SALES", colnames(theta))], 
         is.corr=F, tl.srt=45, title="Lag Effect", mar=c(0,0,1,0))
corrplot(theta[, grep("PZ.REDUCT", colnames(theta))], 
         is.corr=F, tl.srt=45, title="Price Reduction Effect", mar=c(0,0,1,0))
corrplot(theta[, grep("AD.RATE", colnames(theta))], 
         is.corr=F, tl.srt=45, title="Advertising Effect", mar=c(0,0,1,0))
corrplot(theta[, grep("DI.RATE", colnames(theta))], 
         is.corr=F, tl.srt=45, title="Display Effect", mar=c(0,0,1,0))
dev.off()

#########################
# Results of evaluation by MAE/MAPE.
#########################
eval.train=evaluation(theta, Y.block, X.block, 
                      NULL, NULL, 0)
eval.train.2=evaluation(fit$beta.first, Y.block, X.block, 
                        NULL, NULL, 0)
eval.test=evaluation(theta, Y.block.test, X.block.test, 
                     NULL, NULL, 0)
eval.test.2=evaluation(fit$beta.first, Y.block.test, X.block.test, 
                     NULL, NULL, 0)
perf.train=performance(theta.star, theta, K)
eval.train$mae;eval.test$mae;unlist(perf.train)
eval.train.2$mae; eval.test.2$mae
load(file="data/SIM.residual.RData")
# 0.204/0.270

####################
#[2020/02/24 update]
####################
table=data.frame(csv.predict(eval.train, 'train', F), 
                 csv.predict(eval.test, 'test', F))
write.csv(table,"csv/SIM-eval-train-BIC2.csv")


#### (2) CEoptimization from meboot 20 sample  #### 
thetas20=list()
l1s = list()
eval.train = list()
eval.test = list()
perf = list()
me.data = me.dev.large

# Re-estimate 20 times by the same variables X & boostrap Y
for (i in 1:20){
  cat("Sampled model start...", i, "\n")
  n = nrow(me.data$X.block)/20
  start = (1+(i-1)*n)
  end = i*n
  
  Y.block=me.data$Y.vecs[start:end, ]
  Y.vec.new = c(Y.block)
  
  # version 1: Keep right hand side the same, and fluctuate the Response
  fit=var.iter.est(X.big, Y.vec.new, K, 
                   seq(0, 1/K, length.out = 100), 
                   seq(0.001, 1/K, length.out = 100), 
                   T, 'BIC2', 'eigen',
                   10^-2, 10^-1, 10, NA, T, 
                   list(nameleft, nameright), 'standard', 1)
  l1s[[i]] = fit$l1s[['BIC2']]
  
  thetas20[[i]]=theta.fn(nameleft, nameright, fit$beta.est, K)
  eval.train[[i]]=evaluation(thetas20[[i]], Y.block, X.block, NULL, NULL, 0)
  eval.test[[i]]=evaluation(thetas20[[i]], Y.block.test, X.block.test, NULL, NULL, 0)
  perf[[i]]=performance(theta.star, thetas20[[i]], 30)
  
  cat("Train MAE:", eval.train[[i]]$mae, 
      "Test MAE:", eval.test[[i]]$mae,
      "TPR:", perf[[i]]$TPR.en, 
      "TNR:", perf[[i]]$TNR.en, 
      "Sparsity:", perf[[i]]$Sparsity)
  cat("\nSampled model finished...", i, "\n\n")
}

# Perform t test to obtain the less search non-zero theta
# 1.325(80%); 1.725(90%); 2.086(95%); 2.845(99%); n=21, df 20
thetas20[[21]]=theta
tvalue=c(1.325, 1.725, 2.086, 2.845)
ttest.mat = ttest.filter(thetas20, theta)
#ttest.mat = ttest.filter(thetas20, matrix(0, K, M))

search.pos=which(c(t(theta))!=0)
sum(ttest.mat>tvalue[1])
sum(ttest.mat>tvalue[2])
sum(ttest.mat>tvalue[3])
sum(ttest.mat>tvalue[4])

# Predefine searching for the 90% condidence
search.pos=which(c(t(ttest.mat))>tvalue[3])
theta.sd=matrices.fn(thetas20, sd)
theta.sd1=c(t(theta.sd))[search.pos]
theta.new=c(t(theta))
theta.new1=theta.new[search.pos]

p=list()
for(i in 1:100){
  p=c(p, list(c(0.5, 0.25, 0.25)))
}

sum(theta!=0)-sum(c(t(theta))[1:length(theta) %nin% search.pos]!=0)

#### (3) Optimize ceoptim ####
########################
# Optimize obj1. fun.mse.warm.mixed
# optimize obj2. fun.mse.warm.continous
# optimize obj3. fun.mse.warm.discrete
########################
fit.ceoptim=
  CEoptim(f=test, 
          f.arg = list(Y=Y.block, X=X.block, 
                       l1=0, l0=1, k.l0=700, 
                       Theta1.pos=search.pos, 
                       Theta.fixed=theta.new),
          continuous = list(mean=theta.new1, sd=theta.sd1), 
          discrete = list(probs=p),
          maximize = F, verbose = T, N = 10000, 
          iterThr = 100, rho = 0.01, noImproveThr = 5)

fit.ceoptim=
  CEoptim(f=fun.mse.warm.discrete, 
          f.arg = list(Y=Y.block, X=X.block, 
                       l1=0, l0=1, k.l0=700, 
                       Theta1.pos=search.pos, 
                       Theta.fixed=theta.new),
          continuous = list(mean=theta.new1, sd=theta.sd1), 
          discrete = list(probs=p),
          maximize = F, verbose = T, N = 10000, 
          iterThr = 100, rho = 0.01, noImproveThr = 5)
########################
# further optimized result
theta.new=c(t(theta))
theta.new[search.pos] = fit.ceoptim$optimizer$continuous * fit.ceoptim$optimizer$discrete
theta.new[search.pos[fit.ceoptim$optimizer$discrete==2]] = c(t(theta))[search.pos[fit.ceoptim$optimizer$discrete==2]]
theta.new=theta.fn(nameleft, nameright, theta.new, K)
eval.train=evaluation(theta.new, Y.block, X.block, NULL, NULL, 0)
eval.test=evaluation(theta.new, Y.block.test, X.block.test, NULL, NULL, 0)
perf.train=performance(theta.star, theta.new, 30)
eval.train$mae;eval.test$mae;unlist(perf.train)
  
# original iterative estimation results
eval.train=evaluation(theta, Y.block, X.block, NULL, NULL, 0)
eval.test=evaluation(theta, Y.block.test, X.block.test, NULL, NULL, 0)
perf.train=performance(theta.star, theta, 30)
eval.train$mae;eval.test$mae;unlist(perf.train)

# Global goal for the optimization 
theta.goal = c(t(theta.new))
theta.goal[search.pos] = c(t(theta.star))[search.pos]
evaluation(theta.goal, Y.block.test, X.block.test, NULL, NULL, 0)$mae

# plot the result of valued-coef/zeroed-coef
plot.ceoptim(theta.star, theta, theta.new, search.pos)
plot.ceoptim(theta.star, theta, theta.new, search.pos, T)

