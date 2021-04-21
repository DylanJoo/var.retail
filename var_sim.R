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

load(file="data/SIM.data.lag1.RData")
load(file="data/truth/theta.star.RData")
load(file="data/truth/sigma.star.sparse.RData")

#### Prerequisite meta #### 
results=c()
thetas=list()
sigmas=list()

S=10

simulation=function(s, forecast='naive'){
  
  # matrix prepration
  N=nrow(SIM.data[[s]]$Y)/2 #52
  K=ncol(SIM.data[[s]]$Y)   #30
  M=ncol(SIM.data[[s]]$X)
  nameleft=colnames(SIM.data[[1]]$Y)
  nameright=colnames(SIM.data[[1]]$X)
  
  if(forecast=='naive'){
    data=matrixPrep(SIM.data[[s]]$Y, SIM.data[[s]]$X, SIM.data[[s]]$R, 1, N)
    Xtrain=kronecker(diag(K), data$Xtrain)
    Ytrain=c(data$Ytrain)
    
    # Estimation
    fit=var.iter.est(Xtrain, Ytrain, K, 
                     1/seq(1e3, 1/K, length.out = 100), 
                     seq(1e-1, 1/K, length.out = 100)**2, 
                     T, 'BIC1', 'eigen',
                     10^-2, 10^-1, 10, NA, T, 
                     list(nameleft, nameright), 'standard', 1)
    theta=theta.fn(nameleft, nameright, fit$beta.est, K)
    sigma=theta.fn(nameleft, nameleft, fit$sigma.est, K)
    
    # selection/prediction consistent
    perf=performance(theta.star, theta, K)
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
      data=matrixPrep(SIM.data[[s]]$Y, SIM.data[[s]]$X, SIM.data[[s]]$R, i, (i+N-1))
      Xtrain=kronecker(diag(K), data$Xtrain)
      Ytrain=c(data$Ytrain)

      theta=matrix(0, K, M)
      sigma=matrix(0, K, K)
      tpr=rep(0, 2)
      tnr=rep(0, 2)
      maee=rep(0, 2)

      
      # Estimation
      fit=var.iter.est(Xtrain, Ytrain, K, 
                       1/seq(1e3, 1/K, length.out = 100), 
                       seq(1e-1, 1/K, length.out = 100)**2, 
                       T, 'BIC1', 'eigen',
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
  

  cat('##################\nIteration', s, 'Time', i, '\n########################')
  
  return(list(theta=theta/N, 
              sigma=sigma/N,
              tpr=tpr/N,
              tnr=tnr/N,
              maee=maee/N,
              mafe=mafe/N))
}

for(s in 1:S){
  sim=simulation(s, forecast='naive')
  thetas[[s]]=sim$theta
  sigmas[[s]]=sim$sigma
  corrplot(sim$theta, is.corr=F)
  corrplot(sim$sigma, is.corr=F)
  results=cbind(results, c(sim$tpr, sim$tnr, sim$maee, sim$mafe))
}


# Simulation S times with N forecasting
for(s in 1:S){
  theta=rep(0, K*ncol(SIM.data[[s]]$X))
  residual=c()
  ee=c()
  tpr=c()
  tnr=c()
  for(start in 1:N){
    
    #matrix preprocessing
    Ytrain=data$Ytrain
    Xtrain=data$Xtrain
    
    #1-step ahead forecasting
    #Ytest=data$Ynext
    #Xtest=data$Xnext
    
    #Naive forecasting
    Ytest=data$Ytest
    Xteset=data$Xtest
    nameleft=colnames(SIM.data[[1]]$Y)
    nameright=colnames(SIM.data[[1]]$X)
    
    #prerequeisit: mebootstrap setting
    #set.seed(2020)
    #me.dev.large = meboot.oos(SIM.salesG.cate, SIM.exog.cate, 20, "large")
    #meboot.view(X.block, me.dev.large$X.blocks)
    #me.dev.small = meboot.oos(SIM.salesG.cate, SIM.exog.cate, 20, "small")
    #me.dev.similar = meboot.oos(SIM.salesG.cate, SIM.exog.cate, 20, "similar")

    #model training for start to end
    fit=var.iter.est(Xtrain, Ytrain, K, 
                     1/seq(1e3, 1/K, length.out = 100), 
                     seq(1e-1, 1/K, length.out = 100)**2, 
                     T, 'BIC1', 'eigen',
                     10^-2, 10^-1, 10, NA, F, 
                     list(nameleft, nameright), 'standard', 1)
    
    #model evaluation
    #forecast
    residual=rbind(residual, abs(Ytest-fit$beta.est %*% as.matrix(Xtest)))
    ee=rbind(ee, abs(c(t(theta.star))-fit$beta.est))
    theta=theta+fit$beta.est
    
    #consistent
    perf=performance(theta.star, fit$beta.est, K)
    tpr=rbind(tpr, c(perf$TPR.en, perf$TPR.ex))
    tnr=rbind(tnr, c(perf$TNR.en, perf$TNR.ex))
    cat('Iteration', s, 'Time', (start+N), '\n')
    
  }
  result=cbind(result, c('mafe'=mean(apply(residual, 2, mean)), 
                         'maee'=mean(apply(ee, 2, mean)),
                         'tpr.en'=apply(tpr, 2, mean)[1],
                         'tpr.ex'=apply(tpr, 2, mean)[2],
                         'tnr.en'=apply(tnr, 2, mean)[1],
                         'tnr.ex'=apply(tnr, 2, mean)[2]))
  thetas[[s]] = theta.fn(nameleft, nameright, theta/N, K)
}



#### (1) Standard estimation on simulation result  #### 
# Utilize BIC2 for the sparser result 
# Utilize MEBOOT valdiation error for the finest denser result.
Y.block=SIM.data[[1]]$Y[, ]
Y.vec=c(Y.block)
X.big=kronecker(diag(K), SIM.data[[1]]$X)

Y.vec=
X.block.test=cbind(SIM.salesG.cate.test[2:N, ], SIM.exog.cate.test[3:Times, ])
Y.block.test=SIM.salesG.cate.test[3:Times, ]

fit=var.iter.est(X.big, Y.vec, K, 
                 1/seq(1e3, 1/K, length.out = 100), 
                 seq(1e-1, 1/K, length.out = 100)**2, 
                 T, 'BIC1', 'eigen',
                 10^-2, 10^-1, 10, NA, T, 
                 list(nameleft, nameright), 'standard', 1)
theta=theta.fn(nameleft, nameright, fit$beta.est, K)

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
for(i in 1:length(search.pos)){
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

