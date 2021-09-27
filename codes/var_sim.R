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

simulation=function(s, forecast='naive', DATA=SIM.data, estimate='igls', warm=F, thetas.igls=NA){
  
  # matrix prepration
  N=( nrow(DATA[[s]]$Y)-DATA[[s]]$p)/2 # Number of observation
  k.left=ncol(DATA[[s]]$Y)   # Number of categories
  k.right=ncol(DATA[[s]]$Y) + ifelse(is.null(ncol(DATA[[s]]$X)), 0, ncol(DATA[[s]]$X))
  # Number of categories(endogenous variable) + Number of Exogenous variables
  nameleft=colnames(DATA[[s]]$Y)
  nameright=c(rep(nameleft, DATA[[s]]$p), colnames(DATA[[s]]$X))
  
  # step 0: prepare the dataset and spcecify the forecast type

  if(forecast=='naive'){
    # Naive forecast: Train with first half of samples, Test on the other half
    data=matrixPrep(DATA[[s]]$Y, DATA[[s]]$X, DATA[[s]]$R, start=1, end=N, lag=DATA[[s]]$p)
    Xtrain=kronecker(diag(k.left), data$Xtrain)
    Ytrain=c(data$Ytrain)
    
    # Estimation with IGLS
    if(estimate=='igls'){
      # step 1: The naive IGLS estimates.
      fit=var.iter.est(Xtrain, Ytrain, k.left, 
                       1/seq(1e2, k.left, length.out = 100), 
                       1/seq(1e1, k.left, length.out = 100)**2, 
                       T, 'BIC1', c(T, T), 'eigen',
                       10^-2, 10^-1, 10, NA, T, 
                       list(nameleft, nameright), 'standard', 1)
      
      theta=theta.fn(nameleft, nameright, fit$beta.est, k.left)
      sigma=theta.fn(nameleft, nameleft, fit$sigma.est, k.left)
    
    } else if (estimate=='ceoptim'){
      # CEoptim for now, need to specify the IGLS's results
      # step 2a: Specify the IGLS's results (based on same data) 
      #          and pick the considered elements (need further optimization)
      # step 2b: Calculate the empirical probabilities as initial prob (warm=T), 
      #          or default p=0.5. 
      # step 2c: Define the paramenter (objective, init value, optimization setting)
      #          and start the optimization process
      
      if(length(thetas) >= s){
        
        # step 2a 
        theta.igls=c(t(thetas.igls[[s]])) 
        # Selected the coefficients which is non-zero
        # And only optimize based on the picked element
        pick=which(theta.igls!=0)
        
        # step 2b
        if (warm == F |length(pick)==0){
          # 50-50 prob for intialized bernoulli dist.
          pick=1:length(theta.igls)
          prob=rep(list(c(0.5, 0.5)), length(pick))
        } else {
          # Emrpical probabilities for each coefficients
          prob=get_prob(theta.igls[pick])
        }
        
        # step 2c
        # Estimation with CEoptim (detail refer to the CEoptim repo)
        # f.arg: The argument for f
        ## warm: The selected element (to be optimized), record on pick
        ## Theta.init: The initial theta (IGLS's)
        # continuous: gaussian's mu and sigma (of each theta element)
        # discrete: Bernoulli's p (of each theta element)
        fit=CEoptim(f=mse,
                    f.arg=list(Y=data$Ytrain, X=data$Xtrain, warm=pick,
                               Theta.init=theta.igls),
                    continuous=list(mean=theta.igls[pick], sd=rep(0.1, length(pick))),
                    discrete = list(probs=prob),
                    maximize=F, verbose=F, N=10000, iterThr=100, rho=0.5, noImproveThr=5)
        
        theta=mse(fit$optimizer$continuous, fit$optimizer$discrete, warm=pick,
                  Theta.init=theta.igls, ret=T)

      }
      theta=theta.fn(nameleft, nameright, theta, k.left)
      sigma=theta.fn(nameleft, nameleft, cov(data$Ytrain-data$Xtrain%*%t(theta)), k.left)
    }
    
    
    # Theta estimation performances
    perf=performance(theta.star[1:nrow(theta), 1:ncol(theta)], theta, k.left)
    tpr=c(perf$TPR.en, perf$TPR.ex)
    tnr=c(perf$TNR.en, perf$TNR.ex)
    maee=perf$MAEE
    eval.train=evaluation(theta, data$Ytrain, data$Xtrain, NULL, NULL)
    eval.test=evaluation(theta, data$Ytest, data$Xtest, NULL, NULL)
    mafe=c(eval.train$mae, eval.test$mae)
    
    N=1
  } else if(forecast=='rolling'){
    # [Deprecated] The one-step ahead forecasting. 
    # Each simulation run, train on 1:N samples and test on N+1, 
    # then train on 2:N+1 and test on N+2, till the end (2N)
    # Final result take the average of N testing forecast error.
    
    for(i in 1:N){
      data=matrixPrep(DATA[[s]]$Y, DATA[[s]]$X, DATA[[s]]$R, i, (i+N-1))
      Xtrain=kronecker(diag(k.left), data$Xtrain)
      Ytrain=c(data$Ytrain)

      theta=matrix(0, k.left, M)
      sigma=matrix(0, k.left, k.left)
      tpr=rep(0, 2)
      tnr=rep(0, 2)
      maee=rep(0, 2)

      # Estimation
      fit=var.iter.est(Xtrain, Ytrain, k.left, 
                       seq(1e-3, 1/K, length.out = 100), 
                       seq(1e-3, 1/K, length.out = 100), 
                       T, 'BIC2', T, 'eigen',
                       10^-2, 10^-1, 10, NA, F, 
                       list(nameleft, nameright), 'standard', 1)
      
      theta.temp=theta.fn(nameleft, nameright, fit$beta.est, k.left)
      sigma.temp=theta.fn(nameleft, nameleft, fit$sigma.est, k.left)
      theta=theta+theta.temp
      sigma=sigma+sigma.temp
    
      # selection/prediction consistent
      perf=performance(theta.star, theta.temp, k.left)
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

# Simulation Experiments
S=1 # Number of simulation runs
q=10 # Number of category 
sim_type='en' # Type of simulation
load(file=sprintf("data/simulation/sim-%s-%d/data.RData", sim_type, q)) # SIM.data 
load(file=sprintf("data/simulation/sim-%s-%d/theta.star.RData", sim_type , q)) # theta.star

results=matrix(NA, nrow=7, ncol=S)
row.names(results) = c("TPR(en)", "TPR(ex)", "TNR(en)", "TNR(ex)", "MAEE", "MAFE(train)", "MAFE(test)")
thetas=list()
sigmas=list()

results.ce=matrix(NA, nrow=7, ncol=S)
row.names(results.ce) = c("TPR(en)", "TPR(ex)", "TNR(en)", "TNR(ex)", "MAEE", "MAFE(train)", "MAFE(test)")
thetas.ce=list()
sigmas.ce=list()

for(s in 1:S){
  cat('##################\nSimulation', s, '\n##################\n')
  cat('===IGLS===\n')
  go=Sys.time()
  sim=simulation(s, forecast='naive', SIM.data)
  thetas[[s]]=sim$theta
  sigmas[[s]]=sim$sigma
  ret=c(sim$tpr, sim$tnr, sim$maee, sim$mafe)
  results[, s]=ret
  print(Sys.time()-go)
  
  cat('===CEOPTIM===\n')
  go=Sys.time()
  sim.ce=simulation(s, forecast='naive', SIM.data, estimate='ceoptim', warm=T, thetas.igls=thetas)
  thetas.ce[[s]]=sim.ce$theta
  sigmas.ce[[s]]=sim.ce$sigma
  ret.ce=c(sim.ce$tpr, sim.ce$tnr, sim.ce$maee, sim.ce$mafe)
  results.ce[, s]=ret.ce
  print(Sys.time()-go)
  
  cat('===RESULTS===\n') 
  print(results[, s])
  print(results.ce[, s])
  # print(ret.ce.fs)
}

save(thetas, file=sprintf("data/simulation/sim-%s-%d/thetas.RData", sim_type, q))
save(sigmas, file=sprintf("data/simulation/sim-%s-%d/sigmas.RData", sim_type, q))
save(results, file=sprintf("data/simulation/sim-%s-%d/results.RData", sim_type, q))
save(thetas.ce, file=sprintf("data/simulation/sim-%s-%d/thetas.ce.RData", sim_type, q))
save(sigmas.ce, file=sprintf("data/simulation/sim-%s-%d/sigmas.ce.RData", sim_type, q))
save(results.ce, file=sprintf("data/simulation/sim-%s-%d/results.ce.RData", sim_type, q))
