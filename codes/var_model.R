library(glmnet)
library(grplasso)
library(CEoptim)
#######################
#### main function ####
#######################
# /// PARAMETER DESCRIPTION ///
# X.big: Right-hand-sided X input. e.g. SalesG(t-1) & Z(t)
# Y.vec: Left-hand-sided Y input e.g. SalesG(t)
# l1.seq: Lambda strength range LASSO penalty for Theta matrix estimate.
# l2.seq: Lambda strength range LASSO penalty for Inverse covariance matrix estimate.
# std1: (default=F) Select the Theta that within one sigma.
# criteria: (default='BIC2') Lambda selection criteria.
# decomposition: (defualt='eigen') Decomposition method for covaraince matrix. 
# conv.beta: Theta convergence criteria.
# conv.sigma: Sigma convergence criteria.
# max.iter: Maximum iteration.
# verbose: (default=T) Show the temporary result of Theta.
# me.object: Bootstrap sample inputX and inputY
# names: Variables identification for weighted use.
# weighted: (default='standard) Standard/ownfree/ownlight to control the penalty structure of Theta.
# alpha: (default=1) control the LASSO and RIDGE penalty structure.

var.iter.est=function(X.big, Y.vec, k,
                      l1.seq, l2.seq,
                      std1=F, criteria='BIC1', scaling=c(T, T),
                      decomposition='eigen',
                      conv.beta=10^-2, conv.sigma=10^-1,
                      max.iter=10, me.object=NA, verbose=F, 
                      names=NA, weighted='standard', alpha=1){
  # Algorithm:
  
  # step(1): Estimate the theta & select lambda1 @ iteration 1
  # step(2): Calculate the sample covariance matrix 
  # step(3): Estimate the inverse covariance matrix & select lambda2 @ iteration 1
  
  # LOOP-step(4): Decompose the estimated inverse covariance matrix into PP' (P matrix: linear projection)
  # LOOP-step(5): Transform the original X & Y into X~=PX & Y~=PY
  
  # LOOP-step(6): Estimate theta & select lambda1  @ iteration t+1
  # LOOP-step(7): Calculate the sample covariance matrix @ iteration t+1
  # LOOP-step(8): Estimate the inverse covariance matrix & select lambda1 @ iteration t+1
  
  # check
  n=length(Y.vec)/k # number of observation
  m=length(X.big)/k # number of kinds of marketing mix
  iter = 1
  diff.beta=diff.sigma=diff.omega=Inf
  sigma.best = diag(k)
  
  if(anyNA(me.object)){
    me.object=list()
    me.object$X.block=matrix(0, nrow=nrow(X.big)/k, ncol=ncol(X.big)/k)
    me.object$Y.vec=matrix(0, nrow=nrow(X.big)/k, k)
  }

  w=matrix(1, nrow=k, ncol=ncol(X.big)/k)
  # making the weighted matrix that match the x variables
  if(weighted=='ownfree'){
    for(i in 1:k){
      name=sapply(strsplit(names[[1]], '\\.'), function(x) x[2])[i]
      pick=grep(name, names[[2]])
      w[i, ][pick]=0 # set the correponding variable's penalty w to zero
    }
  } else if(weighted=='ownlight'){
    for(i in 1:k){
      name=sapply(strsplit(names[[1]], '\\.'), function(x) x[2])[i]
      pick=grep(name, names[[2]])
      w[i, ][pick]=1/length(pick) # set the correponding variable's penalty w to lighter weight
    }
  }
  
  w=c(t(w))
  
  l1.table=rep(NA, max.iter)


  # step(1)
  fit.beta<-glmnet(x = X.big, y = Y.vec, family = 'gaussian', alpha = alpha,
                   lambda =l1.seq, standardize = T, intercept = T, 
                   penalty.factor = w)
  
  fit.beta$beta = fit.beta$beta[, length(l1.seq):1]

  # Theta selection
  if(criteria=='MEBOOT'){
    l1.choice = var_meboot(me.object$Y.vec, me.object$X.block, fit.beta$beta, k)
    l1.opt = std1.opt(l1.choice, 'right', 1)
    l1.table[iter] = l1.opt
  }else{ # BIC1 or BIC2 or BIC3
    l1.choice = var_theta(Y.vec, X.big, fit.beta$beta, k, criteria, scaling[1])
    l1.opt = which.min(l1.choice)
    l1.table[iter] = l1.opt
  }

  
  # step(2)
  beta.first = beta.init = beta.best = fit.beta$beta[, l1.opt]
  residual.matrix = matrix(Y.vec - X.big %*% beta.best, ncol=k)
  sigma.sample <- (1/n) * t(residual.matrix) %*% (residual.matrix)
  
  # step(3)
  l2.choice=var_omega(sigma.sample, l2.seq, n, scaling[2])
  l2.opt = which.min(l2.choice)
  fit.omega = glasso(s=sigma.sample, nobs=n, rho=l2.seq[l2.opt], penalize.diagonal=F)
  sigma.first = sigma.init = sigma.best = fit.omega$w
  omega.init = omega.best = fit.omega$wi
  
  
  if(verbose){
    cat("iter:", iter, 
        " l1.opt:", l1.seq[l1.opt],
        " l2.opt:", l2.seq[l2.opt], 
        " sparsity:", sum(beta.best!=0)/length(beta.best), 
        " diff.beta:", diff.beta, "\n")
  }
  
  # LOOP-step(4)~LOOP-step(8)
  while((iter<max.iter) & ((diff.beta>conv.beta)|(diff.sigma>conv.sigma))){
    iter = iter + 1
    
    # step(4)
    omega.best <- omega.best
    if(decomposition=='eigen'){
      decomp = eigen(omega.best, symmetric=F)
      C = decomp$vectors
      lambda <- decomp$values
      P = t(C%*% sqrt(diag(lambda)))
    } else if(decomposition=='chol'){
      P = chol(omega.best)
    }
    
    # step(5)
    X.big.tilde <- kronecker(P, diag(n)) %*% X.big
    Y.vec.tilde <- kronecker(P, diag(n)) %*% Y.vec
    
    # step(6)
    fit.beta<-glmnet(x = X.big.tilde, y = Y.vec.tilde, family = 'gaussian', alpha = alpha,
                     lambda =l1.seq, standardize = T, intercept = T, 
                     penalty.factor = w)
    
    fit.beta$beta = fit.beta$beta[, length(l1.seq):1]
    
    # Theta selection
    if(criteria=='MEBOOT'){
      l1.choice = var_meboot(me.object$Y.vec, me.object$X.block, fit.beta$beta, k)
      l1.opt = std1.opt(l1.choice, 'right', 1)
      l1.table[iter] = l1.opt
    }else{ # BIC 1/2/3
      l1.choice = var_theta(Y.vec.tilde, X.big.tilde, fit.beta$beta, k, criteria, scaling[1])
      l1.opt = which.min(l1.choice)
      l1.table[iter] = l1.table[iter]
    }
    
    # step(7)
    beta.best = fit.beta$beta[, l1.opt]
    residual.matrix = matrix(Y.vec - X.big %*% beta.best, ncol=k)
    sigma.sample = (1/n) * t(residual.matrix) %*% (residual.matrix)
    
    # step(8)
    l2.choice=var_omega(sigma.sample, l2.seq, n, scaling[2])
    l2.opt = which.min(l2.choice)
    fit.omega = glasso(s=sigma.sample, nobs=n, rho=l2.seq[l2.opt], penalize.diagonal=F)
    
    sigma.best = fit.omega$w
    omega.best = fit.omega$wi
    
    diff.beta = max(beta.best - beta.init)
    diff.sigma = max(sigma.best - sigma.init)
    
    if(verbose){
      cat("iter:", iter, 
          " l1.opt:", l1.seq[l1.opt],
          " l2.opt:", l2.seq[l2.opt], 
          " sparsity:", sum(beta.best!=0)/length(beta.best), 
          " diff.beta:", diff.beta, "\n")
    }
    beta.init = beta.best
    sigma.init = sigma.best
  }
  out <- list(beta.est=beta.best, beta.first=beta.first,
              sigma.est=sigma.best, sigma.first=sigma.first, 
              sigma.sample=sigma.sample, 
              l1s=l1.table)
}


###################
#### Utilities ####
###################
var_theta=function(Y, X, thetas, k, method, scaling=F){
  # Calculate the BIC value of each theta
  # method: indicates different BIC implementation method
  # Scaling: Standardize the "complexity" and "fitness", to obtain the balanced BIC value
  complexity=rep(0, ncol(thetas))
  fitness=rep(0, ncol(thetas))
  df=c()
  ll=c()
  n=length(Y)/k
  
  ## Bayesian Information Criteria
  if(method=='BIC1'){
    for(i in 1:ncol(thetas)){
      df[i]=sum(thetas[, i]!=0)
      ll[i]= (-1/2) * (n*k*log(2*pi)) + (-1/2) * sum((Y - X %*% as.matrix(thetas[, i]))**2)
    }
  } else if (method=='BIC2'){
    for(i in 1:ncol(thetas)){
      df[i]=sum(thetas[, i]!=0)
      ll[i]= (-1/2) * (n*k) * log( sum((Y - X %*% as.matrix(thetas[, i]))**2) / n*k)
    }
  } else if(method=='BIC3'){
    for(i in 1:ncol(thetas)){
      df[i]=sum(thetas[, i]!=0)
      ll[i]=sum((Y - X %*% as.matrix(thetas[, i]))**2)
    }
  }
  complexity=df*log(n)
  fitness=(-2)*ll
  
  if(scaling){
    complexity=scale(complexity)
    fitness=scale(fitness)
  }
  #plot.ts(complexity + fitness, ylab='theta')
  return(complexity + fitness)
}
var_meboot=function(Y.boot, X.boot, thetas, k){
  # Calculate the rss of each theta
  # [TODO]: Randonmly selected k out-of-samples
  rss = c()
  df = c()
  kl = c()
  for(i in 1:ncol(thetas)){
    B.mat = matrix(thetas[, i], nrow=ncol(X.boot), ncol=ncol(Y.boot), byrow=F)
    resid = as.matrix(Y.boot) - as.matrix(X.boot) %*% B.mat
    rss[i] = sum(resid**2)
  }
  return(rss)
}
var_omega=function(S, l2s, n, scale=F){
  # by definition of BIC | Beta
  bic=c()
  df=c()
  ll=c()
  k=dim(S)[1]
  for(i in 1:length(l2s)){
    fit=glasso(s=S, nobs=n, rho=l2s[i], penalize.diagonal=F)
    df[i]=sum(fit$wi[upper.tri(diag(k))] != 0)
    ll[i] = -(n/2)*log(det(fit$w)) - sum(diag((S*n) %*% fit$wi))
  }
  if(scale){
    term1=scale((-2)*ll)
    term2=scale(df*log(n))
    bic = ifelse(is.nan(term1), rep(0, length(term1)), term1) + 
          ifelse(is.nan(term2), rep(0, length(term2)), term2)
  } else {
    bic = (-2)*ll + df * log(n)
  }
  #plot.ts(bic, ylab='Omega selection')
  return(bic)
}
std1.opt=function(vectors, sparse.side='right', nsd=1){
  # For var_meboot cross-validation, 
  # Select the sparser theta based on validation error (one-standard-deviation)
  std = sd(vectors)
  min.val = min(vectors)
  min.idx = which.min(vectors)

  if(sparse.side == 'left'){
    vectors[(min.idx+1):length(vectors)] = -Inf
  } else if (sparse.side == 'right'){
    vectors[1:(min.idx-1)] = -Inf
  }
  vectors[which(vectors > (min.val + std*nsd))] = -Inf
  smooth.opt = which.max(vectors)
  
  return(smooth.opt)
}
