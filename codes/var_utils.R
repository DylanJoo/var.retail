require(meboot)
require(pls)
library(xlsx)

#### Preprocessing ####
tsTransform=function(series, smooth=F, measurement, offset){
  # Tranform the type of time-series, 
  # series: The times series to be transformed.
  # measurement: Selected a kind of measurement (e.g. logdiff)
  # smooth: Smoothing the zero value (substitute by the minimum)
  # offset: The results need to be removed.
  if(smooth){series[which(series==0)]=min(series[series!=0], 1)}
  series.logdiff = append(0, diff(log(series)))
  series.diff = append(0, diff(series))
  series.vol = append(0, (series[-1] - series[-length(series)]) / series[-length(series)])
  series.bin = as.numeric(series > 0)
  series.score = (series - max(Mode(series)))/max(Mode(series))
  series.full = Reduce(cbind, list(series, series.diff, series.logdiff, series.vol, series.bin, series.score))
  # The growth/volitility cannot be calculated due to truncated start time.
  colnames(series.full) <- c('raw', 'diff', 'log.diff', 'vol', 'bin', 'score')
  
  return(series.full[setdiff(1:nrow(series.full), offset), measurement])
}
tsFixed=function(table){
  # For cleaning the data, uses in the var_preprocess.R
  # Time series debugging function, Remove time series wo/ variation
  zerosd= which(apply(table, 2, sd, na.rm=T)==0)
  nafixed = which(apply(is.na(table), 2, sum) <= nrow(table)/2)
  nadropped = which(apply(is.na(table), 2, sum) > nrow(table)/2)
  for(i in nafixed){
    for(j in which(is.na(table[, i]))){
      table[, i][j] = min(table[, i], na.rm=T) 
    }
  }
  return(table[, setdiff(1:ncol(table), c(nadropped, zerosd))])
}
matrixPrep=function(Y.full, X.full, R.full, start, end, lag=1){
  # Split the train/test set of data
  # [TODO]: The "lag" argument for more than one scenario /line 37 & 40
  K=ncol(Y.full)
  Ytrain=Y.full[lag+(start:end), ]
  Xtrain=cbind( Y.full[(start:end), ], X.full[lag+(start:end), ] )
  
  Ytest=NA
  Xtest=NA
  Ynext=NA
  Xnext=NA
  Resid=NA
  tryCatch({
    Ytest=Y.full[end+(start:end), ]
    Xtest=cbind( Y.full[(end-lag)+(start:end), ], X.full[end+(start:end), ] )
    
    Ynext=Y.full[end+1, ]
    Xnext=X.full[end+1, ]
    Resid=R.full[lag+(start:end), ]
  }, 
  error=function(e){
    print("Invalid testing set and residual.")
  }
  )
  
  
  return(list(
    Ytrain=Ytrain,# 1:N
    Xtrain=Xtrain,
    Ytest=Ytest,  # N+1:2N
    Xtest=Xtest,
    Ynext=Ynext, # Made for "rolling" estimation
    Xnext=Xnext, # Made for "rolling" estimation
    R=Resid
  ))
}

#### Sampling boostrap ####
meboot.oos.raw=function(data, ex, replicates=20, en.std=F, fn=function(x) x, order=1){
  # Bootstrap from raw Dollars-sales data (IRI's), 
  # generate the Dollars sample and then the salesg data
  # data: raw dollars data from IRI
  # ex: preprocessed exogenous data
  Times=nrow(data)
  N=Times-1
  n=N-order
  responses.sales = matrix(NA, nrow=replicates*n, ncol=ncol(data))
  responses.dollars = matrix(NA, nrow=replicates*n, ncol=ncol(data))
  samples.sales = matrix(NA, nrow=replicates*n, ncol=ncol(data))
  samples.dollars = matrix(NA, nrow=replicates*n, ncol=ncol(data))
  
  for(i in 1:ncol(data)){
    # search the largest variation without negative value.
    for(strech in seq(100, 0, -10)){
      dollars=meboot(data[, i], reps=replicates, trim=list(0.1), 
                     scl.adjustment=F, force.clt=F, expand.sd=T, sym=F, 
                     reachbnd=T, fiv=strech)$ensemble
      if(min(dollars)>0){
        break # break if variaion is controlled
      } else {
        data.sd = sd(data[, i])
        data.mu = mean(data[, i])
        z = (min(data[, i])-data.mu)/data.sd
        xadj = list(0.1, xmin=min(data[, i]), xmax=max(data[, i], data.mu+data.sd*z))
        dollars=meboot(data[, i], reps=replicates, trim=xadj, 
                       scl.adjustment=F, force.clt=F, expand.sd=F, sym=F, 
                       reachbnd=F)$ensemble
      }
    }
    
    sales=apply(dollars, 2, fn)
    if(en.std){sales=stdize(sales)}
    
    # Sales sample for counting validation error
    # [TODO] make the larger order to fit the right time. (samples.sales & samples.dollars)
    samples.sales[, i]=c(sales[2:N, ])
    responses.sales[, i] = c(sales[3:Times, ])
    samples.dollars[, i]=c(dollars[2:N, ])
    responses.dollars[, i]=c(dollars[3:Times, ])
  }
  samples.sales = cbind(samples.sales, kronecker(rep(1, replicates), ex[3:Times, ]))
  colnames(samples.sales) <- c(colnames(data), colnames(ex))
  colnames(samples.dollars) <- c(colnames(data))
  colnames(responses.dollars) <- colnames(responses.sales)<- c(colnames(data))

  out <- list(X.block=samples.sales, Y.vec=responses.sales, 
              X.dollars=samples.dollars, Y.dollars=responses.dollars)
}
meboot.oos=function(en, ex, replicates, mode='large', order=1){
  # Bootstrap from sales growth data (simulated), 
  # generate the sales growth samples for cross-validation
  # mode: control the variation of time-series bootstrap trend
  Times=nrow(en)
  N=Times-1
  n=N-order
  responses = matrix(NA, nrow=replicates*n, ncol=ncol(en))
  samples = matrix(NA, nrow=replicates*n, ncol=ncol(en))
  for(i in 1:ncol(en)){
    if(mode=='large'){
      en.new=meboot(en[-1, i], reps=replicates, trim=list(0.1), scl.adjustment=F, force.clt=F, expand.sd=T, sym=F, reachbnd=T, fiv=100)$ensemble
    } else if(mode=='small'){
      en.new=meboot(en[-1, i], reps=replicates, trim=list(0.1), scl.adjustment=F, force.clt=F, expand.sd=T, sym=F, reachbnd=T, fiv=50)$ensemble
    } else if(mode=='similar'){
      en.new=meboot(en[-1, i], reps=replicates, trim=list(0.1), scl.adjustment=F, force.clt=F, expand.sd=F, sym=F, reachbnd=F)$ensemble
    }
    # [TODO] make the larger order to fit the right time. (samples)
    samples[, i] = c(en.new[-N, ])
    responses[, i] = c(en.new[-1, ])
  }
  # bind exogenous
  samples = cbind(samples, kronecker(rep(1, replicates), ex[-(1:(1+order)), ]))
  colnames(samples) <- c(colnames(en), colnames(ex))
  colnames(responses) <- c(colnames(en))
  
  out <- list(X.blocks=samples, Y.vecs=responses)
}

#### Theta visualizing ####
effect.interval=function(thetas, thres=30){
  # Visualize the thetas range of estimated points (standard/ownfree/ownlight)
  # thetas: List of theta results
  catename = sapply(strsplit(rownames(thetas[[1]]), '\\.'), function(x) x[2])
  unpick = grep('SALES', colnames(thetas[[1]]))
  for(i in 1:nrow(thetas[[1]])){
    a = thetas[[1]][i, -(unpick)]
    b = thetas[[2]][i, -(unpick)]
    c = thetas[[3]][i, -(unpick)]
    nonzero = which(a+b+c!=0)
    if(length(nonzero)>thres){nonzero=nonzero[1:thres]}
    mmax = max(a, b, c)
    mmin = min(a, b, c)
    # Visualization
    par(mar = c(15, 5, 2, 2))
    plot.ts(a[nonzero], col=1, type='b', xlab=NA, ylab="estimate", cex=0.8, xaxt='n', 
            ylim=c(mmin, mmax), pch=1, lty=1, main=paste('SalesG', catename[i]))
    lines(b[nonzero], col=2, ylim=c(min, max), type='b', lty=1, cex=0.8)
    lines(c[nonzero], col=4, ylim=c(min, max), type='b', lty=1, cex=0.8)
    legend('topright', c("Standard","Ownlight", 'Ownfree'), cex=1, col=c(1,2,4), lty=1, pch=c(1,1,1), bg='white', pt.cex=0.8)
    axis(1, at=1:length(nonzero), labels=colnames(thetas[[1]])[-(unpick)][nonzero], las=2)
    
  }
}
plot.predict=function(evalobj, re=T){
  # Visualize the line of truth/predicted, includes Dollars sales and sales growth
  # evalobj: The output object of evaluation(), includes the detail evaluation
  catename = sapply(strsplit(colnames(evalobj$y.re), '\\.'), function(x) x[2])
  for(i in 1:nrow(thetas[[1]])){
    if(re){
      a = evalobj$y.re[, i]
      b = evalobj$f.value.re[, i]
      metric='Reconstructed Dollars'
      main.title=paste('Dollars', catename[i])
      plot.ts(a, col=1, ylab=metric, cex=0.8, main=main.title)
      lines(b, col=2, lty=2)
    } else {
      a = evalobj$y[, i]
      b = evalobj$f.value[, i]
      c = evalobj$y
      metric='Sales Growth'
      main.title=paste('SalesG', catename[i])
      plot.ts(a, col=1, ylab=metric, cex=0.8, main=main.title)
      lines(b, col=2, lty=2)
      abline(h=0)
    }
  }
}
plot.ceoptim=function(thetaA, thetaB, thetaC, search.pos, zero=F){
  # Plot the CEoptim estimated, IGLS estimated and the truth value.
  a=c(t(thetaA))[search.pos]
  b=c(t(thetaB))[search.pos]
  c=c(t(thetaC))[search.pos]
  plot(a[a!=0], pch=4, ylab='', xlab='COEF!=0')
  points(b[a!=0], col=2, pch=2)
  points(c[a!=0], col=4, pch=2)
  abline(h=0)
  if(zero){
    plot(a[a==0], pch=4, ylab='', xlab='COEF==0')
    points(b[a==0], col=2, pch=2)
    points(c[a==0], col=4, pch=2)
    abline(h=0)
  }
}
meboot.view=function(X, X.me){
  # Plot and Demo the times series with meboot's time series, and the original as well
  time.original=nrow(X)
  time.me=nrow(X.me) / time.original
  plot.ts(X[, 1], xlab='Weeks', ylab='BEER')
  for(r in 1:time.me){
    lines(X.me[(r-1)*time.original+(1:time.original), 1], col=2)
  }
  lines(X[, 1])
}

csv.predict=function(evalobj, name='train', filename, re=T){
  # Generated the report of forecasting and error results.
  catename = sapply(strsplit(colnames(eval.train$y.re), '\\.'), function(x) x[2])
  table = data.frame()
  K=ncol(eval.train$y)
  Times=nrow(eval.train$y)
  
  # mae
  mae = data.frame(apply(abs(evalobj$err), 2, mean))
  colnames(mae) = 'MAE'
  row.names(mae) = catename
  write.csv(mae, filename)

  # mape
  if (re){
    mape = data.frame(apply(abs(evalobj$err.re) / abs(evalobj$y.re), 2, mean))
    colnames(mape) = 'MAPE'
    write.csv(cbind(mae, mape), filename)
  }
  
  # data
  for(c in 1:ncol(evalobj$y)){
    temp = data.frame(cbind(
        rep(catename[c], Times),
        paste('WEEK', c(1:Times)),
        evalobj$y[, c], evalobj$f.value[, c], evalobj$err[, c],
        evalobj$y.re[, c], evalobj$f.value.re[, c], evalobj$err.re[, c]))
    table = rbind(table, temp)
  }
  colnames(table) <- c('Category', 'Time', 'Y', 'Predicted', 'Error', 'Y.Re', 'Predicted.Re', 'Error.re')
  write.csv(table, filename)
  cat("File created.")
}

#### Evaluate effectiveness ####
performance=function(gt, hat, k){
  # Use in simulation experiment resuls
  # Calculate the TPR/TNR/MAEE b/w theta hat and theta star
  # calculate the MAFE b/w truth and prediction
  if (is.vector(hat)){
    hat = matrix(hat, nrow=nrow(gt), ncol=ncol(gt) , byrow=T)
  }
  
  gt.z <- (gt==0)
  hat.z <- (hat==0)
  
  gt.pos <- (gt>0)
  gt.neg <- (gt<0)
  hat.pos <- (hat>0)
  hat.neg <- (hat<0)
  
  gt.pos.en <- gt.pos[, 1:k]
  gt.neg.en <- gt.neg[, 1:k]
  gt.pos.ex <- gt.pos[, -(1:k)]
  gt.neg.ex <- gt.neg[, -(1:k)]
  gt.z.en <- gt.z[, 1:k]
  gt.z.ex <- gt.z[, -(1:k)]
  
  hat.pos.en <- hat.pos[, 1:k]
  hat.neg.en <- hat.neg[, 1:k]
  hat.pos.ex <- hat.pos[, -(1:k)]
  hat.neg.ex <- hat.neg[, -(1:k)]
  hat.z.en <- hat.z[, 1:k]
  hat.z.ex <- hat.z[, -(1:k)]
  
  out<- list(
    TPR=(sum(gt.pos & hat.pos) + sum(gt.neg & hat.neg))/sum(!gt.z), 
    TNR=sum(gt.z & hat.z)/sum(gt.z),
    TPR.en=(sum(gt.pos.en & hat.pos.en) + sum(gt.neg.en & hat.neg.en))/sum(!gt.z.en),
    TPR.ex=(sum(gt.pos.ex & hat.pos.ex) + sum(gt.neg.ex & hat.neg.ex))/sum(!gt.z.ex),
    TNR.en=sum(gt.z.en & hat.z.en)/sum(gt.z.en),
    TNR.ex=sum(gt.z.ex & hat.z.ex)/sum(gt.z.ex),
    MAEE=mean(abs(gt-hat)),
    Accuracy=(sum(gt.pos&hat.pos)+sum(gt.neg&hat.neg)+sum(gt.z&hat.z))/length(gt), # Zero+pos+neg/all
    Sparsity=(sum(!hat.z)/length(hat))
  )
}
evaluation=function(Theta, Y, X, Y.re, X.re, l1=0){
  # Use in IRI empirical results evaluation
  # Calculate the Y contribution of each endogenmous and exogenous variables
  # training 
  f.value = X %*% t(Theta)
  sales.ctb=X[, grep("SALES", colnames(X))] %*% t(Theta[, grep("SALES", colnames(Theta))])
  pr.ctb= X[, grep("PZ.REDUCT", colnames(X))] %*% t(Theta[, grep("PZ.REDUCT", colnames(Theta))])
  ad.ctb= X[, grep("AD", colnames(X))] %*% t(Theta[, grep("AD", colnames(Theta))])
  di.ctb= X[, grep("DI", colnames(X))] %*% t(Theta[, grep("DI", colnames(Theta))])
  
  err = Y - f.value
  f.value.re = X.re * exp(f.value) #f.value.re2 = X.re * (1+f.value)
  err.re = Y.re - f.value.re
  p.err.re = err.re / Y.re
  
  mae = mean(abs(err))
  mse = mean(err ** 2)
  mape = mean(abs(p.err.re))
  obj = (1/2)*mse + l1*sum(abs(Theta))
  
  out<-list(mae=mae, mape=mape, mse=mse, obj=obj,
            y=Y, y.re=Y.re, 
            f.value=f.value, f.value.re=f.value.re, 
            err=err, err.re=err.re,
            sales.ctb=sales.ctb, pr.ctb=pr.ctb, ad.ctb=ad.ctb, di.ctb=di.ctb)
}
theta.fn=function(row, col, theta, k){
  # Vectorize the coefficient matrix.
  if(is.vector(theta)){
    theta = matrix(theta, k, byrow=T)
  }
  colnames(theta) <- col
  row.names(theta) <- row
  return(theta)
}

#### Testing function ####
unitroot.test=function(phi){
  # ADF statistical test for validate the valid theta
  # Solve the unitroot funciton
  f.unitroot=function(x, p){return(det(diag(ncol(p)) - p * x))}
  tryCatch(
    expr={
      r=uniroot(function(x) f.unitroot(x, phi), interval=c(-1, 1))$root
      return(r)
    },
    error=function(e){return(F)}
  )
}
matrices.fn=function(mat.list, f){
  # Apply the function f for all matrices in the list
  n.row=dim(mat.list[[1]])[1]
  n.col=dim(mat.list[[1]])[2]
  mat=matrix(NA, n.row, n.col)
  for(i in 1:n.row){
    for(j in 1:n.col){
      mat[i, j] = f(sapply(mat.list, function(x) x[i, j]))
    }
  }
  return(mat)
}

#### Optimization ####
get_prob=function(coefs){
  # Get the empirical IGLS probabilities
  # Normalized all the probabilities, based on the estimated theta valud + p_base
  ret = list()
  base = 1-max(abs(coefs))
  for(i in 1:length(coefs)){
    true=base+abs(coefs[i])
    ret[[i]] = c(1-true, true)
  }
  return(ret)
}
ttest.filter=function(theta.samples, theta.tgt){
  # Calculate the whethter the coefs are significantly not equal to zero
  # based on N samples (assume N samples gaussian distributed)
  sig.mat=matrix(NA, nrow(theta), ncol(theta))
  for(i in 1:nrow(theta.tgt)){
    for(j in 1:ncol(theta.tgt)){
      data=sapply(theta.samples, function(x) x[i, j])
      sig.mat[i, j] = abs(t.test(data, mu=theta.tgt[i, j])$statistic)
    }
    sig.mat[is.na(sig.mat)] = 0
  }
  return(sig.mat)
}
mse=function(ThetaC, ThetaZ, Y=NA, X=NA, warm=T, Theta.init=NA, ret=F, pen.f=function(x) 0){
  # Objective function to calculate the MSE
  # ThetaC: continuous value of theta
  # ThetaZ: Discrete valud of theta (1/0)
  if(any(warm)){
    Theta=Theta.init
    Theta[warm]=ThetaC*ThetaZ
  }
  if(ret){return(Theta)}
  Theta = matrix(Theta, nrow=ncol(Y), byrow=T)
  mse=(1/2)*mean( (Y-X %*% t(Theta))**2)
  return(mse + pen.f(Theta))
}

#### Optimization #### 
# [Deprecated] Other Objective functions S for CEoptim 
loglik=function(C, Z, Y=NA, X=NA, sep, warm=T, Theta.init=NA, ret=F){
  cormat = diag(ncol(Y))
  cormat[upper.tri(cormat)] = cormat[lower.tri(cormat)] = C[-(1:sep)]*Z[-(1:sep)]
  cormat[cormat>0.5]=0.5
  cormat[cormat<-0.5]=-0.5
  
  if(any(warm)){
    Theta=Theta.init
    Theta[warm]=C[1:sep] * Z[1:sep]
  }
  if(ret){return(Theta)}
  Theta = matrix(Theta, nrow=ncol(Y), byrow=T)
  R = Y-X %*% t(Theta)
  Sigma = diag(diag(cov(R))) %*% cormat %*% diag(diag(cov(R)))
  
  return( R %*% solve(Sigma) %*% t(R) + nrow(Y) * log(det(Sigma)))
}
mse.c=function(ThetaC, Y=NA, X=NA, warm=T, Theta.init=NA, ret=F, pen.f=function(x) 0){
  if(any(warm)){
    Theta=Theta.init
    Theta[warm]=ThetaC
  }
  if(ret){return(Theta)}
  Theta = matrix(Theta, nrow=ncol(Y), byrow=T)
  mse=(1/2)*mean( (Y-X %*% t(Theta))**2 )
  return(mse + pen.f(Theta))
}
mse.z=function(ThetaZ, Y=NA, X=NA, Theta.init=NA, ret=F, pen.f=function(x) 0){
  Theta=Theta.init * ThetaZ
  if(ret){return(Theta)}
  Theta = matrix(Theta, nrow=ncol(Y), byrow=T)
  mse=(1/2)*mean( (Y-X%*%t(Theta))**2 )
  return(mse + pen.f(Theta))
}
l0=function(X, k, lambda=1){
  return( max((sum(X!=0) - k), 0)*lambda )
}
bic=function(C, Y, X, sep, k.dim){
  Theta = C[1:sep]
  Theta = matrix(Theta, nrow=ncol(Y), byrow=T)
  
  Sigma = diag(C[sep:(sep+k.dim-1)], k.dim)
  Sigma[upper.tri(Sigma)] = Sigma[lower.tri(Sigma)] = C[-(1:(sep+k.dim-1))]
  return( nrow(Y) * log(det(Sigma)) + (Y-X %*% t(Theta)) ** 2 + sum(Theta!=0) * log(nrow(Y)))
}