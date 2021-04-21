#### package to be installed ####
install.packages('CEoptim')
install.packages('corrplot')
install.package('DescTools')
install.packages('glasso')
install.packages('glmnet')
install.packages('meboot')
install.packages('grplasso')
install.packages('pls')
install.packages('xlsx')
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

load("data/salesG.cate.RData")
load("data/exog.cate.RData")
load("data/exog.sku.RData")
load("data/dollars.cate.RData")

#### Prerequisite meta #### 
Times=nrow(salesG.cate) # Observations in data
N=Times-1 # lost one obesrvation due to differecing
P=1
n=N-P # lost p observation due to autoregressive setting
K=ncol(salesG.cate)
thetas=list()
thetas_sku=list()

#### Reformulate into big matrix Y & X #### 
Y.block=salesG.cate[3:Times, ]
Y.vec=c(Y.block)
X.block=cbind(salesG.cate[2:N, ], stdize(exog.cate[3:Times, ]))
X.block.sku=cbind(salesG.cate[2:N, ], stdize(exog.sku[3:Times, ]))
X.big=kronecker(diag(K), X.block)
X.big.sku=kronecker(diag(K), X.block.sku)
nameleft=colnames(salesG.cate)
nameright=colnames(X.block)
nameright.sku=colnames(X.block.sku)

#### Prepare bootstrap X & Y  #### 
fn.logdiff=function(x) return(tsTransform(x, T, 'log.diff', 0))
set.seed(2020)
me.dev = meboot.oos.raw(dollars.cate, stdize(exog.cate), 20, F, fn.logdiff, 1)
# demo thet boostrapped result if neeeded.
meboot.view(dollars.cate[-c(1, nrow(dollars.cate)), ], me.dev$X.dollars)
meboot.view(X.block, me.dev$X.block)

me.dev.sku = meboot.oos.raw(dollars.cate, stdize(exog.sku), 20, F, fn.logdiff, 1)
set.seed(2019)
me.test = meboot.oos.raw(dollars.cate, stdize(exog.cate), 1, F, fn.logdiff, 1)
me.test.sku = meboot.oos.raw(dollars.cate, stdize(exog.sku), 1, F, fn.logdiff, 1)

#### Standard estimation  #### 
# For the sparser result:
fit=var.iter.est(X.big, Y.vec, K, 
                 1/seq(1e3, K, length=100), 
                 seq(1e-2, 1/K, length.out = 100)**2, 
                 T, 'BIC1', 'eigen',
                 10^-2, 10^-1, 10, me.dev, T, 
                 list(nameleft, nameright), 'standard', 1)
thetas[[1]]=theta.fn(nameleft, nameright, fit$beta.est, K)

fit=var.iter.est(X.big, Y.vec, K, 
                 1/seq(1e3, K, length=100), 
                 seq(1e-2, 1/K, length.out = 100)**2, 
                 T, 'MEBOOT', 'eigen',
                 10^-2, 10^-1, 10, me.dev, T, 
                 list(nameleft, nameright), 'standard', 1)
thetas[[1]]=theta.fn(nameleft, nameright, fit$beta.est, K)


#### (1) Visualize the theta effct 
#########################
# Figure of coefficient matrix.
#########################
pdf(file="figure/IRI-theta-meboot.pdf")
corrplot(thetas[[1]][, grep("SALES", colnames(thetas[[1]]))], 
         is.corr=F, tl.srt=45, title="Lag Effect", mar=c(0,0,1,0))
corrplot(thetas[[1]][, grep("PZ.REDUCT", colnames(thetas[[1]]))], 
         is.corr=F, tl.srt=45, title="Price Reduction Effect", mar=c(0,0,1,0))
corrplot(thetas[[1]][, grep("AD.RATE", colnames(thetas[[1]]))], 
         is.corr=F, tl.srt=45, title="Advertising Effect", mar=c(0,0,1,0))
corrplot(thetas[[1]][, grep("DI.RATE", colnames(thetas[[1]]))], 
         is.corr=F, tl.srt=45, title="Display Effect", mar=c(0,0,1,0))
corrplot(fit$sigma.est, is.corr=F, title="Error Covariance")
dev.off()

#### (2) Figure out the performance and evaluation with visulization
#########################
# Results of evaluation by MAE/MAPE.
# [2021/04/19 update]
#########################
# ~18% MAE
eval.train=evaluation(thetas[[1]], 
                      Y.block,
                      X.block,
                      dollars.cate[3:Times, ],
                      dollars.cate[2:N, ], 0)
# ~22% MAE
eval.test=evaluation(thetas[[1]], 
                     me.test$Y.vec, 
                     me.test$X.block, 
                     me.test$Y.dollars, 
                     me.test$X.dollars, 0)
eval.train$mae;eval.train$mape
eval.test$mae;eval.test$mape

####################
#[2021/04/19 update]
####################
csv.predict(eval.train, 'train', "csv/IRI-eval-train-meboot.csv")
csv.predict(eval.test, 'test', 'csv/IRI-eval-test-meboot.csv')

# Estimation visualization
pdf(file="figure/IRI-eval-train-meboot.pdf")
plot.predict(eval.train, F)
dev.off()
pdf(file="figure/IRI-eval-trainre-meboot.pdf")
plot.predict(eval.train, T)
dev.off()

pdf(file="figure/IRI-eval-test-meboot.pdf")
plot.predict(eval.test, F)
dev.off()
pdf(file="figure/IRI-eval-testre-meboot.pdf")
plot.predict(eval.test, T)
dev.off()

#### (2) Category theta effect interval: weighted estimation via LASSO #### 
fit.st=var.iter.est(X.big, Y.vec, K, 
                    seq(0, 1/K, length.out = 100), 
                    seq(0.001, 1/K, length.out = 100), 
                    T, 'BIC2', 'eigen',
                    10^-2, 10^-1, 10, me.dev, T, 
                    list(nameleft, nameright), 'standard', 1)
thetas[[1]]=theta.fn(nameleft, nameright, fit.st$beta.est, K)
fit.ol=var.iter.est(X.big, Y.vec, K, 
                    seq(0, 1/K, length.out = 100), 
                    seq(0.001, 1/K, length.out = 100), 
                    T, 'BIC2', 'eigen',
                    10^-2, 10^-1, 10, me.dev, T, 
                    list(nameleft, nameright), 'ownlight', 1)
thetas[[2]]=theta.fn(nameleft, nameright, fit.ol$beta.est, K)
fit.of=var.iter.est(X.big, Y.vec, K, 
                    seq(0, 1/K, length.out = 100), 
                    seq(0.001, 1/K, length.out = 100), 
                    T, 'BIC2', 'eigen',
                    10^-2, 10^-1, 10, me.dev, T, 
                    list(nameleft, nameright), 'ownfree', 1)
thetas[[3]]=theta.fn(nameleft, nameright, fit.of$beta.est, K)

#########################
# Figure of effect interval on category level.
#########################
pdf(file="figure/IRI-effect-interval-lasso-category-BIC2.pdf")
effect.interval(thetas)
dev.off()


#### (3) SKUs theta effect interval: weighted estimation by LASSO #### 
fit=var.iter.est(X.big.sku, Y.vec, K, 
                 1/seq(1e3, K, length.out = 100), 
                 seq(0.001, 1/K, length.out = 100)**2, 
                 T, 'BIC2', 'eigen',
                 10^-2, 10^-1, 10, NA, T, 
                 list(nameleft, nameright.sku), 'standard', 1)
thetas_sku[[1]]=theta.fn(nameleft, nameright.sku, fit$beta.est, K)
fit.ol=var.iter.est(X.big.sku, Y.vec, K, 
                    seq(0, 1/K, length.out = 100), 
                    seq(0.001, 1/K, length.out = 100), 
                    T, 'BIC2', 'eigen',
                    10^-2, 10^-1, 10, NA, T, 
                    list(nameleft, nameright.sku), 'ownlight', 1)
thetas_sku[[2]]=theta.fn(nameleft, nameright.sku, fit.ol$beta.est, K)
fit.of=var.iter.est(X.big.sku, Y.vec, K, 
                    seq(0, 1/K, length.out = 100), 
                    seq(0.001, 1/K, length.out = 100), 
                    T, 'BIC2', 'eigen',
                    10^-2, 10^-1, 10, NA, T, 
                    list(nameleft, nameright.sku), 'ownfree', 1)
thetas_sku[[3]]=theta.fn(nameleft, nameright.sku, fit.of$beta.est, K)

#########################
# Figure of effect interval on SKUs level.
#########################
pdf(file="figure/IRI-effect-interval-lasso-sku-BIC2.pdf")
effect.interval(thetas_sku)
dev.off()
