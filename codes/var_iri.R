#### package to be installed ####
install.packages('CEoptim')
install.packages('corrplot')
install.package('DescTools')
install.packages('glasso')
install.packages('glmnet')
install.packages('meboot')
install.packages('grplasso')
install.packages('pls')
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
load("data/exog.cate.old.RData")
load("data/exog.sku.RData")
load("data/dollars.cate.RData")

#### Prerequisite meta #### 
Times=nrow(salesG.cate) # Observations in data
N=Times-1 # lost one obesrvation due to differecing
P=1 # Selected Lag order
n=N-P # lost p observation due to autoregressive setting
K.left=ncol(salesG.cate) # Number of category

# Estimated thetas matrix
thetas=list() 
thetas_sku=list()

#### Reformulate into big matrix Y & X ####
# Using the standardized results
data=matrixPrep(stdize(salesG.cate), stdize(exog.cate), NA, 2, N) # data object
Y.vec=c(data$Ytrain) # (53*30)
X.big=kronecker(diag(K.left), data$Xtrain) # (53*30 100*30)
nameleft=colnames(salesG.cate) # Category names
nameright=colnames(data$Xtrain) 

data.sku=matrixPrep(stdize(salesG.cate), stdize(exog.sku), NA, 2, N) # data object
Y.vec.sku=c(data.sku$Ytrain) # (53*30)
X.big.sku=kronecker(diag(K.left), data.sku$Xtrain) # (53*30 422*30)
nameright.sku=colnames(data.sku$Xtrain)


#### Generate MEBoot samples X and Y for cross-validation #### 
## Generated from original IRI data
fn.logdiff=function(x) return(tsTransform(x, T, 'log.diff', 0))
set.seed(2020)
me.dev = meboot.oos.raw(dollars.cate, stdize(exog.cate), 20, F, fn.logdiff, 1)
# includes X/Y of salesG and X/Y of sales Dollar
meboot.view(dollars.cate[-c(1, nrow(dollars.cate)), ], me.dev$X.dollars)
meboot.view(data$Xtrain, me.dev$X.block)

me.dev.sku = meboot.oos.raw(dollars.cate, stdize(exog.sku), 20, F, fn.logdiff, 1)
set.seed(2019)
me.test = meboot.oos.raw(dollars.cate, stdize(exog.cate), 1, F, fn.logdiff, 1)
me.test.sku = meboot.oos.raw(dollars.cate, stdize(exog.sku), 1, F, fn.logdiff, 1)

#### (1) Standard IGLS estimation on IRI#### 
# BIC as criteria, standard weighted
fit=var.iter.est(X.big = X.big, Y.vec = Y.vec, k = K.left, 
                 l1.seq = 1/seq(1e2, K.left, length.out = 100), 
                 l2.seq = 1/seq(1e3, K.left**2, length.out = 100), 
                 std1 = T, criteria = 'BIC1', scaling = c(F, F), decomposition = 'eigen',
                 conv.beta = 10^-2, conv.sigma = 10^-1, max.iter = 10, me.object = NA, 
                 verbose = T, names = list(nameleft, nameright), weighted = 'standard', alpha = 1)
thetas[[1]]=theta.fn(nameleft, nameright, fit$beta.est, K.left)

## Visualize the theta effcts and covariance matrix
pdf(file="figure/IRI-est-bic1.pdf")
corrplot(thetas[[1]][, grep("SALESG", nameright)], is.corr=F, tl.srt=45, mar=c(0,0,1,0))
corrplot(thetas[[1]][, grep("PZ.ADJ|PR.RATE", nameright)], is.corr=F, tl.srt=45, mar=c(0,0,1,0))
corrplot(thetas[[1]][, grep("AD.RATE", nameright)], is.corr=F, tl.srt=45, mar=c(0,0,1,0))
corrplot(thetas[[1]][, grep("DI.RATE", nameright)], is.corr=F, tl.srt=45, mar=c(0,0,1,0))
corrplot(fit$sigma.est, is.corr=F)
dev.off()

# MEBoot CV as criteria
fit=var.iter.est(X.big = X.big, Y.vec = Y.vec, k = K.left, 
                 l1.seq = 1/seq(1e2, K.left, length.out = 100), 
                 l2.seq = 1/seq(1e3, K.left**2, length.out = 100), 
                 std1 = T, criteria = 'MEBOOT', scaling = c(F, F), decomposition = 'eigen',
                 conv.beta = 10^-2, conv.sigma = 10^-1, max.iter = 10, me.object = me.dev, 
                 verbose = T, names = list(nameleft, nameright), weighted = 'standard', alpha = 1)

thetas[[2]]=theta.fn(nameleft, nameright, fit$beta.est, K.left)
pdf(file="figure/IRI-est-MEBOOT.pdf")
corrplot(thetas[[2]][, grep("SALESG", nameright)], is.corr=F, tl.srt=45, mar=c(0,0,1,0))
corrplot(thetas[[2]][, grep("PZ.ADJ|PR.RATE", nameright)], is.corr=F, tl.srt=45, mar=c(0,0,1,0))
corrplot(thetas[[2]][, grep("AD.RATE", nameright)], is.corr=F, tl.srt=45, mar=c(0,0,1,0))
corrplot(thetas[[2]][, grep("DI.RATE", nameright)], is.corr=F, tl.srt=45, mar=c(0,0,1,0))
corrplot(fit$sigma.est, is.corr=F)
dev.off()

## Evaluate the esitmated theta with IRI-data (train) MEBoot samples (test)
eval.train=evaluation(thetas[[1]], 
                      data$Ytrain,
                      data$Xtrain,
                      dollars.cate[3:Times, ],
                      dollars.cate[2:N, ], 0)
eval.test=evaluation(thetas[[1]], 
                     me.test$Y.vec, 
                     me.test$X.block, 
                     me.test$Y.dollars, 
                     me.test$X.dollars, 0)
eval.train$mae;eval.train$mape
eval.test$mae;eval.test$mape

## Output the forecast results, MAE of salesg and MAPE of sales dollars
csv.predict(eval.train, 'train', "csv/IRI-eval-bic1.csv")
csv.predict(eval.test, 'test', 'csv/IRItest-eval-bic1.csv')

# Estimation visualization of evaluation result on MEBoot 
pdf(file="figure/IRI-eval-train-bic1.pdf")
plot.predict(eval.train, F)
dev.off()
pdf(file="figure/IRI-eval-trainre-bic1.pdf")
plot.predict(eval.train, T)
dev.off()

pdf(file="figure/IRI-eval-test-bic1.pdf")
plot.predict(eval.test, F)
dev.off()
pdf(file="figure/IRI-eval-testre-bic1.pdf")
plot.predict(eval.test, T)
dev.off()

#### (2) Category theta effect range #### 
# fit.st: IGLS with standard weighted matrix (sames as (1))
# fit.ol: IGLS with ownlight weighted matrix 
# fit.of: IGLS with ownfree weighted matrix
fit.st=var.iter.est(X.big, Y.vec, K.left, 
                    1/seq(1e2, K.left, length.out = 100), 
                    1/seq(1e3, K.left**2, length.out = 100), 
                    T, 'BIC1', c(F, F), 'eigen',
                    10^-2, 10^-1, 10, NA, T, 
                    list(nameleft, nameright), 'standard', 1)
thetas[[1]]=theta.fn(nameleft, nameright, fit.st$beta.est, K.left)
fit.ol=var.iter.est(X.big, Y.vec, K.left, 
                    1/seq(1e2, K.left, length.out = 100), 
                    1/seq(1e3, K.left**2, length.out = 100), 
                    T, 'BIC1', c(F, F), 'eigen',
                    10^-2, 10^-1, 10, NA, T, 
                    list(nameleft, nameright), 'ownlight', 1)
thetas[[2]]=theta.fn(nameleft, nameright, fit.ol$beta.est, K.left)
fit.of=var.iter.est(X.big, Y.vec, K.left, 
                    1/seq(1e2, K.left, length.out = 100), 
                    1/seq(1e3, K.left**2, length.out = 100), 
                    T, 'BIC1', c(F, F), 'eigen',
                    10^-2, 10^-1, 10, NA, T, 
                    list(nameleft, nameright), 'ownfree', 1)
thetas[[3]]=theta.fn(nameleft, nameright, fit.of$beta.est, K.left)

## Visualize all the effect range of each categories
pdf(file="figure/IRI-effect-range-BIC1.pdf")
effect.interval(thetas)
dev.off()


#### (3) SKUs theta effect interval: weighted estimation by LASSO #### 
# Input X (Salesg @ category + Exogenous @ SKU) 
# fit.st: IGLS with standard weighted matrix
# fit.ol: IGLS with ownlight weighted matrix 
# fit.of: IGLS with ownfree weighted matrix
fit=var.iter.est(X.big.sku, Y.vec.sku, K.left, 
                 1/seq(1e2, K.left, length.out = 100), 
                 1/seq(1e1, K.left, length.out = 100)**2, 
                 T, 'BIC1', c(T, T), 'eigen',
                 10^-2, 10^-1, 10, NA, T, 
                 list(nameleft, nameright.sku), 'standard', 1)
thetas_sku[[1]]=theta.fn(nameleft, nameright.sku, fit$beta.est, K.left)
fit.ol=var.iter.est(X.big.sku, Y.vec.sku, K.left, 
                    1/seq(1e2, K.left, length.out = 100), 
                    1/seq(1e1, K.left, length.out = 100)**2, 
                    T, 'BIC1', c(T, T), 'eigen',
                    10^-2, 10^-1, 10, NA, T, 
                    list(nameleft, nameright.sku), 'ownlight', 1)
thetas_sku[[2]]=theta.fn(nameleft, nameright.sku, fit.ol$beta.est, K.left)
fit.of=var.iter.est(X.big.sku, Y.vec.sku, K.left, 
                    1/seq(1e2, K.left, length.out = 100), 
                    1/seq(1e1, K.left, length.out = 100)**2, 
                    T, 'BIC1',c(T, T), 'eigen',
                    10^-2, 10^-1, 10, NA, T, 
                    list(nameleft, nameright.sku), 'ownfree', 1)
thetas_sku[[3]]=theta.fn(nameleft, nameright.sku, fit.of$beta.est, K.left)

## Visualize all the effect range of top10 SKUs of each categories
pdf(file="figure/IRI-effect-range-sku-BIC1.pdf")
effect.interval(thetas_sku)
dev.off()
