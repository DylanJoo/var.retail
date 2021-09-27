setwd("/Users/jhjoo/Desktop/var/final.codes/")
source("codes/var_utils.R")

load(file="data/IRI.652759.prep.cate.RData")
load(file="data/IRI.652759.prep.sku.RData")

#### category level ####
data = Reduce(function(x, y) merge(x=x, y=y, by="WEEK", all=T), IRI.652759.prep.cate)
data = tsFixed(data)

dollars = sapply(data[, grep('SALES', colnames(data))], 
                function(x) tsTransform(x, T, 'raw', 0))
salesG = sapply(data[, grep('SALES', colnames(data))], 
                    function(x) tsTransform(x, T, 'log.diff', 0))
promotion = sapply(data[, grep('PR.RATE', colnames(data))],
                   function(x) tsTransform(x, F, 'raw', 0))
pz.adjust = sapply(data[, grep('PZ.ADJ', colnames(data))], 
                     function(x) tsTransform(x, F, 'raw', 0))
advertising = sapply(data[, grep('AD.RATE', colnames(data))], 
                     function(x) tsTransform(x, F, 'raw', 0))
display = sapply(data[, grep('DI.RATE', colnames(data))], 
                    function(x) tsTransform(x, F, 'raw', 0))
exog.cate.old=cbind(promotion, advertising, display)
exog.cate=cbind(pz.adjust, advertising, display)

dollars.cate=tsFixed(dollars)
salesG.cate=tsFixed(salesG)
colnames(salesG.cate)=gsub("SALES", "SALESG", colnames(salesG.cate))
exog.cate=tsFixed(exog.cate)
exog.cate.old=tsFixed(exog.cate.old)
##############################################
# SALESG:(2:53) Dollar sales growth of category
# PZ.REDUCT:(2:53) average price reduction ratio of (each SKU/category) 
# AD: (2:53)advertising rate of category
# DI: (2:53)display rate of category
############################################
save(dollars.cate, file="data/dollars.cate.RData")
save(salesG.cate, file="data/salesG.cate.RData")
save(exog.cate, file="data/exog.cate.RData")
save(exog.cate.old, file="data/exog.cate.old.RData")


#### SKU level ####
data.sku = Reduce(function(x, y) merge(x=x, y=y, by="WEEK", all=T), IRI.652759.prep.sku)
data.sku = tsFixed(data.sku)

pr.sku = sapply(data.sku[, grep('PR\\.', colnames(data.sku))], 
                       function(x) tsTransform(x, F, 'raw', 0))
pz.adj.sku = sapply(data.sku[, grep('PZ.ADJ', colnames(data.sku))], 
                   function(x) tsTransform(x, F, 'raw', 0))
advertising.sku = sapply(data.sku[, grep('AD\\.', colnames(data.sku))], 
                     function(x) tsTransform(x, F, 'raw', 0))
display.sku = sapply(data.sku[, grep('DI\\.', colnames(data.sku))], 
                 function(x) tsTransform(x, F, 'raw', 0))
exog.sku=cbind(pz.adj.sku, advertising.sku, display.sku)

exog.sku=tsFixed(exog.sku)
##############################################
# PZ.REDUCT.sku:(2:53) price reduction respect to fixed price
# AD,sku: (2:53) advertising rate of category
# DI.sku: (2:53) display rate of category
############################################
save(exog.sku, file="data/exog.sku.RData")

###########################################
#ADF unit root test
###########################################
load(file='data/dollars.cate.RData')
load(file='data/salesG.cate.RData')
library(tseries)
adf.before=rep(NA, ncol(dollars.cate))

for (c in 1:ncol(dollars.cate)){
  adf.before[c] = adf.test(dollars.cate[, c], 
                    alternative = 'stationary')$p.value
}
plot(adf.before, type = 'h', xlab='Category', ylab='P-Value', lwd=5)
abline(h=0.05, col=2)

adf.after=rep(NA, ncol(salesG.cate))
for (c in 1:ncol(salesG.cate)){
  adf.after[c] = adf.test(salesG.cate[, c], 
                    alternative = 'stationary')$p.value
}
plot(adf.after, type = 'h', xlab='Category', ylab='P-Value', lwd=5)
abline(h=0.05, col=2)
