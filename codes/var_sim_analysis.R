##############################
# Compare the result b/w ceoptim and igls
# Draw the boxplot of S simluation runs
##############################
library(reshape2)

setwd("C://Users/jhjoo/Desktop/var/final.codes/")

sim_type = 'en'
q = 30
theta.star = get(load(file=sprintf("data/simulation/sim-%s-%d/theta.star.RData", sim_type, q)))
thetas = get(load(file=sprintf("data/simulation/sim-%s-%d/thetas.RData", sim_type, q)))
sigmas = get(load(file=sprintf("data/simulation/sim-%s-%d/sigmas.RData", sim_type, q)))
results = get(load(file=sprintf("data/simulation/sim-%s-%d/results.RData", sim_type, q)))
thetas.ce = get(load(file=sprintf("data/simulation/sim-%s-%d/thetas.ce.RData", sim_type, q)))
sigmas.ce = get(load(file=sprintf("data/simulation/sim-%s-%d/sigmas.ce.RData", sim_type, q)))
results.ce = get(load(file=sprintf("data/simulation/sim-%s-%d/results.ce.RData", sim_type, q)))
n_result = ncol(results)
if (n_result != ncol(results.ce)){
  print("Mismatch number of results.")
}
View(results)

data = data.frame(t(cbind(results, results.ce)))
colnames(data) = c("TPR(en)", "TPR(ex)","TNR(en)", "TNR(ex)",
                   "MAEE", "MAFE(tr)", "MAFE(te)")
data$method = c(rep("1_IGLS", n_result), rep("2_CEoptim", n_result))

recover_data=melt(data[, -(5:7)], id.vars="method", variable.name="metric", value.name="value")
maee_data=melt(data[, c(5, 8)], id.vars="method", variable.name="metric", value.name="value")
accuracy_data=melt(data[, -(1:5)], id.vars="method", variable.name="metric", value.name="value")

dev.new(width = 550, height = 300, unit = "px")

# recover (TPR/TNR)
boxplot(recover_data$value~recover_data$method+recover_data$metric, data=recover_data,
        boxwex=0.4, xaxt='n', xlab=NULL, ylab=NULL, col=c(2,4),
        main="Recall of Coefficients")
axis(1, at=seq(1.5 , 20 , 2)[1:4], 
     labels=c("TPR(en)", "TPR(ex)","TNR(en)", "TNR(ex)"), tick=F,
     cex=0.3)

abline(v=0.5,lty=1, col="grey")
abline(v=2.5,lty=1, col="grey")
abline(v=4.5,lty=1, col="grey")
abline(v=6.5,lty=1, col="grey")
abline(v=8.5,lty=1, col="grey")

legend("bottomleft", legend = c("Iterative GLS", "CEoptim"), 
       inset = .01, fill=c(2, 4), bg = 'white',
       col=c(2 , 4), pt.cex = 1, cex = 0.8,  
       horiz = F, title="Method")

# MAEE
boxplot(maee_data$value~maee_data$method+maee_data$metric, data=maee_data,
        boxwex=0.4, xaxt='n', xlab=NULL, ylab=NULL, col=c(2,4), ylim=c(0, 0.02),
        main="Estimation Error")

axis(1, at=seq(1.5 , 20 , 2)[1], 
     labels=c("MAEE"), tick=F,
     cex=0.3)

abline(v=0.5,lty=1, col="grey")
abline(v=2.5,lty=1, col="grey")

legend("bottomleft", legend = c("Iterative GLS", "CEoptim"), 
       inset = .01, fill=c(2, 4), bg = 'white',
       col=c(2 , 4), pt.cex = 1, cex = 0.8,  
       horiz = F, title="Method")

# MAFE (train/test)
boxplot(accuracy_data$value~accuracy_data$method+accuracy_data$metric, data=accuracy_data,
        boxwex=0.4, xaxt='n', xlab=NULL, ylab=NULL, col=c(2,4),
        main="Forecast Error")

axis(1, at=seq(1.5 , 20 , 2)[1:2], 
     labels=c("MAFE(tr)", "MAFE(te)"), tick=F,
     cex=0.3)

abline(v=0.5,lty=1, col="grey")
abline(v=2.5,lty=1, col="grey")
abline(v=4.5,lty=1, col="grey")

legend("topleft", legend = c("Iterative GLS", "CEoptim"), 
       inset = .01, fill=c(2, 4), bg = 'white',
       col=c(2 , 4), pt.cex = 1, cex = 0.8,  
       horiz = F, title="Method")

# # thetas analysis
# [Deprecated] Check the ceoptim theta and igls theta difference

# which.min(results[7, ])
# which.min(results.ce[7, ])
# pick=which(c(theta.star))
# diff=which(c(thetas[[2]])[pick] + c(thetas.ce[[2]])[pick] !=0)
# max_coef=max(c(thetas[[2]][pick], thetas.ce[[2]][pick]))
# min_coef=min(c(thetas[[2]][pick], thetas.ce[[2]][pick]))
# 
# plot(NULL, xlim=c(0,length(pick)), ylim=c(min_coef,max_coef), ylab="value", xlab="Theta elements")
# points(c(theta.star)[pick], type='p', col=1)
# points(c(thetas[[2]])[pick], type='p', col=2)
# points(c(thetas.ce[[2]])[pick], type='p', col=4)

