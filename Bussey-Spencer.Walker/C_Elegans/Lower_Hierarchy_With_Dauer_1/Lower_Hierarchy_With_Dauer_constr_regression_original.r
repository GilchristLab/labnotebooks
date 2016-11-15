#On Git
#Model with noise, with intercept
#Only lower Hierarchy datapoints
#With Dauer
rm(list=ls())
#load libraries
library(rstan)
library(coda)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#the explanatory variables
emp.expr <- read.table(file="Data/Nonoverlapping_Lowest_In_Hierarchy_With_Dauer_Celegans_hermaphrodite_expression.csv", sep=",", as.is=T, header=T)
# data
estm.phi <- read.table(file = "Data/PHI_without_xobs_singlechain.csv", header=T, sep=",")

idx <- emp.expr$Gene.ID %in% estm.phi$ORF_Info
emp.expr <- log(emp.expr[,-c(1,2)])
emp.expr <- emp.expr[idx, ]
row.remove <- NULL
for(i in 1:7300)
  if(any(is.na(emp.expr[i,]))){row.remove <- c(row.remove, i) }
    
emp.expr <- emp.expr[-row.remove,]

#the explanatory variables
dat <- data.frame(emp.expr)
means <- colMeans(dat, na.rm = T)
X <- dat
for(i in 1:length(means))
  X[,i] <- dat[,i] - means[i]
#Adds the Intercept
new.col.names <- c("Intercept", colnames(X))
X <- cbind(rep(1, nrow(X)), X)
colnames(X) <- new.col.names




#the estimated data
y <- log(estm.phi$Phi_Post_Arith_Mean[-row.remove])
y.sd <- estm.phi$logPhi_Post_SD[-row.remove]

#data
celDat <- list(N=nrow(X), K=ncol(X), y=y, X=X, sigmaRoc=y.sd, dummyOnes=rep(1, nrow(X)))

#the model
celeg <- stan(file="constr_regression_original.stan", data = celDat, pars = c("beta", "sigma2", "logLik", "tPhi"), chains=1, iter=2000)

print(celeg)

logLik2 <- as.data.frame(celeg, "logLik")

post <- as.matrix(celeg, "tPhi")
meds <- apply(post,2,mean)
Y <- rapply(X[,-1], c)

sigma2 <- as.data.frame(celeg, "sigma2")
meansigma2 <- colMeans(sigma2)
hdisigma2 <- t(apply(sigma2, 2, quantile, c(0.025, 0.975)))

beta <- as.data.frame(celeg, "beta")
meanBeta <- colMeans(beta)
hdiBeta  <- t(apply(beta, 2, quantile, c(0.025, 0.975)))

res <- rep(NA, nrow(X))
for(i in 1:nrow(X)){
  res[i] <- sum(X[i,-1]*meanBeta[-1], na.rm = T)
}


pdf("lower_hierarchy_with_dauer_original_model_celegans_weighted_expr_vector_beta_noise_fixed_sroc_withIntercept.pdf", height=20)
par(mfrow=c(4,1))
plot(res, y, type="n", xlab = "log(Weighted Mean Emperical Expression)", ylab = "log(ROC Estimated Expression)")
ribModel:::upper.panel.plot(res, y)

plot(meds, Y, type="n", xlab = "log(Estimated 'True' Expression)", ylab = "log(Emperical Expression)")
ribModel:::upper.panel.plot(meds, Y)
#linear regression
t.spent <- c(20,190,60,100,180,690,420,450,570,1440,NA,NA,NA,NA)
rm.idx <- c(13,14,15,16,17)
names(t.spent) <- colnames(dat)
plot(t.spent, meanBeta[-1], type='n', xlab="Time in Lifestage [min]", ylab="Expression Coefficient Omega")
ribModel:::upper.panel.plot(t.spent, meanBeta[-1])

order.occurred <-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
names(order.occurred) <- colnames(dat)
plot(order.occurred,meanBeta[-1],type= 'n', xlab= "Order of Life Stages", ylab="Expression Coefficient Omega")
ribModel:::upper.panel.plot(order.occurred, meanBeta[-1])

#Reads this way if you have the Celegans_mean_expression_confirmed file
#mean.expr <- read.table(file="Data/Celegans_hermaphrodite_mean_expression_confirmed.csv", sep=",", as.is=T, header=T)
#mean.expr <- log(mean.expr[-row.remove,2])

####This is the same graph everytime.
#celeg.mean.expr <- rowMeans(dat, na.rm =T)
#mean.expr <- celeg.mean.expr
#plot(mean.expr, y, type="n", xlab = "log(Mean Emperical Expression)", ylab = "log(ROC Estimated Expression)")
#ribModel:::upper.panel.plot(mean.expr, y,)

dev.off()

pdf("lower_hierarchy_with_dauer_original_model_celegans_weighted_expr_vector_beta_noise_fixed_sroc_singles.pdf", height=20, width=20)
par(mfrow=c(4,5))
cn <- colnames(X)
for(i in 2:14)
{
  plot(as.numeric(X[,i]), y, type="n", xlab = "log(Emperical Expression)", ylab = "log(ROC Estimated Expression)", main = cn[i])
  ribModel:::upper.panel.plot(X[,i], y)
}
par(mfrow=c(4,5))
for(i in 2:14)
{
  plot(as.numeric(X[,i]*meanBeta[i]), y, type="n", xlab = "log(Emperical Expression)", ylab = "log(ROC Estimated Expression)", main = cn[i])
  ribModel:::upper.panel.plot(X[,i]*meanBeta[i], y)
}
dev.off()

#Need to Manually store the mean log likelihood
plot(logLik2[,1], type = "l")
meanLogLik <- mean(logLik2[,1])

plot(density(logLik2[,1]))

#Export Results
#Don't forget to adjust the output file names based on the date.
lifestages <- as.matrix(t.spent)
dnames <- rownames(lifestages)
rownames(lifestages) <- NULL
lifestages <- cbind(dnames, lifestages)
lifestages <- rbind(c("Intercept",NA), lifestages)

betaResults <- as.data.frame(cbind(lifestages, meanBeta, hdiBeta))
write.csv(betaResults, "Lower_Hierarchy_With_Dauer_Results_2Aug2016T1.csv")

meanLogLik <- as.data.frame(mean(logLik2[,1]))
write.csv(meanLogLik, "Lower_Hierarchy_With_Dauer_Mean_Log_Likelihood_2Aug2016T1.csv")