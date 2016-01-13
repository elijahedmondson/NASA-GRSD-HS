# install.packages("survival")
library(splines)
library(survival)
library(KMsurv)
library(OIsurv)


attach(GRSD)
names(GRSD)

time <- days
event <- HCC
X <- cbind(GROUP, sex, family, coat.color)
group <- GROUP

#Log-logistic parametric model coefficients
loglogistic <- survreg(Surv(time,event) ~ X, dist="loglogistic")
summary(loglogistic)

# Cox proportional hazard model - coefficients and hazard rates
coxph <- coxph(Surv(time,event) ~ X, method="breslow")
summary(coxph)








# Descriptive statistics
summary(time)
summary(event)
summary(X)
summary(group)

# Kaplan-Meier non-parametric analysis
kmsurvival <- survfit(Surv(time,event) ~ 1)
summary(kmsurvival)
plot(kmsurvival, mark.time=TRUE, mark=1, col=1, lty=1, lwd=3, cex=1, log=FALSE, xlab="Days Post-Irradiation", ylab="Survival Probability")
title("GRSD Lymphoma")

# Kaplan-Meier non-parametric analysis by group
kmsurvival1 <- survfit(Surv(time, event) ~ group.broad)
summary(kmsurvival1)
plot(kmsurvival1, conf.int="both", mark.time=TRUE, 
     mark=4, col=1, lty=1:5, lwd=2, cex=1, log=FALSE, xscale=1, yscale=1,  
     firstx=0, firsty=1, ymin=0, xlab="Days Post-Irradiation", ylab="Survival Probability")
legend(30, .3, c("Gamma", "HZE", "Unirradiated"), lty = 1:5, lwd=2) 
title("GRSD Lymphoma")

#The previous without 95% CI
kmsurvival1 <- survfit(Surv(time, event) ~ group.broad)
summary(kmsurvival1)
plot(kmsurvival1, conf.int="none", mark.time=TRUE, 
     mark=4, col=1, lty=1:5, lwd=2, cex=1, log=FALSE, xscale=1, yscale=1,  
     firstx=0, firsty=1, ymin=0, xlab="Days Post-Irradiation", ylab="Survival Probability")
legend(30, .3, c("Gamma", "HZE", "Unirradiated"), lty = 1:5, lwd=2) 
title("GRSD Lymphoma")

# Nelson-Aalen non-parametric analysis
nasurvival <- survfit(coxph(Surv(time,event)~1), type="aalen")
summary(nasurvival)
plot(nasurvival, xlab="Time", ylab="Survival Probability")


# Cox proportional hazard model - coefficients and hazard rates
coxph <- coxph(Surv(time,event) ~ X, method="breslow")
summary(coxph)


# Exponential, Weibull, and log-logistic parametric model coefficients
exponential <- survreg(Surv(time,event) ~ X, dist="exponential")
summary(exponential)

weibull <- survreg(Surv(time,event) ~ X, dist="weibull")
summary(weibull)

loglogistic <- survreg(Surv(time,event) ~ X, dist="loglogistic")
summary(loglogistic)

