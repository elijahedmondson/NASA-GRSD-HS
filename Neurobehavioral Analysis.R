library(HZE)
library(ggplot2)
library(doParallel)
library(foreach)
library(Rsamtools)
library(VariantAnnotation)
library(GenomicRanges)
library(regress)
library(MASS)
library(lmtest)
library(car)
library(DOQTL)

############ LOAD FILES ############ 
############ LOAD FILES ############ 
############ LOAD FILES ############ 
############ LOAD FILES ############ 

load(file = "~/Desktop/R/QTL/WD/HS\ HMM\ Rdata/K.Rdata")
load(file = "~/Desktop/R/QTL/WD/HS\ HMM\ Rdata/model.probs.Rdata")
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
outdir = "~/Desktop/"
setwd(outdir)
sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz"
GRSD.pheno <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = GRSD.pheno$row.names, sex = as.numeric(GRSD.pheno$sex == "M"),
                   cohort = as.numeric(GRSD.pheno$Cohort),
                   group = as.character(GRSD.pheno$groups),
                   family = as.character(GRSD.pheno$family),
                   weight = as.numeric(GRSD.pheno$Weight.corrected),
                   unirradiated = as.numeric(GRSD.pheno$Unirradiated),
                   Frz.total = as.numeric(GRSD.pheno$context_pctfrze_total),
                   pig.dis = as.numeric(GRSD.pheno$pigment.dispersion),
                   avgmo.total = as.numeric(GRSD.pheno$context_avgmo_total),
                   tone.frz = as.numeric(GRSD.pheno$cued_tone_pctfrze_total),
                   train.frz = as.numeric(GRSD.pheno$train_deltapctfrze_isi1_isi4),
                   train.shock = as.numeric(GRSD.pheno$train_deltaavgmot_shock1_shock5))

############ SPECIFY COVARIATES ############ 
############ SPECIFY COVARIATES ############ 
############ SPECIFY COVARIATES ############ 
############ SPECIFY COVARIATES ############ 

addcovar = cbind(addcovar, cohort = pheno$cohort)
addcovar = cbind(pheno$row.names, sex = pheno$sex, cohort = pheno$cohort)

HZE <- subset(pheno, group == "HZE")
Gamma <- subset(pheno, group == "Gamma")
Unirradiated <- subset(pheno, group == "Unirradiated")
Allirr <- subset(pheno, unirradiated == 0)



############ MAP ############   
############ MAP ############  
############ MAP ############  
############ MAP ############  

QTL.C1.train.frz = scanone.assoc(pheno = Gamma, pheno.col = "train.frz", probs = model.probs, 
                                 K = K, addcovar = addcovar, markers = MM_snps, sdp.file = sdp.file, ncl = 4)

############ DATA VISUALIZATION ############
############ DATA VISUALIZATION ############  
############ DATA VISUALIZATION ############  
############ DATA VISUALIZATION ############  

layout(matrix(2:1, 2, 1))
boxplot(weight~family,data = pheno, main = "", notch = T, col = c("green", "red", "blue"),
        xlab="", ylab="")

boxplot(train.shock~group,data=cohort1, main="", notch = T, col = c("green", "red", "blue"),
        xlab="", ylab="context_avgmo_total")

layout(matrix(3:1, 3, 1))
hist(Unirradiated$train.shock,  main = "Unirradiated", col = "blue")
hist(HZE$train.shock, main = "HZE", col = "red")
hist(Gamma$train.shock, main = "Gamma", col = "green")





############ MANHATTAN PLOTs ############ 
############ MANHATTAN PLOTs ############ 
############ MANHATTAN PLOTs ############ 
############ MANHATTAN PLOTs ############ 

DOQTL:::plot.scanone.assoc(QTL.C1.train.frz, bin.size = 100, 
                           main = "Gamma: train_deltapctfrze_isi1_isi4")
# 99%
abline(h = 12.73544, col = "green")
# 95%
abline(h = 10.97533, col = "blue")
# 90%
abline(h = 10.24464, col = "yellow")
# 85%
abline(h = 10.02939, col = "orange")
# 80%
abline(h = 9.689125, col = "red")

legend(13,13, title = "1,000 Permutations", 
       c("alpha = 0.01", "alpha = 0.05", "alpha = 0.10", "alpha = 0.15", "alpha = 0.20"),
       lty=c(1,1,1,1,1), lwd=c(2, 2, 2, 2, 2),col=c("green", "blue", "yellow", "orange", "red"))



DOQTL:::plot.scanone.assoc(QTL.C2.frz.total, bin.size = 100, 
                           main = "Cohort 2, HZE (n = 417): train_deltapctfrze_isi1_isi4")
DOQTL:::plot.scanone.assoc(QTL.IRR.avgmo.total, bin.size = 100, 
                           main = "All Irradiated: context_avgmo_total")
DOQTL:::plot.scanone.assoc(QTL.IRR.tone.frz, bin.size = 100, 
                           main = "All Irradiated: cued_tone_pctfrze_total")
DOQTL:::plot.scanone.assoc(QTL.IRR.train.frz, bin.size = 100, 
                           main = "All Irradiated: train_deltapctfrze_isi1_isi4")
DOQTL:::plot.scanone.assoc(QTL.IRR.train.shock, bin.size = 100, 
                           main = "All Irradiated: train_deltaavgmot_shock1_shock5")

save(QTL.IRR.train.shock, QTL.IRR.tone.frz, QTL.IRR.train.frz, QTL.IRR.avgmo.total, QTL.IRR.frz.total, 
     file="~/Desktop/NEURO.unirradiated.QTL.Rdata")

par(mfrow = c(3,1), mar=c(1, 4, 1, 1) + 0.5)
DOQTL:::plot.scanone.assoc(HZE.freeze, bin.size = 100, main = "HZE Ion", ylim=c(0,21))
abline(a = 13, b = 0, col = "red")
DOQTL:::plot.scanone.assoc(gamma.freeze, bin.size = 100, main = "Gamma Ray", ylim=c(0,21))
abline(a = 13, b = 0, col = "red")
DOQTL:::plot.scanone.assoc(un.freeze, bin.size = 100, main = "Unirradiated", ylim=c(0,21))
abline(a = 13, b = 0, col = "red")


layout(matrix(3:1, 3, 1)) 
par(mfrow = c(2,2), mar=c(1, 4, 1, 1) + 0.1)
DOQTL:::plot.scanone.assoc(HZE.freeze, chr=1, bin.size = 100, main = "HZE Ion", ylim=c(0,30))
DOQTL:::plot.scanone.assoc(gamma.freeze, chr=1, bin.size = 100, main = "Gamma ray", ylim=c(0,30))
DOQTL:::plot.scanone.assoc(un.freeze, chr=1, bin.size = 100, main = "Unirradiated", ylim=c(0,30))
DOQTL:::plot.scanone.assoc(total.freeze, chr=1, bin.size = 100, main = "Total Cases", ylim=c(0,30))

