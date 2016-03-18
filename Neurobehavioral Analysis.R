library(DOQTL)
library(doParallel)
library(foreach)
library(Rsamtools)
library(VariantAnnotation)
library(GenomicRanges)
library(regress)
library(MASS)
library(lmtest)
library(HZE)

load(file = "~/Desktop/R/Build/K.Rdata")
load(file = "~/Desktop/R/Build/model.probs.Rdata")
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
outdir = "~/Desktop/"
setwd(outdir)

sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz"

###   Neurobehavioral Analysis   ###

# LOAD PHENO FILES #
GRSD.pheno <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = GRSD.pheno$row.names, sex = as.numeric(GRSD.pheno$sex == "M"),
                   cohort = as.numeric(GRSD.pheno$Cohort),
                   group = as.character(GRSD.pheno$groups),
                   family = as.character(GRSD.pheno$family),
                   weight = as.numeric(GRSD.pheno$weight),
                   unirradiated = as.numeric(GRSD.pheno$Unirradiated),
                   Frz.total = as.numeric(GRSD.pheno$context_pctfrze_total),
                   pig.dis = as.numeric(GRSD.pheno$pigmentdispersion),
                   avgmo.total = as.numeric(GRSD.pheno$context_avgmo_total),
                   tone.frz = as.numeric(GRSD.pheno$cued_tone_pctfrze_total),
                   train.frz = as.numeric(GRSD.pheno$train_deltapctfrze_isi1_isi4),
                   train.shock = as.numeric(GRSD.pheno$train_deltaavgmot_shock1_shock5), 
                   Albino = as.numeric(GRSD.pheno$albino),
                   PulACA = as.numeric(GRSD.pheno$Pulmonary.Adenocarcinoma),
                   HCC = as.numeric(GRSD.pheno$Hepatocellular.Carcinoma),
                   LSA.PreT = as.numeric(GRSD.pheno$PreT))
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))

cohort1 = subset(pheno, cohort == 1)
cohort2 = subset(pheno, cohort == 2)
CO1covar = matrix(cohort1$sex, ncol = 1, dimnames = list(rownames(cohort1), "sex"))
CO2covar = matrix(cohort2$sex, ncol = 1, dimnames = list(rownames(cohort2), "sex"))


HZE.1 <- subset(cohort1, group == "HZE")
HZE.1add = matrix(HZE.1$sex, ncol = 1, dimnames = list(rownames(HZE.1), "sex"))
HZE.2 <- subset(cohort2, group == "HZE")
HZE.2add = matrix(HZE.2$sex, ncol = 1, dimnames = list(rownames(HZE.2), "sex"))
gamma.1 <- subset(pheno, group == "Gamma")
GAMMA.1add = matrix(gamma$sex, ncol = 1, dimnames = list(rownames(gamma), "sex"))
gamma.2 <- subset(pheno, group == "Gamma")
GAMMA.2add = matrix(gamma$sex, ncol = 1, dimnames = list(rownames(gamma), "sex"))
unirradiated.1 <- subset(pheno, group == "Unirradiated")
UN.1add = matrix(unirradiated$sex, ncol = 1, dimnames = list(rownames(unirradiated), "sex"))
unirradiated.2 <- subset(pheno, group == "Unirradiated")
UN.2add = matrix(unirradiated$sex, ncol = 1, dimnames = list(rownames(unirradiated), "sex"))

irradiated <- subset(pheno, unirradiated == 0)
IRRadd = matrix(irradiated$sex, ncol = 1, dimnames = list(rownames(irradiated), "sex"))

### BOXPLOT ###
layout(matrix(2:1, 2, 1))
boxplot(weight~family,data = pheno, main = "", notch = T, col = c("green", "red", "blue"),
        xlab="", ylab="")

boxplot(train.shock~group,data=cohort1, main="", notch = T, col = c("green", "red", "blue"),
        xlab="", ylab="context_avgmo_total")

layout(matrix(3:1, 3, 1))
hist(unirradiated$train.shock, main = "Unirradiated", col = "blue")
hist(HZE$train.shock, main = "HZE", col = "red")
hist(gamma$train.shock, main = "Gamma", col = "green")

### MAP ###
QTL.C1.frz.total = scanone.assoc(pheno = HZE.1, pheno.col = "train.frz", probs = model.probs, K = K, addcovar = HZE.1add, markers = MM_snps, sdp.file = sdp.file, ncl = 4)
QTL.C2.frz.total = scanone.assoc(pheno = HZE.2, pheno.col = "train.frz", probs = model.probs, K = K, addcovar = HZE.2add, markers = MM_snps, sdp.file = sdp.file, ncl = 4)

QTL.IRR.frz.total = scanone.assoc(pheno = CO1covar, pheno.col = "Frz.total", probs = model.probs, K = K, addcovar = IRRadd, markers = MM_snps, sdp.file = sdp.file, ncl = 4)
QTL.IRR.avgmo.total = scanone.assoc(pheno = irradiated, pheno.col = "avgmo.total", probs = model.probs, K = K, addcovar = IRRadd, markers = MM_snps, sdp.file = sdp.file, ncl = 4)
QTL.IRR.tone.frz = scanone.assoc(pheno = irradiated, pheno.col = "tone.frz", probs = model.probs, K = K, addcovar = IRRadd, markers = MM_snps, sdp.file = sdp.file, ncl = 4)
QTL.IRR.train.frz = scanone.assoc(pheno = irradiated, pheno.col = "train.frz", probs = model.probs, K = K, addcovar = IRRadd, markers = MM_snps, sdp.file = sdp.file, ncl = 4)
QTL.IRR.train.shock = scanone.assoc(pheno = irradiated, pheno.col = "train.shock", probs = model.probs, K = K, addcovar = IRRadd, markers = MM_snps, sdp.file = sdp.file, ncl = 4)

### QTL PLOT ###
DOQTL:::plot.scanone.assoc(QTL.C1.frz.total, bin.size = 100, main = "Cohort 1, HZE (n = 190): train_deltapctfrze_isi1_isi4")
DOQTL:::plot.scanone.assoc(QTL.C2.frz.total, bin.size = 100, main = "Cohort 2, HZE (n = 417): train_deltapctfrze_isi1_isi4")
DOQTL:::plot.scanone.assoc(QTL.IRR.avgmo.total, bin.size = 100, main = "All Irradiated: context_avgmo_total")
DOQTL:::plot.scanone.assoc(QTL.IRR.tone.frz, bin.size = 100, main = "All Irradiated: cued_tone_pctfrze_total")
DOQTL:::plot.scanone.assoc(QTL.IRR.train.frz, bin.size = 100, main = "All Irradiated: train_deltapctfrze_isi1_isi4")
DOQTL:::plot.scanone.assoc(QTL.IRR.train.shock, bin.size = 100, main = "All Irradiated: train_deltaavgmot_shock1_shock5")

save(QTL.IRR.train.shock, QTL.IRR.tone.frz, QTL.IRR.train.frz, QTL.IRR.avgmo.total, QTL.IRR.frz.total, file="~/Desktop/NEURO.unirradiated.QTL.Rdata")

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


