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
                   group = as.character(GRSD.pheno$groups),
                   Frz.min1 = as.numeric(GRSD.pheno$context_pctfrze_minute1),
                   Frz.min2 = as.numeric(GRSD.pheno$context_pctfrze_minute2),
                   Frz.min3 = as.numeric(GRSD.pheno$context_pctfrze_minute3),
                   Frz.min4 = as.numeric(GRSD.pheno$context_pctfrze_minute4),
                   Frz.min5 = as.numeric(GRSD.pheno$context_pctfrze_minute5),
                   Frz.total = as.numeric(GRSD.pheno$context_pctfrze_total),
                   pig.dis = as.numeric(GRSD.pheno$pigmentdispersion),
                   Cont.min1 = as.numeric(GRSD.pheno$context_avgmo_minute1),
                   Cont.min2 = as.numeric(GRSD.pheno$context_avgmo_minute2),
                   Cont.min3 = as.numeric(GRSD.pheno$context_avgmo_minute3),
                   Cont.min4 = as.numeric(GRSD.pheno$context_avgmo_minute4),
                   Cont.min5 = as.numeric(GRSD.pheno$context_avgmo_minute5),
                   Cont.total = as.numeric(GRSD.pheno$context_avgmo_total),
                   = as.numeric(GRSD.pheno$),
                   neur1 = as.numeric(GRSD.pheno$context_avgmo_total),
                   neur2 = as.numeric(GRSD.pheno$context_pctfrze_total),
                   neur3 = as.numeric(GRSD.pheno$cued_tone_pctfrze_total),
                   neur4 = as.numeric(GRSD.pheno$train_deltaavgmot_shock1_shock5))

HZE <- subset(pheno, group == "HZE")
HZEadd = matrix(HZE$sex, ncol = 1, dimnames = list(rownames(HZE), "sex"))

gamma <- subset(pheno, group == "Gamma")
GAMMAadd = matrix(gamma$sex, ncol = 1, dimnames = list(rownames(gamma), "sex"))

unirradiated <- subset(pheno, group == "Unirradiated")
UNadd = matrix(unirradiated$sex, ncol = 1, dimnames = list(rownames(unirradiated), "sex"))

addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))

layout(matrix(3:1, 3, 1))
layout(matrix(1:1, 3, 1))
boxplot(neur2~group,data=pheno, main="", notch = T, col = "green",
        xlab="", ylab="context_pctfrze_total")
boxplot(neur1~sex,data=pheno, main="", notch = T, col = "green",
        xlab="", ylab="context_avgmo_total")

hist(HZE$neur2, main = "Total: context_pctfrze_total", col = "green")
lines(density(pheno$neur2), na.rm = T)



qtl = scanone.assoc(pheno = unirradiated, pheno.col = "neur2", probs = model.probs, K = K, 
                    addcovar = UNadd, markers = MM_snps, sdp.file = sdp.file, ncl = 4)

save(qtl, file = "context_pctfrze_total.Rdata")


DOQTL:::plot.scanone.assoc(qtl, bin.size = 100)

DOQTL:::plot.scanone.assoc(qtl, bin.size = 100, main = "Unirradiated: context_pctfrze_total")


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


