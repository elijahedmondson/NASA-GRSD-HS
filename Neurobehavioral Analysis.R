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
###   Neurobehavioral Analysis   ###

# LOAD PHENO FILES #
neuro.pheno600 <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/neuro.pheno600.csv")
neuro.pheno1200 <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/neuro.pheno1200.csv")
neuro.pheno <- rbind(neuro.pheno600, neuro.pheno1200)
library(dplyr)
rename(neuro.pheno, c("X"="row.names"))
names(neuro.pheno)[names(neuro.pheno)=="X"] <- "row.names"

# QC ON PHENO TRAITS #
Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Total-Table 1.csv")
Neuro <- read.csv("~/Desktop/Neurophenotype.EFE/Total-Table 1.csv")
Total1 <- as.data.frame(Total)
Neuro1 <- as.data.frame(Neuro)

 <- merge(x = Total1, y = Neuro1, by.x = "row.names", by.y = "row.names")

test <- left_join(Total, Neuro, by.x = "row.names", by.y = "row.names")

write.csv(test, file = "~/Desktop/test.csv")
save(Neuro, test, Total, file = "Neuro.files.Rdata")


neuro <- read.csv("~/Desktop/R/GRSD.phenotype/Neurophenotype.EFE/Neuro.remove.mismatch/Neuro.Deleted.Mismatch-Table 1.csv")
pheno = data.frame(row.names = neuro$row.names, sex = as.numeric(neuro$sex == "M"),
                   group = as.character(neuro$groups),
                   cohort = as.numeric(neuro$Cohort),
                   batch = as.numeric(neuro$Batch),
                   pig.dis = as.numeric(neuro$pigmentdispersion_bool),
                   days = as.numeric(neuro$days),
                   neur1 = as.numeric(neuro$context_avgmo_total),
                   neur2 = as.numeric(neuro$context_pctfrze_total),
                   neur3 = as.numeric(neuro$cued_tone_pctfrze_total),
                   neur4 = as.numeric(neuro$train_deltaavgmot_shock1_shock5))
HZE <- subset(pheno, group == "HZE")
HZEadd = matrix(HZE$sex, ncol = 1, dimnames = list(rownames(HZE), "sex"))

gamma <- subset(pheno, group == "Gamma")
GAMMAadd = matrix(gamma$sex, ncol = 1, dimnames = list(rownames(gamma), "sex"))

unirradiated <- subset(pheno, group == "Unirradiated")
UNadd = matrix(unirradiated$sex, ncol = 1, dimnames = list(rownames(unirradiated), "sex"))

addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))

layout(matrix(3:1, 3, 1)) 
boxplot(neur2~group,data=pheno, main="", notch = T, col = "green",
        xlab="", ylab="context_pctfrze_total")
boxplot(neur1~sex,data=pheno, main="", notch = T, col = "green",
        xlab="", ylab="context_avgmo_total")

hist(HZE$neur2, main = "Total: context_pctfrze_total", col = "green")
lines(density(pheno$neur2), na.rm = T)


library(DOQTL)

load(file = "~/Desktop/R/Build/K.Rdata")
load(file = "~/Desktop/R/Build/model.probs.Rdata")
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
outdir = "~/Desktop/"
setwd(outdir)

sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz"


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


