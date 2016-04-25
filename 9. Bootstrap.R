### Nonparametric Bootstrap Resampling with Replacement ###

library(BSgenome.Mmusculus.UCSC.mm10)
library(doParallel)
library(VariantAnnotation)
library(GenomicRanges)
library(regress)
library(MASS)
library(DOQTL)
library(lmtest)
library(HZE)
library(dplyr)
library(sm)
options(stringsAsFactors = F)
load(file = "~/Desktop/R/QTL/WD/GRSD.Rdata")
setwd("~/Desktop/files")
outdir = "~/Desktop/files"
Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = Total$row.names, rownames = Total$row.names,
                   sex = as.numeric(Total$sex == "M"),
                   cohort = as.numeric(Total$Cohort),
                   group = as.character(Total$groups),
                   unirradiated = as.numeric(Total$Unirradiated),
                   days = as.numeric(Total$days),
                   PulACA = as.numeric(Total$Pulmonary.Adenocarcinoma),
                   HCC = as.numeric(Total$Hepatocellular.Carcinoma),
                   HSA = as.numeric(Total$Hemangiosarcoma),
                   HS = as.numeric(Total$Histiocytic.Sarcoma),
                   MammACA = as.numeric(Total$Mammary.Gland.Adenocarcinoma),
                   GCT = as.numeric(Total$Granulosa.Cell.Tumor),
                   Thyroid = as.numeric(Total$Thyroid.Tumor),
                   ThyroidAD = as.numeric(Total$Thyroid.Adenoma),
                   STS = as.numeric(Total$Soft.Tissue.Sarcomas),
                   AML = as.numeric(Total$Myeloid.Leukemia),
                   HardACA = as.numeric(Total$Harderian.Gland.Adenocarcinoma),
                   Harderian = as.numeric(Total$Harderian.Tumor),
                   HardAD = as.numeric(Total$Harderian.Gland.Adenoma),
                   LSA.BLL= as.numeric(Total$BLL),
                   LSA.Bmerge= as.numeric(Total$B.merge),
                   LSA.DLBCL= as.numeric(Total$DLBCL),
                   LSA.FBL= as.numeric(Total$FBL),
                   LSA.PreT = as.numeric(Total$PreT),
                   OSA = as.numeric(Total$Osteosarcoma),
                   PitAd = as.numeric(Total$Pituitary.Adenoma),
                   Amyloid = as.numeric(Total$Amyloidosis),
                   NN = as.numeric(Total$non.neoplastic),
                   ectoderm = as.numeric(Total$Ectoderm),
                   endoderm = as.numeric(Total$Endoderm),
                   mesoderm = as.numeric(Total$Mesoderm))
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(row.names(pheno), "sex"))

HZE <- subset(pheno, group == "HZE")
Gamma <- subset(pheno, group == "Gamma")
Un <- subset(pheno, group == "Unirradiated")
All.irr <- subset(pheno, unirradiated == "0")


bootstrap <- HS.assoc.bootstrap(perms = 2, chr = 2, pheno = HZE, pheno.col = "Thyroid",
                                probs, K, addcovar, markers, snp.file, outdir = "~/Desktop/files",
                                tx = "HZE", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/",
                                peakMB = 122584526)


gamma.quant = quantile(gamma$Thyroid, c(0.025,0.975))
print("95% Confidence Interval for QTL:")
print(paste(QUANTILE))



### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
Thy.boot <- read.csv("~/Desktop/QTL Data/Bootstrap.overall-Table 1.csv")

gamma = Thy.boot[which(Thy.boot$TX == "Gamma"),]
hze = Thy.boot[which(Thy.boot$TX == "HZE"),]
un = Thy.boot[which(Thy.boot$TX == "Unirradiated"),]

#HISTOGRAM
layout(matrix(3:1, 3, 1))
hist(un$Thyroid, breaks=150, col="blue", main = "Unirradiated", xlab="Chromosome 2", prob = T, xlim = c(0, 182113224))
lines(density(un$Thyroid), col="black")
hist(gamma$Thyroid, breaks=150, col="green", main = "Gamma", xlab="Chromosome 2", prob = T, xlim = c(0, 182113224))
lines(density(gamma$Thyroid), col="black")
hist(hze$Thyroid, breaks=150, col="red", main = "HZE", xlab="Chromosome 2", prob = T, xlim = c(0, 113282124))
lines(density(hze$Thyroid), col="black", lwd = 1)

#HISTOGRAM WITHOUT DENSITY
layout(matrix(3:1, 3, 1))
hist(un$Thyroid, breaks=25, col="blue", main = "Unirradiated", xlab="Chromosome 2",xlim = c(110000000, 130000000))

hist(gamma$Thyroid, breaks=30, col="green", main = "Gamma", xlab="Chromosome 2",xlim = c(110000000, 130000000))

hist(hze$Thyroid, breaks=20, col="red", main = "HZE", xlab="Chromosome 2",xlim = c(110000000, 130000000))




#KERNAL DENSITY
layout(matrix(3:1, 3, 1))
d = density(hze$Thyroid)
plot(d, col="red", main = "HZE", xlab="Chromosome 2", xlim = c(0, 182113224))
d = density(gamma$Thyroid) 
plot(d, col="green", main = "Gamma", xlab="Chromosome 2", xlim = c(0, 182113224))
d = density(un$Thyroid) 
plot(d, col="blue", main = "Unirradiated", xlab="Chromosome 2", xlim = c(0, 182113224))

hist(hze$Thyroid, prob = T)
lines(density(x))


title(main = "Nonparametric Bootstrap Resampling with Replacement: Distribution of Peak LOD Scores")




thyroid.boot <- factor(Thy.boot$tx, levels = c("Gamma", "HZE", "Unirradiated"),
                       labels = c("Gamma", "HZE", "Unirradiaed"))
polygon(d, col="red", border="blue")


par(mfrow=c(1,1))
sm.density.compare(Thy.boot$Thyroid, Thy.boot$TX, xlab = "Chromosome 2")
title(main = "Thyroid Adenoma: Resample Model Averaging")
colfill = c(2:(2+length(levels(Thy.boot))))
legend("topright", levels(Thy.boot), fill = colfill)

sm.binomial.bootstrap()
