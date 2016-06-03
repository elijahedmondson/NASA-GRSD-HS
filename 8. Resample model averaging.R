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
                   family = as.numeric(Total$family),
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
                   mesoderm = as.numeric(Total$Mesoderm),
                   PSC = as.numeric(Total$Pulmonary.Sarcomatoid.Carcinoma))
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(row.names(pheno), "sex"))
HZE <- subset(pheno, group == "HZE")
Gamma <- subset(pheno, group == "Gamma")
Un <- subset(pheno, group == "Unirradiated")
All.irr <- subset(pheno, unirradiated == "0")



bootstrap <- HS.assoc.bootstrap(perms = 200, chr = 4, pheno = Gamma, pheno.col = "LSA.PreT",
                                probs, K, addcovar, markers, snp.file, outdir = "~/Desktop/files",
                                tx = "Gamma", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/",
                                peakMB = 82878119, window = 8000000)


quant = quantile(thy1$average, c(0.025,0.975))
print("95% Confidence Interval for QTL:")
print(paste(quant))



### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
Thy.boot <- read.csv("~/Desktop/QTL Data/Bootstrap.overall-Table 1.csv")
Thy.boot2 = Thy.boot[which(Thy.boot$AML.LOD > 7),]

gamma = Thy.boot[which(Thy.boot$TX == "gamma"),]
hze = Thy.boot[which(Thy.boot$TX == "HZE"),]
un = Thy.boot[which(Thy.boot$TX == "Unirradiated"),]
allirr = Thy.boot[which(Thy.boot$TX == "All.irradiated"),]

gamma = Thy.boot[which(Thy.boot$TX == "gamma" & Thy.boot$AML.LOD > 5),]
hze = Thy.boot[which(Thy.boot$TX == "HZE" & Thy.boot$AML.LOD > 5),]
un = Thy.boot[which(Thy.boot$TX == "Unirradiated" & Thy.boot$AML.LOD > 5),]
allirr = Thy.boot[which(Thy.boot$TX == "All.irradiated" & Thy.boot$AML.LOD > 5),]



#HISTOGRAM
layout(matrix(3:1, 3, 1))
hist(un$AML.Locus, breaks=150, col="blue", main = "Unirradiated", xlab="Chromosome 2", prob = T, xlim = c(0, 182113224))
lines(density(un$AML.Locus), col="black")
hist(gamma$AML.Locus, breaks=150, col="green", main = "Gamma", xlab="Chromosome 2", prob = T, xlim = c(0, 182113224))
lines(density(gamma$AML.Locus), col="black")
hist(hze$AML.Locus, breaks=150, col="red", main = "HZE", xlab="Chromosome 2", prob = T, xlim = c(0, 182113224))
lines(density(hze$AML.Locus), col="black", lwd = 1)

layout(matrix(4:1, 4, 1))
hist(un$AML.Locus, breaks=150, col="blue", main = "Unirradiated", xlab="Chromosome 2", prob = T, xlim = c(0, 182113224))
lines(density(un$AML.Locus), col="black")
hist(gamma$AML.Locus, breaks=150, col="green", main = "Gamma", xlab="", prob = T, xlim = c(0, 182113224))
lines(density(gamma$AML.Locus), col="black")
hist(hze$AML.Locus, breaks=150, col="red", main = "HZE", xlab="", prob = T, xlim = c(0, 182113224))
lines(density(hze$AML.Locus), col="black", lwd = 1)
hist(allirr$AML.Locus, breaks=150, col="black", main = "All Irradiated", xlab="", prob = T, xlim = c(0, 182113224))
lines(density(All.irradiated$AML.Locus), col="black", lwd = 1)

#HISTOGRAM WITHOUT DENSITY
layout(matrix(3:1, 3, 1))
hist(un$AML.Locus, breaks=150, col="blue", main = "Unirradiated", xlab="Chromosome 2",xlim = c(110000000, 130000000))
hist(gamma$AML.Locus, breaks=150, col="green", main = "Gamma", xlab="Chromosome 2",xlim = c(110000000, 130000000))
hist(hze$AML.Locus, breaks=150, col="red", main = "HZE", xlab="Chromosome 2",xlim = c(110000000, 130000000))




#KERNAL DENSITY
layout(matrix(3:1, 3, 1))
d = density(hze$AML.Locus)
plot(d, col="red", main = "HZE", xlab="Chromosome 2", xlim = c(0, 182113224))
d = density(gamma$AML.Locus)
plot(d, col="green", main = "Gamma", xlab="Chromosome 2", xlim = c(0, 182113224))
d = density(un$AML.Locus)
plot(d, col="blue", main = "Unirradiated", xlab="Chromosome 2", xlim = c(0, 182113224))

layout(matrix(4:1, 4, 1))
d = density(un$AML.Locus)
plot(d, col="blue", main = "Unirradiated", xlab="Chromosome 2", xlim = c(0, 182113224))
d = density(gamma$AML.Locus)
plot(d, col="green", main = "Gamma", xlab="Chromosome 2", xlim = c(0, 182113224))
d = density(hze$AML.Locus)
plot(d, col="red", main = "HZE", xlab="Chromosome 2", xlim = c(0, 182113224))
d = density(allirr$AML.Locus)
plot(d, col="blue", main = "Unirradiated", xlab="Chromosome 2", xlim = c(0, 182113224))

hist(hze$AML.Locus, prob = T)
lines(density(x))


title(main = "Nonparametric Bootstrap Resampling with Replacement: Distribution of Peak LOD Scores")




AML.boot <- factor(Thy.boot$TX, levels = c("All.irradiated", "gamma", "HZE", "Unirradiated"),
                       labels = c("All Irradiated", "Gamma", "HZE", "Unirradiated"))
par(mfrow=c(1,1))
sm.density.compare(Thy.boot$AML.Locus, Thy.boot$TX, xlab = "Chromosome 2", lwd = 2.5)
title(main = "AML Adenoma: Resample Model Averaging")
colfill = c(2:(2+length(levels(AML.boot))))
legend("topright", levels(AML.boot), fill = colfill)


AML.boot <- factor(Thy.boot2$TX, levels = c("All.irradiated", "gamma", "HZE", "Unirradiated"),
                       labels = c("All Irradiated", "Gamma", "HZE", "Unirradiated"))
par(mfrow=c(1,1))
sm.density.compare(Thy.boot2$AML.Locus, Thy.boot2$TX, xlab = "Chromosome 2", lwd = 2.5)
title(main = "AML Adenoma: Resample Model Averaging")
colfill = c(2:(2+length(levels(AML.boot))))
legend("topright", levels(AML.boot), fill = colfill)
#sm.binomial.bootstrap()
