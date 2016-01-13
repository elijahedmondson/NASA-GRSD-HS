#Loading Packages and Data
library(DOQTL)
library(doParallel)
library(foreach)
library(Rsamtools)
library(VariantAnnotation)
library(GenomicRanges)
library(survival)
library(regress)
sample <- read.csv(file = "/Users/elijahedmondson/Desktop/sample.csv")
load(file = "/Users/elijahedmondson/Desktop/GRSD_master.Rdata")
cross = "HS"
options(stringsAsFactors = F)
snp.file = "/Users/elijahedmondson/Desktop/R/Build/mgp.v4.snps.dbSNP.vcf.gz"
sdp.file = "/Users/elijahedmondson/Desktop/R/Build/HS_Sanger_SDPs.txt.bgz"
setwd("/Users/elijahedmondson/Desktop/R/QTL/WD")
pheno = data.frame(row.names = sample$row.names, sex = as.numeric(sample$sex == "M"),  
                   albino = as.numeric(sample$albino),
                   days = as.numeric(sample$days))
covar = data.frame(sex = as.numeric(sample$sex == "M"))
addcovar = covar
rownames(covar) = rownames(pheno)
rownames(addcovar) = rownames(pheno)

#Linkage Mapping
LM.qtl = scanone(pheno = pheno, pheno.col = "days", probs = model.probs, K = K, 
              addcovar = covar, snps = MM_snps)
plot(LM.qtl, main = "test")
save(LM.qtl, file = "")
perms = scanone.perm(pheno = pheno, pheno.col = "days", probs = model.probs, addcovar = covar, 
                     snps = MM_snps, path = "/Users/elijahedmondson/Desktop/R/QTL/perms", 
                     nperm = 10)
thr1 = quantile(perms, probs = 0.90)
thr2 = quantile(perms, probs = 0.95)
thr3 = quantile(perms, probs = 0.99)
plot(qtl, sig.thr = c(thr1, thr2, thr3), main = "")
interval = bayesint(qtl, chr = 7)
interval
mgi = get.mgi.features(chr = interval[1,2], start = interval[1,3], end = interval[3,3], 
                       type = "gene", source = "MGI")
nrow(mgi)
head(mgi)

#Association Mapping
pheno[,2] = as.numeric(pheno[,2]) - 1

AM.qtl = scanone.assoc(pheno = pheno, pheno.col = 2, probs = probs, K = K, 
                    addcovar = covar, markers = markers, cross = "HS", sdp.file = sdp.file, ncl = 2)
plot(AM.qtl, main = "test")
save(AM.qtl, file = "")

png("albino_QTL.png", width = 2000, height = 1600, res = 128)
DOQTL:::plot.scanone.assoc(qtl, bin.size = 100)
dev.off()

png("albino_QTL_chr7.png", width = 2000, height = 1600, res = 128)
DOQTL:::plot.scanone.assoc(qtl, chr = 7, bin.size = 10)
dev.off()







