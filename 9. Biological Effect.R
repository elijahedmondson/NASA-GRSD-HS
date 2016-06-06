
# LOAD PACKAGES #
library(BSgenome.Mmusculus.UCSC.mm10)
library(doParallel)
library(foreach)
library(Rsamtools)
library(VariantAnnotation)
library(DOQTL)
library(GenomicRanges)
library(survival)
library(regress)
library(HZE)
outdir = "~/Desktop/files/"
options(stringsAsFactors = F)
setwd("~/Desktop/files/")
load(file = "~/Desktop/R/QTL/WD/GRSD.Rdata")


# PHENOTYPE #
Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = Total$row.names, sex = as.numeric(Total$sex == "M"),
                   cohort = as.numeric(Total$Cohort),
                   Unirradiated = as.numeric(Total$Unirradiated),
                   family = as.numeric(Total$family),
                   group = as.character(Total$groups),
                   days = as.numeric(Total$days),
                   NN = as.numeric(Total$non.neoplastic),
                   Tumor = as.numeric(Total$Tumor),
                   days2 = as.numeric(Total$Cataract.2.0.Score.Days),
                   cat2 = as.numeric(Total$Cataract.2.0.Score.Event),
                   days3 = as.numeric(Total$Cataract.3.0.Score.Days),
                   cat3 = as.numeric(Total$Cataract.3.0.Score.Event),
                   days4 = as.numeric(Total$Cataract.4.0.Score.Days),
                   cat4 = as.numeric(Total$Cataract.4.0.Score.Event),
                   pigdisp = as.numeric(Total$pigment.dispersion),
                   dilate = as.numeric(Total$Did.Not.Dilate),
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
                   THB.merge = as.numeric(Total$Thyroid.HCC.Bmerge))
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))
pheno["survival"] = rep(1, 1820)
HZE = subset(pheno, group == "HZE")
Gamma = subset(pheno, group == "Gamma")
Unirradiated = subset(pheno, group == "Unirradiated")
All.irr = subset(pheno, Unirradiated == 0)


get.effect.size(pheno, pheno.col = "Thyroid", chr = 2, qtl, probs, , markers,
                sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz",
                dir = "/Users/elijah/Desktop/R/QTL/WD/3.\ CoxPH\ Mapping/"")







