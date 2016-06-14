# LOAD PACKAGES #
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
options(stringsAsFactors = F)
load(file = "~/Desktop/R/QTL/WD/GRSD.Rdata")
load("/Users/elijah/Desktop/R/QTL/WD/hs.colors.Rdata")
setwd("~/Desktop/files")
outdir = "~/Desktop/files"

Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = Total$row.names, rownames = Total$row.names,
                   sex = as.numeric(Total$sex == "M"),
                   LG = as.numeric(Total$light.grey),
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
                   PulMet.transform = as.numeric(Total$PulMet.transform),
                   Metastatic.Tumors = as.numeric(Total$Metastatic.Tumors),
                   Pulmonary.Metastases = as.numeric(Total$Pulmonary.Metastases),
                   HCC.Metastatic.Density = as.numeric(Total$HCC.Metastatic.Density),
                   Tumors.that.could.met = as.numeric(Total$Tumors.that.could.met),
                   Pulmonary.Metastases = as.numeric(Total$Pulmonary.Metastases))
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(row.names(pheno), "sex"))

HZE <- subset(pheno, group == "HZE")
Gamma <- subset(pheno, group == "Gamma")
Un <- subset(pheno, group == "Unirradiated")
All.irr <- subset(pheno, unirradiated == "0")

HCC.met <- Total[which(Total$Hepatocellular.Carcinoma=="1" & Total$HCC.Metastatic.Density>0),]
HCC.met <- Total[which(Total$Hepatocellular.Carcinoma=="1"), ]
pheno = data.frame(row.names = HCC.met$row.names, sex = as.numeric(HCC.met$sex == "M"),
                   HCC.met = as.numeric(HCC.met$HCC.Met),
                   OSA = as.numeric(HCC.met$Osteosarcoma))

PulMET <- pheno[which(pheno$Tumors.that.could.met=="1"), ]

PSC <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/PSC.csv")
pheno = data.frame(row.names = PSC$row.names, rownames = PSC$row.names,
                   sex = as.numeric(PSC$Sex == "M"),
                   PSC = as.numeric(PSC$Sarcomatoid.Score))
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(row.names(pheno), "sex"))



GRSD.assoc(pheno = PulMET, pheno.col = "Pulmonary.Metastases", probs, K, addcovar = addcovar,
           markers, snp.file = "snp.file", outdir = "~/Desktop/files", tx = "All",
           sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")

GRSD.poisson(pheno = pheno, pheno.col = "PSC", probs, K, addcovar = addcovar,
             markers, snp.file = "snp.file", outdir = "~/Desktop/files", tx = "All",
             sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")

perms <- GRSDassoc.perms(perms = 2, chr = 19, pheno = HZE, Xchr = F, addcovar = addcovar,
                         pheno.col = "HCC", probs = probs, K = K, markers = markers,
                         snp.file = snp.file, outdir = "~/Desktop/files", tx = "Test",
                         sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")

bootstrap <- HS.assoc.bootstrap(perms = 2, chr = 2, pheno = HZE, pheno.col = "Thyroid",
                                probs, K, addcovar, markers, snp.file, outdir = "~/Desktop/files",
                                tx = "HZE", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/",
                                peakMB = 122584526)



