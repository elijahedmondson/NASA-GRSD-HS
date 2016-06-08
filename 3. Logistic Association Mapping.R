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
                   NN = as.numeric(Total$non.neoplastic),
                   ectoderm = as.numeric(Total$Ectoderm),
                   endoderm = as.numeric(Total$Endoderm),
                   mesoderm = as.numeric(Total$Mesoderm))
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(row.names(pheno), "sex"))

HZE <- subset(pheno, group == "HZE")
Gamma <- subset(pheno, group == "Gamma")
Un <- subset(pheno, group == "Unirradiated")
All.irr <- subset(pheno, unirradiated == "0")




GRSD.assoc(pheno = Gamma, pheno.col = "LSA.BLL", probs, K, addcovar = addcovar, 
           markers, snp.file = "snp.file", outdir = "~/Desktop/files", tx = "Gamma", 
           sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")



perms <- GRSDassoc.perms(perms = 2, chr = 19, pheno = HZE, Xchr = F, addcovar = addcovar,
                         pheno.col = "HCC", probs = probs, K = K, markers = markers,
                         snp.file = snp.file, outdir = "~/Desktop/files", tx = "Test",
                         sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")

bootstrap <- HS.assoc.bootstrap(perms = 2, chr = 2, pheno = HZE, pheno.col = "Thyroid",
                                probs, K, addcovar, markers, snp.file, outdir = "~/Desktop/files",
                                tx = "HZE", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/",
                                peakMB = 122584526)


addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))


sdp.mat = matrix(as.numeric(intToBits(1:2^8)), nrow = 32)
sdp.mat = sdp.mat[8:1,]
dimnames(sdp.mat) = list(LETTERS[1:8], 1:2^8)

SNP = qtl[[i]]@ranges@start[which.min(qtl[[i]]@elementMetadata@listData$p.value)]
LOD = -log10(qtl[[i]]@elementMetadata@listData$p.value[which(qtl[[i]]@ranges@start == SNP)])


tf = TabixFile(file = sdp.file)
sdps = scanTabix(file = sdp.file, param = GRanges(seqnames = i, ranges = SNP))[[1]]
sdps = strsplit(sdps, split = "\t")
sdps = matrix(unlist(sdps), ncol = 3, byrow = T)
chr  = sdps[1,1]
pos  = as.numeric(sdps[,2])
sdps = as.numeric(sdps[,3])

geno = get.genotype(chr = chr,
                    pos = pos,
                    snp = sdp.mat[,sdps],
                    markers = markers,
                    probs = probs)
