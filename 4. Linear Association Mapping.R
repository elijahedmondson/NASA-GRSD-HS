library(DOQTL)



# 1. GENOTYPE #

load(file = "~/Desktop/R/Build/K.Rdata")
load(file = "~/Desktop/R/Build/model.probs.Rdata")
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
outdir = "~/Desktop/R/QTL/WD/"
setwd(outdir)
sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz"



# 2. PHENOTYPE #

HZE <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/HZE-Table 1.csv")
pheno = data.frame(row.names = HZE$row.names, sex = as.numeric(HZE$sex == "M"),  
                   Albino = as.numeric(HZE$albino),
                   Pulmonary.Adenocarcinoma = as.numeric(HZE$Pulmonary.Adenocarcinoma),
                   HCC = as.numeric(HZE$Hepatocellular.Carcinoma),
                   LSA = as.numeric(HZE$Lymphoma),
                   AML = as.numeric(HZE$Myeloid.Leukemia),
                   Harderian = as.numeric(HZE$Harderian.Tumor),
                   OSA = as.numeric(HZE$Osteosarcoma))

Gamma <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Gamma-Table 1.csv")
pheno = data.frame(row.names = Gamma$row.names, sex = as.numeric(Gamma$sex == "M"),  
                   Albino = as.numeric(Gamma$albino),
                   Pulmonary.Adenocarcinoma = as.numeric(Gamma$Pulmonary.Adenocarcinoma),
                   HCC = as.numeric(Gamma$Hepatocellular.Carcinoma),
                   LSA = as.numeric(Gamma$LSA),
                   AML = as.numeric(Gamma$Myeloid.Leukemia),
                   PreT = as.numeric(Gamma$PreT),
                   FBL = as.numeric(Gamma$FBL))

Unirradiated <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Unirradiated-Table 1.csv")
pheno = data.frame(row.names = Unirradiated$row.names, sex = as.numeric(Unirradiated$sex == "M"),  
                   Albino = as.numeric(Unirradiated$albino),
                   Pulmonary.Adenocarcinoma = as.numeric(Unirradiated$Pulmonary.Adenocarcinoma),
                   HCC = as.numeric(Unirradiated$Hepatocellular.Carcinoma),
                   LSA = as.numeric(Unirradiated$Lymphoma),
                   AML = as.numeric(Unirradiated$Myeloid.Leukemia),
                   Harderian = as.numeric(Unirradiated$Harderian.by.gland))



# 3. COVARIATES #

addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))
samples = intersect(rownames(pheno), rownames(probs))
samples = intersect(samples, rownames(addcovar))
samples = intersect(samples, rownames(K[[1]]))
stopifnot(length(samples) > 0)
print(paste("Found", length(samples), "samples in common."))

pheno = pheno[samples,,drop = FALSE]
addcovar = addcovar[samples,,drop = FALSE]
probs = probs[samples,,,drop = FALSE]



# 4. ASSOCIATION MAPPING #

qtl = scanone.assoc(pheno = pheno, pheno.col = "PreT", probs = model.probs, K = K, 
                        addcovar = addcovar, markers = MM_snps, sdp.file = sdp.file, ncl = 4)

save(qtl, file = "__.Rdata")

DOQTL:::plot.scanone.assoc(qtl, bin.size = 100)

png("PreT.Gamma.LSA_QTL.png", width = 2400, height = 1080, res = 100)
DOQTL:::plot.scanone.assoc(qtl, bin.size = 100)
dev.off()

png("_CHR_4.png", width = 2000, height = 1600, res = 128)
DOQTL:::plot.scanone.assoc(qtl, chr = 4, bin.size = 100)
dev.off()




# 5. LINKAGE MAPPING #

qtl = scanone(pheno = pheno, pheno.col = "sarcomatoid", probs = model.probs, K = K, 
              addcovar = covar, snps = MM_snps)
plot(qtl, main = "PSC")

perms = scanone.perm(pheno = pheno, pheno.col = "PreT", probs = model.probs, addcovar = addcovar, 
                     snps = MM_snps, path = "~/Desktop/", 
                     nperm = 2)

thr1 = quantile(perms, probs = 0.90)
thr2 = quantile(perms, probs = 0.95)
thr3 = quantile(perms, probs = 0.99)

plot(qtl, sig.thr = c(thr1, thr2, thr3), main = "PSC")

setwd("/Users/elijahedmondson/Desktop/R/QTL/WD")
getwd()
save(perms, file = "sarcomatoid.3000perms.Rdata")
save(qtl, file = "sarcomatoid.3000perms.Rdata")
list.files("/Users/elijahedmondson/Desktop/R/QTL/WD")

interval = bayesint(qtlscan, chr = 1)
interval
mgi = get.mgi.features(chr = interval[1,2], start = interval[1,3], 
                       end = interval[3,3], type = "gene", source = "MGI")
nrow(mgi)
head(mgi)

ma = assoc.map(pheno = pheno, pheno.col = "PreT", probs = model.probs, K = K, addcovar = covar, 
               snps = MM_snps, chr = interval[1,2], start = interval[1,3], end = interval[3,3])
coefplot(qtl, chr = 1)
tmp = assoc.plot(ma, thr = 1)
unique(tmp$sdps)




