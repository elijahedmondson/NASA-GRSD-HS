# LOAD PACKAGES #
library(DOQTL)
library(doParallel)
library(foreach)
library(Rsamtools)
library(VariantAnnotation)
library(GenomicRanges)
library(regress)
library(MASS)
library(lmtest)
options(stringsAsFactors = F)
setwd("~/Desktop/R/QTL/WD")
ncl = 4
outdir = "~/Desktop/R/QTL/WD/hq_snps"



# 1. GENOTYPE #

load(file = "~/Desktop/R/QTL/WD/GRSD.Rdata")



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



# 4. DEFINE TRAIT #

rm(HZE, Gamma, Unirradiated)
trait = pheno$PreT
file.prefix = "LSA PreT Gamma"
plot.title = "LSA PreT Gamma, HQ SNPs"
table(trait)
trait
glm(trait ~ addcovar, family = binomial("logit"))
glm(trait ~ addcovar, family = poisson(link = "log"))
glm(trait ~ addcovar, family = gaussian)



# 5. LOGISTIC REGRESSION MODEL #

for(i in 1:length(K)) {
  K[[i]] = K[[i]][samples, samples]
} # for(i)

chrs = c(1:19, "X")
data = vector("list", length(chrs))
names(data) = chrs
for(i in 1:length(chrs)) {
  
  rng = which(markers[,2] == chrs[i])
  data[[i]] = list(probs = probs[,,rng], K = K[[i]],
                   markers = markers[rng,])
  
} # for(i)

rm(probs, K, markers)

setwd(outdir)

# AUTOSOME FUNCTION #
workfxn = function(obj) {
  
  chr = obj$markers[1,2]
  
  setwd(outdir)
  
  strains = sub("/", "_", hs.colors[,2])
  
  hdr = scanVcfHeader(snp.file)
  gr = GRanges(seqnames = chr, range = IRanges(start = 0, 
                                               end = 200e6))
  param = ScanVcfParam(geno = c("GT", "FI"), fixed = "ALT",
                       samples = strains[strains != "C57BL_6J"], which = gr)
  sanger = readVcf(file = snp.file, genome = "mm10", param = param)
  
  # Keep high quality SNPs (quality == 1)
  sanger = sanger[rowSums(geno(sanger)$FI, na.rm = TRUE) == 7]
  
  # Keep polymorphic SNPs.
  keep = which(rowSums(geno(sanger)$GT == "0/0", na.rm = TRUE) < 7)
  sanger = sanger[keep]
  rm(keep)
  
  # We have to do some work to extract the alternate allele.
  alt = CharacterList(fixed(sanger)$ALT)
  alt = unstrsplit(alt, sep = ",")
  
  
  sanger.hdr = data.frame(ID = names(rowRanges(sanger)), CHR = as.character(seqnames(sanger)),
                          POS = start(sanger), REF = as.character(fixed(sanger)$REF),
                          ALT = alt, stringsAsFactors = FALSE)
  rm(alt)
  
  
  sanger = cbind(geno(sanger)$GT[,1:4,drop = FALSE],
                 "C57BL_6J" = "0/0",
                 geno(sanger)$GT[,5:7,drop = FALSE])
  
  sanger = (sanger != "0/0") * 1
  
  # Make the MAF between 1/8 and 4/8.
  flip = which(rowSums(sanger) > 4)
  sanger[flip,] = 1 - sanger[flip,,drop = FALSE]
  rm(flip)
  
  null.mod = glm(trait ~ addcovar, family = binomial(logit))
  #null.mod = glm(trait ~ addcovar, family = poisson(link = "log"))
  null.ll = logLik(null.mod)
  pv = rep(0, nrow(sanger))
  
  glm.fxn = function(snp.rng, local.probs) {
    
    sdp.nums = sanger[snp.rng,] %*% 2^(7:0)
    sdps2keep = which(!duplicated(sdp.nums))
    cur.sdps = sanger[snp.rng,,drop = FALSE][sdps2keep,,drop = FALSE]
    unique.sdp.nums = sdp.nums[sdps2keep]
    m = match(sdp.nums, unique.sdp.nums)
    
    # Multiply the SDPs by the haplotype probabilities.
    cur.alleles = tcrossprod(cur.sdps, local.probs)
    cur.ll = rep(null.ll, nrow(cur.sdps))
    
    # Check for low allele frequencies and remove SDPs with too
    # few samples carrying one allele.
    sdps.to.use = which(rowSums(cur.alleles) > 1.0)
    
    # Run the model at each unique SDP.
    for(j in sdps.to.use) {
      
      
      full.mod = glm(trait ~ addcovar + cur.alleles[j,], family = binomial(logit))
      #full.mod = glm(trait ~ addcovar + cur.alleles[j,], family = poisson(link = "log"))
      cur.ll[j] = logLik(full.mod)
      
    } # for(j)
    
    # This is the LRS.
    cur.ll = cur.ll - null.ll
    
    # Return the results.
    cur.ll[m]
    
  } # glm.fxn()
  
  # SNPs before the first marker.
  snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,1])
    
  } # if(length(snp.rng) > 0)
  
  # SNPs between Markers.
  for(i in 1:(nrow(obj$markers)-1)) {
    
    snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                      sanger.hdr$POS <= obj$markers[i+1,3])
    
    if(length(snp.rng) > 0) {
      
      # Take the mean of the haplotype probs at the surrounding markers.
      pv[snp.rng] = glm.fxn(snp.rng, (obj$probs[,,i] + 
                                          obj$probs[,,i+1]) * 0.5)
      
    } # if(length(snp.rng) > 0)
    
  } # for(i)
  
  # SNPs after the last marker.
  snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])
    
  } # if(length(snp.rng) > 0)
  
  # Convert LRS to p-values using the chi-squared distribution.
  pv = pchisq(2 * pv, df = 1, lower.tail = FALSE)
  pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)
  
  save(pv, file = paste0(file.prefix, "_chr", chr, ".Rdata"))
  
  png(paste0(file.prefix, "_chr", chr,".png"), width = 2000, 
      height = 1600, res = 200)
  plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
  mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
  dev.off()
  
  # Return the positions and p-values.
  return(pv)
  
} # workfxn()


# X FUNCTION #
workfxn.xchr = function(obj) {
  
  chr = obj$markers[1,2]
  
  setwd(outdir)
  
  strains = sub("/", "_", hs.colors[,2])
  
  hdr = scanVcfHeader(snp.file)
  gr = GRanges(seqnames = chr, range = IRanges(start = 0, 
                                               end = 200e6))
  param = ScanVcfParam(geno = c("GT", "FI"), fixed = "ALT",
                       samples = strains[strains != "C57BL_6J"], which = gr)
  sanger = readVcf(file = snp.file, genome = "mm10", param = param)
  
  # Keep high quality SNPs (quality == 1)
  sanger = sanger[rowSums(geno(sanger)$FI, na.rm = TRUE) == 7]
  
  # Keep polymorphic SNPs.
  keep = which(rowSums(geno(sanger)$GT == "0/0", na.rm = TRUE) < 7)
  sanger = sanger[keep]
  rm(keep)
  
  # We have to do some work to extract the alternate allele.
  alt = CharacterList(fixed(sanger)$ALT)
  alt = unstrsplit(alt, sep = ",")
  
  
  sanger.hdr = data.frame(ID = names(rowRanges(sanger)), CHR = as.character(seqnames(sanger)),
                          POS = start(sanger), REF = as.character(fixed(sanger)$REF),
                          ALT = alt, stringsAsFactors = FALSE)
  rm(alt)
  
  
  sanger = cbind(geno(sanger)$GT[,1:4,drop = FALSE],
                 "C57BL_6J" = "0/0",
                 geno(sanger)$GT[,5:7,drop = FALSE])
  
  sanger = (sanger != "0/0") * 1
  
  # Make the MAF between 1/8 and 4/8.
  flip = which(rowSums(sanger) > 4)
  sanger[flip,] = 1 - sanger[flip,,drop = FALSE]
  rm(flip)
  
  null.mod = glm(trait ~ addcovar, family = binomial(logit))
  #null.mod = glm(trait ~ addcovar, family = poisson(link = "log"))
  null.ll = logLik(null.mod)
  pv = rep(0, nrow(sanger))
  
  glm.fxn = function(snp.rng, local.probs) {
    
    sdp.nums = sanger[snp.rng,] %*% 2^(7:0)
    sdps2keep = which(!duplicated(sdp.nums))
    cur.sdps = sanger[snp.rng,,drop = FALSE][sdps2keep,,drop = FALSE]
    unique.sdp.nums = sdp.nums[sdps2keep]
    m = match(sdp.nums, unique.sdp.nums)
    
    # Multiply the SDPs by the haplotype probabilities.
    cur.alleles = tcrossprod(cur.sdps, local.probs)
    cur.ll = rep(null.ll, nrow(cur.sdps))
    
    # Check for low allele frequencies and remove SDPs with too
    # few samples carrying one allele.
    sdps.to.use = which(rowSums(cur.alleles) > 1.0)
    
    sex.col = which(colnames(addcovar) == "sex")
    if(length(sex.col) != 1) {
      stop("One of the columns of addcovar MUST be named 'sex'.")
    } # if(length(sex.col) != 1)
    
    # Run the model at each unique SDP.
    for(j in sdps.to.use) {
      
      
      full.mod = glm(trait ~ addcovar + cur.alleles[j,], family = binomial(logit))
      #full.mod = glm(trait ~ addcovar + cur.alleles[j,], family = poisson(link = "log"))
      cur.ll[j] = logLik(full.mod)
      
    } # for(j)
    
    # This is the LRS.
    cur.ll = cur.ll - null.ll
    
    # Return the results.
    cur.ll[m]
    
  } # glm.fxn()
  
  # SNPs before the first marker.
  snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,1])
    
  } # if(length(snp.rng) > 0)
  
  # SNPs between Markers.
  for(i in 1:(nrow(obj$markers)-1)) {
    
    snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                      sanger.hdr$POS <= obj$markers[i+1,3])
    
    if(length(snp.rng) > 0) {
      
      # Take the mean of the haplotype probs at the surrounding markers.
      pv[snp.rng] = glm.fxn(snp.rng, (obj$probs[,,i] + 
                                        obj$probs[,,i+1]) * 0.5)
      
    } # if(length(snp.rng) > 0)
    
  } # for(i)
  
  # SNPs after the last marker.
  snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])
    
  } # if(length(snp.rng) > 0)
  
  # Convert LRS to p-values using the chi-squared distribution.
  pv = pchisq(2 * pv, df = 1, lower.tail = FALSE)
  pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)
  
  save(pv, file = paste0(file.prefix, "_chr", chr, ".Rdata"))
  
  png(paste0(file.prefix, "_chr", chr,".png"), width = 2600, 
      height = 1200, res = 130)
  plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
  mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
  dev.off()
  
  # Return the positions and p-values.
  return(pv)
  
} # workfxn.xchr()



# 6. MAPPING ANALYSES #

result = vector("list", length(data))
names(result) = names(data)

for(i in 1:19) {
  print(i)
  result[[i]] = workfxn(data[[i]])
} #for(i)

print("X")
result[["X"]] = workfxn.xchr(data[["X"]])

save(result, file = paste0(file.prefix, ".Rdata"))



# 7. FIND THE MAX LOD #

max.LOD = result[[7]]$ID[which(-log10(result[[7]]$pv) > 10)]
max.LOD

max.LOD.position <- result[[7]]$POS[which(-log10(result[[7]]$pv) > 100)]
max.LOD.position

mgi = get.mgi.features(chr = 7, start = 87491804, end = 88689056, type = "gene", source = "MGI")




# 8. PLOTTING #

setwd(outdir)
files = dir(pattern = file.prefix)
files = files[files != paste0(file.prefix, ".Rdata")]
png.files = grep("png$", files)
if(length(png.files) > 0) {
  files = files[-png.files]
}
num = gsub(paste0("^", file.prefix, "_chr|\\.Rdata$"), "", files)
files = files[order(as.numeric(num))]

data = vector("list", length(files))
names(data) = num[order(as.numeric(num))]
for(i in 1:length(files)) {
  
  print(i)
  load(files[i])
  data[[i]] = pv
  data[[i]][,6] = -log10(data[[i]][,6])
  
} # for(i)

num.snps = sapply(data, nrow)
chrs = c(1:19, "X")

xlim = c(0, sum(num.snps))
ylim = c(0, max(sapply(data, function(z) { max(z[,6]) })))



# PLOT ALL CHROMOSOMES #

chrlen = get.chr.lengths()[1:20]
chrsum = cumsum(chrlen)
chrmid = c(1, chrsum[-length(chrsum)]) + chrlen * 0.5
names(chrmid) = names(chrlen)

png(paste0(file.prefix, "_QTL.png"), width = 2600, height = 1200, res = 200)
plot(-1, -1, col = 0, xlim = c(0, max(chrsum)), ylim = ylim, xlab = "",
     ylab = "-log10(p-value)", las = 1, main = plot.title, xaxt = "n")
for(i in 1:length(data)) {
  print(i)
  pos = data[[i]][,3] * 1e-6 + c(0, chrsum)[i]
  points(pos, data[[i]][,6], col = c("black", "grey50")[i %% 2 + 1],
         pch = 20)
} # for(i)
mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.5)
dev.off()


save(result, xlim, ylim, file.prefix, plot.title, file = paste0(file.prefix, "_plotting.Rdata"))




# 9. CONVERT TO GRANGES #

chrs = c(1:19, "X")
qtl = GRangesList(GRanges("list", length(result)))

for(i in 1:length(chrs)) {
  print(i)
  qtl[[i]] <- GRanges(seqnames = Rle(result[[i]]$ID), ranges = IRanges(start = result[[i]]$POS, width = 1), p.value = result[[i]]$pv)
} # for(i)

save(result, file = "~/Desktop/albino.RData")

png("fill .png", width = 2400, height = 1080, res = 100)
DOQTL:::plot.scanone.assoc(qtl, bin.size = 1000)
dev.off()



# 10. SIGNIFICANCE THRESHOLDS #

result = vector("list", length(data))
names(result) = names(data)
females = which(pheno$sex == "0")
males = which(pheno$sex == "1")

perms = matrix(1, nrow = 2, ncol = 2, dimnames = list(1:2, c("A", "X")))

for(p in 1:2) {
  
  new.order = rep(0, length(trait))
  new.order[females] = sample(females)
  new.order[males] = sample(males)
  
  log.perm = trait[new.order]
  trait = log.perm
  
  min.a.pv = 1
  
  for(i in 1:19) {
    print(i)
    result = workfxn(data[[i]])
    min.a.pv = min(min.a.pv, min(result$pv))
  } #for(i)
  
  print("X")
  result = workfxn.xchr(data[["X"]])
  min.x.pv = min(result$pv)
  # Save the minimum p-values.
  perms[p,] = c(-log10(min.a.pv), -log10(min.x.pv))
  
} # for(p)



# 11. TRUNCATED PERMS #

for(p in 1:2) {
  
  new.order = rep(0, length(trait))
  new.order[females] = sample(females)
  new.order[males] = sample(males)
  
  log.perm = trait[new.order]
  trait = log.perm
  
  min.a.pv = 1
  
  
  print(19)
  result = workfxn(data[[19]])
  min.a.pv = min(min.a.pv, min(result$pv))
  
  print("X")
  result = workfxn.xchr(data[["X"]])
  min.x.pv = min(result$pv)
  # Save the minimum p-values.
  perms[p,] = c(-log10(min.a.pv), -log10(min.x.pv))
  
} # for(p)


# 11. SIGNIFICANCE THRESHOLDS # 
get.sig.thr = function(perms, alpha = 0.05, Xchr = TRUE) {
  
  sig.thr = rep(0, length(alpha))
  
  if(Xchr) {
    
    if(!is.matrix(perms)) {
      stop(paste("'perms' is not a matrix. 'perms' must be a matrix",
                 "with 2 columns, named 'A' and 'X'."))
    } # if(!is.matrix(perms))
    
    if(!(all(colnames(perms) %in% c("A", "X")))) {
      stop(paste("The colnames of 'perms' are not equal to 'A' and",
                 "'X'. 'perms' must be a matrix, with 2 columns, named",
                 "'A' and 'X'."))
    } # if(!(all(colnames(perms) %in% c("A", "X"))))
    
    chrlen = get.chr.lengths()
    len.auto = sum(chrlen[1:19])
    len.X = chrlen["X"]
    len.all = len.auto + len.X
    alpha.auto = 1.0 - (1.0 - alpha)^(len.auto / len.all)
    alpha.X    = 1.0 - (1.0 - alpha)^(len.X / len.all)
    
    sig.thr = cbind("A" = quantile(perms[,"A"], probs = 1.0 - alpha.auto, na.rm = TRUE),
                    "X" = quantile(perms[,"X"], probs = 1.0 - alpha.X, na.rm = TRUE))
    rownames(sig.thr) = alpha
    
  } else {
    
    sig.thr = quantile(perms, probs = 1.0 - alpha, na.rm = TRUE)
    names(sig.thr) = alpha
    
  } # else
  
  return(sig.thr)
  
} # get.sig.thr()


#Run function
get.sig.thr(perms, alpha = 0.01, Xchr = TRUE)
