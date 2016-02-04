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
library(HZE)
options(stringsAsFactors = F)
load(file = "~/Desktop/R/QTL/WD/GRSD.Rdata")
setwd("~/Desktop/files")
outdir = "~/Desktop/files"
sapply(pheno, sum)

rm(HZE, Gamma, Unirradiated, Total, addcovar)
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))


GRSD.assoc(pheno, pheno.col = "AML", probs, K, addcovar, markers, snp.file,
           outdir = "~/Desktop/files", tx = "HZE")

GRSDassoc.perms(perms = 50, chr = 1:2, pheno = pheno, Xchr = F, addcovar = addcovar,
                pheno.col = "Thyroid", probs = probs, K = K, markers = markers,
                snp.file = snp.file, outdir = "~/Desktop/files", tx = "Gamma")

maxLOD(file = "/Users/elijah/Desktop/R/QTL/WD/Binary\ Mapping/Total_AML_QTL.Rdata", chr = 2, LODcutoff = 6)

load(file = "/Users/elijah/Desktop/R/QTL/WD/Binary\ Mapping/Total_AML_QTL.Rdata")
max(-log10(qtl[[15]]$p.value))
max.LOD.position <- qtl$end[which(-log10(qtl$p.value) > 6)]

# 1. GENOTYPE #

load(file = "~/Desktop/R/QTL/WD/GRSD.Rdata")



# 2. PHENOTYPE #

Allirr <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Irradiated-Table 1.csv")
pheno = data.frame(row.names = Allirr$row.names, sex = as.numeric(Allirr$sex == "M"),
                   AML = as.numeric(Allirr),
                   OSA = as.numeric(Allirr$Osteosarcoma))

Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Total-Table 1.csv")
pheno = data.frame(row.names = Alir$row.names, sex = as.numeric(Alir$sex == "M"),
                   AML = as.numeric(Alir$AML.transform),
                   OSA = as.numeric(Alir$Osteosarcoma))


HZE <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/HZE-Table 1.csv")
pheno = data.frame(row.names = HZE$row.names, sex = as.numeric(HZE$sex == "M"),
                   Albino = as.numeric(HZE$albino),
                   PulACA = as.numeric(HZE$Pulmonary.Adenocarcinoma),
                   HCC = as.numeric(HZE$Hepatocellular.Carcinoma),
                   HSA = as.numeric(HZE$Hemangiosarcoma),
                   HS = as.numeric(HZE$Histiocytic.Sarcoma),
                   MammACA = as.numeric(HZE$Mammary.Gland.Adenocarcinoma),
                   #GCT = as.numeric(HZE$Granulosa.Cell.Tumor),
                   Thyroid = as.numeric(HZE$Thyroid.Tumor),
                   ThyroidAD = as.numeric(HZE$Thyroid.Adenoma),
                   STS = as.numeric(HZE$Soft.Tissue.Sarcomas),
                   AML = as.numeric(HZE$AML.transform),
                   Harderian = as.numeric(HZE$Harderian.Tumor),
                   HardACA = as.numeric(HZE$Harderian.Gland.Adenocarcinoma),
                   HardAD = as.numeric(HZE$Harderian.Gland.Adenoma),
                   OSA = as.numeric(HZE$Osteosarcoma),
                   LSA.PreT = as.numeric(HZE$PreT),
                   LSA.BLL = as.numeric(HZE$BLL),
                   LSA.FBL = as.numeric(HZE$FBL),
                   LSA.DLBCL = as.numeric(HZE$DLBCL),
                   LSA.Bmerge = as.numeric(HZE$B.merge),
                   Pit = as.numeric(HZE$Pituitary.Adenoma),
                   Black = as.numeric(HZE$black),
                   LG = as.numeric(HZE$light.grey),
                   LB = as.numeric(HZE$light.brown),
                   Grey = as.numeric(HZE$grey),
                   Brown = as.numeric(HZE$creme.brown),
                   PulMets = as.numeric(HZE$Pulmonary.Metastases))

Gamma <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Gamma-Table 1.csv")
pheno = data.frame(row.names = Gamma$row.names, sex = as.numeric(Gamma$sex == "M"),
                   Albino = as.numeric(Gamma$albino),
                   PulACA = as.numeric(Gamma$Pulmonary.Adenocarcinoma),
                   HCC = as.numeric(Gamma$Hepatocellular.Carcinoma),
                   HSA = as.numeric(Gamma$Hemangiosarcoma),
                   HS = as.numeric(Gamma$Histiocytic.Sarcoma),
                   MammACA = as.numeric(Gamma$Mammary.Gland.Adenocarcinoma),
                   GCT = as.numeric(Gamma$Granulosa.Cell.Tumor),
                   Thyroid = as.numeric(Gamma$Thyroid.Tumor),
                   ThyroidAD = as.numeric(Gamma$Thyroid.Adenoma),
                   STS = as.numeric(Gamma$Soft.Tissue.Sarcomas),
                   AML = as.numeric(Gamma$Myeloid.Leukemia),
                   Harderian = as.numeric(Gamma$Harderian.Tumor),
                   HardACA = as.numeric(Gamma$Harderian.Gland.Adenocarcinoma),
                   HardAD = as.numeric(Gamma$Harderian.Gland.Adenoma),
                   OSA = as.numeric(Gamma$Osteosarcoma),
                   LSA.PreT = as.numeric(Gamma$PreT),
                   LSA.BLL = as.numeric(Gamma$BLL),
                   LSA.FBL = as.numeric(Gamma$FBL),
                   LSA.DLBCL = as.numeric(Gamma$DLBCL),
                   LSA.Bmerge = as.numeric(Gamma$B.merge),
                   Pit = as.numeric(Gamma$Pituitary.Adenoma),
                   LG = as.numeric(Gamma$light.grey),
                   LB = as.numeric(Gamma$light.brown),
                   Grey = as.numeric(Gamma$grey),
                   Brown = as.numeric(Gamma$creme.brown))

Unirradiated <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Unirradiated-Table 1.csv")
pheno = data.frame(row.names = Unirradiated$row.names, sex = as.numeric(Unirradiated$sex == "M"),
                   Albino = as.numeric(Unirradiated$albino),
                   PulACA = as.numeric(Unirradiated$Pulmonary.Adenocarcinoma),
                   HCC = as.numeric(Unirradiated$Hepatocellular.Carcinoma),
                   HSA = as.numeric(Unirradiated$Hemangiosarcoma),
                   HS = as.numeric(Unirradiated$Histiocytic.Sarcoma),
                   MammACA = as.numeric(Unirradiated$Mammary.Gland.Adenocarcinoma),
                   #GCT = as.numeric(Unirradiated$Granulosa.Cell.Tumor),
                   Thyroid = as.numeric(Unirradiated$Thyroid.Tumor),
                   ThyroidAD = as.numeric(Unirradiated$Thyroid.Adenoma),
                   STS = as.numeric(Unirradiated$Soft.Tissue.Sarcomas),
                   AML = as.numeric(Unirradiated$Myeloid.Leukemia),
                   Harderian = as.numeric(Unirradiated$Harderian.Tumor),
                   HardACA = as.numeric(Unirradiated$Harderian.Gland.Adenocarcinoma),
                   HardAD = as.numeric(Unirradiated$Harderian.Gland.Adenoma),
                   OSA = as.numeric(Unirradiated$Osteosarcoma),
                   LSA.PreT = as.numeric(Unirradiated$PreT),
                   LSA.BLL = as.numeric(Unirradiated$BLL),
                   LSA.FBL = as.numeric(Unirradiated$FBL),
                   LSA.DLBCL = as.numeric(Unirradiated$DLBCL),
                   LSA.Bmerge = as.numeric(Unirradiated$B.merge),
                   Black = as.numeric(Unirradiated$black),
                   LG = as.numeric(Unirradiated$light.grey),
                   LB = as.numeric(Unirradiated$light.brown),
                   Grey = as.numeric(Unirradiated$grey),
                   Brown = as.numeric(Unirradiated$creme.brown),
                   Pit = as.numeric(Unirradiaed$Pituitary.adenoma))

sapply(pheno, sum)
rm(HZE, Gamma, Unirradiated)

# 3. COVARIATES #

addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))




# FIND THE MAX LOD #

max.LOD = result[[17]]$ID[which(-log10(result[[17]]$pv) > 4)]
max.LOD

max.LOD.position <- result[[16]]$POS[which(-log10(result[[16]]$pv) > 4)]
max.LOD.position

mgi = get.mgi.features(chr = 17, start = 81341699, end = 81433448,
                       type = "gene", source = "MGI")









###### EMBEDDED WITHIN GRSD.assoc() #########


samples = intersect(rownames(pheno), rownames(probs))
samples = intersect(samples, rownames(addcovar))
samples = intersect(samples, rownames(K[[1]]))
stopifnot(length(samples) > 0)
print(paste("Mapping with", length(samples), "samples."))

pheno = pheno[samples,,drop = FALSE]
addcovar = addcovar[samples,,drop = FALSE]
probs = probs[samples,,,drop = FALSE]



# 4. DEFINE TRAIT #

rm(HZE, Gamma, Unirradiated)
#trait <- obj$pheno[,obj$pheno.col]
trait <- pheno$Albino
table(trait)
pheno.col <- "Albino"
file.prefix <- paste(pheno.col, "HZE")
plot.title <- paste(pheno.col, "HZE")
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



# 6. MAPPING ANALYSES #

result = vector("list", length(data))
names(result) = names(data)

for(i in 1:19) {
        print(i)
        result[[19]] = workfxn(data[[19]])
} #for(i)

print("X")
result[["X"]] = GRSDbinom.xchr(data[["X"]])

save(result, file = paste0(file.prefix, ".Rdata"))



# 7. FIND THE MAX LOD #

max.LOD = result[[17]]$ID[which(-log10(result[[17]]$pv) > 4)]
max.LOD

max.LOD.position <- result[[16]]$POS[which(-log10(result[[16]]$pv) > 4)]
max.LOD.position

mgi = get.mgi.features(chr = 17, start = 81341699, end = 81433448,
                       type = "gene", source = "MGI")




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


save(result, xlim, ylim, file.prefix, plot.title,
     file = paste0(file.prefix, "_plotting.Rdata"))




# 9. CONVERT TO GRANGES #

chrs = c(1:19, "X")
qtl = GRangesList(GRanges("list", length(result)))

for(i in 1:length(chrs)) {
        print(i)
        qtl[[i]] <- GRanges(seqnames = Rle(result[[i]]$ID),
                            ranges = IRanges(start = result[[i]]$POS, width = 1),
                            p.value = result[[i]]$pv)
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

perms = matrix(1, nrow = 10, ncol = 2, dimnames = list(1:10, c("A", "X")))

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

                sig.thr = cbind("A" = quantile(perms[,"A"],
                                               probs = 1.0 - alpha.auto, na.rm = TRUE),
                                "X" = quantile(perms[,"X"],
                                               probs = 1.0 - alpha.X, na.rm = TRUE))
                rownames(sig.thr) = alpha

        } else {

                sig.thr = quantile(perms, probs = 1.0 - alpha, na.rm = TRUE)
                names(sig.thr) = alpha

        } # else

        return(sig.thr)

} # get.sig.thr()


#Run function
get.sig.thr(perms, alpha = 0.01, Xchr = TRUE)















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





plot.hs = function(qtl, title, bin.width = 1000, ...) {
        library(GenomicRanges)
        library(BSgenome.Mmusculus.UCSC.mm10)
        new.qtl = NULL
        for(chr in 1:length(qtl)) {

                print(chr)

                # Create SNP bins with given bin.width
                brks = cut(x = 1:length(qtl[[chr]]), breaks = length(qtl[[chr]]) / bin.width)
                # Split up the SNP positions and get the mean.
                pos = split(start(qtl[[chr]]), brks)
                pos = sapply(pos, mean)
                # Split up the p-values and get the max.
                pv = split(mcols(qtl[[chr]])$p.value, brks)
                pv = sapply(pv, min)

                # Make a single new GRanges object to return.
                gr = GRanges(seqnames = seqnames(qtl[[chr]])[1],
                             ranges = IRanges(start = pos, width = 1), p.value = pv)

                if(chr == 1) {
                        new.qtl = gr
                } else {
                        new.qtl = c(new.qtl, gr)
                } # else

        } # for(chr)

        # Get the chromosome lengths.
        chrlen = seqlengths(BSgenome.Mmusculus.UCSC.mm10)
        names(chrlen) = sub("^chr", "", names(chrlen))
        chrlen = chrlen[seqlevels(new.qtl)] * 1e-6

        # Add the chr lengths to the chromosomes for plotting.
        # Switch positions to genome Mb.
        gmb = start(new.qtl) * 1e-6
        for(chr in 2:length(chrlen)) {

                wh = which(seqnames(new.qtl) == names(chrlen)[chr])
                gmb[wh] = gmb[wh] + sum(chrlen[1:(chr - 1)])

        } # for(chr)

        # Get chromosome mid-points for plotting the Chr name.
        chrmid = (chrlen / 2) + cumsum(c(1, chrlen[-length(chrlen)]))

        # Make the plot.
        col = rep(rgb(0,0,0), length(new.qtl))
        even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
        col[even.chr] = rgb(0.7,0.7,0.7)
        plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
             col = col, las = 1, xlab = "", ylab = "-log10(p-value)", main = title)
        mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)

        return(new.qtl)

} #plot.hs
