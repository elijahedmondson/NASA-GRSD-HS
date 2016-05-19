


# Get the SNP at the minimum p-value.
max.snp = markers$SNP_ID[which(markers$Mb_NCBI38 == 3010274)]


# Read in the unique SDPs.

sdp.file = "/Users/elijah/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz"

tf = TabixFile(sdp.file)

sdps = scanTabix(file = sdp.file, param = GRanges(seqnames = obj$markers[1,2], ranges = IRanges(start = 0, end = 200e6)))[[1]]

sdps = strsplit(sdps, split = "\t")
sdps = matrix(unlist(sdps), ncol = 3, byrow = T)
chr = sdps[1,1]
pos = as.numeric(sdps[,2])
sdps = as.numeric(sdps[,3])
# Create a matrix of SDPs.
sdp.mat = matrix(as.numeric(intToBits(1:2^8)), nrow = 32)
sdp.mat = sdp.mat[8:1,]
dimnames(sdp.mat) = list(LETTERS[1:8], 1:2^8)

# get.genotype

geno = get.genotype(chr = 2, pos = max.snp, snp = max.snp,
                    markers, probs = probs)


# Fit the model.
mod = lm(pheno[,11] ~ covar + geno[,1], na.action = na.exclude)

# Get RSS / SST for the genotypes.
anova(mod)[2,2] / sum(anova(mod)[,2])