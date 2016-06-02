
pheno.col = "Thyroid"

# Get the SNP at the minimum p-value.
max.snp = qtl@unlistData[which.min(qtl@unlistData$p.value)]


# Read in the unique SDPs.
sdp.file = "/Users/elijah/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz"
tf = TabixFile(sdp.file)
sdps = scanTabix(file = sdp.file, param = GRanges(seqnames = max.snp@seqnames@values, ranges = max.snp@ranges@start))[[1]]
sdps = strsplit(sdps, split = "\t")
sdps = matrix(unlist(sdps), ncol = 3, byrow = T)
chr  = sdps[1,1]
pos  = as.numeric(sdps[,2])
sdps = as.numeric(sdps[,3])



# Create a matrix of SDPs.
sdp.mat = matrix(as.numeric(intToBits(1:2^8)), nrow = 32)
sdp.mat = sdp.mat[8:1,]
dimnames(sdp.mat) = list(LETTERS[1:8], 1:2^8)

                
sdp.mat[,8]


snp = sdp.mat[,8]
snp = t(snp)
row.names(snp) = hs.colors[,2]



# Get the SNP at the minimum p-value.


max.snp = qtl@unlistData[which.min(qtl@unlistData$p.value)]


geno = get.genotype(chr = chr, 
                    pos = pos, 
                    snp = sdp.mat[,sdps], 
                    markers = markers, 
                    probs = probs)

# Fit the model.
samples = intersect(rownames(pheno), rownames(probs))
samples = intersect(samples, rownames(addcovar))
samples = intersect(samples, rownames(geno))
stopifnot(length(samples) > 0)
pheno = pheno[samples,,drop = FALSE]
geno = geno[samples,,drop = FALSE]
addcovar = addcovar[samples,,drop = FALSE]
probs = probs[samples,,,drop = FALSE]


mod0 = glm(pheno[,pheno.col] ~ addcovar, family = binomial(logit))
summary(anova(mod0))
mod1 = glm(pheno[,pheno.col] ~ addcovar + geno[,1], family = binomial(logit))
summary(anova(mod1))

mod = lm(pheno[,pheno.col] ~ addcovar + geno[,1])

# Get RSS / SST for the genotypes.
anova(mod)[2,2] / sum(anova(mod)[,2])

#perc.var = 100 * (1.0 - (ss / ss.null))




################################################################################
# Given the haplotype probabilities and the location of a single SNP,
# use the haplotype probabilities and the SNP to impute the allele calls in 
# DO mice.
# Daniel Gatti
# Dan.Gatti@jax.org
# Dec. 7, 2015
################################################################################
# Arguments: chr: a single character (or number) that is the chromosome.
#            pos: a positive value that is the SNP location in bp of Mb.
#            snp: The strain pattern of the SNP as 0s and 1s. Must be a named
#                 vector with the strain names.
#            markers: The location of the markers in probs.
#            probs: 3D array containing the haplotype probabilities.


get.genotype = function(chr, pos, snp, markers, probs) {
        
        # Convert the SNP to numbers.
        snp = unlist(snp)
        names(snp) = make.names(sub("_", ".", names(snp)))
        strains = make.names(hs.colors[,2])
        
        #snp = snp[strains]
        #snp = as.numeric(factor(snp)) - 1
        
        
        # Get the slices from the haplotype probs matrix.
        markers = markers[markers[,1] %in% dimnames(probs)[[3]],]
        probs = probs[,,dimnames(probs)[[3]] %in% markers[,1]]
        markers = markers[markers[,2] == chr,]
        probs = probs[,,markers[,1]]
        markers = markers[max(which(markers[,3] < pos)):min(which(markers[,3] > pos)),]
        
        # Get the probs for these markers.
        probs = probs[,,markers[,1], drop = FALSE]
        probs = apply(probs, 1:2, mean)
        
        # Multiply the two matrices and return the result.
        return(probs %*% snp)
        
} # get.genotype()

