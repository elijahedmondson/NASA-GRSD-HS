# Get the SNP at the minimum p-value.
max.snp = assoc[which.min(assoc[,12]),]
geno = get.genotype(chr = max.snp[1,1], pos = max.snp[1,2], snp = max.snp[1,4:11],
                    markers = snps, probs = probs)

# Fit the model.
mod = lm(pheno[,11] ~ covar + geno[,1], na.action = na.exclude)

# Get RSS / SST for the genotypes.
anova(mod)[2,2] / sum(anova(mod)[,2])