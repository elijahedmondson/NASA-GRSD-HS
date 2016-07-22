heat.geno <- lapply(heat.geno, function(x) {
        gsub("0.0", 5, x)
        })
heat.geno <- lapply(heat.geno, function(x) {
        gsub("0.5", 2.5, x)
})
heat.geno <- lapply(heat.geno, function(x) {
        gsub("1.0", 0, x)
})
heat.geno = as.data.frame(heat.geno)

#Yeast
library(qtl)
library(devtools)
library(qtlyeast)
library(qtlhot, quietly = TRUE)

### READ IN Nadia DATA ###
###    P = "A" = AA    ###
###    M = "B" = BB    ###
###  Het = "H" = AB    ###

moms.pops = read.cross(format="csv", dir="~/Desktop/", file = "geno.csv")
plot.map(moms.pops)

geno.image(moms.pops, reorder = 2, main = "Genotype Data: Strains ordered by phenotype severity")

strains = geno.table(moms.pops, scanone.output = F)
plot(strains, lod = 2)

groupclusteredheatmap(strains)

plotPheno(moms.pops, "pheno")

heat.geno <- read.csv("~/Desktop/heat.genot.csv")
heat.geno = na.omit(heat.geno)
heat.geno = as.data.frame(heat.geno)

combined <- cbind(heat.geno$V1,
                  heat.geno$V2,
                  heat.geno$V3)

combined <- cbind(heat.geno[1:77,])
combined = as.matrix(combined)

drops <- c("X")
heat.geno = heat.geno[ , !(names(heat.geno) %in% drops)]
rm(drops)

require(made4)
heatplot(combined, margins = c(5, 13), dend="col", method = "ave", main = "", key = F, labCol=NA)


heat.geno1 = data.matrix(heat.geno)
heatplot(heat.geno[1:70,],dend="row")



library(gplots)
heatmap.2(t(combined), Colv=NA)


### EXAMPLE DATA SET ###

summary(yeast.orf)
plot.map(yeast.orf)

tf.orf <- "YOL084W"
tf.gene <- yeast.annot[yeast.annot$orf == tf.orf, "gene"]
tf.chr <- yeast.annot[yeast.annot$orf == tf.orf, "chr"]

## Invoke QTL hotspot and causal pair library.

cand <- cis.cand.reg$cis.reg[cis.cand.reg$cis.reg$gene == tf.orf, ]
cand.peak <- cand$peak.pos
cand.pos <- cis.cand.reg$cis.reg$phys.pos[cis.cand.reg$cis.reg$gene == tf.orf]

#image of chr 15
cand.mar <- which.min(abs(cand.pos - pull.map(yeast.orf)[[tf.chr]]))
geno.image(yeast.orf)
abline(v = cand.mar, lwd = 4, col = "green")

#1-D scan of genome
yeast.orf <- calc.genoprob(yeast.orf, step = 2)
scan.orf <- scanone(yeast.orf, pheno.col = find.pheno(yeast.orf, tf.orf), method = "hk")
plot(scan.orf)

#1-D scan focused on 15
plot(scan.orf, chr = tf.chr)
abline(v = cand.pos, lwd = 3, col = "green")
lod.thr <- c(summary(perm.orf, 0.001))
abline(h = lod.thr, lwd = 3, col = "red", lty = 2)
