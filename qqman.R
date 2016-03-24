## qqman

library(qqman)
load(file = "/Users/elijah/Desktop/R/QTL/WD/1.\ Linear\ Mapping/NEURO.gamma.QTL.Rdata")

qtl <- QTL.C1.train.frz
        
DOQTL:::plot.scanone.assoc(QTL.C1.train.frz, bin.size = 100, main = "HZE train_deltapctfrze_isi1_isi4: Cohort 1")

qtl1 <- qtl[1]



#x = data.frame("BP," "CHR," "P," and optionally, "SNP)

manhattan(new19, chr = "CHR", bp = "BP", p = "P",
           col = c("gray10", "gray60"), chrlabs = NULL,
           suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
           highlight = NULL, logp = TRUE)
qqman


# Look at a Quantile-Quantile plof of the p-values.
pv = unlist(sapply(qtl, function(z) { z$p.value }))
qqnorm(-log10(pv[sample(1:length(pv), 50000)]))
qqline(-log10(pv[sample(1:length(pv), 50000)]))
# This looks weird. Are there other covariates (like batch) that
# should be in your model? Or is the data full of 0 values with
# only a few non-zero values?


x <- Sys.time()
qq(new6$P)
Sys.time() - x






new1 <- data.frame(BP = qtl$`1`@ranges@start, CHR = as.numeric(1), P = qtl$`1`@elementMetadata$p.value)
new2 <- data.frame(BP = qtl$`2`@ranges@start, CHR = as.numeric(2), P = qtl$`2`@elementMetadata$p.value)
new3 <- data.frame(BP = qtl$`3`@ranges@start, CHR = as.numeric(3), P = qtl$`3`@elementMetadata$p.value)
new4 <- data.frame(BP = qtl$`4`@ranges@start, CHR = as.numeric(4), P = qtl$`4`@elementMetadata$p.value)
new5 <- data.frame(BP = qtl$`5`@ranges@start, CHR = as.numeric(5), P = qtl$`5`@elementMetadata$p.value)
new6 <- data.frame(BP = qtl$`6`@ranges@start, CHR = as.numeric(6), P = qtl$`6`@elementMetadata$p.value)
new7 <- data.frame(BP = qtl$`7`@ranges@start, CHR = as.numeric(7), P = qtl$`7`@elementMetadata$p.value)
new8 <- data.frame(BP = qtl$`8`@ranges@start, CHR = as.numeric(8), P = qtl$`8`@elementMetadata$p.value)
new9 <- data.frame(BP = qtl$`9`@ranges@start, CHR = as.numeric(9), P = qtl$`9`@elementMetadata$p.value)
new10 <- data.frame(BP = qtl$`10`@ranges@start, CHR = as.numeric(10), P = qtl$`10`@elementMetadata$p.value)
new11 <- data.frame(BP = qtl$`11`@ranges@start, CHR = as.numeric(11), P = qtl$`11`@elementMetadata$p.value)
new12 <- data.frame(BP = qtl$`12`@ranges@start, CHR = as.numeric(12), P = qtl$`12`@elementMetadata$p.value)
new13 <- data.frame(BP = qtl$`13`@ranges@start, CHR = as.numeric(13), P = qtl$`13`@elementMetadata$p.value)
new14 <- data.frame(BP = qtl$`14`@ranges@start, CHR = as.numeric(14), P = qtl$`14`@elementMetadata$p.value)
new15 <- data.frame(BP = qtl$`15`@ranges@start, CHR = as.numeric(15), P = qtl$`15`@elementMetadata$p.value)
new16 <- data.frame(BP = qtl$`16`@ranges@start, CHR = as.numeric(16), P = qtl$`16`@elementMetadata$p.value)
new17 <- data.frame(BP = qtl$`17`@ranges@start, CHR = as.numeric(17), P = qtl$`17`@elementMetadata$p.value)
new18 <- data.frame(BP = qtl$`18`@ranges@start, CHR = as.numeric(18), P = qtl$`18`@elementMetadata$p.value)
new19 <- data.frame(BP = qtl$`19`@ranges@start, CHR = as.numeric(19), P = qtl$`19`@elementMetadata$p.value)
new20 <- data.frame(BP = qtl$`X`@ranges@start, CHR = as.numeric(20), P = qtl$`X`@elementMetadata$p.value)

require(plyr)
x <- join_all(list(new1, new2, new3, new4, new5, new6, new7, new8, new9, new10, new11, new12, new13, new14,
           new15, new16, new17, new18, new19, new20), by = "BP", type = "full")
rm(new1, new2, new3, new4, new5, new6, new7, new8, new9, new10, new11, new12, new13, new14,
   new15, new16, new17, new18, new19, new20)


