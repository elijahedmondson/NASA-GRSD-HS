## qqman

library(qqman)

qtl1 <- qtl[1]

qq(qtl1@p.value)

manhattan(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
           col = c("gray10", "gray60"), chrlabs = NULL,
           suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
           highlight = NULL, logp = TRUE, ...)
