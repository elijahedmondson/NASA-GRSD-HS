load("/Users/elijah/Desktop/files/ALLIRR_HCC_QTL.Rdata")
chr = 2

#Max LOD score
top <- max(-log10(qtl[[chr]]$p.value))
top

#Position of Max LOD
max.LOD.position <- qtl[[chr]]@ranges[which(-log10(qtl[[chr]]$p.value) == top)]
max.LOD.position

max.LOD.position <- qtl[[chr]]@ranges[which(-log10(qtl[[chr]]$p.value) > 5)]
max.LOD.position

start = max.LOD.position@start[1]
end = max(max.LOD.position@start[])

mgi = get.mgi.features(chr = chr, start = start, end = end, type = "gene", source = "MGI")
print(mgi$Name)



layout(matrix(3:1, 3, 1))
par(mfrow = c(2,2), mar=c(1, 4, 1, 1) + 0.1)
DOQTL:::plot.scanone.assoc(HZE.days, chr=17, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Gamma.days, chr=17, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Unirradiated.days, chr=17, bin.size = 100, main = "Unirradiated", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(qtl, chr=17, bin.size = 100, main = "Total Cases", ylim=c(0,15))

par(mfrow = c(3,1), mar=c(1, 4, 1, 1) + 0.5)
DOQTL:::plot.scanone.assoc(HZE.days, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Gamma.days, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Unirradiated.days, bin.size = 100, main = "Unirradiated", ylim=c(0,15))

par(mfrow = c(3,1), mar=c(1, 4, 1, 1) + 0.5)
DOQTL:::plot.scanone.assoc(HZE.days, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
abline(a = 13, b = 0, col = "red")
DOQTL:::plot.scanone.assoc(Gamma.days, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
abline(a = 13, b = 0, col = "red")
DOQTL:::plot.scanone.assoc(Unirradiated.days, bin.size = 100, main = "Unirradiated", ylim=c(0,15))
abline(a = 13, b = 0, col = "red")
