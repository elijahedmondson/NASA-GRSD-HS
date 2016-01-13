

result[["5"]] = workfxn(data[["5"]])

qtlscan = scanone.assoc(pheno = pheno, pheno.col = "Harderian", probs = probs, K = K, 
                        addcovar = addcovar, markers = markers, sdp.file = sdp.file, ncl = 1)
DOQTL:::plot.scanone.assoc(AM.qtl, 14, bin.size = 10)

perms = scanone.perm(pheno = pheno, pheno.col = "Harderian", probs = probs, addcovar = addcovar, 
                     snps = markers, path = "/Users/elijah/Desktop/R/QTL/WD/", 
                     nperm = 100)
plot(qtlscan, sig.thr = c(thr1, thr2, thr3), main = "")


load(file = "/Users/elijah/Desktop/R/QTL/WD/Association\ Mapping\ Files/HZE/AMQTL.OSA.Rdata")

#Find the max LOD score#

LOD = -log10(data[[2]]$pv)
LOD

max.LOD.SNP.ID <- data$ID[which(-log10(data[[4]]$pv) > 4)]
max.LOD.SNP.ID

max.LOD.position <- result$POS[which(-log10(result[[4]]$pv) > 4)]
max.LOD.position


# PLOT ONE CHROMOSOME #
load(file = "/Users/elijah/Desktop/R/QTL/WD/hq_snps/Thymic\ LSA\ HZE_plotting.Rdata")

data4 <- data[[4]]

png(paste0("ThyLSA_chr4",".png"), width = 2000, 
    height = 1600, res = 200)
plot(as.numeric(data4$POS) * 1e-6, -log10(data4$pv), pch = 20)
mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr"))
dev.off()




mgi = get.mgi.features(chr = 14, start = 68772000, 
                       end = 68774069, type = "gene", source = "MGI")

png(paste0(file.prefix, "_chr", chr,".png"), width = 2000, 
    height = 1600, res = 200)
plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))


# CONVERT TO GRANGES #

load(file = "~/Desktop/albino.Rdata")
chrs = c(1:19, "X")
qtl = GRangesList(GRanges("list", length(result)))

for(i in 1:length(chrs)) {
  print(i)
  qtl[[i]] <- GRanges(seqnames = Rle(result[[i]]$ID), ranges = IRanges(start = result[[i]]$POS, width = 1), p.value = result[[i]]$pv)
} # for(i)

 save(result, file = "~/Desktop/albino.RData")
 

 
