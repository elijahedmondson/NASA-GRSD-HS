library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(DOQTL)
setwd("~/Desktop/")

# Plot function (w/ binning to average markers and max LOD)
plot.hs.qtl = function(qtl, bin.width = 10000, ...) {

  new.qtl = NULL
  for(chr in 1:length(qtl)) {

    print(chr)

    # Create 100 SNP bins.
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
       col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
  mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)

  return(new.qtl)

} # plot.hs.qtl

setwd("~/Desktop/R/QTL/WD/Heatmap/")

##HZE#############################################
##HZE#############################################
##HZE#############################################


load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_AML_GR_QTL.Rdata")
HZE.AML = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_HardACA_GR_QTL.Rdata")
HZE.HardACA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_HardAD_GR_QTL.Rdata")
HZE.HardAD = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_HCC_GR_QTL.Rdata")
HZE.HCC = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_HS_GR_QTL.Rdata")
HZE.HS = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_HSA_GR_QTL.Rdata")
HZE.HSA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_LSA.BLL_GR_QTL.Rdata")
HZE.LSA.BLL = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_LSA.DLBCL_GR_QTL.Rdata")
HZE.DLBCL = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_LSA.FBL_GR_QTL.Rdata")
HZE.FBL = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_LSA.PreT_GR_QTL.Rdata")
HZE.LSA.PreT = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_MammACA_GR_QTL.Rdata")
HZE.MammACA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_OSA_GR_QTL.Rdata")
HZE.OSA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_Pit_GR_QTL.Rdata")
HZE.pit = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_PulACA_GR_QTL.Rdata")
HZE.PulACA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_STS_GR_QTL.Rdata")
HZE.STS = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_Thyroid_GR_QTL.Rdata")
HZE.Thyroid = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

##Gamma#############################################
##Gamma#############################################
##Gamma#############################################

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_AML_GR_QTL.Rdata")
Gamma.AML = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_GCT_GR_QTL.Rdata")
Gamma.GCT = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_HardACA_GR_QTL.Rdata")
Gamma.HardACA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_HardAD_GR_QTL.Rdata")
Gamma.HardAD = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_HCC_GR_QTL.Rdata")
Gamma.HCC = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_HS_GR_QTL.Rdata")
Gamma.HS = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_HSA_GR_QTL.Rdata")
Gamma.HSA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_LSA.BLL_GR_QTL.Rdata")
Gamma.LSA.BLL = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_LSA.DLBCL_GR_QTL.Rdata")
Gamma.LSA.DLBCL = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_LSA.FBL_GR_QTL.Rdata")
Gamma.LSA.FBL = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_LSA.PreT_GR_QTL.Rdata")
Gamma.LSA.PreT = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_MammACA_GR_QTL.Rdata")
Gamma.MammACA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_OSA_GR_QTL.Rdata")
Gamma.OSA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_Pit_GR_QTL.Rdata")
Gamma.Pit = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_PulACA_GR_QTL.Rdata")
Gamma.STS = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_Thyroid_GR_QTL.Rdata")
Gamma.Thyroid = plot.hs.qtl(qtl)
rm(qtl, file.prefix)


##Background#############################################
##Background#############################################
##Background#############################################
load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_AML_GR_QTL.Rdata")
Unirradiated.AML = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_HardACA_GR_QTL.Rdata")
Unirradiated.HardACA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_HardAD_GR_QTL.Rdata")
Unirradiated.HardAD = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_HCC_GR_QTL.Rdata")
Unirradiated.HCC = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_HS_GR_QTL.Rdata")
Unirradiated.HS = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_HSA_GR_QTL.Rdata")
Unirradiated.HSA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_LSA.BLL_GR_QTL.Rdata")
Unirradiated.LSA.BLL = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_LSA.DLBCL_GR_QTL.Rdata")
Unirradiated.LSA.DLBCL = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_LSA.FBL_GR_QTL.Rdata")
Unirradiated.LSA.FBL = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_LSA.PreT_GR_QTL.Rdata")
Unirradiated.LSA.PreT = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_MammACA_GR_QTL.Rdata")
Unirradiated.MammACA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_OSA_GR_QTL.Rdata")
Unirradiated.OSA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_PulACA_GR_QTL.Rdata")
Unirradiated.PulACA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_STS_GR_QTL.Rdata")
Unirradiated.STS = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_Thyroid_GR_QTL.Rdata")
Unirradiated.Thyroid = plot.hs.qtl(qtl)
rm(qtl, file.prefix)



#Combining the columns#############################################
combined <- cbind((-log10(HZE.AML$p.value)/max-log10(HZE.AML$p.value)),
                  -log10(HZE.GCT$p.value),
                  -log10(HZE.HCC$p.value),
                  -log10(HZE.HGT$p.value),
                  -log10(HZE.HS$p.value),
                  -log10(HZE.HSA$p.value),
                  -log10(HZE.LSA$p.value),
                  -log10(HZE.MalMammary$p.value),
                  -log10(HZE.Metastatic$p.value),
                  -log10(HZE.OSA$p.value),
                  -log10(HZE.PitAd$p.value),
                  -log10(HZE.PulmonaryAdenocarcinoma$p.value),
                  -log10(HZE.PulSarcomatoidCarc$p.value),
                  -log10(HZE.RhSA$p.value),
                  -log10(HZE.STS$p.value),
                  -log10(HZE.ThyTumor$p.value),
                  -log10(Gamma.AML$p.value),
                  -log10(Gamma.GCT$p.value),
                  -log10(Gamma.HCC$p.value),
                  -log10(Gamma.HGT$p.value),
                  -log10(Gamma.HS$p.value),
                  -log10(Gamma.HSA$p.value),
                  -log10(Gamma.LSA$p.value),
                  -log10(Gamma.MalMammary$p.value),
                  -log10(Gamma.Metastatic$p.value),
                  -log10(Gamma.OSA$p.value),
                  -log10(Gamma.PitAd$p.value),
                  -log10(Gamma.PulmonaryAdenocarcinoma$p.value),
                  -log10(Gamma.PulSarcomatoidCarc$p.value),
                  -log10(Gamma.RhSA$p.value),
                  -log10(Gamma.STS$p.value),
                  -log10(Gamma.ThyTumor$p.value),
                  -log10(Background.AML$p.value),
                  -log10(Background.HCC$p.value),
                  -log10(Background.HGT$p.value),
                  -log10(Background.HS$p.value),
                  -log10(Background.HSA$p.value),
                  -log10(Background.LSA$p.value),
                  -log10(Background.MalMammary$p.value),
                  -log10(Background.Metastatic$p.value),
                  -log10(Background.OSA$p.value),
                  -log10(Background.PitAd$p.value),
                  -log10(Background.PulmonaryAdenocarcinoma$p.value),
                  -log10(Background.PulSarcomatoidCarc$p.value),
                  -log10(Background.RhSA$p.value),
                  -log10(Background.STS$p.value),
                  -log10(Background.ThyTumor$p.value))

write.csv(combined, file="~/Desktop/csv.csv")

LOD.fnx <- function(x){
        if(x >= "4")
                return(x / max(combined$))
        if(x < "4")
                return(0)
}
combined$AML.gamma <- sapply(combined$AML.gamma, LOD.fnx)



csv <-read.csv(file="~/Desktop/csv.csv", comment.char="#")
rnames <- csv[,1]
mat_data <- data.matrix(csv[,2:ncol(csv)])

##plot#############################################
mypalette <- colorRampPalette(c("green", "yellow", "red"))(n = 299)

heatmap.2(t(mat_data), Colv=NA, labCol=NA, trace ="row", col=mypalette,
          tracecol = "black",
          RowSideColors = c(
            rep("gray", 0),
            rep("blue", 13),
            rep("black", 13)))

par(lend = 1)
legend(.75, 3, legend = c("Unirradiated", "HZE", "Gamma"),
       col = c("gray", "blue", "black"), lty= 1, lwd = 10)


##Plotting 3 QTL maps for comparison##
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

##plot w perms#############################################

perms <- scanone.perm(pheno, pheno.col = "OSA", probs = model.probs,
                      addcovar = addcovar, snps = MM_snps, nperm=5)
thr1 = quantile(perms, probs = 0.90)
thr2 = quantile(perms, probs = 0.95)
thr3 = quantile(perms, probs = 0.99)

plot(AM.qtl, chr = 14, sig.thr = c(thr1, thr2, thr3), main = "")



##### qqman qqman qqman qqman qqman ######
##### qqman qqman qqman qqman qqman ######
##### FOR QTL FROM HEATMap Function ######
##### qqman qqman qqman qqman qqman ######
##### qqman qqman qqman qqman qqman ######

library(qqman)
qtl <- data.frame(CHR = Thyroid.HZE@seqnames,
                  BP = Thyroid.HZE@ranges@start,
                  P = Thyroid.HZE@elementMetadata@listData$p.value)

qtl <- data.frame(CHR = Thyroid.gamma@seqnames,
                  BP = Thyroid.gamma@ranges@start,
                  P = Thyroid.gamma@elementMetadata@listData$p.value)

chr.function <- function(x){
        if(x == "1")
                return(1)
        if(x == "2")
                return(2)
        if(x == "3")
                return(3)
        if(x == "4")
                return(4)
        if(x == "5")
                return(5)
        if(x == "6")
                return(6)
        if(x == "7")
                return(7)
        if(x == "8")
                return(8)
        if(x == "9")
                return(9)
        if(x == "10")
                return(10)
        if(x == "11")
                return(11)
        if(x == "12")
                return(12)
        if(x == "13")
                return(13)
        if(x == "14")
                return(14)
        if(x == "15")
                return(15)
        if(x == "16")
                return(16)
        if(x == "17")
                return(17)
        if(x == "18")
                return(18)
        if(x == "19")
                return(19)
        else
                return(20)
}
qtl$CHR <- sapply(qtl$CHR, chr.function)

manhattan(qtl, main = 'Gamma Thyroid')
qq(qtl$P)






