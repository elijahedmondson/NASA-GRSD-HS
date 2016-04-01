library(GenomicRanges)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(made4)
library(BSgenome.Mmusculus.UCSC.mm10)
library(DOQTL)
library(HZE)
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
load(file = "~/Desktop/R/QTL/WD/3.\ CoxPH\ Mapping/HZE/Cataract/HZE_Cataract_Latency_QTL.Rdata")
HZE.cataract = plot.hs.qtl(qtl)
rm(qtl)

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

load(file = "~/Desktop/R/QTL/WD/3.\ CoxPH\ Mapping/Gamma/Cataract/Gamma_Cataract_Latency_QTL.Rdata")
Gamma.cataract = plot.hs.qtl(qtl)
rm(qtl)

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
Gamma.PulACA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_STS_GR_QTL.Rdata")
Gamma.STS = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_Thyroid_GR_QTL.Rdata")
Gamma.Thyroid = plot.hs.qtl(qtl)
rm(qtl, file.prefix)


##Unirradiated#############################################
##Unirradiated#############################################
##Unirradiated#############################################

load(file = "~/Desktop/R/QTL/WD/3.\ CoxPH\ Mapping/Unirradiated/Cataract/Unirradiated_Cataract_Latency_QTL.Rdata")
Unirradiated.cataract = plot.hs.qtl(qtl)
rm(qtl)

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



##Total#############################################
##Total#############################################
##Total#############################################

load(file = "~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Total_AML_QTL.Rdata")
Total.AML = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file = "~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Total_HardACA_QTL.Rdata")
Total.HardACA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file = "~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Total_HCC_QTL.Rdata")
Total.HCC = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file = "~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Total_HSA_QTL.Rdata")
Total.HSA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file = "~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Total_LSA_QTL.Rdata")
Total.LSA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file = "~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Total_MammACA_QTL.Rdata")
Total.MammACA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file = "~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Total_OSA_QTL.Rdata")
Total.OSA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file = "~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Total_Pit_QTL.Rdata")
Total.Pit = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file = "~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Total_PulACA_QTL.Rdata")
Total.PulACA = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file = "~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Total_STS_QTL.Rdata")
Total.STS = plot.hs.qtl(qtl)
rm(qtl, file.prefix)

load(file = "~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Total_Thyroid_QTL.Rdata")
Total.Thyroid = plot.hs.qtl(qtl)
rm(qtl, file.prefix)


#Combining the columns#############################################
#Combining the columns#############################################
#Combining the columns#############################################
#Combining the columns#############################################
#Combining the columns#############################################
#Combining the columns#############################################
combined <- data.frame(HZE.AML = -log10(HZE.AML$p.value),
                  HZE.cataract = -log10(HZE.cataract$p.value),
                  HZE.LSA.DLBCL = -log10(HZE.DLBCL$p.value),
                  HZE.LSA.FBL = -log10(HZE.FBL$p.value),
                  HZE.HardACA = -log10(HZE.HardACA$p.value),
                  HZE.HardAD = -log10(HZE.HardAD$p.value),
                  HZE.HCC = -log10(HZE.HCC$p.value),
                  HZE.HS = -log10(HZE.HS$p.value),
                  #HZE.HSA = -log10(HZE.HSA$p.value),
                  HZE.LSA.BLL= -log10(HZE.LSA.BLL$p.value),
                  HZE.LSA.PreT= -log10(HZE.LSA.PreT$p.value),
                  HZE.MammACA= -log10(HZE.MammACA$p.value),
                  HZE.OSA= -log10(HZE.OSA$p.value),
                  HZE.Pit= -log10(HZE.pit$p.value),
                  HZE.PulACA= -log10(HZE.PulACA$p.value),
                  HZE.STS= -log10(HZE.STS$p.value),
                  HZE.Thyroid= -log10(HZE.Thyroid$p.value),
                  Gamma.AML= -log10(Gamma.AML$p.value),
                  Gamma.cataract= -log10(Gamma.cataract$p.value),
                  Gamma.GCT= -log10(Gamma.GCT$p.value),
                  Gamma.HardACA= -log10(Gamma.HardACA$p.value),
                  Gamma.HardAD= -log10(Gamma.HardAD$p.value),
                  Gamma.HCC= -log10(Gamma.HCC$p.value),
                  Gamma.HS= -log10(Gamma.HS$p.value),
                  #Gamma.HSA= -log10(Gamma.HSA$p.value),
                  Gamma.LSA.BLL= -log10(Gamma.LSA.BLL$p.value),
                  Gamma.LSA.DLBCL= -log10(Gamma.LSA.DLBCL$p.value),
                  Gamma.LSA.FBL= -log10(Gamma.LSA.FBL$p.value),
                  Gamma.LSA.PreT= -log10(Gamma.LSA.PreT$p.value),
                  Gamma.MammACA= -log10(Gamma.MammACA$p.value),
                  Gamma.OSA= -log10(Gamma.OSA$p.value),
                  Gamma.Pit= -log10(Gamma.Pit$p.value),
                  Gamma.PulACA = -log10(Gamma.PulACA$p.value),
                  Gamma.STS= -log10(Gamma.STS$p.value),
                  Gamma.Thyroid= -log10(Gamma.Thyroid$p.value),
                  Unirradiated.AML= -log10(Unirradiated.AML$p.value),
                  Unirradiated.cataract= -log10(Unirradiated.cataract$p.value),
                  Unirradiated.HardACA= -log10(Unirradiated.HardACA$p.value),
                  Unirradiated.HardAD= -log10(Unirradiated.HardAD$p.value),
                  Unirradiated.HCC= -log10(Unirradiated.HCC$p.value),
                  Unirradiated.HS= -log10(Unirradiated.HS$p.value),
                  #Unirradiated.HSA= -log10(Unirradiated.HSA$p.value),
                  Unirradiated.LSA.BLL= -log10(Unirradiated.LSA.BLL$p.value),
                  Unirradiated.LSA.DLBCL= -log10(Unirradiated.LSA.DLBCL$p.value),
                  Unirradiated.LSA.FBL= -log10(Unirradiated.LSA.FBL$p.value),
                  Unirradiated.LSA.PreT= -log10(Unirradiated.LSA.PreT$p.value),
                  Unirradiated.MammACA= -log10(Unirradiated.MammACA$p.value),
                  Unirradiated.OSA= -log10(Unirradiated.OSA$p.value),
                  Unirradiated.PulACA= -log10(Unirradiated.PulACA$p.value),
                  Unirradiated.STS= -log10(Unirradiated.STS$p.value),
                  Unirradiated.Thyroid = -log10(Unirradiated.Thyroid$p.value))


LOD.fnx <- function(x){
        if(x >= "3.5")
                return(x)
        if(x < "3.5")
                return(0)
}
combined$HZE.AML <- sapply(combined$HZE.AML, LOD.fnx)
combined$HZE.cataract <- sapply(combined$HZE.cataract, LOD.fnx)
combined$HZE.LSA.DLBCL <- sapply(combined$HZE.LSA.DLBCL, LOD.fnx)
combined$HZE.LSA.FBL <- sapply(combined$HZE.LSA.FBL, LOD.fnx)
combined$HZE.HardACA <- sapply(combined$HZE.HardACA, LOD.fnx)
combined$HZE.HardAD <- sapply(combined$HZE.HardAD, LOD.fnx)
combined$HZE.HCC <- sapply(combined$HZE.HCC, LOD.fnx)
combined$HZE.HS <- sapply(combined$HZE.HS, LOD.fnx)
combined$HZE.LSA.BLL <- sapply(combined$HZE.LSA.BLL, LOD.fnx)
combined$HZE.LSA.PreT <- sapply(combined$HZE.LSA.PreT, LOD.fnx)
combined$HZE.MammACA <- sapply(combined$HZE.MammACA, LOD.fnx)
combined$HZE.OSA <- sapply(combined$HZE.OSA, LOD.fnx)
combined$HZE.Pit <- sapply(combined$HZE.Pit, LOD.fnx)
combined$HZE.PulACA <- sapply(combined$HZE.PulACA, LOD.fnx)
combined$HZE.STS <- sapply(combined$HZE.STS, LOD.fnx)
combined$HZE.Thyroid <- sapply(combined$HZE.Thyroid, LOD.fnx)
combined$Gamma.AML <- sapply(combined$Gamma.AML, LOD.fnx)
combined$Gamma.cataract <- sapply(combined$Gamma.cataract, LOD.fnx)
combined$Gamma.GCT <- sapply(combined$Gamma.GCT, LOD.fnx)
combined$Gamma.HardACA <- sapply(combined$Gamma.HardACA, LOD.fnx)
combined$Gamma.HardAD <- sapply(combined$Gamma.HardAD, LOD.fnx)
combined$Gamma.HCC <- sapply(combined$Gamma.HCC, LOD.fnx)
combined$Gamma.HS <- sapply(combined$Gamma.HS, LOD.fnx)
combined$Gamma.LSA.BLL <- sapply(combined$Gamma.LSA.BLL, LOD.fnx)
combined$Gamma.LSA.DLBCL <- sapply(combined$Gamma.LSA.DLBCL, LOD.fnx)
combined$Gamma.LSA.FBL <- sapply(combined$Gamma.LSA.FBL, LOD.fnx)
combined$Gamma.LSA.PreT <- sapply(combined$Gamma.LSA.PreT, LOD.fnx)
combined$Gamma.MammACA <- sapply(combined$Gamma.MammACA, LOD.fnx)
combined$Gamma.OSA <- sapply(combined$Gamma.OSA, LOD.fnx)
combined$Gamma.Pit <- sapply(combined$Gamma.Pit, LOD.fnx)
combined$Gamma.PulACA <- sapply(combined$Gamma.PulACA, LOD.fnx)
combined$Gamma.STS <- sapply(combined$Gamma.STS, LOD.fnx)
combined$Gamma.Thyroid <- sapply(combined$Gamma.Thyroid, LOD.fnx)
combined$Unirradiated.AML <- sapply(combined$Unirradiated.AML, LOD.fnx)
combined$Unirradiated.cataract <- sapply(combined$Unirradiated.cataract, LOD.fnx)
combined$Unirradiated.HardACA <- sapply(combined$Unirradiated.HardACA, LOD.fnx)
combined$Unirradiated.HardAD <- sapply(combined$Unirradiated.HardAD, LOD.fnx)
combined$Unirradiated.HCC <- sapply(combined$Unirradiated.HCC, LOD.fnx)
combined$Unirradiated.HS <- sapply(combined$Unirradiated.HS, LOD.fnx)
combined$Unirradiated.HSA <- sapply(combined$Unirradiated.HSA, LOD.fnx)
combined$Unirradiated.LSA.BLL <- sapply(combined$Unirradiated.LSA.BLL, LOD.fnx)
combined$Unirradiated.LSA.DLBCL <- sapply(combined$Unirradiated.LSA.DLBCL, LOD.fnx)
combined$Unirradiated.LSA.FBL <- sapply(combined$Unirradiated.LSA.FBL, LOD.fnx)
combined$Unirradiated.LSA.PreT <- sapply(combined$Unirradiated.LSA.PreT, LOD.fnx)
combined$Unirradiated.MammACA <- sapply(combined$Unirradiated.MammACA, LOD.fnx)
combined$Unirradiated.OSA <- sapply(combined$Unirradiated.OSA, LOD.fnx)
combined$Unirradiated.PulACA <- sapply(combined$Unirradiated.PulACA, LOD.fnx)
combined$Unirradiated.STS <- sapply(combined$Unirradiated.STS, LOD.fnx)
combined$Unirradiated.Thyroid <- sapply(combined$Unirradiated.Thyroid, LOD.fnx)



combined.max.divide <- data.frame(HZE.AML = (combined$HZE.AML/max(combined$HZE.AML)),
                  HZE.cataract = (combined$HZE.cataract/max(combined$HZE.cataract)),
                  HZE.LSA.DLBCL = (combined$HZE.LSA.DLBCL)/max(combined$HZE.LSA.DLBCL),
                  HZE.LSA.FBL = (combined$HZE.LSA.FBL/max(combined$HZE.LSA.FBL)),
                  HZE.HardACA = (combined$HZE.HardACA/max(combined$HZE.HardACA)),
                  HZE.HardAD = (combined$HZE.HardAD/max(combined$HZE.HardAD)),
                  HZE.HCC = (combined$HZE.HCC/max(combined$HZE.HCC)),
                  HZE.HS = (combined$HZE.HS/max(combined$HZE.HS)),
                  #HZE.HSA = (combined$HZE.HSA/max(combined$HZE.HSA)),
                  HZE.LSA.BLL= (combined$HZE.LSA.BLL/max(combined$HZE.LSA.BLL)),
                  HZE.LSA.PreT= (combined$HZE.LSA.PreT/max(combined$HZE.LSA.PreT)),
                  HZE.MammACA= (combined$HZE.MammACA/max(combined$HZE.MammACA)),
                  #HZE.OSA= (combined$HZE.OSA/max(combined$HZE.OSA)),
                  HZE.Pit= (combined$HZE.Pit/max(combined$HZE.Pit)),
                  HZE.PulACA= (combined$HZE.PulACA/max(combined$HZE.PulACA)),
                  HZE.STS= (combined$HZE.STS/max(combined$HZE.STS)),
                  HZE.Thyroid= (combined$HZE.Thyroid/max(combined$HZE.Thyroid)),
                  Gamma.AML= (combined$Gamma.AML/max(combined$Gamma.AML)),
                  Gamma.cataract= (combined$Gamma.cataract/max(combined$Gamma.cataract)),
                  Gamma.GCT= (combined$Gamma.GCT/max(combined$Gamma.GCT)),
                  Gamma.HardACA= (combined$Gamma.HardACA/max(combined$Gamma.HardACA)),
                  Gamma.HardAD= (combined$Gamma.HardAD/max(combined$Gamma.HardAD)),
                  Gamma.HCC= (combined$Gamma.HCC/max(combined$Gamma.HCC)),
                  Gamma.HS= (combined$Gamma.HS/max(combined$Gamma.HS)),
                  #Gamma.HSA= (combined$Gamma.HSA/max(combined$Gamma.HSA)),
                  Gamma.LSA.BLL= (combined$Gamma.LSA.BLL/max(combined$Gamma.LSA.BLL)),
                  Gamma.LSA.DLBCL= (combined$Gamma.LSA.DLBCL/max(combined$Gamma.LSA.DLBCL)),
                  Gamma.LSA.FBL= (combined$Gamma.LSA.FBL/max(combined$Gamma.LSA.FBL)),
                  Gamma.LSA.PreT= (combined$Gamma.LSA.PreT/max(combined$Gamma.LSA.PreT)),
                  Gamma.MammACA= (combined$Gamma.MammACA/max(combined$Gamma.MammACA)),
                  #Gamma.OSA= (combined$Gamma.OSA/max(combined$Gamma.OSA)),
                  Gamma.Pit= (combined$Gamma.Pit/max(combined$Gamma.Pit)),
                  Gamma.PulACA= (combined$Gamma.PulACA/max(combined$Gamma.PulACA)),
                  Gamma.STS= (combined$Gamma.STS/max(combined$Gamma.STS)),
                  Gamma.Thyroid= (combined$Gamma.Thyroid/max(combined$Gamma.Thyroid)),
                  Unirradiated.AML= (combined$Unirradiated.AML/max(combined$Unirradiated.AML)),
                  Unirradiated.cataract= (combined$Unirradiated.cataract/max(combined$Unirradiated.cataract)),
                  Unirradiated.HardACA= (combined$Unirradiated.HardACA/max(combined$Unirradiated.HardACA)),
                  Unirradiated.HardAD= (combined$Unirradiated.HardAD/max(combined$Unirradiated.HardAD)),
                  Unirradiated.HCC= (combined$Unirradiated.HCC/max(combined$Unirradiated.HCC)),
                  Unirradiated.HS= (combined$Unirradiated.HS/max(combined$Unirradiated.HS)),
                  #Unirradiated.HSA= (combined$Unirradiated.HSA/max(combined$Unirradiated.HSA)),
                  Unirradiated.LSA.BLL= (combined$Unirradiated.LSA.BLL/max(combined$Unirradiated.LSA.BLL)),
                  Unirradiated.LSA.DLBCL= (combined$Unirradiated.LSA.DLBCL/max(combined$Unirradiated.LSA.DLBCL)),
                  Unirradiated.LSA.FBL= (combined$Unirradiated.LSA.FBL/max(combined$Unirradiated.LSA.FBL)),
                  Unirradiated.LSA.PreT= (combined$Unirradiated.LSA.PreT/max(combined$Unirradiated.LSA.PreT)),
                  Unirradiated.MammACA= (combined$Unirradiated.MammACA/max(combined$Unirradiated.MammACA)),
                  Unirradiated.OSA= (combined$Unirradiated.OSA/max(combined$Unirradiated.OSA)),
                  Unirradiated.PulACA= (combined$Unirradiated.PulACA/max(combined$Unirradiated.PulACA)),
                  Unirradiated.STS= (combined$Unirradiated.STS/max(combined$Unirradiated.STS)),
                  Unirradiated.Thyroid = (combined$Unirradiated.Thyroid/max(combined$Unirradiated.Thyroid)))
colSums(combined.max.divide)





##################################### made4 ##################################### 
##################################### made4 ##################################### 
##################################### made4 ##################################### 
##################################### made4 ##################################### 
require(made4)

heatplot(t(combined.max.divide), margins = c(10, 10), dend="row", method = "complete", main = "", scaleKey = FALSE)



##################################### heatmap.2##################################### 
##################################### heatmap.2##################################### 
##################################### heatmap.2##################################### 
##################################### heatmap.2##################################### 
##################################### heatmap.2##################################### 


mypalette <- colorRampPalette(c("green", "yellow", "red"))(n = 299)

hclust.ave <- function(x) hclust(x, method="median")
heatmap.2(t(combined.max.divide), Colv=NA, labCol=NA, trace ="row", col=mypalette, key = F,
          tracecol = "black", margins = c(2 , 10), hclust = "median", main = "median", scale = "none",
          RowSideColors = c(rep("gray", 15), rep("blue", 16), rep("black", 14)))



par(lend = 1)
legend(.75, 3, legend = c("Unirradiated", "HZE", "Gamma"),
       col = c("gray", "blue", "black"), lty= 1, lwd = 10)


















##Plotting 3 QTL maps for comparison##
##Plotting 3 QTL maps for comparison##
##Plotting 3 QTL maps for comparison##
##Plotting 3 QTL maps for comparison##
##Plotting 3 QTL maps for comparison##
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



