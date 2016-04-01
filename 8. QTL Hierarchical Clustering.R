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


load(file ="~/Desktop/R/QTL/WD/3.\ CoxPH\ Mapping/HZE/Cataract/HZE_Cataract_Latency_CoxPH_plots.Rdata")
qtl.smaller = plot.hs.qtl(data)
save(qtl.smaller, file = "HZE.cataract.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_AML_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.AML.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_HardACA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.HardACA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_HardAD_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.HardAD.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_HCC_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.HCC.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_HS_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.HS.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_HSA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.HSA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_LSA.BLL_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.LSA.BLL.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_LSA.DLBCL_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.DLBCL.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_LSA.FBL_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.FBL.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_LSA.PreT_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.PreT.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_MammACA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.MammACA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_OSA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.OSA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_Pit_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.Pit.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_PulACA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.PulACA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_STS_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.STS.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/HZE_Thyroid_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "HZE.Thyroid.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

##Gamma#############################################
##Gamma#############################################
##Gamma#############################################

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_AML_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.AML.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_GCT_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.GCT.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_HardACA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.HardACA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_HardAD_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.HardAD.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_HCC_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.HCC.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_HS_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.HS.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_HSA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.HSA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_LSA.BLL_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.LSA.BLL.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_LSA.DLBCL_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.LSA.DLBCL.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_LSA.FBL_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.LSA.FBL.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_LSA.PreT_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.LSA.PreT.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_MammACA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.MammACA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_OSA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.OSA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_Pit_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.Pit.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_PulACA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.STS.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Gamma_Thyroid_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Gamma.Thyroid.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)


##Background#############################################
##Background#############################################
##Background#############################################
load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_AML_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.AML.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_HardACA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.HardACA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_HardAD_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.HardAD.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_HCC_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.HCC.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_HS_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.HS.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_HSA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.HSA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_LSA.BLL_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.LSA.BLL.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_LSA.DLBCL_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.LSA.DLBCL.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_LSA.FBL_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.LSA.FBL.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_LSA.PreT_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.LSA.PreT.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_MammACA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.MammACA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_OSA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.OSA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_PulACA_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.PulACA.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_STS_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.STS.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)

load(file ="~/Desktop/R/QTL/WD/2.\ Binary\ Mapping/Unirradiated_Thyroid_GR_QTL.Rdata")
qtl.smaller = plot.hs.qtl(qtl)
save(qtl.smaller, file = "Unirradiated.Thyroid.heatmap.Rdata")
rm(qtl, qtl.smaller, file.prefix)


#Load all files for combination#############################################
#Load all files for combination#############################################
#Load all files for combination#############################################
load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.BHGT.heatmap.Rdata")
HZE.BHGT <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.Ectoderm.heatmap.Rdata")
HZE.Ectoderm <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.Endoderm.heatmap.Rdata")
HZE.Endoderm <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.GCT.heatmap.Rdata")
HZE.GCT <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.HCC.heatmap.Rdata")
HZE.HCC <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.HGAd.heatmap.Rdata")
HZE.HGAd <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.HGT.heatmap.Rdata")
HZE.HGT <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.HS.heatmap.Rdata")
HZE.HS <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.HSA.heatmap.Rdata")
HZE.HSA <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.Intracranial.heatmap.Rdata")
HZE.Intracranial <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.LSA.heatmap.Rdata")
HZE.LSA <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.MalMammAdenocarcinoma.heatmap.Rdata")
HZE.MammAdenocarcinoma <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.MalMammary.heatmap.Rdata")
HZE.MalMammary <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.Mesoderm.heatmap.Rdata")
HZE.Mesoderm <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.Metastatic.heatmap.Rdata")
HZE.Metastatic <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.Myeloid.Leukemia.heatmap.Rdata")
HZE.AML <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.NN.heatmap.Rdata")
HZE.NN <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.OSA.heatmap.Rdata")
HZE.OSA <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.PitAd.heatmap.Rdata")
HZE.PitAd <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.PulmonaryAdenocarcinoma.heatmap.Rdata")
HZE.PulmonaryAdenocarcinoma <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.PulSarcomatoidCarc.heatmap.Rdata")
HZE.PulSarcomatoidCarc <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.RhSA.heatmap.Rdata")
HZE.RhSA <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.STS.heatmap.Rdata")
HZE.STS <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.ThyroidAd.heatmap.Rdata")
HZE.ThyroidAd <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/HZE.ThyTumor.heatmap.Rdata")
HZE.ThyTumor <- qtl.smaller
rm(qtl.smaller)

######################Load GAMMA#######################################
load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.BHGT.heatmap.Rdata")
Gamma.BHGT <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.Ectoderm.heatmap.Rdata")
Gamma.Ectoderm <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.Endoderm.heatmap.Rdata")
Gamma.Endoderm <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.GCT.heatmap.Rdata")
Gamma.GCT <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.HCC.heatmap.Rdata")
Gamma.HCC <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.HGAd.heatmap.Rdata")
Gamma.HGAd <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.HGT.heatmap.Rdata")
Gamma.HGT <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.HS.heatmap.Rdata")
Gamma.HS <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.HSA.heatmap.Rdata")
Gamma.HSA <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.Intracranial.heatmap.Rdata")
Gamma.Intracranial <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.LSA.heatmap.Rdata")
Gamma.LSA <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.MalMammAdenocarcinoma.heatmap.Rdata")
Gamma.MammAdenocarcinoma <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.MalMammary.heatmap.Rdata")
Gamma.MalMammary <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.Mesoderm.heatmap.Rdata")
Gamma.Mesoderm <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.Metastatic.heatmap.Rdata")
Gamma.Metastatic <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.Myeloid.Leukemia.heatmap.Rdata")
Gamma.AML <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.NN.heatmap.Rdata")
Gamma.NN <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.OSA.heatmap.Rdata")
Gamma.OSA <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.PitAd.heatmap.Rdata")
Gamma.PitAd <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.PulmonaryAdenocarcinoma.heatmap.Rdata")
Gamma.PulmonaryAdenocarcinoma <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.PulSarcomatoidCarc.heatmap.Rdata")
Gamma.PulSarcomatoidCarc <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.RhSA.heatmap.Rdata")
Gamma.RhSA <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.STS.heatmap.Rdata")
Gamma.STS <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.ThyroidAd.heatmap.Rdata")
Gamma.ThyroidAd <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Gamma.ThyTumor.heatmap.Rdata")
Gamma.ThyTumor <- qtl.smaller
rm(qtl.smaller)

#########################load BACKGROUND############################
load(file="~/Desktop/R/QTL/WD/Heatmap/Background.BHGT.heatmap.Rdata")
Background.BHGT <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.Ectoderm.heatmap.Rdata")
Background.Ectoderm <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.Endoderm.heatmap.Rdata")
Background.Endoderm <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.HCC.heatmap.Rdata")
Background.HCC <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.HGAd.heatmap.Rdata")
Background.HGAd <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.HGT.heatmap.Rdata")
Background.HGT <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.HS.heatmap.Rdata")
Background.HS <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.HSA.heatmap.Rdata")
Background.HSA <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.Intracranial.heatmap.Rdata")
Background.Intracranial <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.LSA.heatmap.Rdata")
Background.LSA <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.MalMammAdenocarcinoma.heatmap.Rdata")
Background.MammAdenocarcinoma <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.MalMammary.heatmap.Rdata")
Background.MalMammary <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.Mesoderm.heatmap.Rdata")
Background.Mesoderm <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.Metastatic.heatmap.Rdata")
Background.Metastatic <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.Myeloid.Leukemia.heatmap.Rdata")
Background.AML <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.NN.heatmap.Rdata")
Background.NN <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.OSA.heatmap.Rdata")
Background.OSA <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.PitAd.heatmap.Rdata")
Background.PitAd <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.PulmonaryAdenocarcinoma.heatmap.Rdata")
Background.PulmonaryAdenocarcinoma <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.PulSarcomatoidCarc.heatmap.Rdata")
Background.PulSarcomatoidCarc <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.RhSA.heatmap.Rdata")
Background.RhSA <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.STS.heatmap.Rdata")
Background.STS <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.ThyroidAd.heatmap.Rdata")
Background.ThyroidAd <- qtl.smaller
rm(qtl.smaller)

load(file="~/Desktop/R/QTL/WD/Heatmap/Background.ThyTumor.heatmap.Rdata")
Background.ThyTumor <- qtl.smaller
rm(qtl.smaller)

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






