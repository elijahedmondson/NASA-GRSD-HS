setwd("/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE")
load(file = "/Users/elijahedmondson/Desktop/R/QTL/WD/GRSD_master.Rdata")
pheno = data.frame(row.names = HZE$row.names, 
                   sex = as.numeric(HZE$sex == "M"),  
                   OSA = as.numeric(HZE$Osteosarcoma))
covar = data.frame(sex = as.numeric(HZE$sex == "M"))
addcovar = covar
rownames(covar) = rownames(pheno)
rownames(addcovar) = rownames(pheno)

pheno = data.frame(row.names = Untitled$row.names, 
                   sex = as.numeric(Untitled$sex == "M"),  
                   weight = as.numeric(Untitled$weight))

AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HGT", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HGT.Rdata")
rm(AM.qtl)
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Thyroid", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.ThyTumor.Rdata")
rm(AM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "LSA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.LSA.Rdata")
plot(LM.qtl)
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "LSA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.LSA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "NN", probs = model.probs, K = K, 
                addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.NN.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "NN", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.NN.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "PulACA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.PulmonaryAdenocarcinoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "PulACA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.PulmonaryAdenocarcinoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "PulAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.PulAdenoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "PulAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.PulAdenoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "PulSrc", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.PulSarcomatoidCarc.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "PulSrc", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.PulSarcomatoidCarc.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "BHGT", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.BHGT.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "BHGT", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.BHGT.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "CPT", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.CPT.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "CPT", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.CPT.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Dpolyp", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Dpolyp,Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Dpolyp", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Dpolyp.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "EndoPolyp", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.EndoPolyp.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "EndoPolyp", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.EndoPolyp.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Ependymoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Ependymoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Ependymoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Ependymoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "FSA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.FSA.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "FSA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.FSA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "GastSCC", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.GastSCC.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "GastSCC", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.GastSCC.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "GCT", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.GCT.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "GCT", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.GCT.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "HGAca", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.HGAca.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HGAca", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HGAca.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "HGAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.HGAd.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HGAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HGAd.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Hemangioma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Hemangioma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Hemangioma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Hemangioma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "HSA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.HSA.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HSA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HSA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Hepatoblastoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Hepatoblastoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Hepatoblastoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Hepatoblastoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "HCC", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.HCC.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HCC", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HCC.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "HCAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.HCAd.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HCAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HCAd.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "HS", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.HS.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "HS", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.HS.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "IntestinalACA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.IntestinalACA.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "IntestinalACA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.IntestinalACA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "IntestinalNeo", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.IntestinalNeo.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "IntestinalNeo", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.IntestinalNeo.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Intracranial", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Intracranial.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Intracranial", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Intracranial.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Islet", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Islet.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Islet", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Islet.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Myeloid.Leukemia", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Myeloid.Leukemia.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Myeloid.Leukemia", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Myeloid.Leukemia.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Myxosarcoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Myxosarcoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Myxosarcoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Myxosarcoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Odontogenic", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Odontogenic.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Odontogenic", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Odontogenic.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Osteoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Osteoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Osteoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Osteoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "OSA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.OSA.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "OSA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.OSA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Leiomyosarcoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Leiomyosarcoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Leiomyosarcoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Leiomyosarcoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "MalMammary", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.MalMammary.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "MalMammary", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.MalMammary.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "MalOvarian", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.MalOvarian.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "MalOvarian", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.MalOvarian.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "MammAdenoacanthoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.MammAdenoacanthoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "MammAdenoacanthoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.MammAdenoacanthoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "MammAdenocarcinoma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.MammAdenocarcinoma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "MammAdenocarcinoma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.MammAdenocarcinoma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Meningioma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Meningioma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Meningioma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Meningioma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Mesothelioma", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Mesothelioma.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Mesothelioma", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Mesothelioma.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Met", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Met.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Met", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Met.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Metastatic", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Metastatic.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Metastatic", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Metastatic.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Ectoderm", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Ectoderm.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Ectoderm", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Ectoderm.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Endoderm", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Endoderm.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Endoderm", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Endoderm.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Mesoderm", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Mesoderm.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Mesoderm", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Mesoderm.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "OvarianCarc", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.OvarianCarc.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "OvarianCarc", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.OvarianCarc.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Pheo", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Pheo.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Pheo", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Pheo.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "PitAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.PitAd.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "PitAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.PitAd.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "PulAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.PulAd.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "PulAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.PulAd.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "PulMetastasis", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.PulMetastasis.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "PulMetastasis", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.PulMetastasis.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "RCC", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.RCC.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "RCC", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.RCC.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "RhSA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.RhSA.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "RhSA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.RhSA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "STS", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.STS.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "STS", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.STS.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "ThyroidAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.ThyroidAd.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "ThyroidAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.ThyroidAd.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "ThyroidCarc", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.ThyroidCarc.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "ThyroidCarc", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.ThyroidCarc.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "TubulostromalACA", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.TubulostromalACA.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "TubulostromalACA", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.TubulostromalACA.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "TubulostromalAd", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.TubulostromalAd.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "TubulostromalAd", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.TubulostromalAd.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "UndiffSarc", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.UndiffSarc.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "UndiffSarc", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.UndiffSarc.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "UterineStromalSarc", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.UterineStromalSarc.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "UterineStromalSarc", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.UterineStromalSarc.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "UnilHarderian", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.UnilHarderian.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "UnilHarderian", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.UnilHarderian.Rdata")
rm(AM.qtl, LM.qtl)

LM.qtl = scanone(pheno = pheno, pheno.col = "Epidermal", probs = model.probs, K = K, 
                 addcovar = covar, snps = MM_snps)
save(LM.qtl, file = "LMQTL.Epidermal.Rdata")
AM.qtl = scanone.assoc(pheno = pheno, pheno.col = "Epidermal", probs = model.probs, K = K, 
                       addcovar = covar, markers = MM_snps, cross = "HS", 
                       sdp.file = sdp.file, ncl = 2)
save(AM.qtl, file = "AMQTL.Epidermal.Rdata")
rm(AM.qtl, LM.qtl)

setwd("/Users/elijahedmondson/Desktop/R/QTL/WD/HZE.Plots")

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.BHGT.Rdata")
jpeg('Harderian Gland Tumors, Bilateral (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Harderian Gland Tumors, Bilateral (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.CPT.Rdata")
jpeg('Choroid Plexus Tumor (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Choroid Plexus Tumor (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Dpolyp.Rdata")
jpeg('Duodenal Polyp (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Duodenal Polyp (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Ectoderm.Rdata")
jpeg('Ectoderm Derived Tumor (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Ectoderm Derived Tumor (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Endoderm.Rdata")
jpeg('Endoderm Derived Tumor (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Endoderm Derived Tumor (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.EndoPolyp.Rdata")
jpeg('Endometrial Stromal Polyp (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Endometrial Stromal Polyp (Uirradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Ependymoma.Rdata")
jpeg('Ependymoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Ependymoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Epidermal.Rdata")
jpeg('Epidermis Derived Tumors (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Epidermis Derived Tumors (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.FSA.Rdata")
jpeg('Fibrosarcoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Fibrosarcoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.GastSCC.Rdata")
jpeg('Gastric Squamous Cell Carcinoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Gastric Squamous Cell Carcinoma (HZE Ion Irradiated)")
dev.off()

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.GCT.Rdata")
jpeg('Ovarian Granulosa Cell Tumor (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Ovarian Granulosa Cell Tumor (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.HCAd.Rdata")
jpeg('Hepatocellular Adenoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Hepatocellular Adenoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.HCC.Rdata")
jpeg('Hepatocellular Carcinoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Hepatocellular Carcinoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Hemangioma.Rdata")
jpeg('Hemangioma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Hemangioma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Hepatoblastoma.Rdata")
jpeg('Hepatoblastoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Hepatoblastoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.HGAca.Rdata")
jpeg('Harderian Gland Adenocarcinoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Harderian Gland Adenocarcinoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.HGAd.Rdata")
jpeg('Harderian Gland Adenoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Harderian Gland Adenoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.HS.Rdata")
jpeg('Histiocytic Sarcoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Histiocytic Sarcoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.HSA.Rdata")
jpeg('Hemangiosarcoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Hemangiosarcoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.IntestinalACA.Rdata")
jpeg('Intestinal Adenocarcinoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Intestinal Adenocarcinoma (HZE Ion Irradiated)")
dev.off()

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.IntestinalNeo.Rdata")
jpeg('Intestinal Neoplasm (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Intestinal Neoplasm (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Intracranial.Rdata")
jpeg('Intracranial Neoplasm (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Intracranial Neoplasm (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Islet.Rdata")
jpeg('Pancreatic Islet Cell Tumor (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Pancreatic Islet Cell Tumor (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Leiomyosarcoma.Rdata")
jpeg('Leiomyosarcoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Leiomyosarcoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.LSA.Rdata")
jpeg('Lymphoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Lymphoma (HZE Ion Irradiated)")
dev.off()

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.MalMammary.Rdata")
jpeg('Malignant Mammary Tumor (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Malignant Mammary Tumor (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.MalOvarian.Rdata")
jpeg('Malignant Ovarian Tumor (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Malignant Ovarian Tumor (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.MammAdenoacanthoma.Rdata")
jpeg('Mammary Adenoacanthoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Mammary Adenoacanthoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.MammAdenocarcinoma.Rdata")
jpeg('Mammary Adenocarcinoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Mammary Adenocarcinoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Meningioma.Rdata")
jpeg('Meningioma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Meningioma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Mesoderm.Rdata")
jpeg('Mesoderm Derived Tumor (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Mesoderm Derived Tumor (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Mesothelioma.Rdata")
jpeg('Mesothelioma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Mesothelioma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Met.Rdata")
jpeg('Metastatic Disease (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Metastatic Disease (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Myeloid.Leukemia.Rdata")
jpeg('Myeloid Leukemia (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Myeloid Leukemia (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Myxosarcoma.Rdata")
jpeg('Myxosarcoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Myxosarcoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.NN.Rdata")
jpeg('Non-neoplastic (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Non-neoplastic (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Odontogenic.Rdata")
jpeg('Odotogenic Tumors (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Odontogenic Tumors (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.OSA.Rdata")
jpeg('Osteosarcoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Osteosarcoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Osteoma.Rdata")
jpeg('Osteoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Osteoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.OvarianCarc.Rdata")
jpeg('Ovarian Carcinoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Ovarian Carcinoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.Pheo.Rdata")
jpeg('Pheochromocytoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Pheochromocytoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.PitAd.Rdata")
jpeg('Pituitary Adenoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Pituitary Adenoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.PulAd.Rdata")
jpeg('Pulmonary Adenoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Pulmonary Adenoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.PulAdenoma.Rdata")
jpeg('Pulmonary Adenoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Pulmonary Adenoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.PulMetastasis.Rdata")
jpeg('Pulmonary Metastases (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Pulmonary Metastases (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.PulmonaryAdenocarcinoma.Rdata")
jpeg('Pulmonary Adenocarcinoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Pulmonary Adenocarcinoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.PulSarcomatoidCarc.Rdata")
jpeg('Pulmonary Sarcomatoid Carcinoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Pulmonary Sarcomatoid Carcinoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.RCC.Rdata")
jpeg('Renal Cell Carcinoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Renal Cell Carcinoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.RhSA.Rdata")
jpeg('Rhabdomyosarcoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Rhabdomyosarcoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.STS.Rdata")
jpeg('Soft Tissue Sarcoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Soft Tissue Sarcoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.ThyroidAd.Rdata")
jpeg('Thyroid Follicular Adenoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Thyroid Follicular Adenoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.ThyroidCarc.Rdata")
jpeg('Thyroid Follicular Carcinoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Thyroid Follicular Carcinoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.TubulostromalACA.Rdata")
jpeg('Ovarian Tubulostroma Adenocarcinoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Ovarian Tubulostroma Adenocarcinoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.TubulostromalAd.Rdata")
jpeg('Ovarian Tubulostromal Adenoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Ovarian Tubulostromal Adenoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.UndiffSarc.Rdata")
jpeg('Undifferentiated Sarcoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Undifferentiated Sarcoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.UnilHarderian.Rdata")
jpeg('Harderian Gland Tumor, Unilateral (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Harderian Gland Tumor, Unilateral (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)

load(file ="/Users/elijahedmondson/Desktop/R/QTL/WD/Mapping Files/HZE/AMQTL.UterineStromalSarc.Rdata")
jpeg('Uterine Stromal Sarcoma (HZE Ion Irradiated)', width = 1600, height = 800, res = 80)
DOQTL:::plot.scanone.assoc(AM.qtl, bin.size = 100, main = "Uterine Stromal Sarcoma (HZE Ion Irradiated)")
dev.off()
rm(AM.qtl)