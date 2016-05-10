# Function to pull all files from a directory folder and determine QTL of a given LOD score

define.QTL(dir = "/Users/elijah/Desktop/files/completed/", threshold = 5.05)

define.QTL = function(dir = "~/Desktop/R/QTL/WD/2.\ Binomial\ Mapping/Rdata/", threshold = 5.8) {
        
        files <- (Sys.glob(paste0(dir,"*.Rdata")))
        
        for(j in 1:length(files)){
                LODmat = matrix(0, nrow = 20, ncol = 6, dimnames = list(1:20, c("chr", "LOD", "peak", "1.5 min", "1.5 max", "interval range")))
                
                load(file = files[j])
                print(files[j])
                
                for(i in 1:19) {
                        top = max(-log10(qtl[[i]]$p.value))
                        if(top > threshold) {
                                LOD.drop.int = top - 1.5
                                max.LOD.peak <- qtl[[i]]@ranges[which(-log10(qtl[[i]]$p.value) > LOD.drop.int)]
                                max.LOD.position <- qtl[[i]]@ranges[which(-log10(qtl[[i]]$p.value) > LOD.drop.int)]
                                print(paste(i, max.LOD.position@start[1], max(max.LOD.position@start), top))
                                LODmat[i,] = c(i, top, ((min(max.LOD.peak@start) + max(max.LOD.peak@start))/2), 
                                               max.LOD.position@start[1], max(max.LOD.position@start), 
                                               (max(max.LOD.position@start) - max.LOD.position@start[1]))
                        }
                } 
                write.csv(LODmat, file = paste0(files[j], "QTL", ".csv"))
                rm(qtl)
        }
}

define.QTL = function(dir = "~/Desktop/R/QTL/WD/2.\ Binomial\ Mapping/Rdata/", threshold = 5.05) {
        
        files <- (Sys.glob(paste0(dir,"*.Rdata")))
        
        for(j in 1:length(files)){
                LODmat = matrix(0, nrow = 20, ncol = 6, dimnames = list(1:20, c("chr", "LOD", "peak", "1.5 min", "1.5 max", "interval range")))
                
                load(file = files[j])
                print(files[j])
                
                for(i in 1:19) {
                        top = max(-log10(qtl[[i]]$p.value))
                        if(top > threshold & top < 5.73) {
                                LOD.drop.int = top - 1.5
                                max.LOD.peak <- qtl[[i]]@ranges[which(-log10(qtl[[i]]$p.value) > LOD.drop.int)]
                                max.LOD.position <- qtl[[i]]@ranges[which(-log10(qtl[[i]]$p.value) > LOD.drop.int)]
                                print(paste(i, max.LOD.position@start[1], max(max.LOD.position@start), top))
                                LODmat[i,] = c(i, top, ((min(max.LOD.peak@start) + max(max.LOD.peak@start))/2), 
                                               max.LOD.position@start[1], max(max.LOD.position@start), 
                                               (max(max.LOD.position@start) - max.LOD.position@start[1]))
                        }
                } 
                write.csv(LODmat, file = paste0(files[j], "QTL", ".csv"))
                rm(qtl)
        }
}