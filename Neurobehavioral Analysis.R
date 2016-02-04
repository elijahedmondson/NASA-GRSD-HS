###   Neurobehavioral Analysis   ###

# LOAD PHENO FILES #
neuro.pheno600 <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/neuro.pheno600.csv")
neuro.pheno1200 <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/neuro.pheno1200.csv")
neuro.pheno <- rbind(neuro.pheno600, neuro.pheno1200)
library(dplyr)
rename(neuro.pheno, c("X"="row.names"))
names(neuro.pheno)[names(neuro.pheno)=="X"] <- "row.names"

# QC ON PHENO TRAITS #
Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Total-Table 1.csv")

x <- merge(x = Neuro.EFE, y = Total, by.x = "row.names", by.y = "row.names")


total <- left_join(Total, neuro.pheno, by.x = "row.names", by.y = "row.names")

write.csv(neuro.pheno, file = "~/Desktop/neuropheno.csv")
