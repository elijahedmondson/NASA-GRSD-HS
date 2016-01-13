#Working with data 

Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Total-Table 1.csv")
str(Total)

LSA <- Total[ which(Total$Lymphoma == 1), ]

correct <- function(x){
  if(x == "gamma")
    return("Gamma")
  if(x == "HZE-Fe")
    return("HZE")
  if(x == "HZE-Si")
    return("HZE")
  if(x == "sham irradiated")
    return("Unirradiated")
  if(x == "unirradiated")
    return("Unirradiated")
  else
    return(NA)
}

LSA$GROUP <- sapply(LSA$group, correct)

LSA <- Total[ which(Total$Lymphoma == 1), ]
table(LSA$LSA.MI)
table(LSA$LSA.B.T)

hist(LSA$LSA.MI)
hist(LSA$LSA.MI, group = "LSA$LSA.Subtype")

plot(LSA$LSA.MI, LSA$days)
plot(LSA$LSA.Grade, LSA$days)

boxplot(formula = days ~ LSA.Grade, data = LSA)
boxplot(formula = days ~ LSA.B.T, data = LSA)
boxplot(formula = days ~ LSA.Subtype, data = LSA)
boxplot(formula = days ~ LSA.Subtype..simple., data = LSA)
boxplot(formula = days ~ LSA.Notch, data = LSA)
boxplot(formula = days ~ coat.color, data = LSA)

#ggplot2

qplot(data = LSA, LSA.MI, days)
qplot(data = LSA, LSA.MI, days, color = LSA.Subtype)
qplot(data = LSA, LSA.MI, days, geom = c("point", "smooth"))

qplot(data = LSA, LSA.MI, fill = LSA.Subtype)
qplot(data = LSA, log(LSA.MI), fill = group)
qplot(data = LSA, LSA.MI, fill = group)

qplot(data = LSA, LSA.MI, facets = LSA.Subtype..simple. ~ ., binwidth = 5)
qplot(data = LSA, LSA.MI, days, facets = .~LSA.Subtype..simple.)
qplot(data = LSA, LSA.MI, facets = group~ ., binwidth = 5)
qplot(data = LSA, LSA.MI, days, facets = .~group)

begin <- Sys.time()



end <- difftime(Sys.time(), begin, units = 'hours')
end
