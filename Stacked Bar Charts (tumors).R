
Data <- ddply(Data, .(Group), 
              transform, pos = cumsum(Tumor) - (0.5 * Tumor)
)


library(plyr)
library(ggplot2)




HZE <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/HZE-Table 1.csv")
stack = data.frame(Gamma = c("Lymphoma" = sum(Gamma$Lymphoma),
                             "Myeloid Leukemia" = sum(Gamma$Myeloid.Leukemia),
                             "Pulmonary Adenocarcinoma" = sum(Gamma$Pulmonary.Adenocarcinoma),
                             "Hepatocellular Carcinoma" = sum(Gamma$Hepatocellular.Carcinoma),
                             "Hemangiosarcoma" = sum(Gamma$Hemangiosarcoma),
                             "Histiocytic Sarcoma" = sum(Gamma$Histiocytic.Sarcoma),
                             "Mammary Adenocarcinoma" = sum(Gamma$Mammary.Gland.Adenocarcinoma),
                             "Ovarian Granulosa Cell Tumor" = sum(Gamma$Granulosa.Cell.Tumor),
                             "Thyroid Adenoma" = sum(Gamma$Thyroid.Tumor),
                             "Soft Tissue Sarcoma" = sum(Gamma$Soft.Tissue.Sarcomas),
                             "Harderian Gland Tumors" = sum(Gamma$Harderian.Tumor),
                             "Osteosarcoma" = sum(Gamma$Osteosarcoma),
                             "Pituitary Adenoma" = sum(Gamma$Pituitary.Adenoma)),
                   HZE = c("Lymphoma" = sum(HZE$Lymphoma),
                           "Myeloid Leukemia" = sum(HZE$Myeloid.Leukemia),
                           "Pulmonary Adenocarcinoma" = sum(HZE$Pulmonary.Adenocarcinoma),
                           "Hepatocellular Carcinoma" = sum(HZE$Hepatocellular.Carcinoma),
                           "Hemangiosarcoma" = sum(HZE$Hemangiosarcoma),
                           "Histiocytic Sarcoma" = sum(HZE$Histiocytic.Sarcoma),
                           "Mammary Adenocarcinoma" = sum(HZE$Mammary.Gland.Adenocarcinoma),
                           "Ovarian Granulosa Cell Tumor" = sum(HZE$Granulosa.Cell.Tumor),
                           "Thyroid Adenoma" = sum(HZE$Thyroid.Tumor),
                           "Soft Tissue Sarcoma" = sum(HZE$Soft.Tissue.Sarcomas),
                           "Harderian Gland Tumor" = sum(HZE$Harderian.Tumor),
                           "Osteosarcoma" = sum(HZE$Osteosarcoma),
                           "Pituitary Adenoma" = sum(HZE$Pituitary.Adenoma)),
                   Unirradiated = c("Lymphoma" = sum(Unirradiated$Lymphoma),
                                    "Myeloid Leukemia" = sum(Unirradiated$Myeloid.Leukemia),
                                    "Pulmonary Adenocarcinoma" = sum(Unirradiated$Pulmonary.Adenocarcinoma),
                                    "Hepatocellular Carcinoma" = sum(Unirradiated$Hepatocellular.Carcinoma),
                                    "Hemangiosarcoma" = sum(Unirradiated$Hemangiosarcoma),
                                    "Histiocytic Sarcoma" = sum(Unirradiated$Histiocytic.Sarcoma),
                                    "Mammary Adenocarcinoma" = sum(Unirradiated$Mammary.Gland.Adenocarcinoma),
                                    "Ovarian Granulosa Cell Tumor" = sum(Unirradiated$Granulosa.Cell.Tumor),
                                    "Thyroid Adenoma" = sum(Unirradiated$Thyroid.Tumor),
                                    "Soft Tissue Sarcoma" = sum(Unirradiated$Soft.Tissue.Sarcomas),
                                    "Harderian Gland Tumor" = sum(Unirradiated$Harderian.Tumor),
                                    "Osteosarcoma" = sum(Unirradiated$Osteosarcoma),
                                    "Pituitary Adenoma" = sum(Unirradiated$Pituitary.Adenoma)))

library(RColorBrewer)


my.col <- c("#ccf2ff", "#9CB071", "#9CB071",
            "#9CB071", "#9CB071", "#87AFC7",
            "#FFCBA4", "#ff9966", "#E6E600FF",
            "#2B547E", "#E8C034FF", "#404040", 
            "#e6e6e6", "#ff6666", "#617C58", 
            "#87CEEB", "#C48793", "#EDC9AF")
barplot(t(stack), col = my.col, ylab = "Number of Mice", ylim = c(0,800), xlim = c(0,12), width = 2)


legend("bottomright", 
       legend = c(colnames(stack)[18:1]), #in order from top to bottom
       fill = my.col[18:1], # 6:1 reorders so legend order matches graph
       title = "Tumor Histotype")







Group <- c(rep(c("Gamma", "HZE", "Unirradiated"), each = 13))
Histotype <- c(rep(c("Lymphoma", "Myeloid Leukemia", 
                     "Pulmonary Adenocarcinoma", "Hepatocellular Carcinoma", 
                     "Hemangiosarcoma", "Histiocytic Sarcoma", 
                     "Mammary Adenocarcinoma",
                     "Ovarian Granulosa Cell Tumor", "Thyroid Adenoma", 
                     "Soft Tissue Sarcoma", "Harderian Gland Tumor",
                     "Osteosarcoma", "Pituitary Adenoma"), times = 3))
Tumor <- c(stack$Gamma, stack$HZE, stack$Unirradiated)

Data <- data.frame(Group, Histotype, Tumor)
Data
Data <- ddply(Data, .(Group), transform, pos = cumsum(Tumor) - (0.5 * Tumor))


Data$Histotype <- factor(Data$Histotype, levels = c("Lymphoma", 
                                                    "Myeloid Leukemia", "Pulmonary Adenocarcinoma", 
                                                    "Hepatocellular Carcinoma", "Hemangiosarcoma", 
                                                    "Histiocytic Sarcoma", "Mammary Adenocarcinoma", 
                                                    "Ovarian Granulosa Cell Tumor", "Thyroid Adenoma", 
                                                    "Soft Tissue Sarcoma", "Harderian Gland Tumor", 
                                                    "Osteosarcoma", "Pituitary Adenoma"))
Data$Histotype <- factor(Data$Histotype, levels = rev(levels(Data$Histotype)))


wa.col = cbind(values=wes_palette(n=4, name="Royal1"),
               values=wes_palette(n=4, name="Zissou"),
               values=wes_palette(n=4, name="GrandBudapest2"),
               values=wes_palette(n=4, name="Royal2"))

ggplot(Data, aes(x = Group, y= Tumor)) + 
       geom_bar(aes(fill = Histotype), stat="identity", colour = "black") +
       geom_text(aes(label = Histotype, y = pos, size = Tumor)) + 
        scale_radius(range = c(2,10)) +
        scale_fill_manual(values = wa.col) +
        ylab("Number of Mice") + xlab("") +
        theme_bw(base_size = 25)


                             


p + geom_text(aes(label = Frequency), size = 3, hjust = 0.5, vjust = 3, position = "stack") 

