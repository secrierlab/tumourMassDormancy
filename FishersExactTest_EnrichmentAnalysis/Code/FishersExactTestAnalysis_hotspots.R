#########################################################################
###Fishers exact test between samples with and without TMD
###########################################################################

#Load required packages:
library(TCGAbiolinks)
library(stringr)
library(dplyr)
library(maditr)
library(ggplot2)


###############################################
#Download mutation data for cosmic genes only 
###############################################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/GeneLists/")
COSMIC <- read.table("Census_allMon Feb 24 15_27_31 2020.tsv", header = TRUE,sep = "\t")
cosmic_list <- as.character(COSMIC$Gene.Symbol)

##Select solid cancer types:
solid_tumours <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

####Download maf data:
setwd("~/Documents/GitHub/tumourMassDormancy/Data/MutationData/")
genes_to_save <- c("CASP8","HRAS","NRAS","KRAS") #Downloading genes from programmes and from COSMIC

for (tumour in solid_tumours){
  print(paste0(tumour, ": This is ", which(tumour == solid_tumours), " out of ", length(solid_tumours)))    #Show progress
  df.maf <- GDCquery_Maf(tumour, pipelines = "muse")
  df.maf.short <- df.maf[,c("Tumor_Sample_Barcode",
                            "Hugo_Symbol","Chromosome",
                            "Start_Position","End_Position",
                            "Variant_Classification",
                            "Variant_Type",
                            "Reference_Allele",
                            "Tumor_Seq_Allele1",
                            "Tumor_Seq_Allele2",
                            "One_Consequence",
                            "Consequence",
                            "SIFT","PolyPhen","HGVSc","HGVSp","HGVSp_Short")]
  df.maf.short$Cancer <- as.character(tumour)
  df.maf.short <- subset(df.maf.short, Hugo_Symbol %in% genes_to_save)
  
  if (tumour == solid_tumours[1]){
    maf_all <- df.maf.short
  } else {
    maf_all <- rbind(maf_all, df.maf.short)
  }
}

save(maf_all, file = "maf_cosmic_hotspots.RData")



######################################################
########Load the mutation data and process further:
######################################################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/MutationData/")
load("maf_cosmic_hotspots.RData")

maf_all$Patient_short <- sapply(maf_all$Tumor_Sample_Barcode,
                                function(x) paste(str_split(x, "-", simplify = T)[1:4], collapse = "-"))

maf_all.lenient <- maf_all[ which(maf_all$Variant_Classification
                                  %in% c("Missense_Mutation",
                                         "Nonsense_Mutation",
                                         "Nonstop_Mutation",
                                         "Frameshift_Deletion",
                                         "Frameshift_Insertion",
                                         "Inframe_Insertion",
                                         "Inframe_Deletion")),
]
rm(maf_all)


###########################
#CASP8 mutations:
############################
CASP8_maf <- maf_all.lenient[maf_all.lenient$Hugo_Symbol %in% "CASP8",]
#CASP8 hotspots:
#p.E184Q
#p.R233W or p.R233Q
#p.R68* or p.R68Q
CASP8_maf$Mutation <- sapply(CASP8_maf$HGVSp_Short, function(x)
  ifelse(x %in% c("p.E184Q"),"CASP8_E184",
         ifelse(x %in% c("p.R233W","p.R233Q"),"CASP8_R233",
                ifelse(x %in% c("p.R68*","p.R68Q"),"CASP8_R68","Other"))))
table(CASP8_maf$Mutation)

###########################
#HRAS mutations:
###########################
HRAS_maf <- maf_all.lenient[maf_all.lenient$Hugo_Symbol %in% "HRAS",]
#HRAS hotspot mutations:
#p.K117N or p.K117R
#p.Q61R or p.Q61L or p.Q61K or p.Q61H
#p.G60V or p.G60S
#p.A59T
#p.G13V or p.G13S or p.G13R or p.G13D or p.G13C
#p.G12V or p.G12S or p.G12D or p.G12R or p.G12C or p.G12A
#p.A11D
HRAS_maf$Mutation <- sapply(HRAS_maf$HGVSp_Short, function(x)
  ifelse(x %in% c("p.K117N","p.K117R"),"HRAS_K117",
         ifelse(x %in% c("p.Q61R","p.Q61L","p.Q61K","p.Q61H"),"HRAS_Q61",
                ifelse(x %in% c("p.G60V","p.G60S"),"HRAS_G60",
                       ifelse(x %in% c("p.A59T"),"HRAS_A59",
                              ifelse(x %in% c("p.G13V","p.G13S","p.G13R","p.G13D","p.G13C"),"HRAS_G13",
                                     ifelse(x %in% c("p.G12V","p.G12S","p.G12D","p.G12R","p.G12C","p.G12A"),"HRAS_G12",
                                            ifelse(x %in% c("p.A11D"),"HRAS_A11","Other"))))))))

table(HRAS_maf$Mutation)
###########################
#KRAS mutations:
############################
KRAS_maf <- maf_all.lenient[maf_all.lenient$Hugo_Symbol %in% "KRAS",]
#KRAS hotspot mutations:
#p.A146V or p.A146T or p.A146P
#p.K117N
#p.D92Y
#p.R68S or p.R68M 
#p.E63K
#p.E62K
#p.Q61R or p.Q61P or p.Q61L or p.Q61K or p.Q61H or p.Q61E
#p.A59T or p.A59G or p.A59E
#p.I36M
#p.P34L
#p.D33E
#p.Q22K
#p.L19F
#p.G13D or p.G13V or p.G13C
#p.G12V or p.G12S or p.G12R or p.G12D or p.G12C or p.G12A
#p.K5E
KRAS_maf$Mutation <- sapply(KRAS_maf$HGVSp_Short, function(x)
  ifelse(x %in% c("p.A146V","p.A146T","p.A146P"),"KRAS_A149",
         ifelse(x %in% c("p.K117N"),"KRAS_K117",
                ifelse(x %in% c("p.D92Y"),"KRAS_D92",
                       ifelse(x %in% c("p.R68S","p.R68M"),"KRAS_R68",
                              ifelse(x %in% c("p.E63K"),"KRAS_E63",
                                     ifelse(x %in% c("p.E62K"),"KRAS_E62",
                                            ifelse(x %in% c("p.Q61R","p.Q61P","p.Q61L","p.Q61K","p.Q61H","p.Q61E"),"KRAS_Q61",
                                                   ifelse(x %in% c("p.A59T","p.A59G","p.A59E"),"KRAS_A59",
                                                          ifelse(x %in% c("p.I36M"),"KRAS_I36",
                                                                 ifelse(x %in% c("p.P34L"),"KRAS_P34",
                                                                        ifelse(x %in% c("p.D33E"),"KRAS_D33",
                                                                               ifelse(x %in% c("p.Q22K"),"KRAS_Q22",
                                                                                      ifelse(x %in% c("p.L19F"),"KRAS_L19",
                                                                                                      ifelse(x %in% c("p.G13D","p.G13V","p.G13C"),"KRAS_G13",
                                                                                                             ifelse(x %in% c("p.G12V","p.G12S","p.G12R","p.G12D","p.G12C","p.G12A"),"KRAS_G12",
                                                                                                                    ifelse(x %in% c("p.K5E"),"KRAS_K5","Other")))))))))))))))))


table(KRAS_maf$Mutation)
###########################
#NRAS mutations:
NRAS_maf <- maf_all.lenient[maf_all.lenient$Hugo_Symbol %in% "NRAS",]
#p.G12V or p.G12S or p.G12R or p.G12D or p.G12C or p.G12A
#p.G13D or p.G13R
#p.D33H
#p.P34L
#p.Q61R or p.Q61P or p.Q61L or p.Q61K or p.Q61H
#p.A83G
NRAS_maf$Mutation <- sapply(NRAS_maf$HGVSp_Short, function(x)
  ifelse(x %in% c("p.G12V","p.G12S","p.G12R","p.G12D","p.G12C","p.G12A"),"NRAS_G12",
         ifelse(x %in% c("p.G13D","p.G13R"),"NRAS_G13",
                ifelse(x %in% c("p.D33H"),"NRAS_D33",
                       ifelse(x %in% c("p.P34L"),"NRAS_P34",
                              ifelse(x %in% c("p.Q61R","p.Q61P","p.Q61L","p.Q61K","p.Q61H"),"NRAS_Q61",
                                     ifelse(x %in% c("p.A83G"),"NRAS_A83","Other")))))))

table(NRAS_maf$Mutation)



###############################
##Combine hotspot mutation data
################################
maf_hotspot <- rbind(CASP8_maf, NRAS_maf, KRAS_maf, HRAS_maf)



###################################################
##Add hotspot annotation to the original dataframe
##################################################
load("maf_all.lenient.cosmic.RData")
hotspots <- unique(as.character(maf_hotspot$Mutation))
hotspots <- hotspots[!(hotspots %in% "Other")]
for (i in hotspots) {
  
  print(i)
  mutated_samples <- maf_hotspot[maf_hotspot$Mutation %in% i,]
  mutated_samples <- unique(as.character(mutated_samples$Patient_short))
  maf_all.lenient[[i]] <- sapply(maf_all.lenient$Patient_short, function(x)
    ifelse(x %in% mutated_samples,"MT","WT"))
}



########################################################
#Add information about whether a sample has TMD or not
########################################################
#Add information about whether a sample has TMD or not
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("programme_scores_and_TMD_assignments.RData")
prog_expr$Patient_short <- sapply(prog_expr$Barcode,
                                  function(x) paste(str_split(x, "-", simplify = T)[1:4], collapse = "-"))
combined.data <- merge(maf_all.lenient, prog_expr, 
                       by.x = "Patient_short", by.y = "Patient_short")
#####Rename where necessary and remove samples in the middle group:
table(combined.data$TMD_two_categories)
combined.data$TMD <- sapply(combined.data$TMD_two_categories, function(x)
  ifelse(x %in% "YES", "High",
         ifelse(x %in% "NO", "Low","Other")))
table(combined.data$TMD)
combined.data <- combined.data[combined.data$TMD %in% c("High","Low"),]



###########################################################################################
###Select hotspot mutations which have at least 1 mutation in the TMD high and low groups:
###########################################################################################
for (i in hotspots) {
  print(i)
  print(table(combined.data[[i]]))  
}

CASP8_hotspots <- unique(as.character(CASP8_maf$Mutation))
CASP8_hotspots <- CASP8_hotspots[!(CASP8_hotspots %in% "Other")]

KRAS_hotspots <- unique(as.character(KRAS_maf$Mutation))
KRAS_hotspots <- KRAS_hotspots[!(KRAS_hotspots %in% "Other")]

NRAS_hotspots <- unique(as.character(NRAS_maf$Mutation))
NRAS_hotspots <- NRAS_hotspots[!(NRAS_hotspots %in% c("Other","NRAS_P34"))]

HRAS_hotspots <- unique(as.character(HRAS_maf$Mutation))
HRAS_hotspots <- HRAS_hotspots[!(HRAS_hotspots %in% "Other")]


####################################
###KRAS hotspot fisher's exact test
###################################
pvalue <- NULL
OR <- NULL
Upper_confidence_interval <- NULL
Lower_confidence_interval <- NULL
FishersResults <- data.frame(KRAS_hotspots)
for (a in KRAS_hotspots) {
  
  test_table <- table(combined.data$TMD,combined.data[[a]])
  print(test_table)
  test <- fisher.test(test_table)
  p <- test$p.value
  or <- test$estimate
  down <- test$conf.int[1]
  Lower_confidence_interval <- c(Lower_confidence_interval, down)
  up <- test$conf.int[2]
  Upper_confidence_interval <- c(Upper_confidence_interval, up)
  pvalue <- c(pvalue,p)
  OR <- c(OR,or)
  
}
FishersResults$OR <- OR
FishersResults$P_value <- pvalue
p_adj <- p.adjust(pvalue, method = "BH")
FishersResults$P_adjust <- p_adj
FishersResults$LowerConfidenceInterval <- Lower_confidence_interval
FishersResults$UpperConfidenceInterval <- Upper_confidence_interval
FishersResults <- FishersResults[FishersResults$P_adjust <= 0.05,]
FishersResults$yaxis <- 1:2

setwd("~/Documents/GitHub/tumourMassDormancy/FishersExactTest_EnrichmentAnalysis/Figures/")
#Plot the results:
library(ggplot2)
pdf("FishersExactTestResults_KRAS.pdf", width = 4, height = 1)
p <-ggplot(FishersResults, aes(x = OR, y = yaxis))
p + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = UpperConfidenceInterval, xmin = LowerConfidenceInterval), size = .5, height = .2, color = "gray50") +
  geom_point(size = 1, color = "darkblue") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = FishersResults$yaxis, labels = FishersResults$KRAS_hotspots) +
  scale_x_continuous(limits = c(0,1.2))+
  ylab("") +
  xlab("Odds ratio")
dev.off()












###################################
###NRAS hotspot fisher's exact test
##################################
pvalue <- NULL
OR <- NULL
Upper_confidence_interval <- NULL
Lower_confidence_interval <- NULL
FishersResults <- data.frame(NRAS_hotspots)
for (a in NRAS_hotspots) {
  
  test_table <- table(combined.data$TMD,combined.data[[a]])
  print(test_table)
  test <- fisher.test(test_table)
  p <- test$p.value
  or <- test$estimate
  down <- test$conf.int[1]
  Lower_confidence_interval <- c(Lower_confidence_interval, down)
  up <- test$conf.int[2]
  Upper_confidence_interval <- c(Upper_confidence_interval, up)
  pvalue <- c(pvalue,p)
  OR <- c(OR,or)
  
}
FishersResults$OR <- OR
FishersResults$P_value <- pvalue
p_adj <- p.adjust(pvalue, method = "BH")
FishersResults$P_adjust <- p_adj
FishersResults$LowerConfidenceInterval <- Lower_confidence_interval
FishersResults$UpperConfidenceInterval <- Upper_confidence_interval
FishersResults <- FishersResults[FishersResults$P_adjust <= 0.05,]
FishersResults$yaxis <- 1


setwd("~/Documents/GitHub/tumourMassDormancy/FishersExactTest_EnrichmentAnalysis/Figures/")
#Plot the results:
library(ggplot2)
pdf("FishersExactTestResults_NRAS.pdf", width = 4, height = 1)
p <-ggplot(FishersResults, aes(x = OR, y = yaxis))
p + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = UpperConfidenceInterval, xmin = LowerConfidenceInterval), size = .5, height = .2, color = "gray50") +
  geom_point(size = 1, color = "darkblue") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = FishersResults$yaxis, labels = FishersResults$NRAS_hotspots) +
  scale_x_continuous(limits = c(0,1.2))+
  ylab("") +
  xlab("Odds ratio")
dev.off()





##################################
###HRAS hotspot fisher's exact test
##################################
pvalue <- NULL
OR <- NULL
Upper_confidence_interval <- NULL
Lower_confidence_interval <- NULL
FishersResults <- data.frame(HRAS_hotspots)
for (a in HRAS_hotspots) {
  
  test_table <- table(combined.data$TMD,combined.data[[a]])
  print(test_table)
  test <- fisher.test(test_table)
  p <- test$p.value
  or <- test$estimate
  down <- test$conf.int[1]
  Lower_confidence_interval <- c(Lower_confidence_interval, down)
  up <- test$conf.int[2]
  Upper_confidence_interval <- c(Upper_confidence_interval, up)
  pvalue <- c(pvalue,p)
  OR <- c(OR,or)
  
}
FishersResults$OR <- OR
FishersResults$P_value <- pvalue
p_adj <- p.adjust(pvalue, method = "BH")
FishersResults$P_adjust <- p_adj
FishersResults$LowerConfidenceInterval <- Lower_confidence_interval
FishersResults$UpperConfidenceInterval <- Upper_confidence_interval
FishersResults <- FishersResults[FishersResults$P_adjust <= 0.05,]
FishersResults$yaxis <- 1:2


setwd("~/Documents/GitHub/tumourMassDormancy/FishersExactTest_EnrichmentAnalysis/Figures/")
#Plot the results:
library(ggplot2)
pdf("FishersExactTestResults_HRAS.pdf", width = 4, height = 1)
p <-ggplot(FishersResults, aes(x = OR, y = yaxis))
p + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = UpperConfidenceInterval, xmin = LowerConfidenceInterval), size = .5, height = .2, color = "gray50") +
  geom_point(size = 1, color = "darkblue") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = FishersResults$yaxis, labels = FishersResults$HRAS_hotspots) +
  scale_x_continuous(limits = c(0,14))+
  ylab("") +
  xlab("Odds ratio")
dev.off()





####################################
###CASP8 hotspot fishers exact test
####################################
pvalue <- NULL
OR <- NULL
Upper_confidence_interval <- NULL
Lower_confidence_interval <- NULL
FishersResults <- data.frame(CASP8_hotspots)
for (a in CASP8_hotspots) {
  
  test_table <- table(combined.data$TMD,combined.data[[a]])
  print(test_table)
  test <- fisher.test(test_table)
  p <- test$p.value
  or <- test$estimate
  down <- test$conf.int[1]
  Lower_confidence_interval <- c(Lower_confidence_interval, down)
  up <- test$conf.int[2]
  Upper_confidence_interval <- c(Upper_confidence_interval, up)
  pvalue <- c(pvalue,p)
  OR <- c(OR,or)
  
}
FishersResults$OR <- OR
FishersResults$P_value <- pvalue
p_adj <- p.adjust(pvalue, method = "BH")
FishersResults$P_adjust <- p_adj
FishersResults$LowerConfidenceInterval <- Lower_confidence_interval
FishersResults$UpperConfidenceInterval <- Upper_confidence_interval
FishersResults <- FishersResults[FishersResults$P_adjust <= 0.05,]
#No significant results


