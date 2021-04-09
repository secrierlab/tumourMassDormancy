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
genes_to_save <- unique(cosmic_list) #Downloading genes from programmes and from COSMIC

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
                            "SIFT","PolyPhen")]
  df.maf.short$Cancer <- as.character(tumour)
  df.maf.short <- subset(df.maf.short, Hugo_Symbol %in% genes_to_save)
  
  if (tumour == solid_tumours[1]){
    maf_all <- df.maf.short
  } else {
    maf_all <- rbind(maf_all, df.maf.short)
  }
}

save(maf_all, file = "maf_cosmic.RData")



######################################################
########Load the mutation data and process further:
######################################################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/MutationData/")
load("maf_cosmic.RData")

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

#Because aggregation function sums all the values, we give 1 if sum>0,
#and 0 if there are no mutations
maf_all.lenient$Mutation <- 1
maf_all.lenient <- dcast(maf_all.lenient, value.var = "Mutation",
                         formula = Tumor_Sample_Barcode+Patient_short+Cancer~Hugo_Symbol,
                         fun.aggregate = function(x) ifelse(sum(x)>0, yes = 1, no = 0))
save(maf_all.lenient, file = "maf_all.lenient.cosmic.RData")




#####################################################################
#####Fisher's exact test:
#####################################################################

####Load mutation data:
setwd("~/Documents/GitHub/tumourMassDormancy/Data/MutationData/")
load("maf_all.lenient.cosmic.RData")
cosmic_genes <- colnames(maf_all.lenient)
cosmic_genes <- cosmic_genes[-(1:3)]

#Add information about whether a sample has TMD or not
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("programme_scores_and_TMD_assignments.RData")
prog_expr$Patient_short <- sapply(prog_expr$Barcode,
                                  function(x) paste(str_split(x, "-", simplify = T)[1:4], collapse = "-"))
combined.data <- merge(maf_all.lenient, prog_expr, 
                       by.x = "Patient_short", by.y = "Patient_short")

#####Rename where necessary and remove samples in the middle group:
table(combined.data$TMD_three_categories)
combined.data$TMD <- sapply(combined.data$TMD_three_categories, function(x)
  ifelse(x %in% "YES", "High",
         ifelse(x %in% "NO", "Low","Other")))
table(combined.data$TMD)
combined.data <- combined.data[combined.data$TMD %in% c("High","Low"),]
for (a in cosmic_genes) {
  combined.data[[a]] <- sapply(combined.data[[a]], function(x)
    ifelse(x %in% 0,"WT","MT"))
  
}

#remove genes which are mutated at least in 1% of samples (this equates to 19.33 mutations in the dataframe)
mut <- NULL
for (i in cosmic_genes) {
  print(i)
  test <- combined.data[combined.data[[i]] %in% "MT",]
  number_of_mutations <- dim(test)[1]
  mut <- c(mut, number_of_mutations)
}
MutFreq <- data.frame(cosmic_genes, mut)
MutFreq <- MutFreq[MutFreq$mut > 19.33,]
cosmic_genes <- as.character(MutFreq$cosmic_genes)

###Fisher's exact test:
pvalue <- NULL
OR <- NULL
Upper_confidence_interval <- NULL
Lower_confidence_interval <- NULL
FishersResults <- data.frame(cosmic_genes)
for (a in cosmic_genes) {
  
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

##Save the results:
setwd("~/Documents/GitHub/tumourMassDormancy/FishersExactTest_EnrichmentAnalysis/")
save(FishersResults, file = "FishersExactTestResults.RData")

#Save significant results only:
FishersResults <- FishersResults[FishersResults$P_adjust < 0.05,]
write.csv(FishersResults, file = "FishersExactTestResultsSignificant.csv")






###############################################
####Visualise the results:
###############################################
FishersResults <- FishersResults[order(FishersResults$OR),]
FishersResults$yaxis <- 1:15
###Log transform the OR and CI
FishersResults$OR <- log2(FishersResults$OR)
FishersResults$LowerConfidenceInterval <- log2(FishersResults$LowerConfidenceInterval)
FishersResults$UpperConfidenceInterval <- log2(FishersResults$UpperConfidenceInterval)


#Plot the results:
setwd("~/Documents/GitHub/tumourMassDormancy/FishersExactTest_EnrichmentAnalysis/Figures/")
pdf("FishersExactTestResults.pdf", width = 4, height = 3)
p <-ggplot(FishersResults, aes(x = OR, y = yaxis))
p + geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = UpperConfidenceInterval, xmin = LowerConfidenceInterval), size = .5, height = .2, color = "gray50") +
  geom_point(size = 1, color = "darkblue") +
  theme_classic() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = FishersResults$yaxis, labels = FishersResults$cosmic_genes) +
  scale_x_continuous(limits = c(-3.5,4.5))+
  ylab("") +
  xlab("log Odds ratio")
dev.off()



