#########################################################################
###Fishers exact test between samples with and without TMD (Cancerwise)
###########################################################################

#Load required packages:
library(TCGAbiolinks)
library(stringr)
library(dplyr)
library(maditr)
library(ggplot2)

################################################
#Specifiy cancer types of interest: (which have samples in with both low and high TMD)
################################################
solid_tumours <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KIRP","LGG","LIHC","LUAD","LUSC","OV","PAAD","PRAD","READ","SARC","SKCM","STAD","UCEC","UCS","UVM")


###############################################
#Load mutation data for cosmic genes only 
###############################################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/MutationData/")
load("maf_cosmic.RData")



for (a in solid_tumours) {
  
  
  print(a)
  
  ####################################
  ###Further process mutation data
  ####################################
  
  maf_short <- maf_all[maf_all$Cancer %in% a,]
  maf_short$Patient_short <- sapply(maf_short$Tumor_Sample_Barcode,
                                    function(x) paste(str_split(x, "-", simplify = T)[1:4], collapse = "-"))
  maf_short.lenient <- maf_short[ which(maf_short$Variant_Classification
                                        %in% c("Missense_Mutation",
                                               "Nonsense_Mutation",
                                               "Nonstop_Mutation",
                                               "Frameshift_Deletion",
                                               "Frameshift_Insertion",
                                               "Inframe_Insertion",
                                               "Inframe_Deletion")),
  ]
  rm(maf_short)
  #Because aggregation function sums all the values, we give 1 if sum>0,
  #and 0 if there are no mutations
  maf_short.lenient$Mutation <- 1
  maf_short.lenient <- dcast(maf_short.lenient, value.var = "Mutation",
                             formula = Tumor_Sample_Barcode+Patient_short+Cancer~Hugo_Symbol,
                             fun.aggregate = function(x) ifelse(sum(x)>0, yes = 1, no = 0))
  
  
  
  
  
  
  #####################################################################
  #####Fisher's exact test:
  #####################################################################
  
  ##Select cosmic genes of interest:
  cosmic_genes <- colnames(maf_short.lenient)
  cosmic_genes <- cosmic_genes[-(1:3)]
  
  #Add information about whether a sample has TMD or not
  setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
  load("programme_scores_and_TMD_assignments.RData")
  prog_expr$Patient_short <- sapply(prog_expr$Barcode,
                                    function(x) paste(str_split(x, "-", simplify = T)[1:4], collapse = "-"))
  combined.data <- merge(maf_short.lenient, prog_expr, 
                         by.x = "Patient_short", by.y = "Patient_short")
  
  #####Rename where necessary and remove samples in the middle group:
  table(combined.data$TMD_three_categories)
  combined.data$TMD <- sapply(combined.data$TMD_three_categories, function(x)
    ifelse(x %in% "YES", "High",
           ifelse(x %in% "NO", "Low","Other")))
  table(combined.data$TMD)
  combined.data <- combined.data[combined.data$TMD %in% c("High","Low"),]
  print(table(combined.data$TMD))
  
  for (i in cosmic_genes) {
    combined.data[[i]] <- sapply(combined.data[[i]], function(x)
      ifelse(x %in% 0,"WT","MT"))
    
  }
  
  #remove genes which are mutated at least in 5% of samples 
  total_samples <- dim(combined.data)[1]
  threshold <- total_samples / 20
  
  mut <- NULL
  for (i in cosmic_genes) {
    print(i)
    test <- combined.data[combined.data[[i]] %in% "MT",]
    number_of_mutations <- dim(test)[1]
    mut <- c(mut, number_of_mutations)
  }
  
  MutFreq <- data.frame(cosmic_genes, mut)
  MutFreq <- MutFreq[MutFreq$mut > threshold,]
  cosmic_genes <- as.character(MutFreq$cosmic_genes)
  
  
  ###Fisher's exact test:
  pvalue <- NULL
  OR <- NULL
  Upper_confidence_interval <- NULL
  Lower_confidence_interval <- NULL
  FishersResults <- data.frame(cosmic_genes)
  for (i in cosmic_genes) {
    
    test_table <- table(combined.data$TMD,combined.data[[i]])
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
  setwd("~/Documents/GitHub/tumourMassDormancy/FishersExactTest_EnrichmentAnalysis/CancerwiseAnalysisResults/")
  save(FishersResults, file = paste("FishersExactTestResults_",a,".RData",sep = ""))
  
  
}




