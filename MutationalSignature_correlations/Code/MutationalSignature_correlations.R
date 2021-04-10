########################################################################
#####Correlation between programme scores and mutational signatures:
#########################################################################

#Load packages:
library(corrplot)
library(ggpubr)


###############
###Load scores:
###############
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("programme_scores_and_TMD_assignments.RData")



###########################
#Load mutational signatures
###########################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/MutationalSignatures/")
load("sigs.defaultnounknown.ACC.RData")
sigs.all <- sigs.defaultnounknown
CT <- c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS")
for (a in CT) {
  load(paste("sigs.defaultnounknown.",a,".RData",sep = ""))
  sigs.all <- rbind(sigs.all,sigs.defaultnounknown)
  
}
#Remove signatures which are likely to be artifacts (>SBS44)
sigs.all <- sigs.all[,1:49]



#############################
#Merge the two dataframes together:
sigs.all$SampleID <- sapply(rownames(sigs.all), function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
sig_corr <- merge(sigs.all, prog_expr,
                  by.x = "SampleID", by.y = "SampleID")
signatures <- colnames(sigs.all)
sig_corr <- sig_corr[,colnames(sig_corr) %in% c(signatures,"mean_APOBEC","Exhaustion","TumourMassDormancy"),]



##########################################
####Pancancer correlation:
#########################################
pancacer_SBS_matrix <- cor(sig_corr[, 2:53])
pancancer_p_matrix <- cor.mtest(sig_corr[, 2:53])$p
pancacer_SBS_matrix <- pancacer_SBS_matrix[1:49,50:52]
pancancer_p_matrix <- pancancer_p_matrix[1:49,50:52]
pancancer_p_matrix[is.na(pancancer_p_matrix)] <- 1
pancacer_SBS_matrix[is.na(pancacer_SBS_matrix)] <- 0
pancancer_p_matrix <- data.frame(pancancer_p_matrix)
pancacer_SBS_matrix <- data.frame(pancacer_SBS_matrix)
colnames(pancancer_p_matrix) <- colnames(pancacer_SBS_matrix)
rownames(pancancer_p_matrix) <- rownames(pancacer_SBS_matrix)
#Remove non-sigificnat signatures:
pancacer_SBS_matrix <- pancacer_SBS_matrix[!(rownames(pancacer_SBS_matrix) %in% c("SBS7a","SBS7c","SBS7d","SBS9","SBS11","SBS12","SBS17a","SBS19","SBS23","SBS24","SBS25","SBS26","SBS27","SBS28","SBS29","SBS33","SBS34","SBS35","SBS36","SBS37","SBS38","SBS41","SBS42","SBS43")),]
pancancer_p_matrix <- pancancer_p_matrix[!(rownames(pancancer_p_matrix) %in% c("SBS7a","SBS7c","SBS7d","SBS9","SBS11","SBS12","SBS17a","SBS19","SBS23","SBS24","SBS25","SBS26","SBS27","SBS28","SBS29","SBS33","SBS34","SBS35","SBS36","SBS37","SBS38","SBS41","SBS42","SBS43")),]
pancacer_SBS_matrix <- as.matrix(pancacer_SBS_matrix)
pancancer_p_matrix <- as.matrix(pancancer_p_matrix)
setwd("~/Documents/GitHub/tumourMassDormancy/MutationalSignature_correlations/Figures/")
pdf("PancancerMutationalSigantureProgrammeCorrelations.pdf",height = 10,width = 5)
corrplot(pancacer_SBS_matrix,
         method = "color",
         insig = "label_sig",
         p.mat = pancancer_p_matrix,
         cl.pos = "b",
         pch.cex = 1,
         is.corr = F,
         tl.col = "black",
         cl.length = 3,
         cl.cex = 0.8)
dev.off()






##############################################################
##Cancerwise programme correlations with mutational signatures
#############################################################
sig_corr <- merge(sigs.all, prog_expr,
                  by.x = "SampleID", by.y = "SampleID")
sig_corr <- sig_corr[,colnames(sig_corr) %in% c(signatures,"mean_APOBEC","Exhaustion","AngiogenicDormancy","ImmunogenicDormancy","TumourMassDormancy","cancer_type"),]
CT <- unique(sig_corr$cancer_type)
setwd("~/Documents/GitHub/tumourMassDormancy/MutationalSignature_correlations/CancerwiseProgrammeCorrelations/")
CT <- CT[!(CT %in% c("OV","TGCT","PAAD","MESO","CHOL","SKCM","KICH"))]
CT <- CT[order(CT)]
pancacer_SBS_matrix <- data.frame(pancacer_SBS_matrix)
pancancer_p_matrix <- data.frame(pancancer_p_matrix)

################Calculate cancerwise correlations
for (a in CT) {
  
  test.data <- sig_corr[sig_corr$cancer_type %in% a,]
  test.data$cancer_type <- NULL
  
  SBS_matrix <- cor(test.data[, 2:53])
  p_matrix <- cor.mtest(test.data[, 2:53])$p
  SBS_matrix <- SBS_matrix[1:49,50:52]
  p_matrix <- p_matrix[1:49,50:52]
  rownames(p_matrix) <- rownames(SBS_matrix)
  colnames(SBS_matrix) <- c(paste("mean_APOBEC",a,sep = "_"),paste("Exhaustion",a,sep = "_"),paste("TumourMassDormancy",a,sep = "_"))
  colnames(p_matrix) <- c(paste("mean_APOBEC",a,sep = "_"),paste("Exhaustion",a,sep = "_"),paste("TumourMassDormancy",a,sep = "_"))
  SBS_matrix <- SBS_matrix[!(rownames(SBS_matrix) %in% c("SBS7a","SBS7c","SBS7d","SBS9","SBS11","SBS12","SBS17a","SBS19","SBS23","SBS24","SBS25","SBS26","SBS27","SBS28","SBS29","SBS33","SBS34","SBS35","SBS36","SBS37","SBS38","SBS41","SBS42","SBS43")),]
  p_matrix <- p_matrix[!(rownames(p_matrix) %in% c("SBS7a","SBS7c","SBS7d","SBS9","SBS11","SBS12","SBS17a","SBS19","SBS23","SBS24","SBS25","SBS26","SBS27","SBS28","SBS29","SBS33","SBS34","SBS35","SBS36","SBS37","SBS38","SBS41","SBS42","SBS43")),]
  p_matrix[is.na(p_matrix)] <- 1
  SBS_matrix[is.na(SBS_matrix)] <- 0
  p_matrix <- data.frame(p_matrix)
  SBS_matrix <- data.frame(SBS_matrix)
  pancacer_SBS_matrix <- cbind(pancacer_SBS_matrix, SBS_matrix)
  pancancer_p_matrix <- cbind(pancancer_p_matrix, p_matrix)
  
  
}
pancacer_SBS_matrix <- as.matrix(pancacer_SBS_matrix)
pancancer_p_matrix <- as.matrix(pancancer_p_matrix)


######Plot the cancerwise correlations:
setwd("~/Documents/GitHub/tumourMassDormancy/MutationalSignature_correlations/Figures/")
pdf("CancerwiseMutationalSigantureProgrammeCorrelations.pdf",height = 10,width = 17)
corrplot(pancacer_SBS_matrix,
         method = "color",
         insig = "label_sig",
         p.mat = pancancer_p_matrix,
         cl.pos = "b",
         pch.cex = 1,
         is.corr = F,
         tl.col = "black",
         cl.length = 3,
         cl.cex = 0.8)
dev.off()




###Order the columns and replot:
pancacer_SBS_matrix <- data.frame(pancacer_SBS_matrix)
pancancer_p_matrix <- data.frame(pancancer_p_matrix)
pancacer_SBS_matrix <- pancacer_SBS_matrix[,order(colnames(pancacer_SBS_matrix))]
pancancer_p_matrix <- pancancer_p_matrix[,order(colnames(pancancer_p_matrix))]
pancacer_SBS_matrix <- as.matrix(pancacer_SBS_matrix)
pancancer_p_matrix <- as.matrix(pancancer_p_matrix)
pdf("CancerwiseMutationalSigantureProgrammeCorrelationsOrdered.pdf",height = 10,width = 17)
corrplot(pancacer_SBS_matrix,
         method = "color",
         insig = "label_sig",
         p.mat = pancancer_p_matrix,
         cl.pos = "b",
         pch.cex = 1,
         is.corr = F,
         tl.col = "black",
         cl.length = 3,
         cl.cex = 0.8)
dev.off()








