#########################################################################
####TMD Classification:
#########################################################################



##Load the required packages:
library(biomaRt)
library(ggpubr)
library(gtools)
library(corrplot)
library(ggplot2)



#######################
####Define gene lists used to produce the proliferation/apoptosis ratio:
#######################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/GeneLists/")
E2F_targets <- read.table("HallmarkE2FTargetGeneList.txt", header = FALSE, sep = "\t")
E2F_targets <- as.character(E2F_targets$V1)
Apoptosis_markers <- read.table("HallmarkApoptosisGeneList.txt", header = FALSE, sep = "\t")
Apoptosis_markers <- as.character(Apoptosis_markers$V1)
all_genes_list <- unique(c(E2F_targets, Apoptosis_markers))



###############################
#Load the RNA-seq data
###############################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/RNA_seq/")
load("combined_experssion_FPKM.RData") #or load "example_FPKM_data.RData" for example dataframe with 100 samples

##############################
#Load programme scores:
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("program_expression_scores_updated.RData")
rnaseq_data <- combined_data

###############################
###Obtain ENSG ID for all genes
###############################
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
biomart_conversion <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=all_genes_list, mart=human, filters = "hgnc_symbol")
biomart_conversion <- biomart_conversion[biomart_conversion$ensembl_gene_id %in% colnames(rnaseq_data),]

###################################
###Convert HGNC symbols to ENSG ID 
##################################
rnaseq_data  <- rnaseq_data[rownames(rnaseq_data) %in% prog_expr$Barcode,]
for (i in all_genes_list) {
  
  selected_biomart <- biomart_conversion[biomart_conversion$hgnc_symbol %in% i,]
  a <- selected_biomart$ensembl_gene_id
  rnaseq_data[[i]] <- rnaseq_data[[a]]
  
  
}
rnaseq_data <- rnaseq_data[,colnames(rnaseq_data) %in% c(all_genes_list,"cancer_type")]



########################################
#Calculate proliferation/apoptosis ratio:
########################################
#Check that samples are in the same order in the prog_expr and rnaseq_data dataframes:
all(rownames(rnaseq_data) %in% prog_expr$Barcode)
all(rownames(rnaseq_data) == prog_expr$Barcode)
prog_expr$cancer_type <- rnaseq_data$cancer_type
prog_expr$Proliferation <- rowMeans(rnaseq_data[,colnames(rnaseq_data) %in% E2F_targets])
prog_expr$Apoptosis <- rowMeans(rnaseq_data[,colnames(rnaseq_data) %in% Apoptosis_markers])
prog_expr$Proliferation_Apoptosis_ratio <- prog_expr$Proliferation / prog_expr$Apoptosis



############################################################################################
##Correlation plots between different programme scores and the proliferation/apoptosis ratio:
############################################################################################

setwd("~/Documents/GitHub/tumourMassDormancy/TumourMassDormancy_ProliferationApoptosisRatio_relationship/Figures/")

pdf("Dormancy_vs_proliferation_apoptosis_ratio.pdf", width = 5, height = 5)
p <- ggscatter(prog_expr, x = "TumourMassDormancy", y = "Proliferation_Apoptosis_ratio",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
p + labs(x = "Tumour Mass Dormancy", y = "Proliferation/Apoptosis ratio")
dev.off()

pdf("ImmuneDormancy_vs_proliferation_apoptosis_ratio.pdf", width = 5, height = 5)
p <- ggscatter(prog_expr, x = "ImmunogenicDormancy", y = "Proliferation_Apoptosis_ratio",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
p + labs(x = "ImmunogenicDormancy", y = "Prolfieration/Apoptosis ratio")
dev.off()

pdf("AngiogenicDormancy_vs_proliferation_apoptosis_ratio.pdf", width = 5, height = 5)
p <- ggscatter(prog_expr, x = "AngiogenicDormancy", y = "Proliferation_Apoptosis_ratio",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
p + labs(x = "AngiogenicDormancy", y = "Prolifeartion/Apoptosis ratio")
dev.off()


#############################
#Cancer-wise correlation plots:
##############################

pdf("Dormancy_vs_proliferation_apoptosis_ratio_cancerwise.pdf", width = 13, height = 9)
p <- ggscatter(prog_expr, x = "TumourMassDormancy", y = "Proliferation_Apoptosis_ratio",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
print(p + labs(x = "Tumour Mass Dormancy", y = "Proliferation/Apoptosis ratio") + facet_wrap(~cancer_type))
dev.off()

pdf("ImmuneDormancy_vs_proliferation_apoptosis_ratio_cancerwise.pdf", width = 13, height = 9)
p <- ggscatter(prog_expr, x = "ImmunogenicDormancy", y = "Proliferation_Apoptosis_ratio",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
print(p + labs(x = "ImmunogenicDormancy", y = "Proliferation/Apoptosis ratio") + facet_wrap(~cancer_type))
dev.off()

pdf("AngiogenicDormancy_vs_proliferation_apoptosis_ratio_cancerwise.pdf", width = 13, height = 9)
p <- ggscatter(prog_expr, x = "AngiogenicDormancy", y = "Proliferation_Apoptosis_ratio",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
print(p + labs(x = "AngiogenicDormancy", y = "Proliferation/Apoptosis ratio") + facet_wrap(~cancer_type))
dev.off()




###############################################################################################
###Interplay between immunological and angiogenic dormancy
###############################################################################################

setwd("~/Documents/GitHub/tumourMassDormancy/TumourMassDormancy_ProliferationApoptosisRatio_relationship/Figures/")
summary(prog_expr$ImmunogenicDormancy)
summary(prog_expr$AngiogenicDormancy)
plot.data <- prog_expr[!(prog_expr$cancer_type %in% "KICH"),]

###Colour by TMD levels:
pdf("Dormancy_programme_correlations_colour_by_TMD.pdf",height = 9, width = 15)
p <- ggplot(plot.data, aes(x=ImmunogenicDormancy, y=AngiogenicDormancy, color = TumourMassDormancy)) + 
  geom_point(size=3, alpha = 0.5) + 
  scale_color_gradient2(low = "#59ac53", midpoint = mean(prog_expr$TumourMassDormancy),mid = "grey95", high = "#8b5aa8") + theme_classic() +
  facet_wrap(~cancer_type) + geom_hline(yintercept = 0.03183, linetype = "dashed",color = "black") + geom_vline(xintercept = 0.50831, linetype = "dashed",color = "black")
print(p)
dev.off()

###Colour by Proliferation/Apoptosis ratio levels:
pdf("Dormancy_programme_correlations_colour_by_proliferation_apoptosis_ratio.pdf",height = 9, width = 15)
p <- ggplot(plot.data, aes(x=ImmunogenicDormancy, y=AngiogenicDormancy, color = Proliferation_Apoptosis_ratio)) + 
  geom_point(size=3, alpha = 0.5) + 
  scale_color_gradient2(low = "#59ac53", midpoint = mean(prog_expr$Proliferation_Apoptosis_ratio),mid = "grey95", high = "#8b5aa8") + theme_classic() +
  facet_wrap(~cancer_type) + geom_hline(yintercept = 0.03183, linetype = "dashed",color = "black") + geom_vline(xintercept = 0.50831, linetype = "dashed",color = "black")
print(p)
dev.off()


################################################################################################################
###Classification of samples with and without TMD:
############################################################################################################

####Classify patients as having TMD if their TMD score is in the UQ and their Proliferation/apoptosis ratio is < 1
summary(prog_expr$TumourMassDormancy) #UQ = 0.22696
Samples1 <- prog_expr[prog_expr$TumourMassDormancy > 0.22696,]
Samples1 <- as.character(Samples1$Barcode)
Samples2 <- prog_expr[prog_expr$Proliferation_Apoptosis_ratio < 1,]
Samples2 <- as.character(Samples2$Barcode)
TMD_samples <- Samples1[Samples1 %in% Samples2]
prog_expr$TMD_two_categories <- sapply(prog_expr$Barcode, function(x)
  ifelse(x %in% TMD_samples, "YES","NO"))
table(prog_expr$TMD_two_categories) #There are 1586 samples with TMD and 8045 samples without TMD

###Classify samples as having low tumour mass dormancy if their TMD score is in the LQ and their Proliferation/apoptosis ration > 1
summary(prog_expr$TumourMassDormancy) #LQ = -0.14872
Samples1 <- prog_expr[prog_expr$TumourMassDormancy < -0.14872,]
Samples1 <- as.character(Samples1$Barcode)
Samples2 <- prog_expr[prog_expr$Proliferation_Apoptosis_ratio > 1,]
Samples2 <- as.character(Samples2$Barcode)
NO_TMD_samples <- Samples1[Samples1 %in% Samples2]
prog_expr$TMD_three_categories <- sapply(prog_expr$Barcode, function(x)
  ifelse(x %in% TMD_samples, "YES",
         ifelse(x %in% NO_TMD_samples, "NO","MID")))
table(prog_expr$TMD_three_categories) #There are 1586 samples with TMD, 725 samples with very low TMD and 7320 samples in the middle group

###################################################################################################
#######Split samples with TMD into those which angiogenesis dormancy, immunogenic dormancy or both
#################################################################################################
#Samples with angiogenic dormancy will take angiogenic dormancy score in the UQ of the score range
summary(prog_expr$AngiogenicDormancy) #UQ = 0.03183
angiogenesis_high <- prog_expr[prog_expr$AngiogenicDormancy > 0.03183,]
angiogenesis_high <- as.character(angiogenesis_high$Barcode)


#Samples with immunological dormancy will have immunological dormancy score in the UQ of the score range
summary(prog_expr$ImmunogenicDormancy) #UQ = 0.50831
immunogenic_high <- prog_expr[prog_expr$ImmunogenicDormancy > 0.50831,]
immunogenic_high <- as.character(immunogenic_high$Barcode)

angiogenesis_samples <- TMD_samples[TMD_samples %in% angiogenesis_high]
angiogenesis_samples <- angiogenesis_samples[!(angiogenesis_samples %in% immunogenic_high)]
immunogenic_samples <- TMD_samples[TMD_samples %in% immunogenic_high]
immunogenic_samples <- immunogenic_samples[!(immunogenic_samples %in% angiogenesis_high)]
both_angiogenensis_and_angiogenesis <- TMD_samples[!(TMD_samples %in% c(immunogenic_samples, angiogenesis_samples))]
##Add the additional amotation:
prog_expr$TMD_two_categories_detailed <- sapply(prog_expr$Barcode, function(x)
  ifelse(x %in% immunogenic_samples, "Immunogenic Dormancy",
         ifelse(x %in% angiogenesis_samples, "Angiogenic Dormancy",
                ifelse(x %in% both_angiogenensis_and_angiogenesis, "Angiogenic and Immunogenic Dormancy","NO"))))
prog_expr$TMD_three_categories_detailed <- sapply(prog_expr$Barcode, function(x)
  ifelse(x %in% immunogenic_samples, "Immunogenic Dormancy",
         ifelse(x %in% angiogenesis_samples, "Angiogenic Dormancy",
                ifelse(x %in% both_angiogenensis_and_angiogenesis, "Angiogenic and Immunogenic Dormancy",
                       ifelse(x %in% NO_TMD_samples,"NO","MID")))))
prog_expr$TMD <- NULL

#########################
#Save the TMD assignments:
##########################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
save(prog_expr, file = "programme_scores_and_TMD_assignments.RData")


##############################################
####Plot the TMD assignments:
##############################################

setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("programme_scores_and_TMD_assignments.RData")
prog_expr <- prog_expr[!(prog_expr$cancer_type %in% "KICH"),] #Omit the KICH study as it showed low levels of TMD:

setwd("~/Documents/GitHub/tumourMassDormancy/TumourMassDormancy_ProliferationApoptosisRatio_relationship/Figures/")
pdf("Dormancy_programme_correlations_colour_by_TMD_assignments_two_categories.pdf",height = 9, width = 15)
p <- ggplot(prog_expr, aes(x=ImmunogenicDormancy, y=AngiogenicDormancy, color = TMD_two_categories)) + 
  geom_point(size=3, alpha = 0.5) + 
  scale_color_manual(values = c("YES" = "#8b5aa8", "NO" = "grey95")) + theme_classic() +
  facet_wrap(~cancer_type) + geom_hline(yintercept = 0.03183, linetype = "dashed",color = "black") + geom_vline(xintercept = 0.50831, linetype = "dashed",color = "black")
print(p)
dev.off()

pdf("Dormancy_programme_correlations_colour_by_TMD_assignments_three_categories.pdf",height = 9, width = 15)
p <- ggplot(prog_expr, aes(x=ImmunogenicDormancy, y=AngiogenicDormancy, color = TMD_three_categories)) + 
  geom_point(size=3, alpha = 0.5) + 
  scale_color_manual(values = c("YES" = "#8b5aa8", "NO" = "grey40", "MID" = "lightgrey")) + theme_classic() +
  facet_wrap(~cancer_type) + geom_hline(yintercept = 0.03183, linetype = "dashed",color = "black") + geom_vline(xintercept = 0.50831, linetype = "dashed",color = "black")
print(p)
dev.off()

pdf("Dormancy_programme_correlations_colour_by_TMD_assignments_two_categories_detailed.pdf",height = 7, width = 13)
p <- ggplot(prog_expr, aes(x=ImmunogenicDormancy, y=AngiogenicDormancy, color = TMD_two_categories_detailed)) + 
  geom_point(size=2) + 
  scale_color_manual(values = c("Angiogenic and Immunogenic Dormancy" = "#1B9E77", "NO" = "grey95", "Angiogenic Dormancy" = "#D95F02", "Immunogenic Dormancy" = "#7570B3")) + theme_classic() +
  facet_wrap(~cancer_type, ncol = 6) + geom_hline(yintercept = 0.03183, linetype = "dashed",color = "black") + geom_vline(xintercept = 0.50831, linetype = "dashed",color = "black")
print(p + theme(legend.position = "botton"))
dev.off()

pdf("Dormancy_programme_correlations_colour_by_TMD_assignments_three_categories_detailed.pdf",height = 7, width = 13)
p <- ggplot(prog_expr, aes(x=ImmunogenicDormancy, y=AngiogenicDormancy, color = TMD_three_categories_detailed)) + 
  geom_point(size=2) + 
  scale_color_manual(values = c("Angiogenic and Immunogenic Dormancy" = "#1B9E77", "NO" = "grey52","MID" = "grey78", "Angiogenic Dormancy" = "#D95F02", "Immunogenic Dormancy" = "#7570B3")) + theme_classic() +
  facet_wrap(~cancer_type, ncol = 6) + geom_hline(yintercept = 0.03183, linetype = "dashed",color = "black") + geom_vline(xintercept = 0.50831, linetype = "dashed",color = "black")
print(p)
dev.off()





######################################################
###Boxplot programme expression comparison between groups:
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("programme_scores_and_TMD_assignments.RData")
table(prog_expr$TMD_two_categories_detailed)
prog_expr$TMD_two_categories_detailed <- factor(prog_expr$TMD_two_categories_detailed, levels = c("NO","Angiogenic Dormancy","Immunogenic Dormancy","Angiogenic and Immunogenic Dormancy"))
setwd("~/Documents/GitHub/tumourMassDormancy/TumourMassDormancy_ProliferationApoptosisRatio_relationship/Figures/")

my_comparisons <- list(c("NO","Angiogenic and Immunogenic Dormancy"),c("NO","Angiogenic Dormancy"),c("NO","Immunogenic Dormancy"))
pdf("Boxplot_exhaustion_comparison_TMD_groups.pdf",height = 3, width = 6)
p<-ggplot(prog_expr, aes(x=TMD_two_categories_detailed, y=Exhaustion, fill=TMD_two_categories_detailed)) +
  geom_boxplot(outlier.size = 1)
p + scale_color_manual(values = c("Angiogenic and Immunogenic Dormancy" = "#1B9E77", "NO" = "grey95", "Angiogenic Dormancy" = "#D95F02", "Immunogenic Dormancy" = "#7570B3")) + theme_classic() + stat_compare_means(comparisons = my_comparisons, label = "p.signif") # Add pairwise comparisons p-value
dev.off()

my_comparisons <- list(c("NO","Angiogenic and Immunogenic Dormancy"),c("NO","Angiogenic Dormancy"),c("NO","Immunogenic Dormancy"))
pdf("Boxplot_TMD_comparison_TMD_groups.pdf",height = 3, width = 6)
p<-ggplot(prog_expr, aes(x=TMD_two_categories_detailed, y=TumourMassDormancy, fill=TMD_two_categories_detailed)) +
  geom_boxplot(outlier.size = 1)
p + scale_color_manual(values = c("Angiogenic and Immunogenic Dormancy" = "#1B9E77", "NO" = "grey95", "Angiogenic Dormancy" = "#D95F02", "Immunogenic Dormancy" = "#7570B3")) + theme_classic() + stat_compare_means(comparisons = my_comparisons, label = "p.signif")  # Add pairwise comparisons p-value
dev.off()

my_comparisons <- list(c("NO","Angiogenic and Immunogenic Dormancy"),c("NO","Angiogenic Dormancy"),c("NO","Immunogenic Dormancy"))
pdf("Boxplot_APOBEC_comparison_TMD_groups.pdf",height = 3, width = 6)
p<-ggplot(prog_expr, aes(x=TMD_two_categories_detailed, y=mean_APOBEC, fill=TMD_two_categories_detailed)) +
  geom_boxplot(outlier.size = 1)
p + scale_color_manual(values = c("Angiogenic and Immunogenic Dormancy" = "#1B9E77", "NO" = "grey95", "Angiogenic Dormancy" = "#D95F02", "Immunogenic Dormancy" = "#7570B3")) + theme_classic() + stat_compare_means(comparisons = my_comparisons, label = "p.signif")  # Add pairwise comparisons p-value
dev.off()

