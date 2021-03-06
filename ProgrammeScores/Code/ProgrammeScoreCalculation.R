##################################################################################################
####Calculating Scores for APOBEC, Exhaustion, TMD, immunological and angiogenic dormancy programmes:
##################################################################################################

##Load the required packages:
library(biomaRt)
library(ggplot2)
library(stringr)

#######################
####Define gene lists:
#######################

setwd("~/Documents/GitHub/tumourMassDormancy/Data/GeneLists/")

#####Exhaustion gene list:
exhaustion_list <- read.table("exhaustionMarkers.txt", header = FALSE, sep = "\t")
exhaustion_list <- as.character(exhaustion_list$V1)
##TMD gene list:
Gene_list <- read.table("dormancyMarkers_immunologic-angiogenic_updated.txt", header = TRUE,sep = "\t")
dormancy_list_upregulated <- Gene_list[Gene_list$Direction %in% "up",]
dormancy_list_upregulated <- as.character(dormancy_list_upregulated$Gene)
dormancy_list_downregulated <- Gene_list[Gene_list$Direction %in% "down",]
dormancy_list_downregulated <- as.character(dormancy_list_downregulated$Gene)
dormancy_list <- c(dormancy_list_upregulated, dormancy_list_downregulated)
##Immunological dormancy gene list:
immunologic_dormancy <- Gene_list[Gene_list$Type %in% "immunologic",]
immunologic_dormancy_list_upregualted <- immunologic_dormancy[immunologic_dormancy$Direction %in% "up",]
immunologic_dormancy_list_upregualted <- as.character(immunologic_dormancy_list_upregualted$Gene)
immunologic_dormancy_list_downregulated <- immunologic_dormancy[immunologic_dormancy$Direction %in% "down",]
immunologic_dormancy_list_downregulated <- as.character(immunologic_dormancy_list_downregulated$Gene)
immunologic_dormancy_list <- c(immunologic_dormancy_list_downregulated, immunologic_dormancy_list_upregualted)
##Angiogenic dormancy gene list:
angiogenic_dormancy <- Gene_list[Gene_list$Type %in% "angiogenic",]
angiogenic_dormancy_upregulated <- angiogenic_dormancy[angiogenic_dormancy$Direction %in% "up",]
angiogenic_dormancy_upregulated <- as.character(angiogenic_dormancy_upregulated$Gene)
angiogenic_dormancy_downregulated <- angiogenic_dormancy[angiogenic_dormancy$Direction %in% "down",]
angiogenic_dormancy_downregulated <- as.character(angiogenic_dormancy_downregulated$Gene)
angiogenic_dormancy_list <- c(angiogenic_dormancy_downregulated,angiogenic_dormancy_upregulated)
#Apobec programme gene list:
apobec_list <- c("AICDA", "APOBEC1", "APOBEC2", "APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H", "APOBEC4")
#COSMIC genes:
COSMIC <- read.table("Census_allMon Feb 24 15_27_31 2020.tsv", header = TRUE,sep = "\t")
cosmic_list <- as.character(COSMIC$Gene.Symbol)


####################################
###Obtain ENSG ID for programme genes
#####################################
all_genes_list <- c(dormancy_list, exhaustion_list, apobec_list)
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
biomart_conversion <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=all_genes_list, mart=human, filters = "hgnc_symbol")
#Two ENSG numbers are reported for PDCD1 - use ENSG00000188389 NOT ENSG00000276977
#Two ENSG numbers are reported for APOBEC3A - use ENSG00000128383 NOT ENSG00000262156
biomart_conversion <- biomart_conversion[!(biomart_conversion$ensembl_gene_id %in% c("ENSG00000276977","ENSG00000262156")),]


#########################################################
##Load the expression data and change ENSG to HGNC_symbol
#########################################################
#(data is log2 transformed)
setwd("~/Documents/GitHub/tumourMassDormancy/Data/RNA_seq/")
load("combined_experssion_FPKM.RData") # or load "example_FPKM_data.RData" for a shortened data frame with 100 samples
rnaseq_cancer <- combined_data
for (i in all_genes_list) {
  
  selected_biomart <- biomart_conversion[biomart_conversion$hgnc_symbol %in% i,]
  a <- selected_biomart$ensembl_gene_id
  rnaseq_cancer[[i]] <- rnaseq_cancer[[a]]
  
  
}
rnaseq_cancer <- rnaseq_cancer[,colnames(rnaseq_cancer) %in% all_genes_list]
rnaseq_cancer <- data.frame(t(rnaseq_cancer))



#######################################
#Calculate a score for each of the programmes:
#######################################

prog_expr <- data.frame(row.names = colnames(rnaseq_cancer))

prog_expr$mean_APOBEC <- apply(rnaseq_cancer[apobec_list, ], function (x) mean(x, na.rm = T), MARGIN = 2)

prog_expr$Exhaustion <- apply(rnaseq_cancer[exhaustion_list, ], function (x) mean(x), MARGIN = 2)

prog_expr$ImmunogenicDromancyUp <- colSums(rnaseq_cancer[rownames(rnaseq_cancer) %in% immunologic_dormancy_list_upregualted,])
prog_expr$ImmunogenicDromancyDown <- colSums(rnaseq_cancer[rownames(rnaseq_cancer) %in% immunologic_dormancy_list_downregulated,])
prog_expr$ImmunogenicDormancy <- (prog_expr$ImmunogenicDromancyUp - prog_expr$ImmunogenicDromancyDown) / 18
prog_expr$ImmunogenicDromancyDown <- NULL
prog_expr$ImmunogenicDromancyUp <- NULL

prog_expr$AngiogenicDormancyUp <- colSums(rnaseq_cancer[rownames(rnaseq_cancer) %in% angiogenic_dormancy_upregulated,])
prog_expr$AngiogenicDormancyDown <- colSums(rnaseq_cancer[rownames(rnaseq_cancer) %in% angiogenic_dormancy_downregulated,])
prog_expr$AngiogenicDormancy <- (prog_expr$AngiogenicDormancyUp - prog_expr$AngiogenicDormancyDown) / 17
prog_expr$AngiogenicDormancyUp <- NULL
prog_expr$AngiogenicDormancyDown <- NULL

prog_expr$TMDUp <- colSums(rnaseq_cancer[rownames(rnaseq_cancer) %in% dormancy_list_upregulated,])
prog_expr$TMDDown <- colSums(rnaseq_cancer[rownames(rnaseq_cancer) %in% dormancy_list_downregulated,])
prog_expr$TumourMassDormancy <- (prog_expr$TMDUp - prog_expr$TMDDown) / 35
prog_expr$TMDDown <- NULL
prog_expr$TMDUp <- NULL
rownames(prog_expr) <- gsub('\\.', '-', rownames(prog_expr))



#####################################################################
######Randomly select one RNA-seq entry per Patient:
######################################################################
prog_expr$Barcode <- rownames(prog_expr)
prog_expr$SampleID <- sapply(prog_expr$Barcode, function(x)
  paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
prog_expr$PatientID <- sapply(prog_expr$Barcode, function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
Patients <- as.character(prog_expr$PatientID)
Barcodes <- as.character(prog_expr$Barcode)
uniquePatients <- unique(as.character(prog_expr$PatientID))
selection.df <- data.frame(Barcodes, Patients)
Barcodes <- NULL
set.seed(123)
for (i in uniquePatients) {
  print(i)
  selected_data <- selection.df[selection.df$Patients %in% i,]
  selected_data <- as.character(selected_data$Barcodes)
  random_sample <- sample(selected_data, size = 1, replace = FALSE)
  Barcodes <- c(Barcodes, random_sample)
  rm(random_sample)
}
prog_expr <- prog_expr[prog_expr$Barcode %in% Barcodes,]
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
save(prog_expr, file = "program_expression_scores_updated.RData")



################################################
####Plot density distribution of the programme scores
#################################################
#Load programme scores:
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("program_expression_scores_updated.RData")

#########Plots:
setwd("~/Documents/GitHub/tumourMassDormancy/ProgrammeScores/Figures/")

pdf("APOBEC_score_density_plot.pdf",height = 5,width = 5)
ggplot(prog_expr, aes(x=mean_APOBEC))+
  geom_density(color="darkblue", fill="lightblue") + theme_classic()
dev.off()

pdf("TMD_score_density_plot.pdf",height = 5,width = 5)
ggplot(prog_expr, aes(x=TumourMassDormancy))+
  geom_density(color="darkblue", fill="lightblue") + theme_classic()
dev.off()

pdf("Exhaustion_score_density_plot.pdf",height = 5,width = 5)
ggplot(prog_expr, aes(x=Exhaustion))+
  geom_density(color="darkblue", fill="lightblue") + theme_classic()
dev.off()

pdf("Angiogenic_dormancy_score_density_plot.pdf",height = 5,width = 5)
ggplot(prog_expr, aes(x=AngiogenicDormancy))+
  geom_density(color="darkblue", fill="lightblue") + theme_classic()
dev.off()

pdf("Immunogenic_dormancy_score_density_plot.pdf",height = 5,width = 5)
ggplot(prog_expr, aes(x=ImmunogenicDormancy))+
  geom_density(color="darkblue", fill="lightblue") + theme_classic()
dev.off()








