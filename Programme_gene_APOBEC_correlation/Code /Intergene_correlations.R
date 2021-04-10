######################################################
##Correlation between TMD/Exhaustion and APOBEC genes
#######################################################

#Load packages:
library(biomaRt)
library(gtools)
library(corrplot)


######################
####Define gene lists:
######################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/GeneLists/")
exhaustion_list <- read.table("exhaustionMarkers.txt", header = FALSE, sep = "\t")
exhaustion_list <- as.character(exhaustion_list$V1)

Gene_list <- read.table("dormancyMarkers_immunologic-angiogenic_updated.txt", header = TRUE,sep = "\t")

##Create a TMD gene list:
dormancy_list_upregulated <- Gene_list[Gene_list$Direction %in% "up",]
dormancy_list_upregulated <- as.character(dormancy_list_upregulated$Gene)
dormancy_list_downregulated <- Gene_list[Gene_list$Direction %in% "down",]
dormancy_list_downregulated <- as.character(dormancy_list_downregulated$Gene)
dormancy_list <- c(dormancy_list_upregulated, dormancy_list_downregulated)

##Create an immunologic dormancy gene list:
immunologic_dormancy <- Gene_list[Gene_list$Type %in% "immunologic",]
immunologic_dormancy_list_upregualted <- immunologic_dormancy[immunologic_dormancy$Direction %in% "up",]
immunologic_dormancy_list_upregualted <- as.character(immunologic_dormancy_list_upregualted$Gene)
immunologic_dormancy_list_downregulated <- immunologic_dormancy[immunologic_dormancy$Direction %in% "down",]
immunologic_dormancy_list_downregulated <- as.character(immunologic_dormancy_list_downregulated$Gene)
immunologic_dormancy_list <- c(immunologic_dormancy_list_downregulated, immunologic_dormancy_list_upregualted)


##Create an angiogenic dormancy gene list:
angiogenic_dormancy <- Gene_list[Gene_list$Type %in% "angiogenic",]
angiogenic_dormancy_upregulated <- angiogenic_dormancy[angiogenic_dormancy$Direction %in% "up",]
angiogenic_dormancy_upregulated <- as.character(angiogenic_dormancy_upregulated$Gene)
angiogenic_dormancy_downregulated <- angiogenic_dormancy[angiogenic_dormancy$Direction %in% "down",]
angiogenic_dormancy_downregulated <- as.character(angiogenic_dormancy_downregulated$Gene)
angiogenic_dormancy_list <- c(angiogenic_dormancy_downregulated,angiogenic_dormancy_upregulated)



#List of APOBEC genes:
apobec_list <- c("AICDA", "APOBEC1", "APOBEC2", "APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H", "APOBEC4")


all_genes <- c(exhaustion_list, immunologic_dormancy_list_upregualted, immunologic_dormancy_list_downregulated, angiogenic_dormancy_upregulated, angiogenic_dormancy_downregulated)



##############################
###Convert gene lists to ENSG 
##############################
all_genes_list <- c(dormancy_list, exhaustion_list, apobec_list)
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
biomart_conversion <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=all_genes_list, mart=human, filters = "hgnc_symbol")
#Two ENSG numbers are reported for PDCD1 - use ENSG00000188389 NOT ENSG00000276977
#Two ENSG numbers are reported for APOBEC3A - use ENSG00000128383 NOT ENSG00000262156
biomart_conversion <- biomart_conversion[!(biomart_conversion$ensembl_gene_id %in% c("ENSG00000276977","ENSG00000262156")),]



###########################
##Load the expression data:
###########################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/RNA_seq/")
load("combined_experssion_FPKM.RData") #or load "example_FPKM_data.RData" for an example data frame with 100 sample entries
rnaseq_cancer <- combined_data

for (i in all_genes_list) {
  
  selected_biomart <- biomart_conversion[biomart_conversion$hgnc_symbol %in% i,]
  a <- selected_biomart$ensembl_gene_id
  rnaseq_cancer[[i]] <- rnaseq_cancer[[a]]
  
  
}
rnaseq_cancer <- rnaseq_cancer[,colnames(rnaseq_cancer) %in% all_genes_list]
rnaseq_cancer <- data.frame(t(rnaseq_cancer))




##########################################
#Visualise correlation between gene sets:
##########################################

all_markers_list <- list("apobec" = apobec_list, "angiogenic_dormancy" = angiogenic_dormancy_list, "exhaustion" = exhaustion_list, "immunologic_dormancy" = immunologic_dormancy_list, "all_genes" = all_genes)
combs <- combinations(n=5, r=2)
setwd("~/Documents/GitHub/tumourMassDormancy/Programme_gene_APOBEC_correlation/Figures/")
for (n in 1:nrow(combs)){
  
  genes1 <- all_markers_list[[combs[n, 1]]]
  genes2 <- all_markers_list[[combs[n, 2]]]
  name1 <- names(all_markers_list[combs[n, 1]])
  name2 <- names(all_markers_list[combs[n, 2]])
  
  correlation <- cor(t(rnaseq_cancer[c(genes1, genes2), ]), use = "pairwise.complete.obs")
  p_matrix <- cor.mtest(t(rnaseq_cancer[c(genes1, genes2), ]), use = "pairwise.complete.obs")$p
  correlation <- data.frame(correlation)
  p_matrix <- data.frame(p_matrix)
  colnames(p_matrix) <- colnames(correlation)
  rownames(p_matrix) <- rownames(correlation)
  
  correlation <- as.matrix(correlation[genes1, genes2])
  p_matrix <- as.matrix(p_matrix[genes1, genes2])
  
  corrplot(correlation)
  
  
  pdf(file = paste(Sys.Date(), "correlation", name1, name2, ".pdf", sep = "_"), height = 6, width = 12)
  corrplot(correlation,
           method = "color",
           tl.col = "black",
           tl.cex = 1.2,
           cl.cex = 1.2,
           p.mat = p_matrix, 
           insig = "label_sig",
           pch.col = "black",
           pch.cex = 2.5)
  dev.off()
}

