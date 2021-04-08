###################################################################################
#####PHATE analysis using TMD and exhaustion genes
###################################################################################

###############
##Load packages:
library(biomaRt)
library(phateR)
library(data.table)
library(ggplot2)


#######################
####Define gene lists:
#######################
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


####################################################
##Biomart conversion - obtain ENSG IDS for all genes
###################################################
###Convert gene lists to ENSG 
all_genes_list <- c(dormancy_list, exhaustion_list)
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
biomart_conversion <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=all_genes_list, mart=human, filters = "hgnc_symbol")
#Two ENSG numbers are reported for PDCD1 - use ENSG00000188389 NOT ENSG00000276977
biomart_conversion <- biomart_conversion[!(biomart_conversion$ensembl_gene_id %in% c("ENSG00000276977")),]


#################################################
##Load the expression data (no ComBat correction)
#################################################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/RNA_seq/")
load("combined_experssion_FPKM.RData") #or load "example_FPKM_data.RData" for an example dataframe with 100 samples
rnaseq_cancer <- combined_data
for (i in all_genes_list) {
  
  selected_biomart <- biomart_conversion[biomart_conversion$hgnc_symbol %in% i,]
  a <- selected_biomart$ensembl_gene_id
  rnaseq_cancer[[i]] <- rnaseq_cancer[[a]]
  
  
}
rnaseq_cancer <- rnaseq_cancer[,colnames(rnaseq_cancer) %in% c(all_genes_list,"cancer_type")]


#######################################
###Load scores and merge with expr data
######################################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("program_expression_scores_updated.RData")
#Merge dataframe:
rnaseq_cancer$Barcode <- rownames(rnaseq_cancer)
annotation <- merge(prog_expr, rnaseq_cancer,
                    by.x = "Barcode", by.y = "Barcode")
#remove columns with gene expression:
annotation[,10:56] <- NULL
######Make sure expression data has only samples with score annotation
rnaseq_cancer <- rnaseq_cancer[rownames(rnaseq_cancer) %in% annotation$Barcode,]



#########################################
##Run phate using non-combat treated data
########################################
#Select data for TMD and exhaustion genes only 
rnaseq_cancer <- rnaseq_cancer[,colnames(rnaseq_cancer) %in% c(dormancy_list, exhaustion_list)]
setwd("~/Documents/GitHub/tumourMassDormancy/Data/RNA_seq/")
save(rnaseq_cancer, file = "PHATE_expression_data_no_COMBAT.RData")
load("PHATE_expression_data_no_COMBAT.RData") 

#Run PHATE:
setwd("~/Documents/GitHub/tumourMassDormancy/PHATE_analysis/NoComBat_results/")
phate_expr <- as.matrix(rnaseq_cancer)
set.seed(123456)
data_phate <- phate(phate_expr)
phate_coordinates <- data.frame(data_phate$embedding)
save(phate_coordinates, file = "phate_coordinates_no_combat.RData")





######################################
##Load expression data (ComBat treated)
######################################
setwd("~/Documents/TMD_manuscript_data/RNA_seq/")
expr.data.combat <- fread("TCGA_combat_tumor_type_correction.txt", sep = "\t") #Or load "example_FPKM_data_COMBAT.RData" for an example with 100 samples
expr.data.combat <- expr.data.combat[which(expr.data.combat$genes %in% c(dormancy_list, exhaustion_list)),]
all_genes <- expr.data.combat$genes
expr.data.combat$genes <- NULL
expr.data.combat <- as.data.frame(expr.data.combat)
rownames(expr.data.combat) <- all_genes
#Transpose the dataframe
expr.data.combat <- data.frame(t(expr.data.combat))
#Make sure that you select primary tumour samples:
rownames(expr.data.combat) <- gsub('\\.', '-', rownames(expr.data.combat))
expr.data.combat$SampleType <- sapply(rownames(expr.data.combat), function(x)
  strsplit(x,"-")[[1]][5])
expr.data.combat <- expr.data.combat[which(expr.data.combat$SampleType %in% c("01A","01B","01C")),]
expr.data.combat$SampleType <- NULL
#remove cancer type annotation from the barcode:
rownames(expr.data.combat) <- sapply(rownames(expr.data.combat), function(x)
  paste(strsplit(x,"-")[[1]][2:8],collapse="-"))
#Select only samples which have calculated programme scores:
expr.data.combat <- expr.data.combat[rownames(expr.data.combat) %in% annotation$Barcode,]



#########################################
##Run phate using Combat treated data
########################################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/RNA_seq/")
save(expr.data.combat, file = "PHATE_expression_data_COMBAT.RData")
load("PHATE_expression_data_COMBAT.RData") 

#Run PHATE
setwd("~/Documents/GitHub/tumourMassDormancy/PHATE_analysis/ComBat_results/")
phate_expr <- as.matrix(expr.data.combat)
set.seed(123456)
data_phate <- phate(phate_expr)
phate_coordinates_combat <- data.frame(data_phate$embedding)
set.seed(123)
clusters <- kmeans(phate_coordinates_combat, 3)
phate_coordinates_combat$cluster <- clusters$cluster
save(phate_coordinates_combat, file = "phate_coordinates_combat.RData")







###########################################
#Plot PHATE coordinates
###########################################


################################################################
#Plots coloured by cancer type
#Plot PHATE coodinates without COMBAT (coloured by cancer type)
setwd("~/Documents/GitHub/tumourMassDormancy/PHATE_analysis/NoComBat_results/")
load("phate_coordinates_no_combat.RData")
phate_coordinates$Barcode <- rownames(phate_coordinates)
phate_coordinates <- merge(phate_coordinates, annotation,
                           by.x = "Barcode", by.y = "Barcode")
setwd("~/Documents/GitHub/tumourMassDormancy/PHATE_analysis/Figures/")
pdf("PHATE_no_combat_colour_by_cancertype.pdf",height = 5, width = 8)
ggplot(phate_coordinates, aes(x=PHATE1, y=PHATE2, color = cancer_type)) +
  geom_point() +
  theme(legend.title = element_blank(), legend.key.width = unit(1, "cm")) +
  guides(color = guide_legend(override.aes = list(size=5))) + theme_classic()
dev.off()


##Plot PHATE coodinates with COMBAT (coloured by cancer type)
setwd("~/Documents/GitHub/tumourMassDormancy/PHATE_analysis/ComBat_results/")
load("phate_coordinates_combat.RData")
phate_coordinates_combat$Barcode <- rownames(phate_coordinates_combat)
phate_coordinates_combat <- merge(phate_coordinates_combat, annotation,
                                  by.x = "Barcode", by.y = "Barcode")
setwd("~/Documents/GitHub/tumourMassDormancy/PHATE_analysis/Figures/")
pdf("PHATE_combat_colour_by_cancertype.pdf",height = 5, width = 8)
ggplot(phate_coordinates_combat, aes(x=PHATE1, y=PHATE2, color = cancer_type)) +
  geom_point() +
  theme(legend.title = element_blank(), legend.key.width = unit(1, "cm")) +
  guides(color = guide_legend(override.aes = list(size=5))) + theme_classic()
dev.off()
################################################################################



################################################################################
###Color the combat corrected plot by different programme scores:
setwd("~/Documents/GitHub/tumourMassDormancy/PHATE_analysis/Figures/")
pdf("PHATE_combat_colour_by_TMD_programme.pdf",height = 5, width = 8)
ggplot(phate_coordinates_combat, aes(x=PHATE1, y=PHATE2, colour = TumourMassDormancy)) +
  geom_point() + theme_classic() +
  theme(legend.title = element_blank(), legend.key.width = unit(1, "cm")) +
  scale_color_gradient2(low = "green", midpoint = 0.3, mid = "blue", high = "red") +
  theme(legend.position = "right", legend.title = element_text(face = "bold", angle = 90)) +
  guides(colour = guide_colorbar(title.position = "left")) +
  theme(text = element_text(family = "Helvetica", size = 10))
dev.off()

pdf("PHATE_combat_colour_by_exhaustion_programme.pdf",height = 5, width = 8)
ggplot(phate_coordinates_combat, aes(x=PHATE1, y=PHATE2, colour = Exhaustion)) +
  geom_point() + theme_classic() +
  theme(legend.title = element_blank(), legend.key.width = unit(1, "cm")) +
  scale_color_gradient2(low = "green", midpoint = mean(phate_coordinates_combat$Exhaustion), mid = "blue", high = "red") +
  theme(legend.position = "right", legend.title = element_text(face = "bold", angle = 90)) +
  guides(colour = guide_colorbar(title.position = "left")) +
  theme(text = element_text(family = "Helvetica", size = 10))
dev.off()

pdf("PHATE_combat_colour_by_apobec_programme.pdf",height = 5, width = 8)
ggplot(phate_coordinates_combat, aes(x=PHATE1, y=PHATE2, colour = mean_APOBEC)) +
  geom_point() + theme_classic() +
  theme(legend.title = element_blank(), legend.key.width = unit(1, "cm")) +
  scale_color_gradient2(low = "green", midpoint = mean(phate_coordinates_combat$mean_APOBEC), mid = "blue", high = "red") +
  theme(legend.position = "right", legend.title = element_text(face = "bold", angle = 90)) +
  guides(colour = guide_colorbar(title.position = "left")) +
  theme(text = element_text(family = "Helvetica", size = 10))
dev.off()
##############################################################################################




#######################################################
##PlotPHATE coordinates coloured by TMD  assignments:
setwd("~/Documents/GitHub/tumourMassDormancy/PHATE_analysis/ComBat_results/")
load("phate_coordinates_combat.RData")
phate_coordinates_combat$Barcode <- rownames(phate_coordinates_combat)
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("programme_scores_and_TMD_assignments.RData")
phate_coordinates_combat <- merge(phate_coordinates_combat, prog_expr, 
                                  by.x = "Barcode", by.y = "Barcode")
setwd("~/Documents/GitHub/tumourMassDormancy/PHATE_analysis/Figures/")
pdf("PHATE_coordinates_coloured_by_TMD_assignments.pdf",height = 5, width = 8)
ggplot(phate_coordinates_combat, aes(x=PHATE1, y=PHATE2, colour = TMD_two_categories)) +
  geom_point(size = 1) + theme_classic() +
  theme(legend.title = element_blank(), legend.key.width = unit(1, "cm")) + 
  scale_colour_manual(values = c("NO" = "grey90", "YES" = "darkmagenta")) +
  theme(legend.position = "right", legend.title = element_text(face = "bold", angle = 90)) +
  theme(text = element_text(family = "Helvetica", size = 10))
dev.off()

