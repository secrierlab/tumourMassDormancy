######################################################################
######Programme correlation with immune infiltration estimates
########################################################################

#Load packages:
library(biomaRt)
library(ConsensusTME)
library(dplyr)
library(corrplot)
library(ggplot2)
library(ggpubr)
library(reshape)

########################
###Load expression data:
########################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/RNA_seq/")
load("combined_experssion_FPKM.RData")


#Consensus TME requires gene names to be in the hgnc_symbol format instead of ensembl gene id format therefore this needs to be converted:
genes <- colnames(combined_data)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
####Select cancer types of interest:
CT <- unique(as.character(combined_data$cancer_type))




#####################################################
##########Immune infiltration analysis (cancerwise)
#####################################################
#Select reference geneset
bindeaGeneSet <- ConsensusTME::methodSignatures$Bindea
for (a in CT) {
  
  print(a)
  
  expr.data <- combined_data[combined_data$cancer_type %in% a,]
  expr.data$cancer_type <- NULL
  expr.data <- data.frame(t(expr.data))
  expr.data$ensembl_gene_id <- rownames(expr.data)
  expr.data <- merge(expr.data, G_list,
                     by.x = "ensembl_gene_id", by.y = "ensembl_gene_id")
  #For duplicated entries take the mean value:
  expr.data <- expr.data %>% group_by(hgnc_symbol) %>% mutate_each(funs(mean)) %>% distinct
  expr.data <- data.frame(expr.data)
  rownames(expr.data) <- expr.data$hgnc_symbol
  expr.data$ensembl_gene_id <- NULL
  expr.data$hgnc_symbol <- NULL
  expr.data <- as.matrix(expr.data)

  
  ###Immune infiltration analysis:
  immune.infiltration <- ConsensusTME::geneSetEnrichment(expr.data, bindeaGeneSet)
  
  #Save the results:
  setwd("~/Documents/GitHub/tumourMassDormancy/Data/ImmuneInfiltrationEstimates/")
  save(immune.infiltration, file = paste("immune_infiltration_",a,".RData",sep = ""))
  
  
}



##################################################################
########Load immune infiltration estimates for all cancer types:
#################################################################
setwd("~/Documents/Wojciech_manuscript/ImmuneInfiltrationAnalysis")
load("immune_infiltration_ACC.RData")
combined.immune.infiltration <- immune.infiltration
CT <- CT[-1]
for (a in CT) {
  
  print(a)
  
  load(paste("immune_infiltration_",a,".RData",sep = ""))
  combined.immune.infiltration <- cbind(combined.immune.infiltration, immune.infiltration)
  
}

#####Clean up the dataframe
combined.immune.infiltration <- data.frame(t(combined.immune.infiltration))
rownames(combined.immune.infiltration) <- gsub('\\.', '-', rownames(combined.immune.infiltration))
#Select relevant columns and order them 
colnames(combined.immune.infiltration)
combined.immune.infiltration <- combined.immune.infiltration[,colnames(combined.immune.infiltration) %in% c("B_cells","T_gamma_delta_cells","Th1_cells","Th2_cells","Th17_cells","Treg_cells","T_central_memory_cells","T_effector_memory_cells","CD8_T_cells","NK_cells","plasmacytoid_Dendritic_cells","Macrophages","Neutrophils","Eosinophils","Mast_cells","Citotoxic_cells","T_helper_cells","Dendritic_cells","Angiogenesis")]
combined.immune.infiltration <- combined.immune.infiltration[, c("B_cells","T_gamma_delta_cells","Th1_cells","Th2_cells","Th17_cells","Treg_cells","T_central_memory_cells","T_effector_memory_cells","CD8_T_cells","NK_cells","plasmacytoid_Dendritic_cells","Macrophages","Neutrophils","Eosinophils","Mast_cells","Citotoxic_cells","T_helper_cells","Dendritic_cells","Angiogenesis")]
combined.immune.infiltration$Barcode <- rownames(combined.immune.infiltration)




###################################
##Add programme score information:
####################################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("program_expression_scores_updated.RData")
combined_data <- merge(combined.immune.infiltration, prog_expr,
                       by.x = "Barcode", by.y = "Barcode")

#Cleanup the dataframe
combined_data$SampleID <- NULL
combined_data$PatientID <- NULL
combined_data$Barcode <- NULL




##########################################################
##Programme correlation with immune infiltration estimates
#########################################################
corr_matrix <- cor(combined_data, use = "pairwise.complete.obs")
p_matrix <- cor.mtest(combined_data)$p
dimnames(p_matrix) <- dimnames(corr_matrix)
corr_matrix <- corr_matrix[20:24,1:19]
p_matrix <- p_matrix[20:24, 1:19]

setwd("~/Documents/GitHub/tumourMassDormancy/ImmuneInfiltrationAnalysis/Figures/")



col1 <- colorRampPalette(c("#59ac53","grey95","#8b5aa8"))
pdf(file = paste(Sys.Date(), "immune_correlation_prog.pdf", sep = "_"), height = 6, width = 12)
p <- corrplot(corr_matrix,
         tl.col = "black",
         p.mat = p_matrix,
         pch.cex = 0.8,
         pch.col = "firebrick",
         col = col1(100))
dev.off()



#######################################################
####Colour phate plots by immune infiltration estimates
#######################################################
setwd("~/Documents/GitHub/tumourMassDormancy/PHATE_analysis/ComBat_results/")
load("phate_coordinates_combat.RData")
combined.immune.infiltration$Barcode <- rownames(combined.immune.infiltration)
phate_coordinates_combat$Barcode <- rownames(phate_coordinates_combat)
phate_coordinates_combat$cluster <- as.factor(phate_coordinates_combat$cluster)
phate_coordinates_combat$cluster <- factor(recode(phate_coordinates_combat$cluster, "3" = "low", "2" = "high", "1" = "mid"), levels = c("low", "mid", "high"))
phate_coordinates_combat <- merge(phate_coordinates_combat, combined.immune.infiltration,
                                  by.x = "Barcode", by.y = "Barcode")


setwd("~/Documents/GitHub/tumourMassDormancy/ImmuneInfiltrationAnalysis/Figures/")
pdf("PHATE_coloured_by_cytotoxic_cell_infiltration.pdf", height = 5, width = 7)
ggplot(phate_coordinates_combat, aes(x=PHATE1, y=PHATE2, colour = Citotoxic_cells)) +
  geom_point() +
  theme(legend.key.width = unit(1, "cm")) +
  scale_color_gradient2(low = "#59ac53", midpoint = mean(phate_coordinates_combat$Citotoxic_cells),mid = "grey95", high = "#8b5aa8") + theme_classic()
dev.off()


pdf("PHATE_coloured_by_Dendritic_cells_infiltration.pdf", height = 5, width = 7)
ggplot(phate_coordinates_combat, aes(x=PHATE1, y=PHATE2, colour = Dendritic_cells)) +
  geom_point() +
  theme(legend.key.width = unit(1, "cm")) +
  scale_color_gradient2(low = "#59ac53", midpoint = mean(phate_coordinates_combat$Dendritic_cells),mid = "grey95", high = "#8b5aa8") + theme_classic()
dev.off()



pdf("PHATE_coloured_by_T_helper_cells_infiltration.pdf", height = 5, width = 7)
ggplot(phate_coordinates_combat, aes(x=PHATE1, y=PHATE2, colour = T_helper_cells)) +
  geom_point() +
  theme(legend.key.width = unit(1, "cm")) +
  scale_color_gradient2(low = "#59ac53", midpoint = mean(phate_coordinates_combat$T_helper_cells),mid = "grey95", high = "#8b5aa8") + theme_classic()
dev.off()


pdf("PHATE_coloured_by_Angiogenesis_infiltration.pdf", height = 5, width = 7)
ggplot(phate_coordinates_combat, aes(x=PHATE1, y=PHATE2, colour = Angiogenesis)) +
  geom_point() +
  theme(legend.key.width = unit(1, "cm")) +
  scale_color_gradient2(low = "#59ac53", midpoint = mean(phate_coordinates_combat$Angiogenesis),mid = "grey95", high = "#8b5aa8") + theme_classic()
dev.off()





###############################################################
####Boxplot infiltration estimate comparison between TMD groups
###############################################################

setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("programme_scores_and_TMD_assignments.RData")
prog_expr <- prog_expr[,colnames(prog_expr) %in% c("TMD_two_categories_detailed","Barcode")]
prog_expr <- merge(prog_expr, combined.immune.infiltration,
                                  by.x = "Barcode", by.y = "Barcode")
prog_expr <- prog_expr[,colnames(prog_expr) %in% c("Citotoxic_cells","Th1_cells","Th2_cells","Dendritic_cells","TMD_two_categories_detailed","Barcode")]
melt_annot <- melt(prog_expr, id.vars = c("Barcode","TMD_two_categories_detailed"))
melt_annot$TMD_two_categories_detailed <- factor(melt_annot$TMD_two_categories_detailed, levels = c("NO","Angiogenic Dormancy","Immunogenic Dormancy","Angiogenic and Immunogenic Dormancy"))


my_comparisons <- list(c("NO","Immunogenic Dormancy"), c("NO","Angiogenic Dormancy",c("NO","Angiogenic and Immunogenic Dormancy")))
p <- ggboxplot(melt_annot, x = "TMD_two_categories_detailed", y = "value",
               fill = "TMD_two_categories_detailed") + 
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + # Add pairwise comparisons p-value
  stat_compare_means(label.y.npc = "top") + theme_classic()
p <- facet(p, facet.by = "variable", nrow = 1, scales = "fixed") + scale_color_manual(values = c("Angiogenic and Immunogenic Dormancy" = "#1B9E77", "NO" = "grey95", "Angiogenic Dormancy" = "#D95F02", "Immunogenic Dormancy" = "#7570B3")) 
setwd("~/Documents/GitHub/tumourMassDormancy/ImmuneInfiltrationAnalysis/Figures/")
pdf("Boxplot_immune_infiltration_comparison_between_TMD_and_nonTMD_samples.pdf", height = 4,width = 12)
p
dev.off()






