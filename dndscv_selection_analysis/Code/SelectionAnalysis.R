#######################################################
####Signals of selection based on TMD status:
#######################################################
 


#######Load required packages:
library(TCGAbiolinks)
library(stringr)
library(dplyr)
library(maditr)
library(dndscv)
library(ggpubr)



#########################
#Load TMD classification:
#########################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("programme_scores_and_TMD_assignments.RData")
##Select extreme groups
table(prog_expr$TMD_three_categories)
prog_expr <- prog_expr[prog_expr$TMD_three_categories %in% c("YES","NO"),]
samples <- unique(as.character(prog_expr$SampleID))

########################
#Select CT of interest
########################
solid_tumours <- unique(as.character(prog_expr$cancer_type))


######################
####Download maf data:
######################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/MutationData/")

for (tumour in solid_tumours){
  print(paste0(tumour, ": This is ", which(tumour == solid_tumours), " out of ", length(solid_tumours)))    #Show progress
  df.maf <- GDCquery_Maf(tumour, pipelines = "muse")
  df.maf.short <- df.maf[,c("Tumor_Sample_Barcode",
                            "Chromosome",
                            "Start_Position",
                            "Reference_Allele",
                            "Tumor_Seq_Allele2")]
  df.maf.short$Cancer <- as.character(tumour)

  if (tumour == solid_tumours[1]){
    maf_all <- df.maf.short
  } else {
    maf_all <- rbind(maf_all, df.maf.short)
  }
}
save(maf_all, file = "maf_selection_analysis.RData")




###############################################
##Split maf into samples with high and low TMD:
###############################################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/MutationData/")
load("maf_selection_analysis.RData")
maf_all$Tumor_Sample_Barcode <- sapply(maf_all$Tumor_Sample_Barcode,
                                function(x) paste(str_split(x, "-", simplify = T)[1:4], collapse = "-"))

tmd.high <- prog_expr[prog_expr$TMD_three_categories %in% "YES",]
tmd.low <- prog_expr[prog_expr$TMD_three_categories %in% "NO",]
tmd.high <- as.character(tmd.high$SampleID)
tmd.low <- as.character(tmd.low$SampleID)
maf_all[,6] <- NULL
colnames(maf_all) <- c("SampleID", "chr", "pos", "ref", "mut")
maf_high <- maf_all[maf_all$SampleID %in% tmd.high,]
maf_low <- maf_all[maf_all$SampleID %in% tmd.low,]
save(maf_high, file = "maf_high_tmd_samples.RData")
save(maf_low, file = "maf_low_tmd_samples.RData")
rm(maf_all)




##################################################################
###Run dndscv (TMD samples)
##################################################################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/MutationData/")
load("maf_high_tmd_samples.RData")
maf_high <- data.frame(maf_high)
#Run the function
dndsout_high = dndscv(maf_high, refdb="RefCDS_human_GRCh38.p12.rda",cv = NULL)
sel_cv_high = dndsout_high$sel_cv
###Save results:
setwd("~/Documents/GitHub/tumourMassDormancy/dndscv_selection_analysis/")
save(sel_cv_high, file = "selection_results_high_TMD.RData")



##################################################################
###Run the tool (non-tmd samples)
##################################################################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/MutationData/")
load("maf_high_tmd_samples.RData")
maf_low <- data.frame(maf_low)
#Run the function
dndsout_low = dndscv(maf_low, refdb="RefCDS_human_GRCh38.p12.rda", cv = NULL)
sel_cv_low = dndsout_low$sel_cv
setwd("~/Documents/GitHub/tumourMassDormancy/dndscv_selection_analysis/")
###Save results:
save(sel_cv_low, file = "selection_results_low_TMD.RData")




###################################################################
##Select significant genes:
###################################################################

#Load the results:
setwd("~/Documents/GitHub/tumourMassDormancy/dndscv_selection_analysis/")
load("selection_results_high_TMD.RData")
load('selection_results_low_TMD.RData')
sig.high <- sel_cv_high[sel_cv_high$qmis_cv < 0.05,]
sig.high <- as.character(sig.high$gene_name)
sig.low <- sel_cv_low[sel_cv_low$qmis_cv < 0.05,]
sig.low <- as.character(sig.low$gene_name)
all_significant_genes <- unique(c(sig.high, sig.low))



####################################################################
###Combine results dataframes from the 2 runs
####################################################################
sel_cv_high <- sel_cv_high[,colnames(sel_cv_high) %in% c("wmis_cv","gene_name") ]
sel_cv_low <- sel_cv_low[,colnames(sel_cv_low) %in% c("wmis_cv","gene_name") ]
sel_cv_high <- sel_cv_high[sel_cv_high$gene_name %in% all_significant_genes,]
sel_cv_low <- sel_cv_low[sel_cv_low$gene_name %in% all_significant_genes,]
sel_cv_high$dnds_tmd <- sel_cv_high$wmis_cv
sel_cv_low$dnds_no_tmd <- sel_cv_low$wmis_cv
merged.data <- merge(sel_cv_high, sel_cv_low,
                     by.x = "gene_name",by.y = "gene_name")


##################################################################
###Identify genes selected only in samples with/witout TMD
#################################################################
tmd <- sig.high[!(sig.high %in% sig.low)]
notmd <- sig.low[!(sig.low %in% sig.high)]
merged.data$gene_class <- sapply(merged.data$gene_name, function(x)
  ifelse(x %in% tmd,"TMD",
         ifelse(x %in% notmd,"NO_TMD","BOTH")))
table(merged.data$gene_class) #7 genes selected in both conditions, 5 in samples without TMD and 17 in samples with TMD



###########################################
###Plot the results:
############################################
setwd("~/Documents/GitHub/tumourMassDormancy/dndscv_selection_analysis/Figures")
pdf("dnds_selection_results.pdf", width = 7, height = 7)
p <- ggscatter(merged.data, x = "dnds_tmd", y = "dnds_no_tmd",label = merged.data$gene_name,repel = TRUE, color = "gene_class")
p + labs(x = "dndc_tmd", y = "dnds_no_tmd")
dev.off()




