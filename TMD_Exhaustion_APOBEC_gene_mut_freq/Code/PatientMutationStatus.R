#########################################################################################################
#TMD, Exhaustion and APOBEC gene mutation frequency in TCGA solid cancer cohort compared to COSMIC genes
#######################################################################################################

#Load packages
library(TCGAbiolinks)
library(maftools)
library(stringr)
library(maditr)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(svglite)
library(biomaRt)

##Select solid cancer types:
solid_tumours <- c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

#######################
####Define gene lists:
#######################

setwd("~/Documents/GitHub/tumourMassDormancy/Data/GeneLists/")

#Exhaustion gene list
exhaustion_list <- read.table("exhaustionMarkers.txt", header = FALSE, sep = "\t")
exhaustion_list <- as.character(exhaustion_list$V1)

#TMD gene list:
Gene_list <- read.table("dormancyMarkers_immunologic-angiogenic_updated.txt", header = TRUE,sep = "\t")
dormancy_list_upregulated <- Gene_list[Gene_list$Direction %in% "up",]
dormancy_list_upregulated <- as.character(dormancy_list_upregulated$Gene)
dormancy_list_downregulated <- Gene_list[Gene_list$Direction %in% "down",]
dormancy_list_downregulated <- as.character(dormancy_list_downregulated$Gene)
dormancy_list <- c(dormancy_list_upregulated, dormancy_list_downregulated)

##Immunological dormancy
immunologic_dormancy <- Gene_list[Gene_list$Type %in% "immunologic",]
immunologic_dormancy_list_upregualted <- immunologic_dormancy[immunologic_dormancy$Direction %in% "up",]
immunologic_dormancy_list_upregualted <- as.character(immunologic_dormancy_list_upregualted$Gene)
immunologic_dormancy_list_downregulated <- immunologic_dormancy[immunologic_dormancy$Direction %in% "down",]
immunologic_dormancy_list_downregulated <- as.character(immunologic_dormancy_list_downregulated$Gene)
immunologic_dormancy_list <- c(immunologic_dormancy_list_downregulated, immunologic_dormancy_list_upregualted)

##Angiogeneic dormancy
angiogenic_dormancy <- Gene_list[Gene_list$Type %in% "angiogenic",]
angiogenic_dormancy_upregulated <- angiogenic_dormancy[angiogenic_dormancy$Direction %in% "up",]
angiogenic_dormancy_upregulated <- as.character(angiogenic_dormancy_upregulated$Gene)
angiogenic_dormancy_downregulated <- angiogenic_dormancy[angiogenic_dormancy$Direction %in% "down",]
angiogenic_dormancy_downregulated <- as.character(angiogenic_dormancy_downregulated$Gene)
angiogenic_dormancy_list <- c(angiogenic_dormancy_downregulated,angiogenic_dormancy_upregulated)

#APOBEC genes
apobec_list <- c("AICDA", "APOBEC1", "APOBEC2", "APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G", "APOBEC3H", "APOBEC4")

#List of COSMIC genes:
setwd("~/Documents/GitHub/tumourMassDormancy/Data/GeneLists/")
COSMIC <- read.table("Census_allMon Feb 24 15_27_31 2020.tsv", header = TRUE,sep = "\t")
cosmic_list <- as.character(COSMIC$Gene.Symbol)


### Remove genes with missing expression data
setwd("~/Documents/GitHub/tumourMassDormancy/Data/RNA_seq/")
load("combined_experssion_FPKM.RData") #or load "example_FPKM_data.RData" for an example dataframe with 100 samples
all_genes_list <- c(dormancy_list, exhaustion_list, apobec_list)
#Convert HGNC symbols to ENSG
human = useMart("ensembl", dataset="hsapiens_gene_ensembl")
biomart_conversion <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), values=all_genes_list, mart=human, filters = "hgnc_symbol")
#Two ENSG numbers are reported for PDCD1 - use ENSG00000188389 NOT ENSG00000276977
#Two ENSG numbers are reported for APOBEC3A - use ENSG00000128383 NOT ENSG00000262156
biomart_conversion <- biomart_conversion[!(biomart_conversion$ensembl_gene_id %in% c("ENSG00000276977","ENSG00000262156")),]
all_genes_list_ENSG <- biomart_conversion$ensembl_gene_id
setdiff(all_genes_list_ENSG, intersect(all_genes_list_ENSG, colnames(combined_data)))
##All of the genes have reported expression
#



########################
####Download maf data:
#######################
setwd("~/Documents/GitHub/tumourMassDormancy/Data/Mutation_data/")
genes_to_save <- unique(c(all_genes_list, cosmic_list)) #Downloading genes from programmes and from COSMIC
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
save(maf_all, file = "maf_all.RData")




###########################################################
###Summary of mutational data for each patient:
###########################################################

setwd("~/Documents/GitHub/tumourMassDormancy/Data/Mutation_data/")
load("maf_all.RData")

maf_all$Patient_short <- sapply(maf_all$Tumor_Sample_Barcode,
                                function(x) paste(str_split(x, "-", simplify = T)[1:4], collapse = "-"))


maf_all.strict <- maf_all[ which((maf_all$Variant_Classification
                                  %in% c("Missense_Mutation",
                                         "Nonsense_Mutation",
                                         "Nonstop_Mutation",
                                         "Frameshift_Deletion",
                                         "Frameshift_Insertion",
                                         "Inframe_Insertion",
                                         "Inframe_Deletion")) &
                                   (grepl("deleterious",maf_all$SIFT)|
                                      grepl("damaging",maf_all$PolyPhen))),
]


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

maf_all.strict$Mutation <- 1      #Gives every patient-gene mutation value 1
maf_all.strict <- dcast(maf_all.strict, value.var = "Mutation",
                        formula = Tumor_Sample_Barcode+Patient_short+Cancer~Hugo_Symbol,
                        fun.aggregate = function(x) ifelse(sum(x)>0, yes = 1, no = 0))
#Because aggregation function sums all the values, we give 1 if sum>0,
#and 0 if there are no mutations

maf_all.lenient$Mutation <- 1
maf_all.lenient <- dcast(maf_all.lenient, value.var = "Mutation",
                         formula = Tumor_Sample_Barcode+Patient_short+Cancer~Hugo_Symbol,
                         fun.aggregate = function(x) ifelse(sum(x)>0, yes = 1, no = 0))

save(maf_all.lenient, file = "maf_all.lenient.RData")
save(maf_all.strict, file = "maf_all.strict.RData")





################################################
###Summary graph of mutations in cancer patients
################################################

setwd("~/Documents/GitHub/tumourMassDormancy/TMD_Exhaustion_APOBEC_gene_mut_freq/Figures/")
assign.programme <- function(gene){
  all <- c(apobec_list, exhaustion_list, dormancy_list)
  many <- all[duplicated(all)]
  if (gene %in% many) {
    as.factor("Dormancy + \n Exhaustion")
  } else if (gene %in% apobec_list) {
    as.factor("APOBEC")
  } else if (gene %in% exhaustion_list) {
    as.factor("Exhaustion")
  } else if (gene %in% dormancy_list) {
    as.factor("Dormancy")
  } else if (gene %in% cosmic_list) {
    as.factor("COSMIC")
  } else {
    as.factor("None")
  }
}

mafs <- list("Lenient" = maf_all.lenient, "Strict" = maf_all.strict)

for(i in 1:2){
  maf <- mafs[[i]]
  name <- names(mafs)[i]
  
  maf_summary <- maf[, -c(1,2)] %>% group_by(Cancer) %>% summarise_all(funs(sum))
  maf_summary_pan <- sort(colSums(maf_summary[, -1]), decreasing = T)
  
  maf_summary_gg <- melt(maf_summary, variable.name = "Gene", value.name = "Mutations")
  maf_summary_gg$Programme <- sapply(maf_summary_gg$Gene, function(x) assign.programme(x))
  maf_summary_gg$Cancer <- as.factor(maf_summary_gg$Cancer)
  maf_summary_gg$Gene <- as.factor(maf_summary_gg$Gene)
  
  p1 <- ggplot(subset(maf_summary_gg, Gene %in% all_genes_list), 
               aes(x = reorder(Gene, -Mutations), y = Mutations, group = Cancer, fill = Cancer)) +
    geom_bar(stat = "identity", position = "stack") +
    ylab("Number of patients with mutations") +
    xlab("Gene name") +
    ggtitle(label = "Programmes breakdown by cancer type", subtitle = paste0(name, " criteria")) +
    theme_pubr() +
    theme(legend.position = "top", legend.title = element_text(face = "bold", angle = 90)) +
    labs(fill = "TCGA cancer") +
    guides(fill = guide_legend(title.position = "left")) +
    theme(text = element_text(family = "Helvetica")) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
  
  p2 <- ggplot(subset(maf_summary_gg, Gene %in% c(all_genes_list, names(maf_summary_pan)[1:10])), 
               aes(x = reorder(Gene, -Mutations), y = Mutations, group = Cancer, fill = Programme)) +
    geom_bar(stat = "identity") +
    ylab("Number of patients with mutations") +
    xlab("Gene name") +
    ggtitle(label = "Comparison with top 10 COSMIC tier 1 and 2 genes") +
    theme_pubr() +
    theme(legend.position = "right", legend.title = element_text(face = "bold", angle = 90)) +
    labs(fill = "Programme") +
    guides(fill = guide_legend(title.position = "left")) +
    theme(text = element_text(family = "Helvetica")) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
  
  save_plot(p1, filename = paste0("maf_summary", name, "_1.svg"), base_width = 8)
  p1 <- p1 + theme(legend.position = "none")
  save_plot(p1, filename = paste0("maf_summary", name, "_1_nolegend.svg"), base_width = 8)
  
  save_plot(p2, filename = paste0("maf_summary", name, "_2.svg"), base_width = 8)
  p2 <- p2 + theme(legend.position = "none")
  save_plot(p2, filename = paste0("maf_summary", name, "_2_nolegend.svg"), base_width = 8)
}

rm(mafs, maf, maf_summary, maf_summary_gg, i, p1, p2)





