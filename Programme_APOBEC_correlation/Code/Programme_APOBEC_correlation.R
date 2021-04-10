####################################################################################
#####Dormancy and exhaustion programmes correlation with mean APOBEC gene expression
####################################################################################


##Load the required packages:
library(biomaRt)
library(ggpubr)
library(gtools)
library(corrplot)



###############
###Load scores:
###############
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("programme_scores_and_TMD_assignments.RData")




###################################################################
##Correlation plots between  programme scores and APOBEC expression
###################################################################
setwd("~/Documents/GitHub/tumourMassDormancy/Programme_APOBEC_correlation/Figures/")
#Plot correlations
pdf("APOBEC_vs_Exhaustion.pdf", width = 5, height = 5)
p <- ggscatter(prog_expr, x = "mean_APOBEC", y = "Exhaustion",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
p + labs(x = "APOBEC", y = "Exhaustion")
dev.off()

pdf("APOBEC_vs_Dormancy.pdf", width = 5, height = 5)
p <- ggscatter(prog_expr, x = "mean_APOBEC", y = "TumourMassDormancy",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
p + labs(x = "APOBEC", y = "Tumour Mass Dormancy")
dev.off()

pdf("APOBEC_vs_Immune_dormancy.pdf", width = 5, height = 5)
p <- ggscatter(prog_expr, x = "mean_APOBEC", y = "ImmunogenicDormancy",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
p + labs(x = "APOBEC", y = "ImmunogenicDormancy")
dev.off()

pdf("APOBEC_vs_angiogenic_dormancy.pdf", width = 5, height = 5)
p <- ggscatter(prog_expr, x = "mean_APOBEC", y = "AngiogenicDormancy",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
p + labs(x = "APOBEC", y = "AngiogenicDormancy")
dev.off()






##############################
#Cancerwise correlation plots:
##############################

pdf("APOBEC_vs_Exhaustion_cancerwise.pdf", width = 13, height = 9)
p <- ggscatter(prog_expr, x = "mean_APOBEC", y = "Exhaustion",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
print(p + labs(x = "APOBEC", y = "Exhaustion") + facet_wrap(~cancer_type))
dev.off()

pdf("APOBEC_vs_Dormancy_cancerwise.pdf", width = 13, height = 9)
p <- ggscatter(prog_expr, x = "mean_APOBEC", y = "TumourMassDormancy",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
print(p + labs(x = "APOBEC", y = "Tumour Mass Dormancy") + facet_wrap(~cancer_type))
dev.off()


pdf("APOBEC_vs_Immune_dormancy_cancerwise.pdf", width = 13, height = 9)
p <- ggscatter(prog_expr, x = "mean_APOBEC", y = "ImmunogenicDormancy",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
print(p + labs(x = "APOBEC", y = "ImmunogenicDormancy") + facet_wrap(~cancer_type))
dev.off()


pdf("APOBEC_vs_angiogenic_dormancy_cancerwise.pdf", width = 13, height = 9)
p <- ggscatter(prog_expr, x = "mean_APOBEC", y = "AngiogenicDormancy",add = "reg.line", add.params = list(color = "blue", fill = "lightgray"), cor.coef = TRUE,
               conf.int = TRUE, alpha = 0.1,cor.coeff.args = list(method = "pearson", label.sep = "\n"))
print(p + labs(x = "APOBEC", y = "AngiogenicDormancy") + facet_wrap(~cancer_type))
dev.off()



