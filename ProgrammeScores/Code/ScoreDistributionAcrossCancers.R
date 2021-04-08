###############################################################################################
###Visualising variation in TMD, angiogenic and immunogneic dormancy scores across cancer types:
################################################################################################



###Load packages:
library(reshape)
library(ggplot2)


###############
#Load scores:
###############
setwd("~/Documents/GitHub/tumourMassDormancy/Data/ProgrammeScores/")
load("programme_scores_and_TMD_assignments.RData")



####################################################
##Order cancer types according to the average TMD score
#######################################################
CT <- unique(prog_expr$cancer_type)
mean <- NULL
for (i in CT) {
  
  selected.data <- prog_expr[prog_expr$cancer_type %in% i,]
  x <- median(selected.data$TumourMassDormancy)
  mean <- c(mean, x)
}
means <- data.frame(CT, mean)
means <- means[order(means$mean),]
Ordered_cancers <- means$CT

####################################
###Remove unnecessary columns and reshape 
#####################################
prog_expr <- prog_expr[,colnames(prog_expr) %in% c("ImmunogenicDormancy","AngiogenicDormancy","Barcode","cancer_type")]
prog_expr <- melt(prog_expr, id.vars = c("Barcode","cancer_type"))
prog_expr$cancer_type <- factor(prog_expr$cancer_type,
                                levels = Ordered_cancers)


#####
#Plot:
#####
setwd("~/Documents/GitHub/tumourMassDormancy/ProgrammeScores/Figures/")
pdf("immunogenic_and_angiogenic_dormancy_score_across_tissue_types.pdf", width = 10, height = 3)
p <- ggplot(prog_expr, aes(x=cancer_type, y=value, fill=variable)) + 
  geom_boxplot() + scale_color_manual(values = c("AngiogenicDormancy" = "#E7B800", "ImmunogenicDormancy" = "#FC4E07")) + theme_classic() +
  ggtitle("") +
  xlab("Cancer Type") + ylab("Programme Score")
p + rotate_x_text(45) 
dev.off()




