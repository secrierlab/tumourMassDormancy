setwd('~/Documents/GitHub/tumourMassDormancy/')

library(SummarizedExperiment)
library(corrplot)
library(ggplot2)
library(sva)

# Load hypoxia signature scores (hypoxia.scores) and programme scores (prog_expr)
load('Data/ProgrammeScores/Hypoxia_SigScores.Rdata')
load('Data/ProgrammeScores/programme_scores_and_TMD_assignments.RData')

# Merge programme and hypoxia scores
prog_hypoxia.merge <- merge(prog_expr, hypoxia.scores, by = 'Barcode')

# Hypoxia/Dormancy correlations by cancer + save values for volcano plots

dormancy.features <- c('mean_APOBEC','Exhaustion','ImmunogenicDormancy','AngiogenicDormancy','TumourMassDormancy')
hypoxia.features <- c('buffa','west','winter')

solid_tumours <- sort(names(table(prog_hypoxia.merge$Study)))
cor.stats <- matrix(NA, nrow = 0, ncol = 10)
for (tumour in solid_tumours) {
  print(paste0('Running hypoxia/dormancy correlations for tumour: ', tumour))
  pcc.tumour <- prog_hypoxia.merge[prog_hypoxia.merge$Study == tumour, ]
  SBS_matrix <- cor(pcc.tumour[,c(dormancy.features, hypoxia.features)])
  p_matrix <- cor.mtest(pcc.tumour[,c(dormancy.features, hypoxia.features)])$p
  dimnames(p_matrix) <- dimnames(SBS_matrix)
  
  SBSmat <- SBS_matrix[hypoxia.features, dormancy.features]
  pmat <- p_matrix[hypoxia.features, dormancy.features]
  
  # Row to add
  row.to.add <- c(as.numeric(SBSmat['buffa',]), as.numeric(pmat['buffa',]))
  cor.stats <- rbind(cor.stats, row.to.add)
  
}
rownames(cor.stats) <- solid_tumours
colnames(cor.stats) <- c(dormancy.features, paste0(dormancy.features,'.pval'))

# Adjust p-values via a BH correction to account for multiple testing
for (col in paste0(dormancy.features,'.pval')) {
  cor.stats[,col] <- p.adjust(cor.stats[, col], method = 'BH')
}

cor.stats <- data.frame(cor.stats)
cor.stats$Study <- rownames(cor.stats)

# Volcano plots for TumourMassDormancy, ImmunogenicDormancy, and AngiogenicDormancy
p_tmd <- ggplot(cor.stats, aes(x = TumourMassDormancy, y = -log10(TumourMassDormancy.pval))) + geom_point() + 
  geom_text(aes(label = ifelse(cor.stats$TumourMassDormancy.pval < .05, Study, NA)), nudge_y = -0.2) + 
  labs(x = 'Tumour Mass Dormancy Correlation', y = '-log10(Adjusted Correlation Significance)') +
  xlim(-.6, .6) + ylim(0, 8) + theme_bw() +
  geom_hline(yintercept = -log10(0.05), color = 'red', linetype = 'dashed') + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed')
print(p_tmd)
ggsave(filename = 'Hypoxia/Figures/Corr_by_Study_TMD.pdf', p_tmd)

p_immuno <- ggplot(cor.stats, aes(x = ImmunogenicDormancy, y = -log10(ImmunogenicDormancy.pval))) + geom_point() + 
  geom_text(aes(label = ifelse(cor.stats$ImmunogenicDormancy.pval < .05, Study, NA)), nudge_y = -0.2) + 
  labs(x = 'Immunogenic Dormancy Correlation', y = '-log10(Adjusted Correlation Significance)') +
  xlim(-.6, .6) + ylim(0, 8) + theme_bw() +
  geom_hline(yintercept = -log10(0.05), color = 'red', linetype = 'dashed') + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed')
print(p_immuno)
ggsave(filename = 'Hypoxia/Figures/Corr_by_Study_Immuno.pdf', p_immuno)

p_angio <- ggplot(cor.stats, aes(x = AngiogenicDormancy, y = -log10(AngiogenicDormancy.pval))) + geom_point() + 
  geom_text(aes(label = ifelse(cor.stats$AngiogenicDormancy.pval < .05, Study, NA)), nudge_y = -0.2) + 
  labs(x = 'Angiogenic Dormancy Correlation', y = '-log10(Adjusted Correlation Significance)') +
  xlim(-.6, .6) + ylim(0, 8) + theme_bw() +
  geom_hline(yintercept = -log10(0.05), color = 'red', linetype = 'dashed') + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed')
print(p_angio)
ggsave(filename = 'Hypoxia/Figures/Corr_by_Study_Angio.pdf', p_angio)


