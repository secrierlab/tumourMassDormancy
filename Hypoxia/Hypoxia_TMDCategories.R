library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(gridExtra)

load('Data/Hypoxia_SigScores.Rdata') # Hypoxia signature scores
load('Data/programme_scores_and_TMD_assignments.RData') # TMD programme scores

# Merge programmes, and re-label TMD_three_categories_detailed for ease
df.full <- cbind(prog_expr, hypoxia.scores[,c('buffa','west','winter')])
df.full$TMD_three_categories_detailed <- sapply(df.full$TMD_three_categories_detailed,
                                                function(x) ifelse(x == 'Angiogenic and Immunogenic Dormancy','A/ID',
                                                                   ifelse(x == 'Angiogenic Dormancy', 'AD',
                                                                          ifelse(x == 'Immunogenic Dormancy', 'ID',x))))

# Boxplot comparing hypoxia distributions across the five TMD_three_categories_detailed groups
p_3d <- ggplot(df.full, aes(x = TMD_three_categories_detailed, y = buffa, color = TMD_three_categories_detailed)) + 
  geom_boxplot() + labs(x = 'Tumour Mass Dormancy (Category)', y = 'Hypoxia Score (Buffa)')

# Add significance to p_3d
p_3d.comparisons <- list(c('MID','NO'), c('ID','NO'), c('AD','NO'), c('A/ID','NO'))
p_3d <- p_3d + stat_compare_means(comparisons = p_3d.comparisons, aes(label = ..p.signif..)) +
  scale_color_manual(values = c('#1B9E77','#D95F02','#7570B3','grey','grey50')) +
  theme_bw() + theme(legend.position = 'none')
print(p_3d)

ggsave(filename = 'Figures/Hypoxia_by_TMD.pdf', plot = p_3d)

# Alter Buffa thresholds to cluster to the nearest 10
df.full$Buffa.Group <- sapply(df.full$buffa, function(x) 10 * round(x/10))
df.full$Buffa.Group <- as.factor(df.full$Buffa.Group)

# Determine ranking of TCGA studies by TMD prevalence
tmd.by.cancer <- table(df.full$cancer_type, df.full$TMD_three_categories)
tmd.by.cancer <- data.frame(YES = tmd.by.cancer[,'YES'],
                            NO = tmd.by.cancer[,'NO'] + tmd.by.cancer[,'MID'])
tmd.by.cancer$TMD_Perc <- tmd.by.cancer$YES/(tmd.by.cancer$YES + tmd.by.cancer$NO)
tmd.by.cancer <- tmd.by.cancer[order(tmd.by.cancer$TMD_Perc, decreasing = TRUE),]

df.full.plot <- transform(df.full,
                          cancer_type = factor(cancer_type, levels = rownames(tmd.by.cancer)))

# Full plot (excluding KICH)
# KICH is removed on account of 64/65 samples being labelled as 'MID'
p_full <- ggplot(df.full.plot[df.full.plot$cancer_type != 'KICH',], aes(x = Buffa.Group, fill = TMD_three_categories_detailed)) +
  geom_bar(position = 'fill') +
  labs(x = 'Hypoxia Group', y = 'TMD Group') + 
  scale_x_discrete(drop = FALSE, breaks = c(-50, 0, 50)) +
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_manual('TMD Group', values = c('A/ID' = '#1B9E77',
                                            'AD' = '#D95F02',
                                            'ID' = '#7570B3',
                                            'MID' = 'grey',
                                            'NO' = 'grey50')) +
  theme(legend.position = 'none')

print(p_full)
ggsave(filename = 'Figures/TMD_by_HypoxiaGroup_full.pdf', plot = p_full)

# Facet plot (excluding KICH)
p_facet <- p_full + facet_wrap(~ cancer_type, nrow = 5, ncol = 6) + theme(legend.position = 'top')
print(p_facet)
ggsave(filename = 'Figures/TMD_by_HypoxiaGroup_facet.pdf', plot = p_facet)
