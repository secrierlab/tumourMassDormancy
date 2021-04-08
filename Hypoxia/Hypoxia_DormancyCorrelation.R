library(SummarizedExperiment)
library(corrplot)
library(ggplot2)
library(sva)

# Load immunological dormancy gene list
dormancy.gene.list <- read.table('Data/dormancyMarkers_immunologic-angiogenic_updated2.txt', header = TRUE)
genes.dormancy <- dormancy.gene.list$Gene

# Load exhaustion gene list
exhaustion.gene.list <- read.table('Data/exhaustionMarkers.txt')
genes.exhaustion <- exhaustion.gene.list$V1

# Load hypoxia signature scores (hypoxia.scores) and programme scores (prog_expr)
load('Data/Hypoxia_SigScores.Rdata')
load('Data/programme_scores_and_TMD_assignments.RData')

# Combine TMD programme genes
all_genes_list <- unique(c(genes.dormancy, genes.exhaustion))

# Load FPKM normalized expression data
se.list <- list()
col.req <- c('barcode','patient','sample','Study','definition')
se.files <- list.files(path = '../../Data/GeneExpression/FPKM_Normalized/',
                       pattern = 'Rdata',
                       full.names = TRUE)

for (file in se.files) {
  study <- strsplit(file, '[.]')[[1]][6]
  print(paste0('Loading expression data from study: ',study,'...'))
  load(file)
  se$Study <- study
  colData(se) <- colData(se)[,col.req]
  se <- se[,se$definition == 'Primary solid Tumor']
  se <- se[rowData(se)$external_gene_name %in% c(genes.dormancy,genes.exhaustion),]
  
  se.list[[study]] <- se
}

se.comb <- do.call(cbind, unlist(se.list))
se.comb <- se.comb[, prog_expr$Barcode]

se.comb.assay <- assay(se.comb)
rownames(se.comb.assay) <- rowData(se.comb)$external_gene_name

se.comb.log2 <- log2(se.comb.assay + 1) # Log2 normalization

# COMBAT correction
combat_cancer <- ComBat(dat = se.comb.log2,
                        batch = se.comb$Study,
                        mod = NULL,
                        par.prior = TRUE,
                        prior.plots = FALSE)

# Merge programme and hypoxia scores
prog_hypoxia.merge <- merge(prog_expr, hypoxia.scores, by = 'Barcode')

# Pan-cancer continuous TMD/hypoxia score plots
p_apobec <- ggplot(prog_hypoxia.merge, aes(x = buffa, y = mean_APOBEC)) + 
  geom_point(alpha = .05) + geom_smooth(method = 'lm') + theme_bw() +
  ggtitle(paste0('APOBEC programme: R^2 = ', round(cor(x = prog_hypoxia.merge$buffa, y = prog_hypoxia.merge$mean_APOBEC),2)))
ggsave(filename = 'Figures/HypoxiaCor_Apobec.pdf', plot = p_apobec)

p_id <- ggplot(prog_hypoxia.merge, aes(x = buffa, y = ImmunogenicDormancy)) + 
  geom_point(alpha = .05) + geom_smooth(method = 'lm') + theme_bw() +
  ggtitle(paste0('Immunological Dormancy programme: R^2 = ', round(cor(x = prog_hypoxia.merge$buffa, y = prog_hypoxia.merge$ImmunogenicDormancy),2)))
ggsave(filename = 'Figures/HypoxiaCor_Immuno.pdf', plot = p_id)

p_ang <- ggplot(prog_hypoxia.merge, aes(x = buffa, y = AngiogenicDormancy)) + 
  geom_point(alpha = .05) + geom_smooth(method = 'lm') + theme_bw() +
    ggtitle(paste0('Angiogenic Dormancy programme: R^2 = ', round(cor(x = prog_hypoxia.merge$buffa, y = prog_hypoxia.merge$AngiogenicDormancy),2)))
ggsave(filename = 'Figures/HypoxiaCor_Angio.pdf', plot = p_ang)

p_tmd <- ggplot(prog_hypoxia.merge, aes(x = buffa, y = TumourMassDormancy)) + 
  geom_point(alpha = .05) + geom_smooth(method = 'lm') + theme_bw() +
  ggtitle(paste0('Tumour Mass Dormancy programme: R^2 = ', round(cor(x = prog_hypoxia.merge$buffa, y = prog_hypoxia.merge$TumourMassDormancy),2)))
ggsave(filename = 'Figures/HypoxiaCor_TMD.pdf', plot = p_tmd)

p_exh <- ggplot(prog_hypoxia.merge, aes(x = buffa, y = Exhaustion)) + 
  geom_point(alpha = .05) + geom_smooth(method = 'lm') + theme_bw() +
  ggtitle(paste0('Exhaustion programme: R^2 = ', round(cor(x = prog_hypoxia.merge$buffa, y = prog_hypoxia.merge$Exhaustion),2)))
ggsave(filename = 'Figures/HypoxiaCor_Exhaustion.pdf', plot = p_exh)

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
ggsave(filename = 'Figures/Corr_by_Study_TMD.pdf', p_tmd)

p_immuno <- ggplot(cor.stats, aes(x = ImmunogenicDormancy, y = -log10(ImmunogenicDormancy.pval))) + geom_point() + 
  geom_text(aes(label = ifelse(cor.stats$ImmunogenicDormancy.pval < .05, Study, NA)), nudge_y = -0.2) + 
  labs(x = 'Immunogenic Dormancy Correlation', y = '-log10(Adjusted Correlation Significance)') +
  xlim(-.6, .6) + ylim(0, 8) + theme_bw() +
  geom_hline(yintercept = -log10(0.05), color = 'red', linetype = 'dashed') + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed')
print(p_immuno)
ggsave(filename = 'Figures/Corr_by_Study_Immuno.pdf', p_immuno)

p_angio <- ggplot(cor.stats, aes(x = AngiogenicDormancy, y = -log10(AngiogenicDormancy.pval))) + geom_point() + 
  geom_text(aes(label = ifelse(cor.stats$AngiogenicDormancy.pval < .05, Study, NA)), nudge_y = -0.2) + 
  labs(x = 'Angiogenic Dormancy Correlation', y = '-log10(Adjusted Correlation Significance)') +
  xlim(-.6, .6) + ylim(0, 8) + theme_bw() +
  geom_hline(yintercept = -log10(0.05), color = 'red', linetype = 'dashed') + 
  geom_vline(xintercept = 0, color = 'red', linetype = 'dashed')
print(p_angio)
ggsave(filename = 'Figures/Corr_by_Study_Angio.pdf', p_angio)


