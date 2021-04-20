################################
# Single Random forest classifier of APOBEC enrichment by TMD and exhaustion programme gene expression for demonstrative and analytical purposes
################################

library(randomForest)
library(pROC)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)
library(ggpubr)

setwd('~/Documents/GitHub/tumourMassDormancy')

# Load TMD gene list (required for direction labelling, so exhaustion list not required)
dormancy.gene.list <- read.table('Data/GeneLists/dormancyMarkers_immunologic-angiogenic_updated.txt', header = TRUE)

# Load and attach APOBEC enrichment scores
load('Data/ProgrammeScores/apobec.group.replic.Rdata')
apobec.group.replic$SampleCode <- sapply(apobec.group.replic$samples,
                                         function(x) substr(x,1,16))

# Load in Z-normalised data
se.list <- list()
col.req <- c('barcode','patient','sample','Study','definition')
se.files <- list.files(path = '~/Documents/PhD/Projects/Apobec_Dormancy/Data/GeneExpression/Z_Score',
                       pattern = 'Rdata',
                       full.names = TRUE)

for (file in se.files) {
  study <- strsplit(file, '[.]')[[1]][2]
  print(paste0('Loading expression data from study: ',study,'...'))
  load(file)
  se$Study <- study
  colData(se) <- colData(se)[,col.req]
  
  # Add APOBEC information
  se$Apobec.Total <- NA
  shared.samples <- intersect(se$sample, apobec.group.replic$SampleCode)
  
  se.sample.index <- match(shared.samples, se$sample)
  apobec.sample.index <- match(shared.samples, apobec.group.replic$SampleCode)
  se$Apobec.Total[se.sample.index] <- apobec.group.replic$Apobec.Total[apobec.sample.index]
  se <- se[, !is.na(se$Apobec.Total)]
  se$Apobec.Enriched <- se$Apobec.Total > 50
  
  se.list[[study]] <- se
  
}

# Combine data
se.comb <- do.call(cbind, unlist(se.list))
se.comb.assay <- assay(se.comb)
rownames(se.comb.assay) <- rowData(se)$external_gene_name # Replace Ensembl ID

# Create RF input
set.seed(1234)

rf.input.full <- data.frame(t(se.comb.assay))
rf.input.full$Apobec.Enriched <- se.comb$Apobec.Enriched

# Extract APOBEC enrichment indices
index.apobecTRUE <- which(rf.input.full$Apobec.Enriched)
index.apobecFALSE <- which(!rf.input.full$Apobec.Enriched)

# Extract samples
#   All samples labelled as TRUE for APOBEC enrichment
#   A randomly-selected equal number of samples labelled as FALSE
#   This will enable a balanced sample for random forest classification
samples.include <- c(index.apobecTRUE,
                     sample(index.apobecFALSE, length(index.apobecTRUE), replace = FALSE))

rf.input <- rf.input.full[samples.include, ]
rf.input$Apobec.Enriched <- as.factor(rf.input.full$Apobec.Enriched[samples.include])

# Apply model
model <- randomForest(Apobec.Enriched ~ .,
                      data = rf.input,
                      do.trace = 100,
                      ntree = 1000,
                      proximity = TRUE,
                      importance = TRUE)
print(model)

# MDS plot to separate high- and low-APOBEC-enriched samples
distance.matrix <- dist(1-model$proximity)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE)
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100,1)
mds.values <- mds.stuff$points
mds.data <- data.frame(Sample = rownames(mds.values),
                       X = mds.values[,1],
                       Y = mds.values[,2])
mds.data <- cbind(mds.data, rf.input)

mds.data <- mds.data[sample(1:nrow(mds.data), nrow(mds.data), replace = FALSE),]
p_mds <- ggplot(data = mds.data, aes(x = X, y = Y, label = Sample)) + 
  geom_point(aes(color = Apobec.Enriched), size = .5) +
  theme_bw() + scale_color_manual(values = c('lightskyblue', 'navy')) +
  xlab(paste('MDS1 - ', mds.var.per[1],'%', sep='')) +
  ylab(paste('MDS2 - ', mds.var.per[2],'%', sep=''))
print(p_mds)
ggsave(filename = 'ApobecEnrichmentAnalysis/Figures/MDS_ApobecEnrich.pdf', plot = p_mds)

# Waterfall plot displaying expression of PLG
#   Importance analysis dictates that PLG is the main driver behind the classifier
#   Almost all samples in the breakout cluster, defined by a specific PLG expression value,
#     are labelled as non-APOBEC-enriched
mds.data.reorder <- mds.data
mds.data.reorder <- mds.data[order(mds.data$PLG, decreasing = TRUE),]
mds.data.reorder$row.index <- 1:nrow(mds.data.reorder)
p_plg <- ggplot(data = mds.data.reorder, aes(Sample, fill = Apobec.Enriched)) + 
  geom_rect(aes(x = Sample, xmin = row.index - .45, xmax = row.index + .45,
                ymin = 0, ymax = PLG)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position = c(0.1, 0.05)) +
  scale_fill_manual(values = c('lightskyblue','navy')) +
  ylab('Normalised PLG Expression')
print(p_plg)
ggsave(filename = 'ApobecEnrichmentAnalysis/Figures/PLG_Waterfall.pdf', plot = p_plg)

# Gene-specific differential expression analysis dependent on APOBEC-enrichment
# For each gene included in se.comb.assay (all samples):
#   1. Separate samples into Apobec-enriched and non-Apobec-enriched
#   2. Calculate difference in median Z-score between the two groups
#   3. Conduct a Wilcox test to determine difference in medians

wilcox.pvals = median.diffs <- vector(length = ncol(rf.input.full) - 1)
for (i in 1:length(wilcox.pvals)) {
  
  # Calculate differences in medians
  med.apobec <- median(rf.input.full[rf.input.full$Apobec.Enriched, i])
  med.noapobec <- median(rf.input.full[!rf.input.full$Apobec.Enriched, i])
  median.diffs[i] <- med.apobec - med.noapobec
  
  # Wilcox test
  w <- wilcox.test(rf.input.full[, i] ~ rf.input.full$Apobec.Enriched)
  wilcox.pvals[i] <- w$p.value
  
}

# Account for multiple testing using a Benjamini-Hochberg correction
wilcox.pvals <- p.adjust(wilcox.pvals, method = 'BH')

# Collate results and present as a volcano plot
median.diffs.df <- data.frame(Gene = colnames(rf.input.full)[1:ncol(rf.input.full)-1],
                              Diff.Median = median.diffs,
                              Wilcox.pval = wilcox.pvals)
median.diffs.df$Direction <- 'Upregulated'
median.diffs.df$Direction[which(median.diffs.df$Gene %in% dormancy.gene.list$Gene[dormancy.gene.list$Direction == 'down'])] <- 'Downregulated'
median.diffs.df$Label <- ifelse(median.diffs.df$Wilcox.pval < 1e-2, median.diffs.df$Gene, NA)

p_wilcox <- ggplot(median.diffs.df, aes(x = Diff.Median, y = -log10(Wilcox.pval), color = Direction, label = Gene)) +
  geom_point() + theme_bw() + scale_color_manual(values = c('orange', 'purple')) +
  geom_text(label = median.diffs.df$Label, color = 'black', nudge_y = -.3, size = 3) +
  geom_vline(xintercept = 0, linetype = 'dashed') + geom_hline(yintercept = -log10(0.01), linetype = 'dashed') +
  labs(x = 'Difference in Medians', y = '-log10 Wilcox test (adjusted)')
print(p_wilcox)
ggsave(filename = 'ApobecEnrichmentAnalysis/Figures/ApobecGroups_Volcano.pdf', plot = p_wilcox)





