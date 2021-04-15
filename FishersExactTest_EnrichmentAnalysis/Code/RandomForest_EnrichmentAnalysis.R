################################
# Random Forest analysis of mutation enrichment in TMD vs NO TMD samples
################################

setwd('~/Documents/GitHub/tumourMassDormancy')

library(dplyr)
library(randomForest)
library(ggplot2)
library(pROC)

# Load strict and lenient mutation data
load('Data/MutationData/maf_all.strict.RData')
load('Data/MutationData/maf_all.lenient.RData')

# Convert to data frames
maf_all.strict <- as.data.frame(maf_all.strict)
maf_all.lenient <- as.data.frame(maf_all.lenient)

rownames(maf_all.strict) <- maf_all.strict$Patient_short
rownames(maf_all.lenient) <- maf_all.lenient$Patient_short

# Extract genes with mutation rate > threshold
threshold <- 0.01

maf_strict.mut <- maf_all.strict[, 4:ncol(maf_all.strict)] # Extract only mutation data
maf_strict.mut <- maf_strict.mut[, apply(maf_strict.mut, 2, mean) > threshold]
maf_strict.thres <- cbind(maf_all.strict[,1:3], maf_strict.mut)

maf_lenient.mut <- maf_all.lenient[, 4:ncol(maf_all.lenient)]
maf_lenient.mut <- maf_lenient.mut[, apply(maf_lenient.mut, 2, mean) > threshold]
maf_lenient.thres <- cbind(maf_all.lenient[,1:3], maf_lenient.mut)

# Load in programme scores
load('Data/ProgrammeScores/programme_scores_and_TMD_assignments.RData')

# Attach expression programmes, and delete non-matchers
maf_strict.complete <- merge(x = maf_strict.thres, y = prog_expr,
                             by.x = 'Patient_short', by.y = 'SampleID')
maf_lenient.complete <- merge(x = maf_lenient.thres, y = prog_expr,
                              by.x = 'Patient_short', by.y = 'SampleID')

# RUN FOR STRICT MUTATION ASSIGNMENT

# Extract columns of interest (genes and TMD_three_categories)
genes.strict <- colnames(maf_strict.mut)
rf_strict.full <- maf_strict.complete[, c(genes.strict, 'TMD_three_categories')]

# Extract samples indices of interest (TMD_three_categories == 'YES' or 'NO')
sample_index.TMD <- which(rf_strict.full$TMD_three_categories == 'YES')
sample_index.noTMD <- which(rf_strict.full$TMD_three_categories == 'NO')

# Set seed and randomly selected noTMD samples to generate balanced datasets
set.seed(1234)
samples_index.include <- c(sample_index.noTMD,
                           sample(sample_index.TMD, length(sample_index.noTMD), replace = FALSE))
rf_strict.data <- rf_strict.full[samples_index.include, ]
rf_strict.data$TMD_three_categories <- as.factor(rf_strict.data$TMD_three_categories)

# Run random forest analysis
Sys.time(); model <- randomForest(TMD_three_categories ~ .,
                                  data = rf_strict.data,
                                  do.trace = 10,
                                  ntree = 500,
                                  proximity = TRUE,
                                  importance = TRUE); Sys.time()
print(model)

# Save importance analysis and identify mutation-enriched genes
imp.df <- data.frame(model$importance[order(model$importance[,'MeanDecreaseGini'], decreasing = TRUE),])
imp.df$Gene <- rownames(imp.df)
p_strict <- ggplot(head(imp.df, 10), aes(x = MeanDecreaseGini, y = reorder(Gene, MeanDecreaseGini))) + 
  geom_bar(stat = 'identity') +
  labs(x = 'Mean Decrease in Gini Index', y = 'Gene') +
  ggtitle('TMD Mutation Enrichment: Strict Assignment')
print(p_strict)

save(imp.df, file = 'FishersExactTest_EnrichmentAnalysis/RandomForest_Importance_Strict.RData')
ggsave(filename = 'FishersExactTest_EnrichmentAnalysis/Figures/RandomForest_ImportanceTop10_Strict.pdf', plot = p_strict)


# RUN FOR LENIENT MUTATION ASSIGNMENT

# Extract columns of interest (genes and TMD_three_categories)
genes.lenient <- colnames(maf_lenient.mut)
rf_lenient.full <- maf_lenient.complete[, c(genes.lenient, 'TMD_three_categories')]

# Extract samples indices of interest (TMD_three_categories == 'YES' or 'NO')
sample_index.TMD <- which(rf_lenient.full$TMD_three_categories == 'YES')
sample_index.noTMD <- which(rf_lenient.full$TMD_three_categories == 'NO')

# Set seed and randomly selected noTMD samples to generate balanced datasets
set.seed(1234)
samples_index.include <- c(sample_index.noTMD,
                           sample(sample_index.TMD, length(sample_index.noTMD), replace = FALSE))
rf_lenient.data <- rf_lenient.full[samples_index.include, ]
rf_lenient.data$TMD_three_categories <- as.factor(rf_lenient.data$TMD_three_categories)

# Run random forest analysis
Sys.time(); model <- randomForest(TMD_three_categories ~ .,
                                  data = rf_lenient.data,
                                  do.trace = 10,
                                  ntree = 500,
                                  proximity = TRUE,
                                  importance = TRUE); Sys.time()
print(model)

# Save importance analysis and identify mutation-enriched genes
imp.df <- data.frame(model$importance[order(model$importance[,'MeanDecreaseGini'], decreasing = TRUE),])
imp.df$Gene <- rownames(imp.df)
p_lenient <- ggplot(head(imp.df, 10), aes(x = MeanDecreaseGini, y = reorder(Gene, MeanDecreaseGini))) + 
  geom_bar(stat = 'identity') +
  labs(x = 'Mean Decrease in Gini Index', y = 'Gene') +
  ggtitle('TMD Mutation Enrichment: Lenient Assignment')
print(p_lenient)

save(imp.df, file = 'FishersExactTest_EnrichmentAnalysis/RandomForest_Importance_Lenient.RData')
ggsave(filename = 'FishersExactTest_EnrichmentAnalysis/Figures/RandomForest_ImportanceTop10_Lenient.pdf', plot = p_strict)






