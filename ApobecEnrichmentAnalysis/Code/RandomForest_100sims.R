################################
# 100 Random forest classifiers of APOBEC enrichment by TMD and exhaustion programme gene expression
################################

library(randomForest)
library(pROC)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)
library(ggpubr)

setwd('~/Documents/GitHub/tumourMassDormancy')

# Load and attach APOBEC enrichment scores
load('Data/ProgrammeScores/apobec.group.replic.Rdata')
apobec.group.replic$SampleCode <- sapply(apobec.group.replic$samples,
                                         function(x) substr(x,1,16))

# Load in Z-normalised data
se.list <- list()
col.req <- c('barcode','patient','sample','Study','definition')
se.files <- list.files(path = 'Data/RNA_seq/Z_Scores',
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
rf.input.full <- data.frame(t(se.comb.assay))
rf.input.full$Apobec.Enriched <- se.comb$Apobec.Enriched

# Extract APOBEC enrichment indices
index.apobecTRUE <- which(rf.input.full$Apobec.Enriched)
index.apobecFALSE <- which(!rf.input.full$Apobec.Enriched)

# Create objects for saving results
n.runs <- 100 # Each run takes around 40 seconds
auc.data.list <- list()
auc.vals <- vector(length = n.runs)

importance.gini.decrease <- data.frame(
  Gene = colnames(rf.input.full)[1:ncol(rf.input.full)-1])

# Run model n.runs times and save AUC and importance results
set.seed(1234)
for (i in 1:n.runs) {
  print(paste0('Applying random forest model: run ', i, ' of ', n.runs, '...'))
  
  # Select samples
  samples.include <- c(index.apobecTRUE,
                       sample(index.apobecFALSE, length(index.apobecTRUE), replace = FALSE))
  rf.input <- rf.input.full[samples.include, ]
  rf.input$Apobec.Enriched <- as.factor(rf.input.full$Apobec.Enriched[samples.include])
  
  # Apply model
  model <- randomForest(Apobec.Enriched ~ .,
                        data = rf.input,
                        do.trace = 100,
                        ntree = 1000,
                        importance = TRUE)
  
  # Generate, plot, and save ROC data
  roc <- roc(rf.input$Apobec.Enriched, model$votes[,2])
  roc.df <- data.frame(Sensitivity = roc$sensitivities,
                       Specificity = 1 - roc$specificities)
  auc.data.list[[i]] <- roc.df
  
  auc.vals[i] <- auc(roc)
  print(auc(roc))
  
  # Save importance values
  importance.gini.decrease[[paste0('run',i)]] <- model$importance[, 'MeanDecreaseGini']
  
}

# Save output
ggsave('ApobecEnrichmentAnalysis/Figures/AUC_100sims.pdf', device = 'pdf')
save(auc.data.list, file = 'ApobecEnrichmentAnalysis/Results/AUC_100sims_data.Rdata')
save(auc.vals, file = 'ApobecEnrichmentAnalysis/Results/AUC_100sims.Rdata')

# Collate AUC results into a single plot
#   Round specificity values to 2 decimal points
#   Extract the mean, minimum, and maximum Sensitivity values for each Specificity values
auc.all <- do.call(rbind, auc.data.list)
auc.all$Spec.Bin <- round(auc.all$Specificity, 2)

auc.summary2 <- auc.all %>% 
  group_by(Spec.Bin) %>%
  summarise(
    N = n(),
    Sensitivity.Min = min(Sensitivity),
    Sensitivity.Max = max(Sensitivity),
    Sensitivity = mean(Sensitivity)
  )

p_auc <- ggplot(auc.summary2, aes(x = Spec.Bin)) +
  geom_ribbon(data = auc.summary2, aes(ymin = Sensitivity.Min, ymax = Sensitivity.Max), fill = 'plum') + 
  geom_line(aes(y = Sensitivity)) +
  geom_line(aes(y = Sensitivity.Max), color = 'purple4') +
  geom_line(aes(y = Sensitivity.Min), color = 'purple4') +
  theme_bw() +
  labs(x = 'Specificity', y = 'Sensitivity')
print(p_auc)
ggsave(filename = 'ApobecEnrichmentAnalysis/Figures/AUC_100sims.pdf', plot = p_auc)

# Generate histogram of AUC values
pdf(file = 'ApobecEnrichmentAnalysis/Figures/AUC_values_hist.pdf')
hist(auc.vals, breaks = seq(0.81, 0.86, by = 0.005),
     xlab = 'AUC', main = paste0('Mean AUC = ', round(mean(auc.vals), 4)),
     col = 'plum')
abline(v = mean(auc.vals), col = 'purple4')
dev.off()

# Collate importance analysis
rownames(importance.gini.decrease) <- importance.gini.decrease$Gene
importance.gini.decrease <- importance.gini.decrease[, -1]
save(importance.gini.decrease, file = 'ApobecEnrichmentAnalysis/Results/MeanDecreaseGini_100sims.Rdata')

# Print importance analysis
importance.summary <- data.frame(Gene = rownames(importance.gini.decrease),
                                 MeanDecreaseGini = apply(importance.gini.decrease, 1, mean),
                                 sd = apply(importance.gini.decrease, 1, sd))
importance.summary <- importance.summary[order(importance.summary$MeanDecreaseGini, decreasing = TRUE),]
importance.summary <- head(importance.summary, 10)

p_importance <- ggplot(importance.summary, aes(x = MeanDecreaseGini, y = reorder(Gene, MeanDecreaseGini))) +
  geom_bar(stat = 'identity') +
  geom_errorbar(aes(xmin = MeanDecreaseGini - sd, xmax = MeanDecreaseGini + sd), width = .5) +
  ylab('Gene') + xlab('Mean Decrease in Gini Index') + theme_bw()
print(p_importance)
ggsave('ApobecEnrichmentAnalysis/Figures/Importance_Top10.pdf', plot = p_importance, width = 5, height = 5)




