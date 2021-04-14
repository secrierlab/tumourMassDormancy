################################
# Application of tSNE and EM clustering for APOBEC-enrichment definition from mutational signature profiles
################################

library('Rtsne')
library('ggplot2')
library('factoextra')
library('NbClust')
library('stats')
library('EMCluster')
library('dplyr')
library('ggpubr')

# Obtain mutational signature profiles
setwd('~/Documents/GitHub/tumourMassDormancy')
files = list.files('Data/MutationalSignatures', full.names = TRUE)

sigs.complete <- NULL
for (file in files) {
  study = strsplit(file, split = '[.]')[[1]][3]
  load(file)
  sigs.defaultnounknown$Study = study
  sigs.complete <- rbind(sigs.complete, sigs.defaultnounknown)
  print(paste0('Download complete: ', study, ' , Total samples: ', nrow(sigs.complete)))
}

# Remove zero columns
sigs.input <- sigs.complete[,apply(sigs.complete[,1:ncol(sigs.complete)-1], 2, sum) > 0]

# tSNE Demonstration
set.seed(1234)
print('Running tSNE dimensionality reduction...')

sigs.tsne.full = Rtsne(sigs.input[,1:ncol(sigs.input)-1], perplexity = 100, max_iter = 5000, pca = TRUE, check_duplicates = FALSE)
sigs.result.tsne <- cbind(sigs.input, as.data.frame(sigs.tsne.full$Y))
colnames(sigs.result.tsne)[33:34] <- c('TSNE1', 'TSNE2')

p_tsne.study <- ggplot(sigs.result.tsne, aes(x=TSNE1, y=TSNE2, colour=Study)) +
  geom_point(size = .5) + theme_bw()
print(p_tsne.study)
ggsave(filename = 'ApobecEnrichmentAnalysis/Figures/tSNE_Study.pdf', plot = p_tsne.study)


# Expectation-Maximization clustering to determine optimal number of clusters
emclust.optnum <- data.frame(n.clust = 1:20, likelihood = NA)
for (n.clust in emclust.optnum$n.clust) {
  print(paste0('Testing... n.clust = ',n.clust))
  sigs.tsne.g <- emgroup(sigs.result.tsne[,c('TSNE1','TSNE2')], nclass = n.clust)
  emclust.optnum$likelihood[n.clust] <- sigs.tsne.g$llhdval
}

pdf(file = 'ApobecEnrichmentAnalysis/Figures/EMClustering_Likelihoods.pdf', height = 5)
plot(x = emclust.optnum$n.clust, y = emclust.optnum$likelihood,
     pch = 16, type = 'b', xlab = 'Num. clusters', ylab = 'Likelihood', main = 'Optimal Num. Clusters (likelihood)')
abline(v = 10, col = 'red'); abline(v = 17, col = 'red')
dev.off()

# Demonstrate result of cluster analysis with n.clust = 10 and 17, as well as colouring by Apobec
sigs.tsne.g1 <- emgroup(sigs.result.tsne[,c('TSNE1','TSNE2')], nclass = 10)
sigs.tsne.g2 <- emgroup(sigs.result.tsne[,c('TSNE1','TSNE2')], nclass = 17)
sigs.tsne.ex <- sigs.result.tsne
sigs.tsne.ex$Cluster10 <- as.factor(sigs.tsne.g1$class)
sigs.tsne.ex$Cluster17 <- as.factor(sigs.tsne.g2$class)
sigs.tsne.ex$Apobec <- sigs.tsne.ex$SBS2 + sigs.tsne.ex$SBS13

p_tsne.clust10 <- ggplot(sigs.tsne.ex, aes(x = TSNE1, y = TSNE2, colour = Cluster10)) + 
  geom_point(size = .5) + theme_bw()
print(p_tsne.clust10)
ggsave(filename = 'ApobecEnrichmentAnalysis/Figures/tSNE_10clust.pdf', plot =  p_tsne.clust10)

p_tsne.clust17 <- ggplot(sigs.tsne.ex, aes(x = TSNE1, y = TSNE2, colour = Cluster17)) + 
  geom_point(size = .5) + theme_bw()
print(p_tsne.clust17)
ggsave(filename = 'ApobecEnrichmentAnalysis/Figures/tSNE_17clust.pdf', plot =  p_tsne.clust17)

p_tsne.apobec <- ggplot(sigs.tsne.ex, aes(x = TSNE1, y = TSNE2, colour = Apobec)) + 
  geom_point(size = .5) + theme_bw() +
  scale_colour_gradient(low = 'lightskyblue', high = 'navy')
print(p_tsne.apobec)
ggsave(filename = 'ApobecEnrichmentAnalysis/Figures/tSNE_Apobec.pdf', plot =  p_tsne.apobec)


# Procedure for defining and testing replicability of APOBEC-enrichment
set.seed(1234)
apobec.group.replic <- data.frame(samples = rownames(sigs.input)) # store Apobec cluster inclusion here
num.runs = 100 # Each run takes approximately 3 minutes

# For each run:
#       1. Run tSNE with perplexity = 100 and max_iter = 5000
#       2. Apply EM clustering with 10 classes
#       3. Define APOBEC cluster as that with the highest median SBS2+SBS13 enrichment
#       4. Add APOBEC cluster label to apobec.group.replic

for (run in 1:num.runs) {
  print(paste0('Run ',run,' of ',num.runs,'... ',Sys.time()))
  
  # Run tSNE
  sigs.tsne.full = Rtsne(sigs.input[,1:ncol(sigs.input)-1], perplexity = 100, max_iter = 5000, pca = TRUE, check_duplicates = FALSE)
  sigs.tsne <- cbind(sigs.input, as.data.frame(sigs.tsne.full$Y))
  colnames(sigs.tsne)[33:34] <- c('TSNE1', 'TSNE2')
  
  # EM clustering, with 10 clusters
  sigs.tsne$Apobec <- sigs.tsne$SBS2 + sigs.tsne$SBS13
  sigs.tsne.g <- emgroup(sigs.tsne[,c('TSNE1','TSNE2')], nclass = 10)
  sigs.tsne$Cluster10 <- as.factor(sigs.tsne.g$class)
  
  # Find APOBEC cluster (cluster with highest median Apobec enrichment)
  apobec.med <- sigs.tsne %>%
    group_by(Cluster10) %>%
    summarise(Med.Apobec = median(Apobec, na.rm = TRUE))
  apobec.clust <- apobec.med$Cluster10[apobec.med$Med.Apobec == max(apobec.med$Med.Apobec)]
  
  # Find samples in APOBEC cluster and add to replicability test
  sigs.tsne$Apobec.Cluster <- sigs.tsne$Cluster10 == apobec.clust
  apobec.group.replic[[paste0('run',run)]] <- sigs.tsne$Apobec.Cluster
}
apobec.group.replic$Apobec.Total <- apply(apobec.group.replic[,-1], 1, sum)
apobec.clust.tot <- data.frame(table(apobec.group.replic$Apobec.Total))
colnames(apobec.clust.tot)[1] <- 'InApobecCluster'

# Write APOBEC cluster assignments
save(apobec.group.replic,
     file = 'Data/ProgrammeScores/apobec.group.replic.Rdata')

# If required:
load('Data/ProgrammeScores/apobec.group.replic.Rdata')

# Load in APOBEC enrichment scores calculated by Sigminer and add them to apobec.group.replic
apobec_enrich.files <- list.files(path = 'Data/Sigminer_APOBECEnrich',
                                  pattern = 'Rdata',full.names = TRUE)
apobec.group.replic$APOBEC_Enriched = apobec.group.replic$APOBEC_Enrichment <- NA

# For each file: load, add $APOBEC_Enrichment and $APOBEC_Enriched to apobec.group.replic
for (file in apobec_enrich.files) {
  load(file) # filename: apobec
  # Remove samples which don't have signature profiles
  apobec <- apobec[apobec$Tumor_Sample_Barcode %in% apobec.group.replic$samples]
  samples.index <- match(apobec$Tumor_Sample_Barcode, apobec.group.replic$samples)
  apobec.group.replic$APOBEC_Enrichment[samples.index] <- apobec$APOBEC_Enrichment
  apobec.group.replic$APOBEC_Enriched[samples.index] <- apobec$APOBEC_Enriched
}

# Comparisons
table(apobec.group.replic$APOBEC_Enriched, apobec.group.replic$Apobec.Total, dnn = c('APOBEC_Enriched','APOBEC Cluster Assignment'))

apobec.group.replic_tab <- apobec.group.replic %>%
  group_by(Apobec.Total, APOBEC_Enriched) %>% summarise(n())
colnames(apobec.group.replic_tab) <- c('InApobecCluster','APOBEC_Enriched','Freq')
p_score <- ggplot(data = apobec.group.replic_tab, aes(x=InApobecCluster, y=Freq, fill=APOBEC_Enriched)) +
  geom_bar(stat='identity')
print(p_score)
ggsave(filename = 'ApobecEnrichmentAnalysis/Figures/ApobecEnrich_Barplot.pdf', plot = p_score)

sigs.tsne.ex$APOBEC_Enriched <- apobec.group.replic$APOBEC_Enriched[match(apobec.group.replic$samples, rownames(sigs.tsne))]
sigs.tsne.ex$APOBEC_Enrichment <- apobec.group.replic$APOBEC_Enrichment[match(apobec.group.replic$samples, rownames(sigs.tsne))]

p_tsne.enrich <- ggplot(sigs.tsne.ex, aes(x=TSNE1, y=TSNE2, colour=APOBEC_Enriched)) +
  geom_point(size = .5) +
  ggtitle(paste0('Apobec Cluster: Enrichment'))
plot(p_tsne.enrich)
ggsave(filename = 'ApobecEnrichmentAnalysis/Figures/ApobecEnrich_tSNE.pdf', plot = p_tsne.enrich)




