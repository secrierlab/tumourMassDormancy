# PROCEDURE FOR CALCULATING Z SCORES FROM FPKM NORMALIZED DATA

library(org.Hs.eg.db)

# Only dormancy and exhaustion genes need to be extracted
setwd('~/Documents/GitHub/tumourMassDormancy/')
dormancy_gene_list <- read.table('Data/GeneLists/dormancyMarkers_immunologic-angiogenic_updated.txt', h = T)
genes.dormancy <- dormancy_gene_list$Gene

exhaustion_gene_list <- read.table('Data/GeneLists/exhaustionMarkers.txt', h = F)
genes.exhaustion <- exhaustion_gene_list$V1

# Load in programme scores
load('Data/ProgrammeScores/programme_scores_and_TMD_assignments.RData')

# Load in example data, log2/scale, and save
load('Data/RNA_seq/example_FPKM_data.RData')

# Define gene names and extract matched genes
gene.names <- mapIds(org.Hs.eg.db,
                     keys = colnames(combined_data),
                     keytype = 'ENSEMBL',
                     column = 'SYMBOL')
tmd.index <- gene.names %in% c(genes.dormancy, genes.exhaustion)

combined_data.sub <- combined_data[, tmd.index]
colnames(combined_data.sub) <- gene.names[tmd.index]

# Calculate Z scores via log2, followed by scale()
combined_data.z <- t(scale(log2(combined_data.sub + 1)))

# Save output
save(combined_data.z, file = 'Data/RNA_seq/example_Znorm_data.RData')








