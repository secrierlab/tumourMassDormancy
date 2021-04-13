setwd('~/Documents/GitHub/tumourMassDormancy/')

library(ggplot2)

# Load in hypoxia signatures and programme scores
load('Data/GeneLists/HypoxiaSignatures.Rdata')
load('Data/ProgrammeScores/programme_scores_and_TMD_assignments.RData')

# Load in dormancy and exhaustion markers (will be removed from hypoxia scores)
dormancy_markers <- read.table('Data/GeneLists/dormancyMarkers_immunologic-angiogenic_updated.txt', header = TRUE)
exhaustion_markers <- read.table('Data/GeneLists/exhaustionMarkers.txt')

genes.dormancy <- dormancy_markers$Gene
genes.exhaustion <- exhaustion_markers$V1

# Remove genes from hypoxia signatures which appear in dormancy signature
hypoxia.sigs$buffa <- hypoxia.sigs$buffa[!(hypoxia.sigs$buffa %in% genes.dormancy)]
hypoxia.sigs$west <- hypoxia.sigs$west[!(hypoxia.sigs$west %in% genes.dormancy)]
hypoxia.sigs$winter <- hypoxia.sigs$winter[!(hypoxia.sigs$winter) %in% genes.dormancy]

# Load expression data (all genes included in the three hypoxia signatures) - se.df
load('Data/RNA_seq/Hypoxia_Score_FPKM_data.RData')

# Verification that the hypoxia and programme score data frames are ordered identically
table(rownames(se.df) == prog_expr$Barcode)

# Define a function to convert expression scores into binary hypoxia profiles
hypoxia.convert <- function(expr.vec) {
  median.expr <- median(expr.vec)
  hypoxia.vec <- sapply(expr.vec, function(x) ifelse(x > median.expr, 1, -1))
  return(hypoxia.vec)
}

# Create data frame for collating hypoxia scores
# NB, Samples are in the same order for se.df and prog_expr
hypoxia.scores <- data.frame(Barcode = rownames(se.df),
                             Study = prog_expr$cancer_type)

# For each hypoxia signature:
#   - extract the relevant genes, ensuring they are present in the dataset
#   - convert into the binary hypoxia format
#   - sum across each sample
#   - add the results to the dataframe 'hypoxia.scores'

for (sig in c('buffa','west','winter')) {
  hypoxia.sig.genes <- hypoxia.sigs[[sig]]
  hypoxia.sig.genes <- hypoxia.sig.genes[hypoxia.sig.genes %in% colnames(se.df)]
  se.sig <- se.df[,hypoxia.sig.genes]
  se.hypoxia <- apply(se.sig, 2, hypoxia.convert)
  se.hypoxia.scores <- apply(se.hypoxia, 1, sum)
  hypoxia.scores[[sig]] <- se.hypoxia.scores
}

# Save hypoxia scores
save(hypoxia.scores, file = 'Data/ProgrammeScores/Hypoxia_SigScores.Rdata')

# Compare with Bhandari scores (pre-calculated, included genes within TMD programmes)
bhandari <- read.table('Data/ProgrammeScores/BhandariHypoxiaScores.txt', h=T)

hypoxia.scores$patient <- sapply(hypoxia.scores$Barcode,
                                 function(x) paste(strsplit(x,'-')[[1]][1:3],collapse = '.'))

# Match samples within our hypoxia scores and those calculated by Bhandari et al.
hypoxia.merge <- merge(hypoxia.scores, bhandari, by.x = 'patient', by.y = 'patient_id')

# Buffa
p_buffa <- ggplot(hypoxia.merge, aes(x = buffa, y = Buffa_hypoxia_score_pan_cancer)) + 
  geom_point(alpha = .05) + theme_bw() +
  labs(x = 'Buffa (calculated)', y = 'Buffa (Bhandari et al)')
print(p_buffa)
ggsave(filename = 'Hypoxia/Figures/Buffa_BhandariComparison.pdf', plot = p_buffa)


# West
p_west <- ggplot(hypoxia.merge, aes(x = west, y = West_hypoxia_score_pan_cancer)) +
  geom_point(alpha = .05) + theme_bw() +
  labs(x = 'West (calculated)', y = 'West (Bhandari et al)')
print(p_west)
ggsave(filename = 'Hypoxia/Figures/West_BhandariComparison.pdf', plot = p_west)

# Winter
p_winter <- ggplot(hypoxia.merge, aes(x = winter, y = Winter_hypoxia_score_pan_cancer)) +
  geom_point(alpha = .05) + theme_bw() +
  labs(x = 'Winter (calculated)', y = 'Winter (Bhandari et al)')
print(p_winter)
ggsave(filename = 'Hypoxia/Figures/Winter_BhandariComparison.pdf', plot = p_winter)




