# SCRIPT FOR DEFINING HYPOXIA SCORES

# Procedure:
#   Collate FPKM-normalized data for all samples
#   Extract genes involved in established mRNA hypoxia scores
#     (not including genes included in dormancy/exhaustion scores)
#   Extract samples where dormancy scores have been defined
#   Across the pan-cancer sample, for each hypoxia signature calculate score as follows:
#     - For each gene, assign +1 to samples with the top 50% of mRNA abundance values (value > median)
#     -  and assign -1 to samples with the bottom 50% of mRNA abundance values (value < median)
#     - Sum across all genes in a signature to obtain the score
library(SummarizedExperiment)

load('Data/HypoxiaSignatures.Rdata') # Genes constituting defined hypoxia signatures (excluding TMD programme)

# Load program scores to show which samples are included
load('Data/programme_scores_and_TMD_assignments.RData')

dormancy_markers <- read.table('Data/dormancyMarkers_immunologic-angiogenic_updated2.txt', header = TRUE)
exhaustion_markers <- read.table('Data/exhaustionMarkers.txt')

genes.dormancy <- dormancy_markers$Gene
genes.exhaustion <- exhaustion_markers$V1

# Remove genes from hypoxia signatures which appear in dormancy signature
hypoxia.sigs$buffa <- hypoxia.sigs$buffa[!(hypoxia.sigs$buffa %in% genes.dormancy)]
hypoxia.sigs$west <- hypoxia.sigs$west[!(hypoxia.sigs$west %in% genes.dormancy)]
hypoxia.sigs$winter <- hypoxia.sigs$winter[!(hypoxia.sigs$winter) %in% genes.dormancy]

# Load FPKM-normalized expression data
#   and extract all genes included in the 3 hypoxia signatures
se.list <- list()
col.req <- c('barcode','patient','sample','Study','definition')
# Specify path to FPKM_Normalized values
se.files <- list.files(path = '../../Data/GeneExpression/FPKM_Normalized/',
                       pattern = 'Rdata',
                       full.names = TRUE)

for (file in se.files) {
  study <- strsplit(file, '[.]')[[1]][6] # Alter to match study code in file name
  print(paste0('Loading expression data from study: ',study,'...'))
  load(file)
  se$Study <- study
  colData(se) <- colData(se)[,col.req]
  se <- se[,se$definition == 'Primary solid Tumor'] # Extract primary tumours
  se <- se[rowData(se)$external_gene_name %in% unique(unlist(hypoxia.sigs)),] # Extract genes contained across all hypoxia signatures

  se.list[[study]] <- se
}

se.comb <- do.call(cbind, unlist(se.list))
se.comb <- se.comb[, prog_expr$Barcode]

se.comb.assay <- assay(se.comb)
rownames(se.comb.assay) <- rowData(se.comb)$external_gene_name

# Define function to convert expression scores into binary hypoxia profiles
hypoxia.convert <- function(expr.vec) {
  median.expr <- median(expr.vec)
  hypoxia.vec <- sapply(expr.vec, function(x) ifelse(x > median.expr, 1, -1))
  return(hypoxia.vec)
}

# Create data frame for collating hypoxia scores
hypoxia.scores <- data.frame(Barcode = se.comb$barcode,
                             Study = se.comb$Study)

# For each hypoxia signature:
#   - extract the relevant genes, ensuring they are present in the dataset
#   - convert into the binary hypoxia format
#   - sum across each sample
#   - add the results to the dataframe 'hypoxia.scores'

for (sig in c('buffa', 'west', 'winter')) {
  hypoxia.sig.genes <- hypoxia.sigs[[sig]]
  hypoxia.sig.genes <- hypoxia.sig.genes[hypoxia.sig.genes %in% rownames(se.comb.assay)]
  se.sig <- se.comb.assay[hypoxia.sig.genes,]
  se.hypoxia <- apply(se.sig, 1, hypoxia.convert)
  se.hypoxia.scores <- apply(se.hypoxia, 1, sum)
  hypoxia.scores[[sig]] <- se.hypoxia.scores
}

# Save hypoxia scores
save(hypoxia.scores, file = 'Data/Hypoxia_SigScores.Rdata')

# Compare with Bhandari scores (pre-calculated, included genes within TMD programmes)
bhandari <- read.table('Data/BhandariHypoxiaScores.txt', h=T)

hypoxia.scores$patient <- sapply(hypoxia.scores$Barcode,
                                 function(x) paste(strsplit(x,'-')[[1]][1:3],collapse = '.'))

hypoxia.merge <- merge(hypoxia.scores, bhandari, by.x = 'patient', by.y = 'patient_id')
p_bhandari <- ggplot(hypoxia.merge, aes(x = buffa, y = Buffa_hypoxia_score_pan_cancer)) + 
  geom_point(alpha = .05) + theme_bw() +
  labs(x = 'Buffa (calculated)', y = 'Buffa (Bhandari et al)')
print(p_bhandari)
ggsave(filename = 'Figures/Buffa_BhandariComparison.pdf', plot = p_bhandari)

