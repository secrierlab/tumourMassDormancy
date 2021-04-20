# TCGA Mutation Data Download and APOBEC scoring
library('TCGAbiolinks')
library('sigminer')
library('maftools')

setwd('~/Documents/PhD/Data/')

studies = unique(sigs.input$Study)

for (study in studies[4:length(studies)]) {
  mut <- GDCquery_Maf(tumor = study, pipelines = 'mutect2')
  mut_maf <- read.maf(maf = mut, useAll = FALSE)
  mt_tally_sbs <- sig_tally(
    mut_maf,
    ref_genome = 'BSgenome.Hsapiens.UCSC.hg38',
    useSyn = TRUE,
    mode = 'SBS'
  )
  
  apobec <- mt_tally_sbs$APOBEC_scores
  save(apobec,
       file = paste0('Sigminer_APOBECEnrich/TCGA.Sigminer.APOBEC.Enrich.',study,'.Rdata'))
  
}

