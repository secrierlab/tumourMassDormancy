# Tumour mass dormancy quantification

This repository contains the workflow for tumour mass dormancy quantification across TCGA cancers.

![alt text](PHATE_analysis/Figures/Figure1_top.jpg)

# Table of contents

## ApobecEnrichmentAnalysis

This folder contains the code which establishes how the tumour mass dormancy and exhaustion programmes are predictive of enrichment of APOBEC mutational signatures.

- ApobecEnrich_define.R: Runs iterations of tSNE and expectation-maximization clustering to create Apobec mutational signature enrichment labels, and compare them to established Apobec enrichment scores.
- RandomForest_100sims.R: Creates 100 random forests from balanced APOBEC-enriched/non-enriched samples, evaluates the overall predictive capacity of the models, and details gene importance.
- RandomForest_exemplary.R: Creates a single random forest using the method as detailed above, and runs further visualisation, including association with normalised PLG expression (Fig 4F) and the corresponding multidimensional scaling plot (Fig 4G). Finally, runs Wilcoxon tests to determine differences in median expressions between Apobec-enriched/non-enriched samples.
- Zscore_normalization.R: Runs standard Z-score normalisation on exemplary FPKM-normalised data.

## Hypoxia

This folder contains the code for calculating hypoxia signature scores and associating these with tumour mass dormancy categories.

- Hypoxia_ScoreCalc.R: Calculates hypoxia scores using the signatures defined by Buffa, West, and Winter, following the procedure detailed by Bhandari et al. Compares the scores calculated here with those determined by Bhandari et al. from a slightly smaller TCGA cohort.
- Hypoxia_TMDCategories.R: Displays the varying distributions of Buffa signature hypoxia scores associated with the detailed tumour mass dormancy categories, as well as the shifting proportions of tumour mass dormancy categorisation associated with hypoxia scores at both the pan-cancer and cancer-specific levels.
- Hypoxia_DormancyCorrelation.R: Displays correlations between the hypoxia and various tumour mass dormancy-associated scores, as well as creating volcano plots visualising correlations on a cancer-specific basis.

## Survival

This folder contains the code for conducting survival analysis associated with tumour mass dormancy categorisation.

- Survival_TMD.R: Visualises survival curves using the pan-cancer cohort, and as a comparison of early/late stage cancers. Runs Cox proportional hazard models to account for features including age at diagnosis, gender, cancer type, and tumour stage in both binary (early/late) and full (Stage I-IV) modes. Conducts a Fisher's exact test to determine significant differences in tumour mass dormancy between early/late stage cancers.

# Copyright

This code is free and is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY. See the GNU General Public License for more details.
