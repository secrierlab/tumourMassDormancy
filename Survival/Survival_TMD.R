library(TCGAbiolinks)
library(survminer)
library(survival)
library(finalfit)

# Download survival data into your current directory, and save within tcga.clinical list
# NOTE: This will download the data to your current directory. 
# Would recommend downloading the data whilst in your preferred directory
studies <- list.dirs(path = 'GDCdata', full.names = FALSE, recursive = FALSE)

tcga.clinical <- list()
for (study in studies) {
  print(paste0('Download data: ', study, '...'))
  query <- GDCquery_clinic(study, type = 'clinical')
  
  query <- query[, c('submitter_id',
                     'tumor_stage',
                     'age_at_diagnosis',
                     'days_to_last_follow_up',
                     'gender',
                     'vital_status',
                     'days_to_death')]
  query$Study <- strsplit(study, '-')[[1]][2]
  
  tcga.clinical[[strsplit(study, '-')[[1]][2]]] <- query
}
surv.data <- do.call(rbind, tcga.clinical)

# Add TMD categories
load('programme_scores_and_TMD_assignments.RData')
prog_expr <- prog_expr[, c('PatientID',
                           'TMD_two_categories', 'TMD_three_categories',
                           'TMD_two_categories_detailed', 'TMD_three_categories_detailed')]
surv.tmd <- merge(x = surv.data, y = prog_expr,
                  by.x = 'submitter_id', by.y = 'PatientID')

# Basic survival analysis
surv.tmd$time <- ifelse(surv.tmd$vital_status == 'Dead',
                        surv.tmd$days_to_death, surv.tmd$days_to_last_follow_up)
surv.tmd$status <- ifelse(surv.tmd$vital_status == 'Dead', 1, 0)
surv.tmd2 <- surv.tmd[!is.na(surv.tmd$time),]

# Fit survival curves
surv.tmd2$TMD <- sapply(surv.tmd2$TMD_three_categories_detailed,
                        function(x) ifelse(x == 'Angiogenic and Immunogenic Dormancy', 'A/ID',
                                           ifelse(x == 'Angiogenic Dormancy', 'AD',
                                                  ifelse(x == 'Immunogenic Dormancy', 'ID',x))))
surv.tmd2$TMD <- factor(surv.tmd2$TMD, levels = c('NO','MID','AD','ID','A/ID'))
surv.tmd2$TMD_three_categories <- factor(surv.tmd2$TMD_three_categories, levels = c('NO','MID','YES'))

# Kaplan-Meier plots
fit_tmd <- survfit(Surv(time, status) ~ TMD_three_categories, data = surv.tmd2) # standard categories
p_tmd <- ggsurvplot(fit_tmd, pval = TRUE, conf.int = TRUE,
                    risk.table = TRUE, risk.table.y.text.col = TRUE)
print(p_tmd)
ggsave(filename = 'Figures/KaplanMeier_TMD.pdf', plot = p_tmd$plot)

fit_tmd2 <- survfit(Surv(time, status) ~ TMD, data = surv.tmd2) # detailed categories
p_tmd2 <- ggsurvplot(fit_tmd2, pval = TRUE, conf.int = TRUE,
                     risk.table = TRUE, risk.table.y.text.col = TRUE,
                     palette = c('grey','grey95','#D95F02','#7570B3','#1B9E77'))
print(p_tmd2)
ggsave(filename = 'Figures/KaplanMeier_TMDdetailed.pdf', plot = p_tmd2$plot)

# Kaplan-Meier plots by Tumour Stage
surv.tmd2$TumorStage <- sapply(surv.tmd2$tumor_stage,
                               function(x) ifelse(x %in% c('stage i','stage ia','stage ib'),'Stage 1',
                                                  ifelse(x %in% c('stage ii','stage iia','stage iib','stage iic'),'Stage 2',
                                                         ifelse(x %in% c('stage iii','stage iiia','stage iiib','stage iiic'),'Stage 3',
                                                                ifelse(x %in% c('stage iv','stage iva','stage ivb','stage ivc'),'Stage 4',NA)))))

surv.early <- surv.tmd2[surv.tmd2$TumorStage %in% c('Stage 1','Stage 2'), ]
surv.late <- surv.tmd2[surv.tmd2$TumorStage %in% c('Stage 3','Stage 4'), ]

p_early <- ggsurvplot(survfit(Surv(time, status) ~ TMD_three_categories, data = surv.early), 
                    pval = TRUE); print(p_early); ggsave(filename = 'Figures/KaplanMeier_EarlyStage.pdf', plot = p_early$plot)
p_late <- ggsurvplot(survfit(Surv(time, status) ~ TMD_three_categories, data = surv.late), 
                    pval = TRUE); print(p_late); ggsave(filename = 'Figures/KaplanMeier_LateStage.pdf', plot = p_late$plot)

# Hazard Ratio comparisons: TMD_three_categories
form1 <- as.formula('Surv(time, status) ~ TMD_three_categories + age_at_diagnosis + gender')
form2 <- as.formula('Surv(time, status) ~ TMD_three_categories + age_at_diagnosis + gender + Study')
form3 <- as.formula('Surv(time, status) ~ TMD_three_categories + age_at_diagnosis + gender + TumorStage')
form4 <- as.formula('Surv(time, status) ~ TMD_three_categories + age_at_diagnosis + gender + Study + TumorStage')
formulas <- c(form1, form2, form3, form4)

for (i in 1:length(formulas)) {
  form.i <- formulas[[i]]
  print(paste0('Plotting: ',form.i,'...'))
  g <- ggforest(model = coxph(form.i, data = surv.tmd2))
  ggsave(filename = paste0('Figures/HazardRatio_exp',i,'.pdf'), plot = g,
         width = 8, height = 10)
}

# Hazard Ratio comparisons: TMD_three_categories
hr.collate <- data.frame(formula = c('TMD/Age/Gender','+ Study','+ Stage','+ Study + Stage'),
                         HR = NA, L = NA, U = NA)
for (i in 1:length(formulas)) {
  r <- coxph(formulas[i][[1]], data = surv.tmd2)
  hr.collate$HR[i] <- coef(r)[2]
  hr.collate$L[i] <- confint(r)[2,1]
  hr.collate$U[i] <- confint(r)[2,2]
}
hr.collate$formula <- factor(hr.collate$formula, levels = hr.collate$formula)

p_hr <- ggplot(hr.collate, aes(x = 1, y = HR)) + theme_bw() +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L)) +
  labs(x = element_blank(), y = 'log Hazard Ratio (95% CI)') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  facet_wrap(~ formula, nrow = 1); print(p_hr); ggsave(filename = 'Figures/CoxAnalysis_TMD.pdf', plot = p_hr)

# Hazard Ratio comparisons: TMD_three_categories_detailed
form1 <- as.formula('Surv(time, status) ~ TMD + age_at_diagnosis + gender')
form2 <- as.formula('Surv(time, status) ~ TMD + age_at_diagnosis + gender + Study')
form3 <- as.formula('Surv(time, status) ~ TMD + age_at_diagnosis + gender + TumorStage')
form4 <- as.formula('Surv(time, status) ~ TMD + age_at_diagnosis + gender + Study + TumorStage')
formulas <- c(form1, form2, form3, form4)

hr.collate <- data.frame()
for (form in formulas) {
  r <- coxph(form, data = surv.tmd2)
  r.res <- data.frame(TMD = c('AD','ID','A/ID'),
                      HR = coef(r)[2:4],
                      L = confint(r)[2:4,1],
                      U = confint(r)[2:4,2],
                      formula = ifelse(form == form1, 'TMD/Age/Gender',
                                       ifelse(form == form2, '+ Study',
                                              ifelse(form == form3, '+ Stage',
                                                     ifelse(form == form4, '+ Study + Stage',NA)))))
  hr.collate <- rbind(hr.collate, r.res)
}
hr.collate$formula <- factor(hr.collate$formula, levels = c('TMD/Age/Gender','+ Study','+ Stage','+ Study + Stage'))
hr.collate$TMD <- factor(hr.collate$TMD, levels = c('ID','AD','A/ID'))

p_hr <- ggplot(hr.collate, aes(x = TMD, y = HR)) + theme_bw() +
  geom_point(aes(colour = TMD), size = 4) + geom_errorbar(aes(ymax = U, ymin = L)) +
  labs(x = 'Tumour Mass Dormancy Group', y = 'log Hazard ratio (95% CI)') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red') +
  scale_colour_manual(values = c('A/ID' = '#1B9E77',
                                 'AD' = '#D95F02',
                                 'ID' = '#7570B3')) +
  facet_wrap(~ formula, nrow = 1); print(p_hr); ggsave(filename = 'Figures/CoxAnalysis_TMD_detailed.pdf', plot = p_hr)

# TMD categorisation by tumour stage
p_tmdStage <- ggplot(surv.tmd2[!is.na(surv.tmd2$TumorStage),], aes(x = TumorStage, fill = TMD_three_categories)) +
  geom_bar(position = 'fill') + scale_y_continuous(labels = scales::percent)
print(p_tmdStage)
ggsave(filename = 'Figures/Barplot_TMD_by_Stage.pdf', plot = p_tmdStage)

p_tmdStage <- ggplot(surv.tmd2[!is.na(surv.tmd2$TumorStage),], aes(x = TumorStage, fill = TMD)) +
  geom_bar(position = 'fill') + scale_y_continuous(labels = scales::percent)
print(p_tmdStage)
ggsave(filename = 'Figures/Barplot_TMDdet_by_Stage.pdf', plot = p_tmdStage)

