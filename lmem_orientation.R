library(lme4)
library(car)

input_file <- snakemake@input[["data"]]
output_file <- snakemake@output[[1]]
#input_file <- '//exports.hps.kau.science.roche.com/pred/rpmda/BN40423-v2/jira/RPMDA-10324-effect-phone-location-balance-test/results/k_means/clusters_combined_digital_features.csv'

data <- read.table(input_file, sep=",", header=T)

digital_features = unique(data$feature_digital)

results <- NULL
for (d in digital_features) {
  
  data_subset = data[(data$feature_digital==d),]
  
  lmem <- lmer(value_digital ~ cluster + (1|subject_id), data=data_subset, REML=F)
  summary_lmem <- summary(lmem)
  results_lmem <- Anova(lmem)
  
  results_lmem$feature_digital <- d
  
  results = rbind(results, results_lmem)
}

write.csv(results, output_file, row.names=TRUE)
