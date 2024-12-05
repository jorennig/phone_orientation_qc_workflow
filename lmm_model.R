library(tidyverse)
library(lme4)

data_file <- snakemake@input[["data"]]
output_file <- snakemake@output[["model_object"]]
lmm_formula_string <- snakemake@params[["lmm_formula"]]
lmm_formula <- as.formula(lmm_formula_string)

possibly_lmm = possibly(.f = lmer, otherwise = NA)
lmm_model <- function(df) {
  possibly_lmm(lmm_formula, 
  data = df, 
  REML = TRUE)
}

d <- read_tsv(data_file, col_types = list(
  subject_id = col_character()
))

d_grouped <- d %>% group_by(feature) %>% nest()

d_grouped_model <- d_grouped %>%
  mutate(model = map(data, lmm_model))

write_rds(d_grouped_model, output_file)











