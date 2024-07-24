# This script builds E2G "modules". The modules are BED files which contain E2G
# predictions for a given tissue and given type of E2G model. The tissues and
# model types are defined below. These BED files will eventually be used for
# benchmarking various models of E2G predictions.

library(data.table)
library(R.utils)

# source helper functions for building module annotations
source('workflow/scripts/e2g_helper_functions.R')

# create output directory for E2G modules
output_dir <- 'results/e2g_modules/'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# define path to gene score file and read in gene scores
genescore_path = snakemake@input[[1]]
gene_scores = read.delim(genescore_path, header = F)
scores = gene_scores[,2]
names(scores) = gene_scores[,1]

# define list of model types
models <- c(
  'ENCODE_E2G_DNase_only',
  'ENCODE_E2G_DNase_HiC',
  'ABC_DNase_only',
  'ENCODE_E2G_ext',
  'ABC_ext',
  'scE2G'
)

for (i in 1:length(models)) {

  # define current model
  model <- models[i]

  # define tissues and build E2G modules for threshold-based models
  if (model %in% c('ENCODE_E2G_DNase_only', 'ENCODE_E2G_DNase_HiC', 'ABC_DNase_only')) {
    tissues <- c(
      'ALL', # all biosamples combined
      'BLD', # blood-related biosamples
      'BRN', # brain-related biosamples
      'K562', # K562-related biosamples
      'GM12878', # GM12878-related biosamples
      'Mono' # monocyte-related biosamples
    )
    for (j in 1:length(tissues)) {
      tissue <- tissues[j]
      build_e2g_module_thresholded(
        scores = scores,
        model = model,
        tissue = tissue,
        output_dir = output_dir
      )
    }
  }

  # define tissues and build E2G modules for extended models
  if (model %in% c('ENCODE_E2G_ext', 'ABC_ext')) {
    tissues < c(
      'K562', # K562-related biosamples
      'GM12878' # GM12878-related biosamples
    )
    for (j in 1:length(tissues)) {
      tissue <- tissues[j]
      build_e2g_module_extended(
        scores = scores,
        model = model,
        tissue = tissue,
        output_dir = output_dir
      )
    }
  }

  # build E2G modules for scE2G model in K562 cells
  if (model == 'scE2G') {
    build_e2g_module_sce2g(
      scores = scores,
      output_dir = output_dir
    )
  }
}

