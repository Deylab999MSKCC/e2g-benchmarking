# This script creates a demo E2G prediction file for input to the snakemake 
# pipeline.
#
# Author: Karthik Guruvayurappan

library(data.table)
library(R.utils)

source('workflow/scripts/e2g_helper_functions.R')

# read in E2G predictions from ENCODE E2G DNase models
model <- 'ENCODE_E2G_DNase_only'
tissue <- 'K562'

# define directory to predictions and get list of files
# the files contain E2G predictions for different ENCODE biosamples
pred_dir <- paste0(
    'resources/',
    model,
    '/thresholded_predictions'
)
pred_files <- list.files(pred_dir, pattern = 'thresholded_predictions')
pred_files <- subset_pred_files(pred_files, tissue)

# create dataframe to hold combined E2G predictions
pooled_e2g <- c()

for (i in 1:length(pred_files)) {

    # read in current biosample predictions
    current_pred <- data.frame(fread(paste0(pred_dir, '/', pred_files[i])))

    # reformat current biosample predictions
    current_pred <- cbind.data.frame(
      current_pred$X.chr, 
      current_pred$start,
      current_pred$end,
      current_pred$TargetGene,
      current_pred$Score
    )
    colnames(current_pred) <- c('chr', 'start', 'end', 'TargetGene', 'Score')

    # add to pooled E2G links
    pooled_e2g <- rbind(pooled_e2g, current_pred)
}

# read in list of genes to keep
genescores <- read.table('resources/gene_scores.txt')
pooled_e2g <- pooled_e2g[pooled_e2g$TargetGene %in% genescores$V1, ]

# write to output BED file
write.table(
    pooled_e2g,
    paste0('resources', '/test_', model, '_', tissue, '.bed'),
    sep = '\t', 
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

