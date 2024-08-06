library(data.table)

# source helper functions for building module annotations
source('workflow/scripts/e2g_helper_functions.R')

build_e2g_predictions <- function(model, tissue) {

  # define models with thresholded predictions
  thresholded_models <- c(
    'ENCODE_E2G_DNase_only',
    'ABC_DNase_only',
    'ENCODE_E2G_DNase_HiC'
  )

  # define baseline models
  baseline_models <- c(
    'dist_to_tss',
    'nearest_expressed_gene'
  )

  # subset enhancer-gene prediction files to use based on tissue + model
  if (model %in% thresholded_models) {
    input_path <- paste0('resources/', model, '/')
    input_path <- paste0(input_path, 'thresholded_predictions/')
    pred_files <- list.files(input_path, pattern = "thresholded_predictions")
    pred_files <- subset_pred_files(pred_files, tissue)
  }

  if (model %in% baseline_models) {
    input_path <- 'resources/baseline_models/'
    pred_files <- read.table(paste0(input_path, 'preds.txt'), header = F)[ , 1]
    pred_files <- subset_pred_files(pred_files, tissue)
    pred_files <- paste0(
      pred_files,
      '/',
      model, 
      ".tsv.gz"
    )
  }

  # create list to hold combined e2g predictions
  pooled_e2g <- list()

  # read through each E2G prediction set and add to pooled list
  for (i in 1:length(pred_files)) {

    # read in E2G predictions for each cell line
    e2g_preds <- data.frame(fread(paste0(input_path, "/", pred_files[i])))

    # additional filtering based on model type and processing
    if (model %in% thresholded_models) {
      e2g_preds <- cbind.data.frame(
        e2g_preds$X.chr, 
        e2g_preds$start, 
        e2g_preds$end, 
        e2g_preds$TargetGene, 
        e2g_preds$Score
      )
    }
    if (model == 'dist_to_tss') {
      e2g_preds <- e2g_preds[which(e2g_preds$Score < 53976), ]
      e2g_preds <- cbind.data.frame(
        e2g_preds$X.chr, 
        e2g_preds$start, 
        e2g_preds$end, 
        e2g_preds$TargetGene,
        1
      )
    }
    if (model == 'nearest_expressed_gene') {
      e2g_preds <- e2g_preds[which(e2g_preds$Score > 0), ]
      e2g_preds <- cbind.data.frame(
        e2g_preds$X.chr, 
        e2g_preds$start, 
        e2g_preds$end, 
        e2g_preds$TargetGene,
        1
      )
    }

    colnames(e2g_preds) <- c(
      "chr", 
      "start", 
      "end", 
      "TargetGene", 
      "Score"
    )
    pooled_e2g[[i]] <- e2g_preds
  }

  # concatenate all of the predictions together
  pooled_e2g <- do.call(rbind, pooled_e2g)

  # define output file
  output_dir <- 'results/e2g_predictions/'
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  output_file <- paste0(
    output_dir,
    model,
    '_',
    tissue,
    '_',
    'e2g_predictions.txt'
  )

  # write pooled E2G links to output file
  fwrite(
    pooled_e2g,
    output_file,
    sep = "\t", 
    quote=FALSE, 
    row.names=FALSE, 
    col.names=TRUE
  )
}