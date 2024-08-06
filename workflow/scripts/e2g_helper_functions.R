subset_pred_files <- function(pred_files, tissue_name) {

  # keep all files if tissue name is ALL
  if (tissue_name == 'ALL') {
    keep_pred_files = pred_files
    return (keep_pred_files)
  }
  
  # keep blood-related biosamples if tissue is BLD
  if (tissue_name == 'BLD') {

    # define search terms
    search_terms <- c(
      'T-', 'T_cell', 'T_follicular_helper_cell', 'hemato', 'B_cell', '_CD', 
      'B-', 'monocyte', 'killer', 'GM1', 'K562', 'myeloid', 'lymphoid'
    )

    # keep all predictions that match search terms
    keep_pred_files <- c()
    for (i in 1:length(search_terms)) {
      search_pred_files <- pred_files[grep(search_terms[i], pred_files)]
      keep_pred_files <- c(keep_pred_files, search_pred_files)
    }

    # remove duplicates and return
    keep_pred_files <- unique(keep_pred_files)
    return (keep_pred_files)
  }

  # keep brain-related biosamples if tissue is BRN
  if (tissue_name == 'BRN') {

    # define search terms
    search_terms <- c(
      '_astro', 'neur', 'cerebell', 'cortex', 'brain', 'glia'
    )

    # keep all predictions that match search terms
    keep_pred_files <- c()
    for (i in 1:length(search_terms)) {
      search_pred_files <- pred_files[grep(search_terms[i], pred_files)]
      keep_pred_files <- c(keep_pred_files, search_pred_files)
    }

    # remove duplicates and return
    keep_pred_files <- unique(keep_pred_files)
    return (keep_pred_files)
  }

  # keep K562 biosamples if tissue is K562
  if (tissue_name == 'K562') {
    keep_pred_files <- pred_files[grep('K562', pred_files)]
    return (keep_pred_files)
  }

  # keep Mono biosamples if tissue is Mono
  if (tissue_name == 'Mono') {
    keep_pred_files <- pred_files[grep('mono', pred_files)]
    return (keep_pred_files)
  }

  # keep GM12878 biosamples if tissue is Mono
  if (tissue_name == 'GM12878') {
    keep_pred_files <- pred_files[grep('GM12878', pred_files)]
    return (keep_pred_files)
  }
}


build_e2g_module_thresholded <- function(scores, model, tissue, output_dir) {
  
  # define directory to predictions and get list of files
  # the files contain E2G predictions for different ENCODE biosamples
  pred_dir <- paste0(
    'resources/',
    model,
    '/thresholded_predictions'
  )
  pred_files <- list.files(pred_dir, pattern = 'thresholded_predictions')

  # filter prediction files to match tissue of interest
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
      current_pred$TargetGene
    )
    colnames(current_pred) <- c('chr', 'start', 'end', 'TargetGene')

    # add to pooled E2G links
    pooled_e2g <- rbind(pooled_e2g, current_pred)
  }

  # add gene scores to BED file and write output file
  target_gene_scores <- match(pooled_e2g$TargetGene, names(scores))
  target_gene_scores <- as.numeric(scores)[target_gene_scores]
  target_gene_scores[is.na(target_gene_scores)] <- 0
  pooled_e2g <- cbind(pooled_e2g[, c(1:3)], target_gene_scores)
  write.table(
    pooled_e2g,
    paste0(output_dir, '/', model, '_', tissue, '.bed'),
    sep = '\t', 
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}


build_e2g_module_extended <- function(scores, model, tissue, output_dir) {
  
  # read in predictions using model + tissue name combination
  pooled_e2g <- c()
  if (model == 'ENCODE_E2G_ext') {
    pooled_e2g = data.frame(fread(paste0(
      'resources/',
      model,
      '/predictions_ext/',
      tissue,
      '/encode_e2g_predictions_threshold0.336.tsv.gz'
    )))
    pooled_e2g = cbind.data.frame(
      pooled_e2g[,1], pooled_e2g[,2], pooled_e2g[,3], pooled_e2g[,6]
    )
  }
  
  if (model == 'ABC_ext') {
    pooled_e2g = data.frame(fread(paste0(
      'resources/',
      model,
      '/predictions_ext/',
      tissue,
      '/Predictions/',
      'ABC_pred_filtered_0.027_selfPromoters.tsv.gz'
    )))
    pooled_e2g = cbind.data.frame(
      pooled_e2g[,1], pooled_e2g[,2], pooled_e2g[,3], pooled_e2g[,11]
    )
  }

  colnames(pooled_e2g) = c("chr", "start", "end", "TargetGene")
  target_gene_scores = match(pooled_e2g$TargetGene, names(scores))
  target_gene_scores = as.numeric(scores)[target_gene_scores]
  target_gene_scores[is.na(target_gene_scores)] = 0
  pooled_e2g = cbind(pooled_e2g[,c(1:3)], target_gene_scores)

  # write to output
  write.table(
    pooled_e2g,
    paste0(output_dir, '/', model, '_', tissue, '.bed'),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}


build_e2g_module_sce2g <- function(scores, output_dir) {

  # define path to predictions
  pred_path <- 'resources/scE2G/IGVF_scE2G/'

  # define prediction files
  pred_files <- c(
    'K562_Xu_multiome_7features_genomewide_predictions.tsv.gz',
    'K562_Xu_multiome_7features_intx_genomewide_predictions.tsv.gz',
    'K562_Xu_scATAC_6features_genomewide_predictions.tsv.gz',
    'K562_Xu_scATAC_6features_intx_genomewide_predictions.tsv.gz',
    'K562_Xu_scATAC_6features_genomewide_predictions.tsv.gz',
    'Pairs.Kendall.tsv.gz'
  )

  # define score types and thresholds for each prediction file type
  pred_types <- c(
    rep('ENCODE_rE2G', 4),
    'ABC',
    'Kendall'
  )
  pred_thresholds <- c(
    0.171,
    0.176,
    0.154,
    0.22,
    0.015,
    0.011
  )

  for (i in 1:length(pred_files)) {

    # read in E2G predictions
    pooled_e2g = data.frame(fread(paste0(
      pred_path, 
      pred_files[i]
    )))

    # threshold based on prediction type
    if (pred_types[i] == 'ENCODE_rE2G') {
      pooled_e2g <- pooled_e2g[
        which(pooled_e2g$ENCODE.rE2G.Score > pred_thresholds[i]), 
      ]
      pooled_e2g = cbind.data.frame(
        pooled_e2g[,1], pooled_e2g[,2], pooled_e2g[,3], pooled_e2g[,6]
      )
    }

    if (pred_types[i] == 'ABC') {
      pooled_e2g <- pooled_e2g[
        which(pooled_e2g$ABC.Score > pred_thresholds[i]), 
      ]
      pooled_e2g = cbind.data.frame(
        pooled_e2g[,1], pooled_e2g[,2], pooled_e2g[,3], pooled_e2g[,6]
      )
    }

    if (pred_types[i] == 'Kendall') {
      pooled_e2g <- pooled_e2g[
        which(pooled_e2g$Kendall > pred_thresholds[i]), 
      ]
      pooled_e2g = cbind.data.frame(
        pooled_e2g[,1], pooled_e2g[,2], pooled_e2g[,3], pooled_e2g[,4]
      )
    }

    # reformat table
    colnames(pooled_e2g) <- c('chr', 'start', 'end', 'TargetGene')
    target_gene_scores = match(pooled_e2g$TargetGene, names(scores))
    target_gene_scores = as.numeric(scores)[target_gene_scores]
    target_gene_scores[is.na(target_gene_scores)] = 0
    pooled_e2g = cbind(pooled_e2g[,c(1:3)], target_gene_scores)

    # write output file
    model_name <- strsplit(pred_files[i], '.', fixed = TRUE)[[1]][1]
    write.table(
      pooled_e2g,
      paste0(
        output_dir, '/', 'scE2G', '_', model_name, '_', pred_types[i], '.bed'
      ),
      sep = '\t',
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
  }
}

