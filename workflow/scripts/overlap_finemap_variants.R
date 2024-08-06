library(data.table)
library(R.utils)

# read in command line arguments
args <- commandArgs(trailingOnly = TRUE)
annot_variants_dir <- args[1]
output_dir <- args[2]
finemap_var_dir <- args[3]

# create list of finemap variant files to use
finemap_var_files <- c(
  'FINEMAP_Ulirsch_ALL_RBC_10', 
  'FINEMAP_Ulirsch_ALL_RBC_50',
  'FINEMAP_Ulirsch_ALL_MCH_10',
  'FINEMAP_Ulirsch_ALL_MCH_50',
  'FINEMAP_Ulirsch_ALL_MCV_10',
  'FINEMAP_Ulirsch_ALL_MCV_50',
  'FINEMAP_Ulirsch_ALL_Hb_10',
  'FINEMAP_Ulirsch_ALL_Hb_50',
  'FINEMAP_Ulirsch_ALL_HbA1c_10',
  'FINEMAP_Ulirsch_ALL_HbA1c_50',
  'FINEMAP_Omer_blood_RBC_DISTRIB_WIDTH_10',
  'FINEMAP_Omer_blood_RBC_DISTRIB_WIDTH_50'
)

# iterate through finemap variant lists and read files in
finemap_var_lists <- list()

for (i in 1:length(finemap_var_files)) {

  # get current file
  file <- finemap_var_files[i]
  var_list <- read.table(paste0(finemap_var_dir, file))[, 1]

  # add to lists of finemapped variants
  finemap_var_lists[[i]] <- var_list
}

finemap_var_e2g_overlap <- matrix(0, 22, 12) # matrix of finemapped variants + e2g overlaps
finemap_vars <- matrix(0, 22, 12) # matrix of fine-mapped variants per chrom per list
var_e2g_overlaps <- rep(0, 22)
num_common_variants <- 0

for (chrom in 1:22) {

  # read in annotated variants for current chromosome
  annot_variants <- data.frame(fread(paste0(
    annot_variants_dir,
    '/',
    strsplit(annot_variants_dir, '/')[[1]][3],
    '.',
    chrom,
    '.annot.gz'
  )))

  # add to number of total common variants
  num_common_variants <- num_common_variants + nrow(annot_variants)

  # iterate through lists of finemapped variants
  chrom_finemap_var_e2g_overlap <- rep(0, 12)
  chrom_finemap_vars <- rep(0, 12)

  # whether variants overlap an ENCODE E2G annotation
  chrom_variant_e2g_overlap <- annot_variants[, 5]

  for (i in 1:length(finemap_var_lists)) {

    # create empty vector for each annotated variant
    finemap_annot <- rep(0, nrow(annot_variants))

    # annotate chromosome variants overlapping finemapped variants
    finemap_annot[match(intersect(annot_variants$SNP, finemap_var_lists[[i]]), annot_variants$SNP)] = 1

    # sum number of finemapped variants that overlap E2G annotation
    chrom_finemap_var_e2g_overlap[i] <- sum(chrom_variant_e2g_overlap[finemap_annot == 1])

    # record number of fine-mapped variants per list
    chrom_finemap_vars[i] <- sum(finemap_annot)
  }

  # add back to matrix containing info across all chromosomes
  finemap_var_e2g_overlap[chrom, ] <- chrom_finemap_var_e2g_overlap
  finemap_vars[chrom, ] <- chrom_finemap_vars

  # sum number of common variants overlapping E2G
  var_e2g_overlaps[chrom] <- sum(chrom_variant_e2g_overlap)
}

# sum number of overlaps per fine-mapped variant list
sum_finemap_var_e2g_overlap <- colSums(finemap_var_e2g_overlap)

# sum number of overlaps between common variants and E2G
sum_var_e2g_overlaps <- sum(var_e2g_overlaps)

# sum number of finemapped variants per variant list
sum_finemap_vars <- colSums(finemap_vars)

# perform bootstrapping of E2G annotations to get confidence intervals
boot_finemap_var_e2g_overlap <- matrix(0, 100, 12)

for (chrom in 1:22) {
  
  # record bootstraps per chromosome
  chrom_finemap_var_e2g_overlap <- c()

  # read in annotated variants for current chromosome
  annot_variants <- data.frame(fread(paste0(
    annot_variants_dir,
    '/',
    strsplit(annot_variants_dir, '/')[[1]][3],
    '.',
    chrom,
    '.annot.gz'
  )))

  # whether variants overlap an ENCODE E2G annotation
  chrom_variant_e2g_overlap <- annot_variants[, 5]

  # compute lists of finemapped variant indices
  finemap_var_indices <- list()

  for (i in 1:length(finemap_var_lists)) {
    finemap_var_indices[[i]] <- match(intersect(annot_variants$SNP, finemap_var_lists[[i]]), annot_variants$SNP)
  }

  for (i in 1:100) {

    # resample whether variants overlap E2G annotation
    boot_chrom_variant_e2g_overlap <- sample(
      chrom_variant_e2g_overlap,
      length(chrom_variant_e2g_overlap),
      replace = F
    )

    boot_chrom_finemap_var_e2g_overlap <- rep(0, length(finemap_var_lists))

    for (j in 1:length(finemap_var_lists)) {

      # create empty vector for each annotated variant
      finemap_annot <- rep(0, nrow(annot_variants))

      # annotate chromosome variants overlapping finemapped variants
      finemap_annot[finemap_var_indices[[j]]] = 1

      # sum number of finemapped variants that overlap E2G annotation
      boot_chrom_finemap_var_e2g_overlap[j] <- sum(boot_chrom_variant_e2g_overlap[finemap_annot == 1])
    }

    # add to bootstrap iterations
    chrom_finemap_var_e2g_overlap <- rbind(
      chrom_finemap_var_e2g_overlap,
      boot_chrom_finemap_var_e2g_overlap
    )

  }

  # add to bootstrap across all chromosomes
  boot_finemap_var_e2g_overlap <- boot_finemap_var_e2g_overlap + chrom_finemap_var_e2g_overlap
}

# create matrices to hold bootstrap enrichments, precision, and recall
enr_bootmat <- matrix(0, 100, 12)
prec_bootmat <- matrix(0, 100, 12)
recall_bootmat <- matrix(0, 100, 12)


for(nboot in 1:100){

  # bootstrap chromosomes
  idx = sample(1:22, 22, replace=TRUE)

  # compute bootstrap enrichment value
  enr_bootmat[nboot, ] = (
    colSums(finemap_var_e2g_overlap[idx, ]) / colSums(finemap_vars[idx, ])
  ) / (sum(var_e2g_overlaps[idx]) / num_common_variants)
  prec_bootmat[nboot, ] = (
    colSums(finemap_var_e2g_overlap[idx, ]) / sum(var_e2g_overlaps[idx])
  )
  recall_bootmat[nboot, ] = (colSums(finemap_var_e2g_overlap[idx, ]) / colSums(finemap_vars[idx, ]))
}

# compute point estimates for true enrichment, precision, and recall (per finemap variant set)
enr = (sum_finemap_var_e2g_overlap/sum_var_e2g_overlaps) / (sum_finemap_vars / num_common_variants)
precision = (sum_finemap_var_e2g_overlap / sum_var_e2g_overlaps)
recall = (sum_finemap_var_e2g_overlap / sum_finemap_vars)

# compute standard deviations from bootstrap estimates (bootstrap chromosomes)
senr = apply(enr_bootmat, 2, sd)
sprec = apply(prec_bootmat, 2, sd)
srecall = apply(recall_bootmat, 2, sd)

# compute bootstrap enrichments, precision, and recall (bootstrap labels)
boot_enr = matrix(0, 100, 12)
boot_prec = matrix(0, 100, 12)
boot_recall = matrix(0, 100, 12)

for(mm in 1:100){
  boot_enr[mm, ] = (boot_finemap_var_e2g_overlap[mm, ] / sum_var_e2g_overlaps) / (sum_finemap_vars / num_common_variants)
  boot_prec[mm, ] = (boot_finemap_var_e2g_overlap[mm, ] / sum_var_e2g_overlaps)
  boot_recall[mm, ] = (boot_finemap_var_e2g_overlap[mm, ] / sum_finemap_vars)
}

# compute p-values for enrichment, precision, and recall using bootstrap (resampling)
penr = c()
pprec = c()
precall = c()

for(mm in 1:12){
  penr = c(penr, pnorm(enr[mm], mean(boot_enr[,mm]), sd(boot_enr[,mm]), lower.tail = F))
  pprec = c(pprec, pnorm(precision[mm], mean(boot_prec[,mm]), sd(boot_prec[,mm]), lower.tail = F))
  precall = c(precall, pnorm(recall[mm], mean(boot_recall[,mm]), sd(boot_recall[,mm]), lower.tail = F))
}

# write to output data frame
outdf = cbind(enr, senr, penr, precision, sprec, pprec, recall, srecall, precall)
colnames(outdf) = c("ENR", "sENR", "pENR", "Precision", "s.Precision", "p.Precision", "Recall", "s.Recall", "p.Recall")
rownames(outdf) = c("Ulirsch_ALL_RBC_10", "Ulirsch_ALL_RBC_50",
                    "Ulirsch_ALL_MCH_10", "Ulirsch_ALL_MCH_50",
                    "Ulirsch_ALL_MCV_10", "Ulirsch_ALL_MCV_50",
                    "Ulirsch_ALL_Hb_10", "Ulirsch_ALL_Hb_50",
                    "Ulirsch_ALL_HbA1c_10", "Ulirsch_ALL_HbA1c_50",
                    "Ulirsch_ALL_RBCdistwidth_10", "Ulirsch_ALL_RBCdistwidth_50")

module_name <- strsplit(annot_variants_dir, '/')[[1]][3]
write.table(outdf, file = paste0(output_dir, "/", module_name, ".FINEMAP.RES.txt"),
            col.names = T, row.names = T, sep = "\t", quote=F)