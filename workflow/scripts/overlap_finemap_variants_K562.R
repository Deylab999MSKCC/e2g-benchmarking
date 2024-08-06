library(data.table)
library(R.utils)

###############################################################################
# Part 1: Compute values for enrichment of fine-mapped variants in e2g preds
###############################################################################

# read in lists of finemapped variants
finemapped_variant_lists <- list()

for (i in 23:34) {
    
    # read in finemapped variants and add to list
    finemapped_variant_file <- snakemake@input[[i]]
    finemapped_variant_list <- read.table(finemapped_variant_file)[, 1]
    finemapped_variant_lists[[i - 22]] <- finemapped_variant_list
}

# create matrices to hold various values to track
total_common_variants <- 0      # number of common variants across all chromosomes
common_variants_e2g_overlap_by_chrom <- rep(0, 22) # number of common variants overlapping an enhancer per chromosome
finemapped_variants_by_chrom_by_list <- matrix(0, 22, 12) # number of finemapped variants per chrom (row) and per list (column)
finemapped_variants_e2g_overlap_by_chrom_by_list <- matrix(0, 22, 12) # number of finemapped variants with E2G overlap per chrom (row) and per list (column)

for (chrom in 1:22) {

    # read in annotated variants for current chromosome
    chrom_annot_variants <- data.frame(fread(
        snakemake@input[[chrom]]
    ))

    # add number of common variants to total common variant count
    total_common_variants <- total_common_variants + nrow(chrom_annot_variants)

    # add number of common variants overlapping an enhancer per chromosome
    chrom_common_variants_e2g_overlap <- chrom_annot_variants[, 5]
    common_variants_e2g_overlap_by_chrom[chrom] <- sum(chrom_common_variants_e2g_overlap)

    # create vectors to hold information across lists
    chrom_finemapped_variants_by_list <- rep(0, length(finemapped_variant_lists)) # for this chrom, finemapped variants per list
    chrom_finemapped_variants_e2g_overlap_by_list <- rep(0, length(finemapped_variant_lists)) # for this chrom, finemapped variants overlapping enhancer per list

    for (i in 1:length(finemapped_variant_lists)) {

        # create empty vector
        chrom_finemapped_variants <- rep(0, nrow(chrom_annot_variants))

        # annotate chromosome variants overlapping finemapped variants
        chrom_finemapped_variants[match(intersect(chrom_annot_variants$SNP, finemapped_variant_lists[[i]]), chrom_annot_variants$SNP)] = 1

        # add total number of finemapped variants for this list
        chrom_finemapped_variants_by_list[i] <- sum(chrom_finemapped_variants)

        # add total number of finemapped variants overlapping E2G annotations for this list
        chrom_finemapped_variants_e2g_overlap_by_list[i] <- sum(chrom_common_variants_e2g_overlap[chrom_finemapped_variants == 1])
    }

    # add back overlap lists (finemap and E2G + finemap) to matrices
    finemapped_variants_by_chrom_by_list[chrom, ] <- chrom_finemapped_variants_by_list
    finemapped_variants_e2g_overlap_by_chrom_by_list[chrom, ] <- chrom_finemapped_variants_e2g_overlap_by_list
}

# sum number of finemapped variants per variant list
finemapped_variants_by_list <- colSums(finemapped_variants_by_chrom_by_list)

# sum number of finemap + e2g overlaps per fine-mapped variant list
finemapped_variants_e2g_overlap_by_list <- colSums(finemapped_variants_e2g_overlap_by_chrom_by_list)
print(common_variants_e2g_overlap_by_chrom)

# sum number of overlaps between common variants and E2G
common_variants_e2g_overlap <- sum(common_variants_e2g_overlap_by_chrom)

# compute point estimates for enrichment, precision, and recall (per finemap variant set)
enr_by_list = (finemapped_variants_e2g_overlap_by_list / common_variants_e2g_overlap) / (finemapped_variants_by_list / total_common_variants)
precision_by_list = (finemapped_variants_e2g_overlap_by_list / common_variants_e2g_overlap)
recall_by_list = (finemapped_variants_e2g_overlap_by_list / finemapped_variants_by_list)


###############################################################################
# Part 2: Bootstrapping enrichment, precision, and recall for standard deviation
###############################################################################

# create matrices to hold bootstrap enrichments, precision, and recall
enr_by_boot_by_list <- matrix(0, 100, 12)
precision_by_boot_by_list <- matrix(0, 100, 12)
recall_by_boot_by_list <- matrix(0, 100, 12)

for (nboot in 1:100) {

    # bootstrap chromosomes
    idx = sample(1:22, 22, replace=TRUE)

    # compute bootstrap enrichment value
    enr_by_boot_by_list[nboot, ] = (
        colSums(finemapped_variants_e2g_overlap_by_chrom_by_list[idx, ]) / colSums(finemapped_variants_by_chrom_by_list[idx, ])
    ) / (sum(common_variants_e2g_overlap_by_chrom[idx]) / total_common_variants)

    # compute bootstrap precision value
    precision_by_boot_by_list[nboot, ] = (
        colSums(finemapped_variants_e2g_overlap_by_chrom_by_list[idx, ]) / sum(common_variants_e2g_overlap_by_chrom[idx])
    )

    # compute bootstrap recall value
    recall_by_boot_by_list[nboot, ] = (colSums(finemapped_variants_e2g_overlap_by_chrom_by_list[idx, ]) / colSums(finemapped_variants_by_chrom_by_list[idx, ]))
}

# compute standard deviation from bootstrap estimates (bootstrap chromosomes)
sd_enr_by_list = apply(enr_by_boot_by_list, 2, sd)
sd_precision_by_list = apply(precision_by_boot_by_list, 2, sd)
sd_recall_by_list = apply(recall_by_boot_by_list, 2, sd)

###############################################################################
# Part 3: Permutation testing/bootstrapping for p-values of enrichments
###############################################################################

# perform bootstrapping of E2G annotations to get confidence intervals
finemapped_variants_e2g_overlap_by_boot_by_list <- matrix(0, 100, 12) # row is bootstrap iteration, column is list

for (chrom in 1:22) {
  
    # read in annotated variants for current chromosome
    chrom_annot_variants <- data.frame(fread(
        snakemake@input[[chrom]]
    ))

    # whether variants overlap an ENCODE E2G annotation
    chrom_common_variants_e2g_overlap <- chrom_annot_variants[, 5]

    # compute indices of fine-mapped variants for chromosome (for faster computation)
    chrom_finemapped_variant_idxs_by_list <- list()

    for (i in 1:length(finemapped_variant_lists)) {
        chrom_finemapped_variant_idxs_by_list[[i]] <- match(intersect(chrom_annot_variants$SNP, finemapped_variant_lists[[i]]), chrom_annot_variants$SNP)
    }

    # create variable for chromosome-level results
    chrom_finemapped_variants_e2g_overlap_by_boot_by_list <- c()

    # perform bootstrap iterations computing finemap variant + e2g by list for this chromosome
    for (i in 1:100) {

        # resample whether variants overlap E2G annotation
        boot_chrom_common_variants_e2g_overlap <- sample(
            chrom_common_variants_e2g_overlap,
            length(chrom_common_variants_e2g_overlap),
            replace = F
        )

        boot_chrom_finemapped_variants_e2g_overlap_by_list <- rep(0, length(finemapped_variant_lists))

        for (j in 1:length(finemapped_variant_lists)) {

            # create empty vector for each annotated variant
            chrom_finemapped_variants <- rep(0, nrow(chrom_annot_variants))

            # annotate chromosome variants overlapping finemapped variants
            chrom_finemapped_variants[chrom_finemapped_variant_idxs_by_list[[j]]] = 1

            # sum number of finemapped variants that overlap E2G annotation
            boot_chrom_finemapped_variants_e2g_overlap_by_list[j] <- sum(boot_chrom_common_variants_e2g_overlap[chrom_finemapped_variants == 1])
        }

        # add to bootstrap iterations
        chrom_finemapped_variants_e2g_overlap_by_boot_by_list <- rbind(
            chrom_finemapped_variants_e2g_overlap_by_boot_by_list,
            boot_chrom_finemapped_variants_e2g_overlap_by_list
        )
    }

    # add to bootstrap across all chromosomes
    finemapped_variants_e2g_overlap_by_boot_by_list <- finemapped_variants_e2g_overlap_by_boot_by_list + chrom_finemapped_variants_e2g_overlap_by_boot_by_list
}

# create matrices to hold bootstrap enrichments, precision, and recall
enr_by_boot_by_list <- matrix(0, 100, 12)
precision_by_boot_by_list <- matrix(0, 100, 12)
recall_by_boot_by_list <- matrix(0, 100, 12)

for (nboot in 1:100) {
    enr_by_boot_by_list[nboot, ] = (finemapped_variants_e2g_overlap_by_boot_by_list[nboot, ] / common_variants_e2g_overlap) / (finemapped_variants_by_list / total_common_variants)
    precision_by_boot_by_list[nboot, ] = (finemapped_variants_e2g_overlap_by_boot_by_list[nboot, ] / common_variants_e2g_overlap)
    recall_by_boot_by_list[nboot, ] = (finemapped_variants_e2g_overlap_by_boot_by_list[nboot, ] / finemapped_variants_by_list)
}

# compute p-values for enrichment, precision, and recall using bootstrap (resampling)
p_enr_by_list = c()
p_precision_by_list = c()
p_recall_by_list = c()

for (i in 1:12) {
    p_enr_by_list = c(p_enr_by_list, pnorm(enr_by_list[i], mean(enr_by_boot_by_list[,i]), sd(enr_by_boot_by_list[,i]), lower.tail = F))
    p_precision_by_list = c(p_precision_by_list, pnorm(precision_by_list[i], mean(precision_by_boot_by_list[,i]), sd(precision_by_boot_by_list[,i]), lower.tail = F))
    p_recall_by_list = c(p_recall_by_list, pnorm(recall_by_list[i], mean(recall_by_boot_by_list[,i]), sd(recall_by_boot_by_list[,i]), lower.tail = F))
}

# write to output data frame
outdf = cbind(enr_by_list, sd_enr_by_list, p_enr_by_list, precision_by_list, sd_precision_by_list, p_precision_by_list, recall_by_list, sd_recall_by_list, p_recall_by_list)
colnames(outdf) = c("ENR", "sENR", "pENR", "Precision", "s.Precision", "p.Precision", "Recall", "s.Recall", "p.Recall")
rownames(outdf) = c("Ulirsch_ALL_RBC_10", "Ulirsch_ALL_RBC_50",
                    "Ulirsch_ALL_MCH_10", "Ulirsch_ALL_MCH_50",
                    "Ulirsch_ALL_MCV_10", "Ulirsch_ALL_MCV_50",
                    "Ulirsch_ALL_Hb_10", "Ulirsch_ALL_Hb_50",
                    "Ulirsch_ALL_HbA1c_10", "Ulirsch_ALL_HbA1c_50",
                    "Ulirsch_ALL_RBCdistwidth_10", "Ulirsch_ALL_RBCdistwidth_50")

write.table(outdf, file = snakemake@output[[1]],
            col.names = T, row.names = T, sep = "\t", quote=F)

