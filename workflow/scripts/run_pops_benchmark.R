library(GenomicRanges)
library(dplyr)
library(data.table)

# configurable parameters
tissue <- 'ALL'
use_pops <- FALSE
num_pops <- 2
num_eg_genes <- 1

# read in enhancer-gene links
e2g_preds <- data.frame(fread(snakemake@input[[1]]))

# retrieve necessary columns for the benchmark (matching IGVF format)
e2g_preds <- e2g_preds[, c(
    'ElementChr',
    'ElementStart',
    'ElementEnd',
    'GeneSymbol',
    'Score'
)]
colnames(e2g_preds) <- c('chr', 'start', 'end', 'TargetGene', 'Score')

# create lists of credible set variants per trait
credible_set_variant_lists <- list()
traits <- list.files('resources/CredibleSets/')
for (i in 1:length(traits)) {
    credible_set_variants <- read.table(
        paste0("resources/CredibleSets/", traits[i], "/", "variant.list.txt"), 
        header = T
    )
    credible_set_variant_lists[[i]] <- credible_set_variants
}
names(credible_set_variant_lists) <- traits

# create list of all common variants across all chromosomes
common_variants_lists <- list()
for (chrom in 1:22) {
    chrom_common_variants = data.frame(fread(paste0(
        'resources/BIMS_hg38/1000G.EUR.QC.', chrom, '.bim'
    )))
    colnames(chrom_common_variants) = c("CHR", "SNP", "CM", "BP", "REF", "ALT")
    common_variants_lists[[chrom]] = chrom_common_variants
}
names(common_variants_lists) <- paste0('chr', 1:22)

# read in list of causal genes from PoPS
pops_genes <- data.frame(fread("resources/pops_genes.tsv"))
if (tissue == "BLD") {
    blood_traits <- c("Eosino", "Lym", "Mono", "Neutro", "Plt", "RBC", "MCV")
    pops_genes = pops_genes[which(pops_genes$Disease %in% blood_traits == T), ]
}

# create genomic ranges for enhancers
enhancer_granges = GRanges(
    seqnames = e2g_preds$chr,
    ranges = IRanges(start=e2g_preds$start, end = e2g_preds$end)
)

# add chrom, start, and end column to PoPS-gene credible sets
pops_genes$chr <- as.character(sapply(
    pops_genes$CredibleSet, function(x) return(strsplit(x, ":")[[1]][1])
))
pops_genes$cred_start <- as.numeric(sapply(
    pops_genes$CredibleSet, function(x) return(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])
))
pops_genes$cred_end <- as.numeric(sapply(
    pops_genes$CredibleSet, function(x) return(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2])
))

# get list of unique credible sets
unique_credible_sets <- unique(pops_genes$CredibleSet)

# create vectors to store output information
true_positives <- c()
positives <- c()
causals <- c()
creds_vec <- c()

# read in SNPs in promoters (for downstream filtering)
coding_promoter_snps = read.table("resources/Promoter_Coding_variants.txt", header=F)[,1]

for (i in 1:length(unique_credible_sets)) {
    
    # subset PoPS credible set-gene pairs to only include current credible set
    pops_cred_set_gene_pairs <- pops_genes[which(pops_genes$CredibleSet == unique_credible_sets[i]), ]

    # get causal genes from the credible set
    cred_set_causal_genes = pops_cred_set_gene_pairs$TargetGene[which(pops_cred_set_gene_pairs$truth == 1)]

    if (length(cred_set_causal_genes) > 1) {
        print(i)
    }

    # get unique diseases from the credible set
    disease_names = unique(pops_cred_set_gene_pairs$Disease)

    # select top N PoPS-implicated genes for credible set
    if (use_pops) {
      cred_set_pops_genes <- pops_cred_set_gene_pairs$TargetGene[order(pops_cred_set_gene_pairs$POPS.Rank, decreasing = F)[1:num_pops]]
    }

    # get credible set chromosome
    cred_set_chrom <- as.character(strsplit(unique_credible_sets[i], ":")[[1]][1])

    # create vectors to hold SNP IDs and PIPs
    rsids = c()
    pips = c()

    # get rsIDs and PIPs for all SNPs in current credible set for diseases (implicated by PoPS?)
    for (num_trait in 1:length(disease_names)) {
      idx =  which(credible_set_variant_lists[[disease_names[num_trait]]]$CredibleSet == unique_credible_sets[i])
      rsids = c(rsids, credible_set_variant_lists[[disease_names[num_trait]]]$rsid[idx])
      pips = c(pips, credible_set_variant_lists[[disease_names[num_trait]]]$pip[idx])
    }

    # filter for PIP > 0.1
    rsids <- rsids[which(pips > 0.10)]

    # filter SNP IDs that overlap promoters
    rsids <- setdiff(rsids, coding_promoter_snps)

    # exit for loop if there are no variants to look at
    if(length(intersect(rsids, common_variants_lists[[cred_set_chrom]]$SNP)) == 0){
      next
    }

    # get BP positions of rsIDs for this credible set
    rsids_bps <- common_variants_lists[[cred_set_chrom]]$BP[match(intersect(rsids, common_variants_lists[[cred_set_chrom]]$SNP), common_variants_lists[[cred_set_chrom]]$SNP)]
    rsids_starts <- rsids_bps - 1
    rsids_ends = rsids_bps + 1

    # create genomic ranges for credible set variants
    rsids_gr = GRanges(
        seqnames = cred_set_chrom,
        ranges = IRanges(start = rsids_starts, end = rsids_ends)
    )

    # overlap credible set variants with enhancer-gene predictios
    cred_set_eg_overlaps = findOverlaps(rsids_gr, enhancer_granges, type = "any", select = "all")

    if (length(cred_set_eg_overlaps) > 0) {

        # print(i)

        # get e2g preds genes and scores that overlap credible set variants
        cred_set_e2g_preds = e2g_preds[unique(subjectHits(cred_set_eg_overlaps)), c("TargetGene", "Score")]

        # take maximum score per gene
        cred_set_e2g_gene_scores = tapply(cred_set_e2g_preds$Score, cred_set_e2g_preds$TargetGene, max)

        # get gene scores for each gene in PoPS cred set-gene pairs
        pops_cred_set_gene_scores = rep(0, length(pops_cred_set_gene_pairs$TargetGene)) # all possible genes the credible set is evaluated against? what are the rows on the PoPS table?
        names(pops_cred_set_gene_scores) = pops_cred_set_gene_pairs$TargetGene
        pops_cred_set_gene_scores[match(intersect(names(cred_set_e2g_gene_scores), pops_cred_set_gene_pairs$TargetGene), pops_cred_set_gene_pairs$TargetGene)] = cred_set_e2g_gene_scores[match(intersect(names(cred_set_e2g_gene_scores), pops_cred_set_gene_pairs$TargetGene), names(cred_set_e2g_gene_scores))]

        # if there are no genes implicated by E2G linking scores
        if (max(pops_cred_set_gene_scores) == 0) {
            true_positives = c(true_positives, 0)
            positives = c(positives, 0)
        } 
        else {

            # determine top N genes implicated by E2G scores
            cred_set_pred_genes = names(pops_cred_set_gene_scores)[order(pops_cred_set_gene_scores, decreasing = T)[1:pmin(num_eg_genes, length(pops_cred_set_gene_scores))]]

            # intersect with PoPS (if specified)
            if (use_pops) {
                cred_set_pred_genes = intersect(cred_set_pred_genes, cred_set_pops_genes)
            }
            
            # record number of causal genes correctly predicted and number of genes predicted
            true_positives = c(true_positives, length(intersect(cred_set_pred_genes, cred_set_causal_genes)))
            positives = c(positives, length(cred_set_pred_genes))
        }
    }
    else {
        true_positives <- c(true_positives, 0)
        positives <- c(positives, 0)
    }

    # record numbers of causal genes per credible set
    causals = c(causals, length(cred_set_causal_genes))
}

print(true_positives)
print(positives)
print(causals)

# compute precision and recall of true causal genes
precision_length = sum(true_positives[which(positives !=  0)]) / length(positives[which(positives !=  0)])
recall_length = sum(true_positives[which(causals !=  0)]) / length(causals[which(causals !=  0)])

precision_sum = sum(true_positives[which(positives !=  0)]) / sum(positives[which(positives !=  0)])
recall_sum = sum(true_positives[which(causals !=  0)]) / sum(causals[which(causals !=  0)])

precision_length_boot = c()
precision_sum_boot = c()
recall_length_boot = c()
recall_sum_boot = c()

# bootstrap precision and recall measurements
for(nboot in 1:100){
   idx = sample(1:length(positives), length(positives), replace = T)
   true_positives_boot = true_positives[idx]
   positives_boot = positives[idx]
   relevant_boot = causals[idx]

   precision_length_boot = c(precision_length_boot, sum(true_positives_boot[which(positives_boot !=  0)]) / length(positives_boot[which(positives_boot !=  0)]))
   recall_length_boot = c(recall_length_boot, sum(true_positives_boot[which(relevant_boot !=  0)]) / length(relevant_boot[which(relevant_boot !=  0)]))

   precision_sum_boot = c(precision_sum_boot, sum(true_positives_boot[which(positives_boot !=  0)])/sum(positives_boot[which(positives_boot !=  0)]))
   recall_sum_boot = c(recall_sum_boot, sum(true_positives_boot[which(relevant_boot !=  0)])/sum(relevant_boot[which(relevant_boot !=  0)]))
}

# compute standard deviations for both types of precisions and recalls
sd_precision_length = sd(precision_length_boot)
sd_recall_length = sd(recall_length_boot)
sd_precision_sum = sd(precision_sum_boot)
sd_recall_sum = sd(recall_sum_boot)

# write output df
outdf <- data.frame(cbind(precision_length, sd_precision_length, recall_length, sd_recall_length, precision_sum, sd_precision_sum, recall_sum, sd_recall_sum))
colnames(outdf) <- c('Length.Precision', 'SD.Length.Precision', 'Length.Recall', 'SD.Length.Recall', 'Sum.Precision', 'SD.Sum.Precision', 'Sum.Recall', 'SD.Sum.Recall')
write.table(
    outdf,
    snakemake@output[[1]],
    col.names = T, 
    row.names = F, 
    sep = "\t", 
    quote=F
)

