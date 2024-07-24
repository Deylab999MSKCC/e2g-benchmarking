# E2G Benchmarking Pipeline

This repository contains a pipeline to benchmark various enhancer-gene (E2G) predictions using fine-mapped GWAS variants.

## Setup

There are several required input files to run this pipeline. [ADD DATA DOWNLOAD SCRIPT ONCE DATA IS AVAILABLE]

1. Gene scores: [INSERT DESCRIPTION OF GENE SCORES HERE]. These scores file needs to be downloaded to `resources/`. [INSERT A LOCATION TO DOWNLOAD GENE SCORES]
2. E2G predictions: There are six types of E2G predictions that are benchmarked (ABC_DNase_only, ABC_ext, ENCODE_E2G_DNase_HiC, ENCODE_E2G_DNase_only, ENCODE_E2G_ext, scE2G). Each of these models requires a subfolder within the `resources/` folder. For ABC_DNase_only, ENCODE_E2G_DNase_HiC, and ENCODE_E2G_DNase_only, the predictions for each file need to be placed in a subfolder called `thresholded_predictions`. One example path for these files would be: `resources/ABC_DNase_only/thresholded_predictions`. [INSERT INSTRUCTIONS FOR OTHER FILES LATER]
3. Common variants: These are common variants for each chromosome from the 1000 Genomes Project. These files need to be downloaded to `resources/BIMS_hg38`.
4. Finemapped variants: These are lists of fine-mapped variants per trait at different PIP thresholds. These files need to be downloaded to `resources/finemapped_variant_lists`.

## Step 1: Build E2G "modules"

The first step is to build E2G "modules" for each model type. The modules are BED files with the chromosome, start, and end of the enhancer, and the gene score for the enhancer's target gene. There are separate modules for each tissue. The tissue will determine the cell lines from which enhancer-gene links will be derived. 

## Step 2: Clean BED files

One caveat from the previous step is that there may be duplicate/overlapping enhancer regions since the same enhancer may be present in multiple cell lines or for multiple genes. The next step is to sort and merge the enhancer-gene predictions. The maximum value for the score column is taken. The columns of the BED files will stay the same.

## Step 3: Annotate variants

The next step is to determine whether common variants overlap with enhancers from the E2G links. The columns of the .annot.gz files are the chromosome, base pair (position), SNP ID, position in centimorgans, and whether the variant overlaps an ENCODE-E2G link (the value transferred here is the score from the E2G modules).

## Step 4: Overlap E2G predictions with fine-mapped GWAS variants

## Step 5: Compare E2G, PoPS, and E2G + PoPS

