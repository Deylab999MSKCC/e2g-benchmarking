# E2G GWAS Benchmarking

This repository contains a snakemake pipeline to benchmark enhancer-gene (E2G) predictions using information from genome-wide association studies (GWAS).

There are two benchmarking methods:
1. Overlapping fine-mapped GWAS variants with predicted enhancers
2. Using E2G predictions to link GWAS credible sets to causal genes

This README will explain the steps required to run this benchmarking pipeline.

## Step 1: Prepare input enhancer-gene prediction file

As input, this pipeline takes enhancer-gene predictions generated from your method of choice. The required input file is a BED file containing the following columns: 
1. chromosome (ie: chr1, chr2, etc.)
2. start coordinate
3. end coordinate
4. gene name (ie: SAMD11, NOC2L, etc.)
5. score (how confident is the model in the enhancer-gene prediction)

## Step 2: Merge enhancer regions

For overlapping fine-mapped GWAS variants with predicted enhancers, we are only interested in the enhancer regions themselves, and not their gene links.
First, we sort the enhancer regions by genomic position. Then, since there may be overlapping enhancer regions (ie: a region linked to 2 genes, two slightly 
offset regions from different biosamples, etc.), we merge overlapping enhancer regions together.

## Step 3: Annotate common variants overlapping enhancers

The next step is to take every common variant in hg38 (from the 1000 Genomes Project), and mark whether that variant overlaps a predicted enhancer region.
For runtime, this is done separately for each chromosome.

## Step 4: Compute enrichment, precision, and recall of finemapped variants within enhancer regions.

The next step is to compute enrichment, precision, and recall values for the presence of fine-mapped variants within enhancer regions. The definitions are as follows:
1. Enrichment: (# finemapped variants overlapping enhancers / # common variants overlapping enhancers) / (# finemapped variants / # common variants). In other words, do the enhancer annotations contain a greater proportion of finemapped variants than the genome?
2. Precision: (# finemapped variants overlapping enhancers / # common variants overlapping enhancers). In other words, what fraction of the variants overlapping enhancers are also finemapped for a GWAS trait?
3. Recall: (# finemapped variants overlapping enhancers / # finemapped variants). In other words, what fraction of finemapped variants are recovered by the enhancer annotations?

This can be done for different sets of traits. In this pipeline, we have currently implemented this for K562-related traits (red blood cell-related traits) and GM1287-related traits (white blood cell-related). For each of these, there are different lists of finemapped variants per trait, and the enrichment, precision, and recall are computed on the level of each list. For each value, standard deviations are computed by performing 100 bootstrap iterations, bootstrapping across chromosomes. P-values for each are computed by using permutation testing, shuffling the overlaps between common variants and enhancer regions. The 100 permutations are used to compute a normal null distribution, from which p-values are computed.





