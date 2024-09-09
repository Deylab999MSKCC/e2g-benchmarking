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

## Step 4: Compute enrichment, precision, and recall of finemapped variants within enhancer regions

The next step is to compute enrichment, precision, and recall values for the presence of fine-mapped variants within enhancer regions. The definitions are as follows:
1. Enrichment: (# finemapped variants overlapping enhancers / # common variants overlapping enhancers) / (# finemapped variants / # common variants). In other words, do the enhancer annotations contain a greater proportion of finemapped variants than the genome?
2. Precision: (# finemapped variants overlapping enhancers / # common variants overlapping enhancers). In other words, what fraction of the variants overlapping enhancers are also finemapped for a GWAS trait?
3. Recall: (# finemapped variants overlapping enhancers / # finemapped variants). In other words, what fraction of finemapped variants are recovered by the enhancer annotations?

This can be done for different sets of traits. In this pipeline, we have currently implemented this for K562-related traits (red blood cell-related traits) and GM1287-related traits (white blood cell-related). For each of these, there are different lists of finemapped variants per trait, and the enrichment, precision, and recall are computed on the level of each list. For each value, standard deviations are computed by performing 100 bootstrap iterations, bootstrapping across chromosomes. P-values for each are computed by using permutation testing, shuffling the overlaps between common variants and enhancer regions. The 100 permutations are used to compute a normal null distribution, from which p-values are computed.

## Step 5: Run credible set-gene benchmarking

Another way to benchmark enhancer-gene predictions using GWAS is to benchmark how well enhancer-gene predictions can link GWAS credible sets to their target genes. In the ENCODE-E2G paper, one approach presented to do this is to intersect predictions from the ENCODE-rE2G enhancer-gene model with PoPS scores. PoPS scores specify whether a gene is likely to be relevant for a disease. One way to link credible sets to genes with PoPS would then be to take the gene with the highest PoPS score in the vicinity (within 2 MB) of a credible set. Here, be default, we use an approach that intersects the top E2G gene with the top 2 PoPS genes to establish credible set-gene links. We then compute the precision and recall against a silver-standard dataset curated by the ENCODE consortium (see ENCODE-rE2G paper).

# To reproduce ENCODE-rE2G GWAS Benchmarking

## Step 1: Download ENCODE-rE2G data from Synapse.

All of the enhancer-gene predictions from various methods are avaiable on [Synapse](https://www.synapse.org/Synapse:syn52234275/files/) using the link provided in the paper. First, this data needs to be downloaded from Synapse to the `resources/` folder.

## Step 2:







