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

## Step 4: 



