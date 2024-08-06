# This file annotates common variants from the 1KGP with whether or not they
# overlap an E2G link.
#
# Author: Karthik Guruvayurappan

import numpy as np
import pandas as pd
from pybedtools import BedTool
import gzip

ENHANCERS_FILE = snakemake.input[0]
BIMFILE = snakemake.input[1]
ANNOT_FILE = snakemake.output[0]

# read in unique enhancer regions from E2G predictions
enhancers = BedTool(ENHANCERS_FILE).sort()

# read in BIM file containing common variants
df_bim = pd.read_csv(
    BIMFILE,
    sep = '\t', 
    usecols = [0,1,2,3], 
    names = ['CHR','SNP','CM','BP']
)

# convert BIM file into BED file
bimbed = [['chr'+str(x1), x2, x2, 1] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
bimbed = BedTool(bimbed)

# intersect common variants with E2G annotations
overlap_variants = bimbed.intersect(enhancers, wa = True)
overlap_pos = [variant.start for variant in overlap_variants]

# create overlaps column for all common variants
df_bim['overlaps_enhancer'] = df_bim['BP'].isin(overlap_pos)

# process and write to output file
df_bim['overlaps_enhancer'] = df_bim['overlaps_enhancer'].astype(int)
df_bim = df_bim[[
    'CHR',
    'BP',
    'SNP',
    'CM',
    'overlaps_enhancer'
]]

df_bim.columns = [
    'CHR',
    'BP',
    'SNP',
    'CM',
    'ANNOT'
]

with gzip.open(ANNOT_FILE, 'wb') as f:
    df_bim.to_csv(f, sep = "\t", index = False)

