import pandas as pd
import numpy as np
from pybedtools import BedTool
import gzip
import os
import sys

# parse input arguments to get input files, output file names, and BIM dir
input_bed = sys.argv[1]
bimfile_dir = sys.argv[2]
output_dir = sys.argv[3]

# parse arguments to make output file name
output_file = input_bed.split('/')[-1]
output_file = output_file.split('.')[0]

# create directory for module
module_name = output_file
output_dir = output_dir + module_name + '/'
os.makedirs(output_dir)
output_file = output_dir + output_file


def make_annot_files(bed_for_annot, bimfile, annot_file):
    df_bim = pd.read_csv(bimfile,
            delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    iter_bim = [['chr'+str(x1), x2, x2, 1] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    bimbed = BedTool(iter_bim)
    annotbed = bimbed.intersect(bed_for_annot, wb=True)
    bp = [x.start for x in annotbed]
    score = [float(x.fields[7]) for x in annotbed]
    df_int = pd.DataFrame({'BP': bp, 'ANNOT':score})
    df_annot = pd.merge(df_bim, df_int, how='left', on='BP')
    df_annot.fillna(0, inplace=True)
    temp = df_annot[['ANNOT']].astype(float)
    df_annot = pd.concat([df_bim.iloc[:,[0,3,1,2]], temp], axis = 1)
    if annot_file.endswith('.gz'):
        with gzip.open(annot_file, 'wb') as f:
            df_annot.to_csv(f, sep = "\t", index = False)
    else:
        df_annot.to_csv(annot_file, sep="\t", index=False)


for numchr in range(1, 23, 1):
    bimfile = bimfile_dir + "/" + "1000G.EUR.QC." + str(numchr) + ".bim"
    annot_file = output_file + "." + str(numchr) + ".annot.gz"
    bed_for_annot = BedTool(input_bed).sort()
    make_annot_files(bed_for_annot, bimfile, annot_file)
