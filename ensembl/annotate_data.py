"""
(c) MGH Center for Integrated Diagnostics
"""

from __future__ import print_function
from __future__ import absolute_import
import sys
import os
import pandas as pd
from ensembl import Ensembl

__author__ = 'Allison MacLeay'

def create_loc(chrom, pos, read_pos):
    pos = int(pos) + int(read_pos)
    return '{}:{}-{}'.format(chrom, pos, int(pos) + 1)

def annotate_single(inf, ouf):
    df = pd.read_csv(inf, sep='\t')
    genes = []
    ens = Ensembl()
    for i, row in df.iterrows():
        try:
            gene = ens.annotate_one_gene(create_loc(row['chrom'], row['pos'], row['read_pos']))[0]
        except IndexError:
            gene = 'N/A'
        genes.append(gene)
    df['gene'] = genes
    df.to_csv(ouf, sep='\t', index=False)

if __name__ == '__main__':
    dfile = sys.argv[1]
    ofile = os.path.join(os.path.dirname(dfile), '{}_annotated.txt'.format(os.path.basename(dfile).split('.')[0]))
    annotate_single(dfile, ofile)