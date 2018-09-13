"""
(c) MGH Center for Integrated Diagnostics
"""

from __future__ import print_function
from __future__ import absolute_import
import argparse
import matplotlib.pyplot as plt
import pandas as pd

__author__ = 'Allison MacLeay'


def plot_data(ifile, data_col, label_col, sort=True):
    df = pd.read_csv(ifile, sep='\t')
    if sort:
        df = df.sort(data_col, ascending=False)
    plt.bar(range(df.shape[0]), df[data_col])
    plt.show()
    pass

if __name__ == '__main__':
    args = argparse.ArgumentParser()
    args.add_argument('--input', required=True)
    args.add_argument('--data', default='average_coverage')
    args.add_argument('--labels', default='sample')
    nargs = args.parse_args()
    plot_data(nargs.input, nargs.data, nargs.labels)