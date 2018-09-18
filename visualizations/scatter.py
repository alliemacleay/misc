"""
(c) MGH Center for Integrated Diagnostics
"""

from __future__ import print_function
from __future__ import absolute_import
import argparse
import matplotlib.pyplot as plt
import pandas as pd


__author__ = 'Allison MacLeay'


def plot_data(ifile, x, y, label_col):
    df = pd.read_csv(ifile, sep='\t')
    plt.scatter(df[x], df[y])
    plt.show()
    pass

if __name__ == '__main__':
    args = argparse.ArgumentParser()
    args.add_argument('--input', required=True)
    args.add_argument('--x', default='background_alk_average_coverage')
    args.add_argument('--y', default='ati_average_coverage')
    args.add_argument('--labels', default='sample')
    nargs = args.parse_args()
    plot_data(nargs.input, nargs.x, nargs.y, nargs.labels)