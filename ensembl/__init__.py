"""
(c) MGH Center for Integrated Diagnostics

to execute run:
pyensembl install --release 75 --species homo_sapiens

"""

from __future__ import print_function
from __future__ import absolute_import
from pyensembl import EnsemblRelease

__author__ = 'Allison MacLeay'

class Ensembl(object):
    def __init__(self):
        self.db = EnsemblRelease(75)

    def annotate_one_gene(self, location):
        chrom, start, stop = self.parse_location(location)
        return self.db.gene_names_at_locus(chrom, start, stop)

    @staticmethod
    def parse_location(loc):
        start, stop = loc.split('-')
        chrom, start = start.split(':')
        start, stop = int(start), int(stop)
        return chrom, start, stop