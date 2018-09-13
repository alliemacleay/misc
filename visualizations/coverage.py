from __future__ import print_function
from __future__ import absolute_import
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
import matplotlib.collections as collections
import pandas as pd
from pyensembl import EnsemblRelease

__author__ = 'Allison MacLeay'


class FusionViz(object):
    def __init__(self, chrom, start, stop, hist):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.dict = self._to_dict(hist)  # hist is [x, y]
        self.hist = self.xy()

    def _to_dict(self, hist):
        """ pad coverage data to match length of region """
        dct = {}
        for i, pos in enumerate(hist[0]):
            dct[pos] = hist[1][i]
        return dct

    def xy(self):
        """ return array of buffered x and y"""
        x, y = [], []
        for i in range(self.start, self.stop):
            x.append(i)
            y.append(self.dict.get(i, 0))
        return [x, y]

    def plot(self, saveto=False):
        lhist = len(self.hist[0])

        plt.subplot(2, 1, 1)
        n, fracs, patches = plt.hist(self.invert(), bins=lhist/10)
        plt.title('Coverage')
        plt.ylabel('colored histogram')

        # We'll color code by height, but you could use any scalar
        fracs = n / n.max()

        # we need to normalize the data to 0..1 for the full range of the colormap
        norm = colors.Normalize(fracs.min(), fracs.max())

        # Now, we'll loop through our objects and set the color of each accordingly
        for thisfrac, thispatch in zip(fracs, patches):
            color = plt.cm.viridis(norm(thisfrac))
            thispatch.set_facecolor(color)
        barplot_sp = plt.subplot(2, 1, 2)
        _ = plt.bar(self.hist[0], self.hist[1])
        exon_array = self.parse_ref_exons()
        coll = collections.BrokenBarHCollection.span_where(self.hist[0], ymin=0, ymax=350, where=exon_array > 0, color='g')
        barplot_sp.add_collection(coll)
        plt.xticks(np.arange(self.start, self.stop, 1000))
        plt.xlabel('position')
        plt.ylabel('coverage')
        plt.show()
        if saveto:
            plt.savefig(saveto)
            print('Plot saved to {}'.format(saveto))

    def invert(self):
        iv = []
        for i, a in enumerate(self.hist[0]):
            for x in range(self.hist[1][i]):
                iv.append(a + self.start)
        return iv

    def d3(self):
        pass

    def load_ensembl_ref(self, ens_db, version, rid=None):
        """ Download, load, and index ensembl data """
        ens_db.download(version)
        ens_db.index()
        if rid is not None:
            return ens_db.transcript_by_id(rid)
        else:
            return None

    def get_exon_numbers(self, ens_db, gene):
        dct = {'start': [], 'id': []}
        e = ens_db.exon_ids_of_gene_name(gene)
        for exid in e:
            exon = ens_db.exon_by_id(exid)
            dct['start'].append(exon.start)
            dct['id'].append(exid)
        df = pd.DataFrame(dct)
        df = df.sort('start')
        df['number'] = range(1, df.shape[0] + 1)
        return df

    def parse_ref_exons(self):
        """ Return fasta reference with only the sequences needed"""
        ens_db = EnsemblRelease(75)
        try:
            exons = ens_db.exons_at_locus(self.chrom, self.start, self.stop)
        except ValueError as e:
            # Load pyensembl db
            raise e
        exon_array = np.zeros(self.stop - self.start)
        exon_numbers = self.get_exon_numbers(ens_db, exons[0].gene_name)
        for exobj in exons:
            ex_number = exon_numbers[exon_numbers['id'] == exobj.id]['number'].values[0]
            start = exobj.start - self.start
            stop = exobj.end - self.start
            i = start
            while i < stop:
                exon_array[i] = 1
                i += 1
        # 2:29,448,326-29,448,432 exon 19
        return exon_array




