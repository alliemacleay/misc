from __future__ import print_function
import os
import re
import gzip
import itertools
import argparse
import time

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# Edited Martin Aryee's function to demultiplex to reduce
# the amount of memory required to run.
# Martin's original demultiplex function can be found
# here:
# https://github.com/aryeelab/umi/wiki
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__ = 'Allison MacLeay'


parser = argparse.ArgumentParser()
parser.add_argument('--read1', required=True)
parser.add_argument('--read2', required=True)
parser.add_argument('--index1', required=True)
parser.add_argument('--index2', required=True)
parser.add_argument('--min_reads', type=int, default=10000)
parser.add_argument('--sample_barcodes')
parser.add_argument('--out_dir', default='.')
args = vars(parser.parse_args())
out_dir = args['out_dir']

#args = {'out_dir':'/PHShome/ma695/tmp', 'min_reads':10}
#base = '/data/joung/sequencing_bcl/131007_M01326_0075_000000000-A6B33/Data/Intensities/BaseCalls'
#args['read1'] = os.path.join(base, 'Undetermined_S0_L001_R1_001.fastq.gz')
#args['read2'] = os.path.join(base, 'Undetermined_S0_L001_R2_001.fastq.gz')
#args['index1'] = os.path.join(base, 'Undetermined_S0_L001_I1_001.fastq.gz')
#args['index2'] = os.path.join(base, 'Undetermined_S0_L001_I2_001.fastq.gz')

def fq(file):
    if re.search('.gz$', file):
        fastq = gzip.open(file, 'rb')
    else:
        fastq = open(file, 'r')
    with fastq as f:
        while True:
            l1 = f.readline()
            if not l1:
                break
            l2 = f.readline()
            l3 = f.readline()
            l4 = f.readline()
            yield [l1, l2, l3, l4]

def get_sample_id(i1, i2):
    sample_barcode = get_seq(i1, i2)
    if sample_names.has_key(sample_barcode):
        return sample_names[sample_barcode]
    else:
        return sample_barcode

def get_seq(i1, i2):
    seq1=i1[1]
    seq2=i2[1]
    return seq1[1:8] + seq2[1:8]

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

sample_names = {}
if not args['sample_barcodes']==None:
    for line in open(args['sample_barcodes'], 'r'):
        fields = line.strip().split('\t')
        if len(fields)==2:
            sampleid, barcode = fields
            sample_names[barcode] = sampleid

outfiles_r1 = {}
outfiles_r2 = {}
outfiles_i1 = {}
outfiles_i2 = {}

total_count = 0
count = {}

# Create count dictionary first
start = time.time()
for i1,i2 in itertools.izip(fq(args['index1']), fq(args['index2'])):
    sample_id = get_sample_id(i1, i2)

    # Increment read count and create output buffers if this is a new sample barcode
    if not count.has_key(sample_id):
        count[sample_id] = 0
    count[sample_id] += 1

for r1,r2,i1,i2 in itertools.izip(fq(args['read1']), fq(args['read2']), fq(args['index1']), fq(args['index2'])):
    # the original demultiplex stored sequences in a buffer to execute in 1N instead of 2N
    # this version minimizes the memory requirement by running in 2N
    total_count += 1
    if total_count % 1000000 == 0:
        print ("Processed %d reads in %.1f minutes." % (total_count, (time.time()-start)/60))
    sample_id=get_sample_id(i1,i2)
    if count[sample_id] < args['min_reads']:
	# Write remaining buffered reads to a single fastq.
	# (These reads correspond to barcodes that were seen less than min_reads times)
        undetermined_r1 = open(os.path.join(out_dir, 'undetermined.r1.fastq'), 'a')
        undetermined_r2 = open(os.path.join(out_dir, 'undetermined.r2.fastq'), 'a')
        undetermined_i1 = open(os.path.join(out_dir, 'undetermined.i1.fastq'), 'a')
        undetermined_i2 = open(os.path.join(out_dir, 'undetermined.i2.fastq'), 'a')
        for line in r1:
            print (line, file=undetermined_r1, end="")
        for line in r2:
            print (line, file=undetermined_r2, end="")
        for line in i1:
            print (line, file=undetermined_i1, end="")
        for line in i2:
            print (line, file=undetermined_i2, end="")
        undetermined_r1.close()
        undetermined_r2.close()
        undetermined_i1.close()
        undetermined_i2.close()
    else:
        outfiles_r1[sample_id] = open(os.path.join(out_dir, '%s.r1.fastq' % sample_id), 'a')
        outfiles_r2[sample_id] = open(os.path.join(out_dir, '%s.r2.fastq' % sample_id), 'a')
        outfiles_i1[sample_id] = open(os.path.join(out_dir, '%s.i1.fastq' % sample_id), 'a')
        outfiles_i2[sample_id] = open(os.path.join(out_dir, '%s.i2.fastq' % sample_id), 'a')
        for line in r1:
            print (line, file=outfiles_r1[sample_id], end="")
        for line in r2:
            print (line, file=outfiles_r2[sample_id], end="")
        for line in i1:
            print (line, file=outfiles_i1[sample_id], end="")
        for line in i2:
            print (line, file=outfiles_i2[sample_id], end="")
        outfiles_r1[sample_id].close()
        outfiles_r2[sample_id].close()
        outfiles_i1[sample_id].close()
        outfiles_i2[sample_id].close()

num_fastqs = len([v for k,v in count.iteritems() if v>=args['min_reads']])
print('Wrote FASTQs for the %d sample barcodes out of %d with at least %d reads.' % (num_fastqs, len(count), args['min_reads']))
