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
parser.add_argument('--p5_barcodes')
parser.add_argument('--p7_barcodes')
parser.add_argument('--out_dir', default='.')
args = vars(parser.parse_args())
out_dir = args['out_dir']

#args = {'out_dir':'/PHShome/ma695/tmp', 'min_reads':10}
#base = '/data/joung/sequencing_bcl/131007_M01326_0075_000000000-A6B33/Data/Intensities/BaseCalls'
#args['read1'] = os.path.join(base, 'Undetermined_S0_L001_R1_001.fastq.gz')
#args['read2'] = os.path.join(base, 'Undetermined_S0_L001_R2_001.fastq.gz')
#args['index1'] = os.path.join(base, 'Undetermined_S0_L001_I1_001.fastq.gz')
#args['index2'] = os.path.join(base, 'Undetermined_S0_L001_I2_001.fastq.gz')

swap={}
do_swap=1
fargs=['read1','read2','index1','index2']
for f in fargs:
    name=args[f]
    if name.find('_R1_')>0:
        swap['read1']=name
    elif name.find('_R2_')>0:
        swap['read2']=name
    elif name.find('_I1_')>0:
        swap['index1']=name
    elif name.find('_I2_')>0:
        swap['index2']=name
    else:
        do_swap=0 # one or more files do not adhere to schema.  Can not confidently swap

if (do_swap==1) & (len(swap)==4):
    # swapping files for names that are passed in in the wrong order
    # but follow the schema of containing _R1_
    for f in fargs:
        args[f]=swap[f]

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a combined barcode key file from 2 files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def create_key(f1,f2):
    bc_dict={}
    bc_dictA={}
    bc_dictP={}
    add_file_to_dict(f1,bc_dictA)
    add_file_to_dict(f2,bc_dictP)

    for idA in bc_dictA.keys():
        for idP in bc_dictP.keys():
            id=idA+'_'+idP
            bc=bc_dictA[idA][1:]+bc_dictP[idP][1:]
            bc_dict[id]=bc
    return bc_dict

def add_file_to_dict(fname,d):
    """
    helper function - add a file to the dictionary
    :param fname: an array of strings with filenames to parse
    :param d: a dictionary (dict)
    :return:
    """
    HEADER=1 # skip first line if equal to 1
    fh=open(fname,'r')
    for line in fh:
        if HEADER == 1:
            HEADER=HEADER-1
            continue
        [id,seq]=line.strip().split('\t')
        d[id]=seq
    fh.close()
    return

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# helper functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#          MAIN
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# parse barcode file to tag with A and P ids
sample_names = {}
if not args['sample_barcodes']==None:
    for line in open(args['sample_barcodes'], 'r'):
        fields = line.strip().split('\t')
        if len(fields)==2:
            sampleid, barcode = fields
            sample_names[barcode] = sampleid

if ('p5_barcodes' in args.keys()) & ('p7_barcodes' in args.keys()):
    sample_names = create_key(args['p5_barcodes'], args['p7_barcodes'])

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

print ("Read count complete in %.1f minutes." % ( (time.time()-start)/60))

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
	if 'undetermined_r1' not in vars():
		undetermined_r1 = open(os.path.join(out_dir, 'undetermined.r1.fastq'), 'w')
	if 'undetermined_r2' not in vars():
		undetermined_r2 = open(os.path.join(out_dir, 'undetermined.r2.fastq'), 'w')
	if 'undetermined_i1' not in vars():
		undetermined_i1 = open(os.path.join(out_dir, 'undetermined.i1.fastq'), 'w')
	if 'undetermined_i2' not in vars():
		undetermined_i2 = open(os.path.join(out_dir, 'undetermined.i2.fastq'), 'w')
        for line in r1:
            print (line, file=undetermined_r1, end="")
        for line in r2:
            print (line, file=undetermined_r2, end="")
        for line in i1:
            print (line, file=undetermined_i1, end="")
        for line in i2:
            print (line, file=undetermined_i2, end="")
    else:
	if sample_id not in outfiles_r1.keys():
		outfiles_r1[sample_id] = open(os.path.join(out_dir, '%s.r1.fastq' % sample_id), 'w')
		outfiles_r2[sample_id] = open(os.path.join(out_dir, '%s.r2.fastq' % sample_id), 'w')
		outfiles_i1[sample_id] = open(os.path.join(out_dir, '%s.i1.fastq' % sample_id), 'w')
		outfiles_i2[sample_id] = open(os.path.join(out_dir, '%s.i2.fastq' % sample_id), 'w')
        for line in r1:
            print (line, file=outfiles_r1[sample_id], end="")
        for line in r2:
            print (line, file=outfiles_r2[sample_id], end="")
        for line in i1:
            print (line, file=outfiles_i1[sample_id], end="")
        for line in i2:
            print (line, file=outfiles_i2[sample_id], end="")

undetermined_r1.close()
undetermined_r2.close()
undetermined_i1.close()
undetermined_i2.close()

for sample_id in outfiles_r1.keys():
	outfiles_r1[sample_id].close()
for sample_id in outfiles_r2.keys():
	outfiles_r2[sample_id].close()
for sample_id in outfiles_i1.keys():
	outfiles_i1[sample_id].close()
for sample_id in outfiles_i2.keys():
	outfiles_i2[sample_id].close()

num_fastqs = len([v for k,v in count.iteritems() if v>=args['min_reads']])
print('Wrote FASTQs for the %d sample barcodes out of %d with at least %d reads in %.1f minutes.' % (num_fastqs, len(count), args['min_reads'], (time.time()-start)/60 ))
