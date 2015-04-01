#!/usr/bin/python

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Saving unused lego
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__= 'Allison MacLeay'

import sys
import os
import argparse


class umitag(Tool):
    """
    : Unused - must be done in batch
    """
    name = "Umi tag"
    inputs = ['umi_fastq_dir']
    outputs = ['1.fastq','2.fastq']
    mem_req = 6 * 1024
    cpu_req = 4
    time_req = 12 * 60

    #python ../umitag.py --read1_in mysample.r1.fastq --read2_in mysample.r2.fastq --read1_out mysample.r1.umitagged.fastq \
    # --read2_out mysample.r2.umitagged.fastq --index1 mysample.i1.fastq --index2 mysample.i2.fastq

    def cmd(self, i, s, p):
        return """
        python {s[umitag]} --read1_in {fastq_r1}
        --read2_in {fastq_r2}
        --read1_out $OUT.1.fastq
        --read2_out $OUT.2.fastq
        --index1 {fastq_i1}
        --index2 {fastq_i2}
        """, \
            {
                'fastq_r1': os.path.join(str(i['umi_fastq_dir'][0]), p['barcode_index'] + '.r1.fastq'),
                'fastq_r2': os.path.join(str(i['umi_fastq_dir'][0]), p['barcode_index'] + '.r2.fastq'),
                'fastq_i1': os.path.join(str(i['umi_fastq_dir'][0]), p['barcode_index'] + '.i1.fastq'),
                'fastq_i2': os.path.join(str(i['umi_fastq_dir'][0]), p['barcode_index'] + '.i2.fastq'),
            }

class get_file_names(Tool):
    name = "Get file names"
    inputs = ['umi_fastq_dir']
    outputs = ['txt']
    mem_req = 6 * 1024
    cpu_req = 4
    time_req = 12 * 60

    def cmd(self,i,s,p):
        return """
        {s[get_files]} {input_dir} >> $OUT.txt
        """, \
            {
                'input_dir':    i['umi_fastq_dir'][0],
            }

#-----------------------------------------
#	MAIN
# run umitag.py for all files in a directory
# that have the same prefix
#-----------------------------------------
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Run umitag utility in batches of similarly prefixed names.")
	parser.add_argument('--dir', default='.', help='directory containing output of umi demultiplex')
	parser.parse_args()
