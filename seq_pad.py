#!/usr/bin/python

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__= 'Allison MacLeay'

import sys
import os
import argparse
import glob
import gzip

#-----------------------------------------
#	MAIN
# run umitag.py for all files in a directory
# that have the same prefix
#-----------------------------------------
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Run umitag utility in batches of similarly prefixed names.")
	parser.add_argument('--dir', default='.', help='directory containing fastq output')
	parser.add_argument('--out', default='seq_pad_out', help='directory for output')
	parser.add_argument('--len', default=147, help='length to trim and pad to')
	args=parser.parse_args()

	l=int(args.len)
	names={}
	os.system("mkdir -p " + args.out)
	files = glob.glob(os.path.join(args.dir,"*.fastq.gz"))
	for f in files:
		pfx = '.'.join(f.split('.')[:-2]).split('/')[-1]
		if pfx in names.keys():
			sys.stderr.write( pfx + " is duplicated!\nExiting\n")
			sys.exit(1)
		names[pfx]=1
		fh = gzip.open(f,'r')
		out = gzip.open(os.path.join(args.out, pfx + "_padded.fastq.gz"),'wb')

		ct=0
		for line in fh:
			line = line.strip()
			ct+=1
			if ct%4 == 2:
				#sequence
				if len(line) < l:
					line = line + ('N'* (l-len(line)))
				#print line[:l]
				out.write(line[:l] + '\n')
			elif ct%4 == 0:
				#quality
				if len(line) < l:
					line = line + ('#'* (l-len(line)))
				#print line[:l]
				out.write(line[:l] + '\n')
			else:
				out.write(line + '\n')

		fh.close()
		out.close()
