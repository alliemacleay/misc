#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# changing author of this commit
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__= 'Allison MacLeay'

import sys
import os
import os.path
import argparse

#-----------------------------------------
#	MAIN
#-----------------------------------------
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Generate barcode files from one file")
	parser.add_argument('--file', required=True, help='file containing I5 and I7 barcodes.')
	parser.add_argument('--dir', default='touched', help='folder to write all files to.  Defaults to touched')
	parser.add_argument('--name', default="manflds", help='name of output.')
	args=parser.parse_args()
	if hasattr(args,'file'):
		file=args.file
	else:
		exit("No file provided")

	fh=open(file,'r')
	cmd= 'touch '+args.dir+ '/'
	for line in fh:
		os.system(cmd + line)

