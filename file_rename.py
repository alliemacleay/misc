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
	parser=argparse.ArgumentParser(description="Rename files based on barcode ids from an input file")
	parser.add_argument('--string', required=True, help='Text to add to filename')
	parser.add_argument('--dir', required=True, help='directory containing files to rename.')
	parser.add_argument('--start', default=0, help='index of start position for barcode search. first position is 0')
	args=parser.parse_args()
	data=False
	if hasattr(args,'string'):
		txt=args.string
	else:
		exit("No text provided")
	files=next(os.walk(args.dir))[2]
	ct=0
	ct_renamed=0
	for fqfile in files:
		ct+=1
		fnew=fqfile[:int(args.start)] + txt + fqfile[int(args.start):]
		print "rename " + fqfile + " to " + fnew
		ct_renamed+=1
		os.rename(os.path.join(args.dir,fqfile),os.path.join(args.dir,fnew))
	print "Renamed " + str(ct_renamed) + " of " + str(ct) + " files found in " + args.dir
