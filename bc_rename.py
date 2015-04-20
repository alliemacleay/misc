#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# changing author of this commit
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__= 'Allison MacLeay'

import sys
import os
import os.path
import argparse


def add_file_to_dict(fname,length,d):
    """
    helper function - add a file to the dictionary
    :param fname: an array of strings with filenames to parse
    :param length: the length of the sequence index. first 
    :    	   characters are dropped if the length is smaller
    :param d: a dictionary (dict)
    :return:
    """
    HEADER=1 # skip first line if equal to 1
    fh=open(fname,'r')
    for line in fh:
        if HEADER == 1:
            HEADER=HEADER-1
            continue
	
        line=line.strip().split('\t')
	if len(line)<2:
	    continue
        [id,seq]=line
        d[seq[-length:]]=id
    fh.close()
    return
#-----------------------------------------
#	MAIN
#-----------------------------------------
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Rename files based on barcode ids from an input file")
	parser.add_argument('--file', required=True, help='file containing I5 and I7 barcodes.')
	parser.add_argument('--dir', required=True, help='directory containing files to rename.')
	parser.add_argument('--start', default=0, help='index of start position for barcode search. first position is 0')
	parser.add_argument('--len', default=7, help='length of barcode.  Default is 7.')
	args=parser.parse_args()
	data=False
	dbc={}
	dbc['header']=''
	if hasattr(args,'file'):
		file=args.file
	else:
		exit("No file provided")
	barcodes={}
	add_file_to_dict(args.file,args.len,barcodes)
	print barcodes
	files=next(os.walk(args.dir))[2]
	ct=0
	ct_renamed=0
	for fqfile in files:
		ct+=1
		if len(fqfile)>args.start+args.len:
			seq=fqfile[args.start:args.start+args.len]
			if seq in barcodes.keys():
				fnew=fqfile.replace(seq,barcodes[seq])
				print "rename " + fqfile + " to " + fnew
				ct_renamed+=1
				os.rename(os.path.join(args.dir,fqfile),os.path.join(args.dir,fnew))
	print "Renamed " + str(ct_renamed) + " of " + str(ct) + " files found in " + args.dir
