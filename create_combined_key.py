#!/usr/bin/python

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create a combined barcode key file from 2 files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__= 'Allison MacLeay'

import sys
import os
import argparse

def create_key(f1,f2,out):
    bc_dict={}
    bc_dictA={}
    bc_dictP={}
    add_file_to_dict(f1,bc_dictA)
    add_file_to_dict(f2,bc_dictP)

    # create new file with from both inputs
    # -- this may not be necessary --
    fh=open(out,'w')
    for idA in bc_dictA.keys():
        for idP in bc_dictP.keys():
            id=idA+'_'+idP
            bc=bc_dictA[idA]+bc_dictP[idP]
            bc_dict[id]=bc
            print >> fh, id + "\t" + bc
    fh.close()
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

#-----------------------------------------
#	MAIN
# 
# 
#-----------------------------------------
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Create one barcode index file from 2")
	parser.add_argument('--file1', required='true', help='first barcode file')
	parser.add_argument('--file2', required='true', help='first barcode file')
	parser.add_argument('--out', default='combined_key.txt', help='full path to output file')
	args=parser.parse_args()
	create_key(args.file1,args.file2,args.out)
