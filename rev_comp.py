#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# changing author of this commit
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__= 'Allison MacLeay'

import sys
import os
import os.path
import argparse

def read_file():
	txt=''
	fh=open(file,'r')
	for line in fh:
		if line[0]!='>':
			txt+=line
	return txt

def get_rev_dna():
	return ["ACGTU","TGCAA"]

def get_rev_rna():
	return ["ACGTU","UGCAA"]

def seq_replace(seq,trans):
	seq_map = []
	seq_new = ''
	for i in range(len(seq)):
		pos=(trans[0].find(seq[i]))
		if pos<0:
			nt='N'
		else:
			nt=trans[1][pos]
		seq_new+=nt
	return seq_new	

def answer(txt,out):
	if out != 'stdout':
		sys.stdout=open(out,'w')
	print txt
	sys.stdout.close()
	return

#-----------------------------------------
#	MAIN
#-----------------------------------------
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="return the reverse compliment of a dna sequence")
	parser.add_argument('--input', help='sequence to generate reverse compliment')
	parser.add_argument('--file', help='file containing sequence input.')
	parser.add_argument('--out', default='stdout', help='output file')
	parser.add_argument('--rna', action='store_false', help='returns rna')
	args=parser.parse_args()
	input = ''
	if hasattr(args, 'input'):
		input=args.input
	elif hasattr(args,'file'):
		file=args.file
		input=read_file()
	else:
		help()
		exit("No input provided")
	trans_dict=[]
	if(args.rna):
		trans_dict=get_rev_rna()
	else:		
		trans_dict=get_rev_dna()
	input = seq_replace(input,trans_dict)
	output = input[::-1] # reverses string
	answer(output,args.out)

