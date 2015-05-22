#!/usr/bin/python

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__= 'Allison MacLeay'

import sys
import os
import argparse
import gzip

def read_stats(s,l):
	rlen=len(l)
	s['read_count']+=1
	if rlen not in s['length'].keys():
		s['length'][rlen]=0
	Nct=l.count('N')
	Npt=Nct/float(rlen)*100
	if Npt > s['max_Ns']:
		s['max_Ns']=Npt
	s['N_count']+=Nct
	s['bases']+=rlen
	s['length'][rlen]+=1
	return
def qual_stats(s,l):
	r=range(len(l))
	t=0
	i=0
	for i in r:
		t=t+ord(l[i])-33
	i+=1
	ave_float = t/float(i)
	average = float("{0:.2f}".format(ave_float))
	if average not in s['average'].keys():
		s['average'][average]=0
	s['average'][average]+=1
	s['ave_total']+=ave_float
	if ave_float<s['ave_min']:
		s['ave_min']=ave_float
	if ave_float>s['ave_max']:
		s['ave_max']=ave_float
	return
def bin_qual_scores(s):
	BIN_NUM=10
	brange=(s['ave_max']-s['ave_min'])/BIN_NUM
	low=s['ave_min']
	s['ave_bin']={}
	for i in range(BIN_NUM):
		s['ave_bin'][low]=0
		low=low+brange
	for k in s['average'].keys():
		for low in s['ave_bin'].keys():
			if abs(k-low)<brange:
				s['ave_bin'][low]+=s['average'][k]
	return
def print_stats(s):
	total = s['read_count']
	last = 0
	keyarr=['bases','read_count','N_count']
	#keyarr.extend(['max_Ns','length','ave_max','ave_min','ave_total','ave_bin'])]
	for i in keyarr:
		print str(i) + ' ' + str(s[i])
	print "\nLength Distribution"
	for k in sorted(s['length'].keys()):
		print "{0}: {1}\t{2:.2f}%".format(k,s['length'][k],s['length'][k]/float(total)*100)
	print "\nQuality Score Distribution"
	print "Max qscore:{0}\nMin qscore:{1}\nMean qscore:{2}".format(s['ave_max'],s['ave_min'],s['ave_total']/float(total))
	for k in sorted(s['ave_bin'].keys()):
		if last != 0:
			print "{0:.2f} to {1:.2f}: {2:.2f}%".format(last,k,s['ave_bin'][k]/float(total)*100)
		last=k
	print "\nMaximum percentage of Ns in one read: {0:.2f}%".format(s['max_Ns'])
	print "Percent of Ns in reads: {0:.2f}%".format(s['N_count']/float(s['bases'])*100)
	return
#-----------------------------------------
#	MAIN
# run umitag.py for all files in a directory
# that have the same prefix
#-----------------------------------------
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Compute stats about a fastq file")
	parser.add_argument('--file', required=True, help='fastq file to analyze')
	args=parser.parse_args()

	fh=''
	ext=args.file.split('.')[-1]
	if ext == 'gz':
		fh=gzip.open(args.file,'r')
	else:
		fh=open(args.file,'r')
	stats={}
	stats['read_count']=0
	stats['length']={}
	stats['average']={}
	stats['ave_min']=100
	stats['ave_max']=0
	stats['ave_total']=0
	stats['N_count']=0
	stats['max_Ns']=0
	stats['bases']=0
	ct=0
	print "Processing "+args.file
	for line in fh:
		ct+=1
		if ct%4 == 2:
			read_stats(stats,line.strip())
		if ct%4 == 0:
			qual_stats(stats,line.strip())

	fh.close()
	
	stats['mean_qual']=stats['ave_total']/ct
	bin_qual_scores(stats)
	stats['lines']=ct
	print_stats(stats)
		
