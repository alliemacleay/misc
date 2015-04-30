#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# changing author of this commit
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__= 'Allison MacLeay'

import sys
import os
import os.path
import argparse
import glob
import gzip

def debug_here():
	# break at line 12 or 14
	return
def get_files(dir,name):
	flist=[]
	fpath=os.path.join(dir,name)
	fpath=glob.glob(fpath)
	for file in fpath:
		flist.append(file)	
	for f in next(os.walk(dir))[1]:
		folder=os.path.join(dir,f)
		fpath=os.path.join(folder,name)
		fpath=glob.glob(fpath)
		for file in fpath:
			flist.append(file)	
	return flist

def get_fid(fname):
	fid=''
	id=['','']
	for p in fname.split('_'):
		if p[0]=='A':
			id[0]=p
		elif p[0]=='P':
			id[1]=p
	fid='_'.join(id)
	return fid

def get_key(h):
	debug_here()
	info=h.split(':')
	return info[6]

def get_stats(d,flist):
	ct_same=0
	ct_diff=0
	ct_miss_left=0
	ct_miss_right=0

	if 'count_files' not in d.keys():
		d['count_files']=0	
	if 'totals' not in d.keys():
		d['totals']={}
	if 'pass' not in d.keys():
		d['pass']=0
	else:
		ct_miss_right=d['totals'][d['pass']]

	d['pass']=d['pass']+1
	d['totals'][d['pass']]=0

	d['fq']={}

	for file in flist:
		fh=''
		fid=get_fid(file.split('/')[-1])
		if fid not in d['fq'].keys():
				d['fq'][fid]={}
				d['fq'][fid]['fname'+str(d['pass'])]=file
				d['count_files']=d['count_files']+1
		if file.split('.')[-1]=='gz':
			fh=gzip.open(file,'rb')
		else:
			fh=open(file,'r')
		t=0
		ct_line=0
		entry=[]
		for line in fh:
			entry.append(line)
			t+=1
			if t==4:
				# parse data
				k=get_key(entry[0])
				d['totals'][d['pass']]=d['totals'][d['pass']]+1
				if (d['totals'][d['pass']] % 10000 ==0):
					print str(d['totals'][d['pass']]) +' reads in pass '+str(d['pass'])
				if d['pass']>1:
					if k in d['fq'][fid].keys():
						ct_miss_right-=1
						# diff
						cmp_res='same'
						diff=[]
						cmp=d['fq'][fid][k]['data']
						for l in range(len(cmp)):
							if cmp[l] != entry[l+1]:
								cmp_res='diff'
								diff.append(str(cmp[l]) + '<>' + str(entry[l+1]))
						d['fq'][fid][k]['diff']=cmp_res
						d['fq'][fid][k]['diffa']=diff
						del d['fq'][fid][k]['data']
					else:
						ct_miss_left+=1
						# add
						d['fq'][fid][k]={}
						d['fq'][fid][k]['data']=entry
						d['fq'][fid][k]['diff']='missing on left'
						d['fq'][fid][k]['right']=ct_line
				else:
					# add
					d['fq'][fid][k]={}
					d['fq'][fid][k]['data']=entry
					d['fq'][fid][k]['diff']='missing on right'
					d['fq'][fid][k]['left']=ct_line
				t=0
				entry=[]
					
	d['totals']['miss_l']=ct_miss_left
	d['totals']['miss_r']=ct_miss_right
	return d

def write_stats(d,out,quiet):
	if out != 'stdout':
		sys.stdout=open(out,'w')
	file=['f1','f2']
	
	print header_row()
	for fid in d['fq'].keys():
		if type(d['fq'][fid]) is not dict:
			continue #break
		for kval in d['bc'][fid].keys():
			if type(d['bc'][fid][kval]) is not dict:
				continue
			if 'diff' in d['fq'][fid][kval].keys():
					if d['fq'][fid][kval]['diff'] != 'same':
						if quiet & (d['fq'][fid][kval]['diff'].find('miss') > -1):
							pass
							# don't print
						else:
								#print diff_row_to_string(fid,kval,d['bc'][fid][kval])
								print row_to_string(fid,kval,d['fq'][fid][kval])
	print totals_to_string(d['totals'])
				
	sys.stdout.close()
	return

#-----------------------------------------
#	MAIN
#-----------------------------------------
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Compare unordered fastq files")
	parser.add_argument('--dir1', required=True, help='First directory containing fusion caller output. Could be substituted for a filename.')
	parser.add_argument('--dir2', required=True, help='Second directory containing fusion caller output.  Could be substituted for a filename.')
	parser.add_argument('--out', default="stdout", help='name of output.')
	parser.add_argument('--miss_off', action='store_true', help='skip missing file output.')
	args=parser.parse_args()

	stats={}
	stats['f1']={}
	stats['f2']={}
	stats['f1']['filename']=args.dir1
	stats['f2']['filename']=args.dir2
	if (os.path.isdir(args.dir1)) and (os.path.isdir(args.dir2)):
		files1=get_files(args.dir1,'*.fastq.gz')
		files2=get_files(args.dir2,'*.fastq.gz')
	elif (os.path.isfile(args.dir1)) and (os.path.isfile(args.dir1)):
		files1=[args.dir1]
		files2=[args.dir2]
	else:
		exit("Error with dir1 or dir2.  Both must be valid files or both must be directories")
	get_stats(stats,files1)
	get_stats(stats,files2)
	write_stats(stats,args.out,args.miss_off)
