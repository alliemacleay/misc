#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# changing author of this commit
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__= 'Allison MacLeay'

import sys
import os
import os.path
import argparse
import glob

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

	if 'bc' not in d.keys():
		d['bc']={}
	columns=[1,2,5,6]
	key=0
	for file in flist:
		HEADER=True
		fid=get_fid(file.split('/')[-1])
		if fid not in d['bc'].keys():
				d['bc'][fid]={}
				d['bc'][fid]['fname'+str(d['pass'])]=file
				d['count_files']=d['count_files']+1
		fh=open(file,'r')
		for line in fh:
			t=0
			data=line.strip().split('\t')
			if HEADER:
				d['header']=data
				HEADER=False
				continue
			d['totals'][d['pass']]+=1
			kval=data[key]
			cmp=''
			if d['pass']>1:
				# compare
				cmp='same'
				if kval in d['bc'][fid].keys():
					if ct_miss_right > 0:
						ct_miss_right-=1
				else:
					d['bc'][fid][kval]={}
					cmp='missing on left'
					ct_miss_left+=1

			for t in columns:
				colname=d['header'][int(t)]

				#compare
				if d['pass']>1:
					if colname not in d['bc'][fid][kval].keys():
						if cmp=='same':
							cmp='missing column'
						#d['bc'][fid][kval][colname]=''
						d['bc'][fid][kval][colname]=data[t]

					# record new column values
					if d['bc'][fid][kval][colname] != data[t]:
						if 'diff' not in d['bc'][fid][kval].keys():
							d['bc'][fid][kval]['diff']={}
						if cmp=='same':
							cmp='diff'
						d['bc'][fid][kval]['diff'][colname]=str(d['bc'][fid][kval][colname]) + ' <> ' + str(data[t])
				else:
					# load first pass
					if kval not in d['bc'][fid].keys():
						d['bc'][fid][kval]={}
					d['bc'][fid][kval]['cmp']='missing on right'
					d['bc'][fid][kval][colname]=data[t]
			if cmp != '':
				d['bc'][fid][kval]['cmp']=cmp
				if cmp=='same':
					ct_same+=1
				elif cmp=='diff':
					ct_diff+=1
	if 'totals' not in d.keys():
		d['totals']={}
	if d['pass']>1:
		d['totals']['same']=ct_same
		d['totals']['diff']=ct_diff
		d['totals']['miss_l']=ct_miss_left
		d['totals']['miss_r']=ct_miss_right
	return

def diff_row_to_string(bc,fusion,row):
	txt=''
	vals=[]
	txt+=fusion
	cmp=row['cmp']
	txt+='\t'+cmp
	if cmp == 'diff':
			for k in row['diff'].keys():
				if type(row['diff'][k]) is str:
						vals.append(row[k])
			txt+='\t'.join(vals)
	
	return txt

def diff_to_string(dcol):
	txt=''
	for k in dcol.keys():
		dbp=''
		if (k.find('reakPoint')>-1) or (k.find('overage')>-1):
			bp=dcol[k].split(' <> ')
			if len(bp[0].split(':'))>1:
				for i in range(len(bp)):
					bp[i]=bp[i].split(':')[1]
			dbp=int(bp[1])-int(bp[0])
			dbp=' (change: ' + str(dbp) + ')'
		if txt!='':
			txt+='\t'
		txt+= '~~col: '+k + ' *** ' + dcol[k] +'' + dbp + '  '

	return txt

def row_to_string(bc,fusion,row):
	txt=''+bc
	vals=[]
	txt+='\t'+fusion
	for col in row.keys():
		if (type(row[col]) is not dict) and (col != 'diff'):
			txt+='\t'+str(row[col])
		if col=='diff':
			txt=txt+ '\tdifferent\t'+diff_to_string(row['diff'])
	return txt

def totals_to_string(totals):
	nm=''
	txt=''
	fn=['left','right']
	for c in totals.keys():
		if type(c) is int:
			nm=fn[c-1] + ' file'
		else:
			nm=str(c)
		txt=txt + nm + '\t' + str(totals[c]) + '\t'
	return txt


def header_row():
	txt=''
	columns=['Barcode id','Fusion_isoform','Fusion_suggested','Left_BreakPoint','Coverage','Right_Breakpoint','Comparison','Differences']
	txt='\t'.join(columns)
	return txt


def write_stats(d,out,quiet):
	if out != 'stdout':
		sys.stdout=open(out,'w')
	file=['f1','f2']
	
	print header_row()
	for fid in d['bc'].keys():
		if type(d['bc'][fid]) is not dict:
			continue #break
		for kval in d['bc'][fid].keys():
			if type(d['bc'][fid][kval]) is not dict:
				continue
			if 'cmp' in d['bc'][fid][kval].keys():
					if d['bc'][fid][kval]['cmp'] != 'same':
						if quiet & (d['bc'][fid][kval]['cmp'].find('miss') > -1):
							pass
							# don't print
						else:
								#print diff_row_to_string(fid,kval,d['bc'][fid][kval])
								print row_to_string(fid,kval,d['bc'][fid][kval])
	print totals_to_string(d['totals'])
				
	sys.stdout.close()
	return

#-----------------------------------------
#	MAIN
#-----------------------------------------
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Generate barcode files from one file")
	parser.add_argument('--dir1', required=True, help='First directory containing fusion caller output.')
	parser.add_argument('--dir2', required=True, help='Second directory containing fusion caller output.')
	parser.add_argument('--out', default="stdout", help='name of output.')
	parser.add_argument('--miss_off', action='store_true', help='skip missing file output.')
	args=parser.parse_args()

	stats={}
	stats['f1']={}
	stats['f2']={}
	stats['f1']['filename']=args.dir1
	stats['f2']['filename']=args.dir2
	files1=get_files(args.dir1,'*.fusions.detango.txt')
	files2=get_files(args.dir2,'*.fusions.detango.txt')
	get_stats(stats,files1)
	get_stats(stats,files2)
	write_stats(stats,args.out,args.miss_off)

