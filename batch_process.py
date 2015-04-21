#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	There's no good way in cosmos to create a tag
#	based on output or to conditionally run a process
#	This is a wrapper to run the umi utilities.
#
# Run umitag.py for a directory in batches
# This is a helper script for Martin Aryee's
# scripts to demultiplex Illumina sequencing
# reads with sample specific and molecule
# specific tags.
# https://github.com/aryeelab/umi/wiki 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__= 'Allison MacLeay'

import sys
import os
import argparse
import time

#-----------------------------------------
# Get umitag.py bash command for specified
# file name prefix
#-----------------------------------------
def get_cmd(group, umitag, params):
	path='.'
	outdir=''
	cmd=''
	if 'path' in params.keys():
		path=params['path']
	if 'out' in params.keys():
		outdir=params['out']
		if(outdir[0]=='/'):
			outdir=outdir[1:]
		else:
			outdir=os.path.abspath(outdir)[1:]
		if(outdir[-1] != '/'):
			outdir=outdir+'/'
	fastq_r1=path + '/' + group +'.r1.fastq'
	fastq_r2=path + '/' + group+'.r2.fastq'
	out1= '' + group+'.out1.fastq'
	out2= '' + group+'.out2.fastq'
	fastq_i1=path + '/' + group +'.i1.fastq'
	fastq_i2=path + '/' + group+'.i2.fastq'
	cmd= (cmd+""" \
	python {umitag} --read1_in {fastq_r1} \
	--read2_in {fastq_r2} \
        --read1_out {outdir}{out1} \
        --read2_out {outdir}{out2} \
        --index1 {fastq_i1} \
        --index2 {fastq_i2} \
        """ ).format(
		umitag  =	umitag,
		fastq_r1= 	fastq_r1,
		fastq_r2= 	fastq_r2,
		outdir=		'/' + outdir,
		out1=		out1,
		out2=		out2,
		fastq_i1=	fastq_i1,
		fastq_i2=	fastq_i2,
	)
	return cmd
	
#-----------------------------------------
# Return all unique file prefixes
#-----------------------------------------
def get_names(dir):
	
	files=next(os.walk(dir))[2]
	files=map(lambda x: x.split('.')[0], files)
	dfiles=set(files)
	return dfiles

#-----------------------------------------
# Delay completion of script until all
#  files are written
#-----------------------------------------
def check_done(file_num,path):
	start=time.time()
	timeout=(10*60) # 10 minutes
	done=0
	while done==0:
		files=next(os.walk(path))[2]
		all_closed=0
		if len(files)==file_num:
			all_closed=1
			for f in files:
				if f.find('tmp')>0:
					all_closed=0
					print 'found tmp file ' + f + '.  Waiting...'
					continue
		if all_closed==1:
			done=1
		elif (time.time()-start)>timeout:
			print 'Job timed out'
		else:
			time.sleep(5)
			print 'checking for job completion after waiting %d seconds' % (time.time()-start)
			print 'searching for ' + str(file_num) + ', found ' + str(len(files)) + ' files in ' + path

#-----------------------------------------
#	MAIN
# run umitag.py for all files in a directory
# that have the same prefix
#-----------------------------------------
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Run umitag utility in batches of similarly prefixed names.")
	parser.add_argument('--dir', default='.', help='directory containing output of umi demultiplex')
	parser.add_argument('--umi_tool',default='umitag.py', help='umi tag tool absolute path. default is umitag.py in current directory')
	parser.add_argument('--out', default='tagout', help='directory to deposit output files')
	parser.add_argument('--log', default='batch_log', help='directory to deposit bsub log files')
	parser.add_argument('--bsub_off',  action='store_true', help='turn bsub off to test on systems without lsf')
	#parser.add_argument('--undet',  action='store_false', help='include reads less than parameter set my min reads.  Default will skip files named undetermined')
	args=parser.parse_args()

	p={}
	if hasattr(args,'dir'):
		p['path']=args.dir
	if hasattr(args, 'out'):
		p['out']=args.out
		os.system('mkdir -p '+args.out)
	if hasattr(args, 'log'):
		os.system('mkdir -p '+args.log)
		os.system('ls ' + p['path'] + ' >> ' + args.log + '/ls_inputdir.txt')
	f=get_names(args.dir)
	if len(f)<1:
		print "Error: No file prefixes were found in "+args.dir+"."
	count_lsf=0
	for tag in f:
		if (tag.find('undetermined') > -1 ):
			# skip undeterminded for now
			cmd='echo skipping undetermined files'
		elif(args.bsub_off):
			cmd=get_cmd(tag, args.umi_tool, p)
		else:
			cmd='bsub -q medium -u am282 -o ' + os.path.join(args.log,'lsf_out.log') + ' -e ' + os.path.join(args.log,'lsf_err.log') + ' ' + get_cmd(tag, args.umi_tool, p)
			# Keep track of lsf job for listener
			count_lsf=count_lsf+2
		
		print 'batch process running command:\n' + cmd
		os.system(cmd)
		
	if (count_lsf>0):
		check_done(count_lsf,args.out)
	print 'batch_process done'
			
