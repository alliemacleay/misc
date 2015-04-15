#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate a manifest file for cosmos based on the
# the contents of a directory
# prints to STDOUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__= 'Allison MacLeay'

import sys
import os
import argparse

#-----------------------------------------
#	MAIN
# run umitag.py for all files in a directory
# that have the same prefix
#-----------------------------------------
if __name__ == '__main__':
	parser=argparse.ArgumentParser(description="Generate a manifest file for cosmos output")
	parser.add_argument('--dir', default='.', help='directory containing fastq output.')
	parser.add_argument('--batch_manifest', help='manifest file for multiplexed run.')
	parser.add_argument('--batch_id', default='demux', help='id for the multiplexed fastq.  This will be the first part of the name of the outputted manifest files')
	parser.add_argument('--header', required=True, help='Header for manifest file')
	parser.add_argument('--out', default='manfiles', help='directory to deposit output files')
	args=parser.parse_args()

	p={}
	batch_manifest=''
	header=''
	text=''
	if hasattr(args,'dir'):
		p['path']=args.dir
	if hasattr(args, 'out'):
		p['out']=args.out
		os.system('mkdir -p '+args.out)
	if hasattr(args, 'batch_manifest'):
		batch_manifest=args.batch_manifest
	if hasattr(args, 'header'):
		header=args.header
	if hasattr(args, 'batch_id'):
		batch_id=args.batch_id
	files=next(os.walk(p['path']))[2]
	if len(files)<1:
		print "Error: No files were found in "+args.dir+"."
	f = open(batch_manifest, 'r')
	for line in f:
		if not line.startswith("#"):  # This should really only happen once
			line = line.strip().split("\t")
			cid_id, nae_id, p5, p7 = line[2], line[3], line[6], line[7]
			ia=cid_id.find('umi_demux,')
			if ia>-1:
				cid_id=cid_id[ia+9:]
			queue={}
			for i in files:
				prefix=i.split('.')[0]
				if prefix in queue.keys():
					# remove files from queue and write both to manifest
					fastq1 = queue[prefix]
					fastq2 = i
					if(fastq1.find('1')<0):
						# switch the 2 files
						tmp=fastq1
						fastq1=fastq2
						fastq2=tmp
					fastq1 = os.path.join(p['path'], fastq1)
					fastq2 = os.path.join(p['path'], fastq2)
					if not header == '':
						text = header + "\n"
					elif not batch_header == '':
						text = batch_header
					else:
						text = "# HEADER"
					text += "\t".join(line + [fastq1, fastq2, "", ""]) + "\n"
					manfname=os.path.join(p['out'],batch_id + '_' + prefix + '.manifest')
					print str(manfname)
					manfile=open(manfname,'w')
					manfile.write(text)
					manfile.close()
					del queue[prefix]
				else:
					# add to queue
					queue[prefix]=i
		else:
			batch_header=line
	print 'manifests generator done'
