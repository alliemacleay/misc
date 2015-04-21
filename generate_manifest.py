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
# Subroutines
#-----------------------------------------

# return barcode prefix
# assumes it is the last 2 fields of the filename
def get_bc_prefix(filename):
	prefix=filename.split('.')[0]
	if prefix.find('_')>-1:
		prefix='_'.join(prefix.split('_')[-2:])
	return prefix

# - create a dictionary of manifest files grouped by 
# prefix name.
# - Zip files if needed
#
# ex.  dict['A104_P109']['files']	= [ 	'A104_P109.r1.fastq.gz','A104_P109.r2.fastq.gz',
#						'A104_P109.i1.fastq.gz','A104_P109.i2.fastq.gz']
#      dict['A104_P109']['unmapped']	= {} // ideally this would be empty otherwise a dictionary
#						of unmapped columns
#      dict['A104_P109']['A104_P109.r1.fastq.gz'] = 'full/path/name/A104_P109.r1.fastq.gz'
#      dict['A104_P109']['r1']		= 'A104_P109.r1.fastq.gz'
#	...etc.
def get_grouped_manifest_files(directory, drmaa_cmd):
	# default drmaa_cmd is empty string
	# preferably bsub -u <user_name> -o <out_log> -e <err_log> -q <medium>
	if not 'drmaa_cmd' in vars():
		drmaa_cmd = ''
	files=next(os.walk(directory))[2]
	inputs={}
	if len(files)<1:
		print "Error: No files were found in "+directory+"."
		return inputs

	needs_zipping = False

	# create a mapping for the file names
	rcols={'index1'	:'r3',
		'i1'	:'r3',
		'index2':'r2',
		'i2'	:'r2',
		'read1'	:'r1',
		'r1'	:'r1',
		'read2'	:'r4',
		'r2'	:'r4',
		'out1'	:'r1',
		'out2'	:'r4',
		}
	# create a reverse dictionary to verify mappings
	rcols_rev=list(set(rcols.values()))
	rcols_rev_dict=dict(zip(rcols_rev,[0]*len(rcols_rev)))
	for i in files:
		prefix=get_bc_prefix(i)
		# initialize dictionary for this prefix
		if prefix not in inputs.keys():
			inputs[prefix]={}
			inputs[prefix]['files']=[]
			inputs[prefix]['unmapped']=rcols_rev_dict.copy()
			inputs[prefix]['unmapped_files']=[]
		# add file information for this prefix
		fullname=os.path.join(directory,i)
		if(fullname.split('.')[-1] != 'gz'):
			fullname=fullname + '.gz'   # will be zipped if not already
			needs_zipping=True
		inputs[prefix]['files'].append(i)
		inputs[prefix][i]=fullname
		fastq_name=i.lower().replace('.','_').split('_')
		for section in fastq_name:
			if section in rcols:
				if rcols[section] in inputs[prefix]['unmapped'].keys():
					inputs[prefix][rcols[section]]=i
					del inputs[prefix]['unmapped'][rcols[section]]
					break
		else:
			inputs[prefix]['unmapped_files'].append(i)
	if(needs_zipping == True):
		cmd = 'echo zipping_files_in_'+ directory + ' && gzip ' + directory + '/*.fastq'
		if(drmaa_cmd != ''):
			cmd = drmaa_cmd + ' ' + cmd
			os.system(cmd)
	return inputs

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
	parser.add_argument('--header', required=False, help='Header for manifest file')
	parser.add_argument('--out', default='manfiles', help='directory to deposit output files')
	parser.add_argument('--drmaa_cmd', default='', help='drmaa to execute zipping files if needed')
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
	fdict=get_grouped_manifest_files(p['path'],args.drmaa_cmd)
	#print fdict
	#exit()

	not_found={}
	f = open(batch_manifest, 'r')
	manfname=os.path.join(p['out'],batch_manifest.split('/')[-1].split('.')[0] + '_DEMULTIPLEXED' + '.manifest')
	print str(manfname)
	manfile=open(manfname,'w')
	for line in f:
		if not line.startswith("#"):  # This should really only happen once
			line = line.strip().split("\t")
			if len(line)<7:		# line is blank or incorrect format
				continue
			cid_id, nae_id, p5, p7 = line[2], line[3], line[6], line[7]
			ia=cid_id.find('umi_demux,')
			if ia>-1:
				cid_id=cid_id[ia+9:]
			prefix=p5 + '_' + p7
			if prefix in fdict.keys():
				if len(line)<8:
					for ln in range(len(line)-1,8):
						line.append('')
				line[8]='' # Clear run_folder
				for fqty in ['r1','r2','r3','r4']:
					if fqty in fdict[prefix].keys():
						line.append(fdict[prefix][fqty])
					else:
						line.append('')
				if len(fdict[prefix]['unmapped_files'])>0:
					print 'Warning: Some files did not adhere to naming conventions and were not assigned written to manifest'
					print '\tFiles: '+'\n'.join(fdict[prefix]['unmapped_files'])
				text += "\t".join(line) + "\n"
				manfile.write(text)
				del fdict[prefix]
			else:
				not_found[prefix]=line
		else:
			manfile.write(line)	
	manfile.close()
	f.close()
	print 'manifests generator done'
