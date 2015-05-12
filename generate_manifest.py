#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Generate a manifest file for cosmos based on the
# the contents of a directory
# prints to STDOUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__= 'Allison MacLeay'

import sys
import os
import argparse
import time

#-----------------------------------------
# Subroutines
#-----------------------------------------

# return file prefix
def get_bc_prefix(filename):
	prefix=filename.split('.')[0]
	return prefix

# - create a dictionary of manifest files grouped by 
# prefix name.
# - Zip files if needed
#
# ex.  dict['input']['A104_P109']['files']	= [ 	'A104_P109.r1.fastq.gz','A104_P109.r2.fastq.gz',
#						'A104_P109.i1.fastq.gz','A104_P109.i2.fastq.gz']
#      dict['input]['A104_P109']['unmapped']	= {} // ideally this would be empty otherwise a dictionary
#						of unmapped columns
#      dict['input']['A104_P109']['A104_P109.r1.fastq.gz'] = 'full/path/name/A104_P109.r1.fastq.gz'
#      dict['input']['A104_P109']['r1']		= 'A104_P109.r1.fastq.gz'
#	...etc.
def get_grouped_manifest_files(directory, drmaa_cfg):
	# default drmaa_cmd is empty string
	# preferably bsub -u <user_name> -o <out_log> -e <err_log> -q <medium>
	drmaa_cmd=''
	ret_dict={}
	ret_dict['zip_info']={}
	ret_dict['zip_info']['job_type']=''
	if 'drmaa_cfg' in vars():
		ret_dict['zip_info']['job_type']=drmaa_cfg
		if drmaa_cfg=='lsf':
			drmaa_cmd = 'bsub -u am282 -o o_zip.log -e e_zip.log -q medium'
	files=next(os.walk(directory))[2]
	inputs={}
	if len(files)<1:
		print "Error: No files were found in "+directory+"."
		return inputs

	needs_zipping = False

	# create a mapping for the file names
	rcols={#'index1'	:'r3',
		#'i1'	:'r3',
		#'index2':'r2',
		#'i2'	:'r2',
		'read1'	:'r1',
		'r1'	:'r1',
		'read2'	:'r4',
		'r2'	:'r4',
		'out1'	:'r1',
		'out2'	:'r2',
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
					inputs[prefix][rcols[section]]=fullname
					del inputs[prefix]['unmapped'][rcols[section]]
					break
		else:
			inputs[prefix]['unmapped_files'].append(i)
	if(needs_zipping == True):
		ret_dict['zip_info']['job_id'] = zip_fastq_files(directory,drmaa_cmd)
	else:
		print 'Files already zipped'
	ret_dict['inputs']=inputs
	return ret_dict

def zip_fastq_files(folder,drmaa_cmd):
	job_id=''
	cmd = 'echo zipping_files_in_'+ folder + ' && gzip ' + os.path.join(folder , '*.fastq')
	if(drmaa_cmd != ''):
		cmd = drmaa_cmd + ' "' + cmd + '"'
	out = os.popen(cmd).read()
	if(out.find('<') > -1 ) & (out.find('>') > -1 ):
		# bsub output found
		job_id=out.split('<')[1].split('>')[0]
	return job_id

def check_done(p):
	start=time.time()
	timeout=(30*60)	# 30 minutes
	done=False
	if 'job_type' in p.keys():
		if p['job_type'] != 'lsf':
			print 'No support for tracking jobs that are not lsf'
			done=True
	if 'job_id' not in p.keys():
		print 'Error: Missing job id'
		done=True
	else:
		job_id=p['job_id']
	while (done == False):
		if (time.time()-start) > timeout:
			print 'Job timed out after 30 minutes.'
			done=True
		status = get_job_status(job_id)
		if status == 'PEND':
			pass
		elif status == 'RUN':
			pass
		else:
			done=True
		if not done:
			run_dur=round((time.time()-start),0)
			if (run_dur % 30) < 6:
				print 'waiting for job ' + job_id + ' status to be DONE. Status is ' + status + ' at ' + str(time.time()-start) + ' seconds'
			time.sleep(5)
	return

def get_job_status(job_id):
	status=''
	field_num=3
	ret=os.popen('bjobs ' + job_id).read()
	ret=ret.split('\n')[1]
	for n in range(field_num-1):
		sn=ret.find(' ')
		ret=ret[sn:]
		while(ret.find(' ') == 0):
			ret=ret[1:]
	status=ret[:ret.find(' ')]
	return status

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
	parser.add_argument('--drmaa_config', default='', help='drmaa configuration to execute zipping files if needed')
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
	freturn=get_grouped_manifest_files(p['path'],args.drmaa_config)
	fdict=freturn['inputs']
	zip_info=freturn['zip_info']
	#print fdict
	#exit()

	not_found={}
	f = open(batch_manifest, 'r')
	man_name=batch_manifest.split('/')[-1].split('.')[0] + '_DEMULTIPLEXED' + '.manifest'
	manfname=os.path.join(p['out'],man_name)
	mantmp='.tmp_' + man_name
	print str(manfname)
	manfile=open(mantmp,'w')
	for line in f:
		if not line.startswith("#"):  # This should really only happen once
			found_prefix = False
			line = line.strip().split("\t")
			if len(line)<7:		# line is blank or incorrect format
				continue
			cid_id, nae_id, p5, p7 = line[2], line[3], line[6], line[7]
			ia=cid_id.find('umi_demux,')
			if ia>-1:
				cid_id=cid_id[ia+9:]
			prefix=p5 + '_' + p7
			prefix_key=''
			for rec_name in fdict.keys():
				if rec_name.find(prefix)>-1:
					found_prefix = True
					prefix_key=rec_name
			if found_prefix:
				if len(line)<8:
					for ln in range(len(line)-1,8):
						line.append('')
				line[8]='' # Clear run_folder
				for fqty in ['r1','r2','r3','r4']:
					if fqty in fdict[prefix_key].keys():
						line.append(fdict[prefix_key][fqty])
					else:
						line.append('')
				if len(fdict[prefix_key]['unmapped_files'])>0:
					print 'Warning: Some files did not adhere to naming conventions and were not assigned written to manifest'
					print '\tFiles: '+'\n'.join(fdict[prefix_key]['unmapped_files'])
				text = "\t".join(line) + "\n"
				manfile.write(text)
				del fdict[prefix_key]
			else:
				not_found[prefix]=line
		else:
			manfile.write(line)	
	manfile.close()
	f.close()
	if len(not_found) > 0:
		print 'Warning: the following keys were in the batch manifest and did not correspond with any fastq output files:'
		for i in not_found.keys():
			print i + ' '.join(not_found[i])
	print 'wrote to ' + manfname
	if 'job_id' in zip_info.keys():
		check_done(zip_info)
	cp_cmd='mv '+mantmp+ ' ' + manfname
	os.system(cp_cmd) 
	print 'manifests generator done'
