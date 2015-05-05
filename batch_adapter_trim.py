# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# There's no good way in cosmos to create a tag
# based on output or to conditionally run a process
# This is a wrapper to run the umi utilities.
#
# Run a program (SeqPrep) for a directory in batches
# This is a helper script for Martin Aryee's
# scripts to demultiplex Illumina sequencing
# reads with sample specific and molecule
# specific tags.
# https://github.com/aryeelab/umi/wiki 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__author__ = 'Allison MacLeay'

import sys
import os
import argparse
import time

# -----------------------------------------
# Get bash command for specified
# file name prefix
# -----------------------------------------
def get_cmd(group, prog, params):
    path = '.'
    outdir = ''
    cmd = ''
    job_group = ''
    a1 = ''
    a2 = ''
    if 'path' in params.keys():
        path = params['path']
    if 'out' in params.keys():
        outdir = params['out']
        if (outdir[0] == '/'):
            outdir = outdir[1:]
        else:
            outdir = os.path.abspath(outdir)[1:]
        if (outdir[-1] != '/'):
            outdir = outdir + '/'
    if 'job_group' in params.keys():
        job_group = params['job_group']
    if 'a1' in params.keys():
        a1 = params['a1']
    if 'a1' in params.keys():
        a2 = params['a2']
    fastq_r1 = os.path.join(path, group + '.r1.fastq.gz')
    fastq_r2 = os.path.join(path ,group + '.r2.fastq.gz')
    out1 = '' + group + '_trim.r1.fastq.gz'
    out2 = '' + group + '_trim.r2.fastq.gz'
    cmd = (cmd + """ \
    {prog} -f {fastq_r1} \
-r {fastq_r2} \
-1 {outdir}{out1} \
-2 {outdir}{out2} \
-A {adapter1} \
-B {adapter2} \
        """ ).format(
        prog=prog,
        fastq_r1=fastq_r1,
        fastq_r2=fastq_r2,
        outdir='/' + outdir,
        out1=out1,
        out2=out2,
        adapter1=a1,
        adapter2=a2,
    )
    return cmd


# -----------------------------------------
# Return all unique file prefixes
# -----------------------------------------
def get_names(dir):
    files = next(os.walk(dir))[2]
    files = map(lambda x: x.split('.')[0], files)
    dfiles = set(files)
    return dfiles


# -----------------------------------------
# Delay completion of script until all
# files are written
#-----------------------------------------
def check_done(group_id,ct):
    start = time.time()
    timeout = (24 * 60 * 60)  # 24 hours
    done = 0

    while done==0:
        if (time.time()-start)>timeout:
            print 'Job timed out'
            done=1
        else:
            if are_jobs_done(group_id,ct):
                print 'No running jobs found in group ' + group_id
                done=1
            time.sleep(10)
            print 'checking for job completion after waiting %d seconds' % (time.time()-start)

def get_group_id(group):
    i = os.popen("bjgroup | grep "+ group).read()
    i = i.split('\n')
    num=''
    max=0
    for line in i:
        if are_jobs_done(line.split()[0],0):
            # clean old jobs	
            os.system("bgdel " + line)
        num=line.strip().split('_')[0].split('/')[-1]
        if str(num)=='':
            num=0
        if type(num) is str:
                if not num.isdigit():
                    num=0
        if int(num) > max:
            max = int(num)
    print 'last id was ' + str(max)
    return (group + '/' + str(max+1) + "_trim")

def are_jobs_done(group,lsf_ct):
    status=os.popen("bjobs -g " + group).read().split('\n')
    ct=0
    group_status=True
    for line in status:
        ct+=1
        if ct==1:
            continue
        status_arr=line.strip().split()
        if len(status_arr)>2:
            status=status_arr[2]
        else:
            status=''
        if status == 'RUN':
            group_status=False
        if( (ct-1) != lsf_ct):
            print 'Different number of jobs passed in ('+str(lsf_ct)+') and recovered (' + str(ct-1) + ')'
        return group_status


#-----------------------------------------
#	MAIN
# run a program (SeqPrep) for all files in a directory
# that have the same prefix
#-----------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
            description="Run command for removing adapter sequecnes in batches of similarly prefixed names.")
    parser.add_argument('--dir', default='.', help='directory containing output of umi demultiplex')
    parser.add_argument('--script', default='./SeqPrep', help='SeqPrep absolute path.  default is SeqPrep in current directory')
    parser.add_argument('--a1', required=True, help='Adapter 1')
    parser.add_argument('--a2', required=True, help='Adapter 2')
    parser.add_argument('--out', default='tagout', help='directory to deposit output files')
    parser.add_argument('--log', default='batch_log', help='directory to deposit bsub log files')
    parser.add_argument('--bsub_off', action='store_true', help='turn bsub off to test on systems without lsf')
    #parser.add_argument('--undet',  action='store_false', help='include reads less than parameter set my min reads.  Default will skip files named undetermined')
    args = parser.parse_args()

    p = {}
    lsf_group=''
    lsf_group_cmd=''
    if hasattr(args, 'dir'):
        p['path'] = args.dir
    if hasattr(args, 'out'):
        p['out'] = args.out
        os.system('mkdir -p ' + args.out)
    if hasattr(args, 'log'):
        os.system('mkdir -p ' + args.log)
        os.system('ls ' + p['path'] + ' >> ' + args.log + '/ls_inputdir.txt')
    if hasattr(args, 'a1'):
        p['a1']=args.a1
    if hasattr(args, 'a2'):
        p['a2']=args.a1
    f = get_names(args.dir)
    if len(f) < 1:
        print "Error: No file prefixes were found in " + args.dir + "."
    count_lsf = 0
    if not args.bsub_off:
        lsf_group = get_group_id("/demux")
        lsf_group_cmd=' -g ' + lsf_group
    for tag in f:
        if (tag.find('undetermined') > -1 ):
            # skip undeterminded for now
            cmd = 'echo skipping undetermined files'
        elif (args.bsub_off):
            cmd = get_cmd(tag, args.script, p)
        else:
            cmd = 'bsub -q medium -u am282 -o ' + os.path.join(args.log, 'lsf_out.log') + ' -e ' + os.path.join(
                args.log, 'lsf_err.log') + lsf_group_cmd + ' ' + get_cmd(tag, args.script, p)
            # Keep track of lsf job for listener
            count_lsf = count_lsf + 1

        print 'batch process running command:\n' + cmd
        os.system(cmd)

    if (count_lsf > 0):
        if lsf_group != '':
            check_done(lsf_group, count_lsf)
    print 'batch_process done'
			
