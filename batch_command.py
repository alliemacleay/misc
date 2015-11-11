__author__ = 'Allison MacLeay'

import sys
import os
import argparse
import time
from glob import glob
from tempfile import mkdtemp

"""
Run a command on every file in a directory
"""

def get_cmd(fname, cmd, params):
    """ return command """

    out_name = os.path.join(params['out'], params['name'] + fname)
    cmd = cmd.replace('$INPUT', fname)
    cmd = cmd.replace('$OUT', out_name)

    return cmd


def bsub_cmd(user, log, flags):
    """ return bsub command"""
    if flags != '':  # make sure it is buffered by spaces
        if flags[0] != ' ':
            flags = ' ' + flags
        if flags[-1] != ' ':
            flags = flags + ' '
    return 'bsub -q medium -u ' + user + ' -o ' + flags + os.path.join(log, 'lsf_out.log') + ' -e ' + os.path.join(log,
                                                                                                                   'lsf_err.log') + ' '


def get_names(path, pattern):
    """ Return all unique file prefixes """
    return glob(os.path.join(path, pattern))


def check_done(file_num, path):
    """
    Delay completion of script until all
    files are written
    """
    start = time.time()
    timeout = (24 * 60 * 60)  # 24 hours
    done = 0
    while done == 0:
        files = next(os.walk(path))[2]
        all_closed = 0
        if len(files) == file_num:
            all_closed = 1
            for f in files:
                if f.find('tmp') > 0:
                    all_closed = 0
                    print 'found tmp file ' + f + '.  Waiting...'
                    continue
        if all_closed == 1:
            done = 1
        elif (time.time() - start) > timeout:
            done = 1
            print 'Job timed out'
        else:
            time.sleep(5)
            print 'checking for job completion after waiting %d seconds' % (time.time() - start)
            print 'searching for ' + str(file_num) + ', found ' + str(len(files)) + ' files in ' + path


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="Batch command helper",
                                     description="Run umitag utility in batches of similarly prefixed names.  "
                                                 "Pass in a command in the format ls $INPUT >> $OUTPUT")
    parser.add_argument('--dir', default='.', help='directory containing input files')
    parser.add_argument('--cmd', default='ls $INPUT',
                        help='command to run on each file in the specified directory')
    parser.add_argument('--out', default='tmp', help='directory to deposit output files')
    parser.add_argument('--batch_name', default='batch', help='name to prepend to processed files')
    parser.add_argument('--pattern', default='*.*', help='match pattern for input files.  default=*.*')
    parser.add_argument('--verbose', default=False, help='verbosity.  default=False')

    # optional bsub options
    parser.add_argument('--bsub', action='store_true', help='use bsub')
    parser.add_argument('--log', help='directory to deposit bsub log files')
    parser.add_argument('--bsub_user', default='am282', help='user name for bsub command. default=am282')
    parser.add_argument('--bsub_mod', default='', help='extra parameters for bsub command')
    parser.add_argument('--output_count', help='number of expected output files per input.  Process will wait to '
                                               'complete until all files are created.  Leave this flag out to avoid '
                                               'the process waiting')
    args = parser.parse_args()

    if args.out == 'tmp':
        print 'WARNING: Output directory was not supplied.  Temporary directory will be used instead.'
    p = {}
    p['path'] = args.dir
    p['name'] = args.batch_name
    if args.out != 'tmp':
        os.system('mkdir -p {}'.format(p['out']))
    else:
        p['out'] = mkdtemp()

    # bsub options
    if args.bsub:
        if not hasattr(args, 'log'):
            raise Exception("ERROR: If using bsub you must specify a directory for the log files.")
        os.system('mkdir -p {}'.format(args.log))
        os.system('ls {} >> {}'.format(p['path'], os.path.join(args.log, 'ls_inputdir.txt')))
        expected_output = int(args.output_count) if hasattr(args, 'output_count') else 0

    files = get_names(args.dir, args.pattern)
    if len(files) < 1:
        print "Error: No file prefixes were found in {}.".format(args.dir)

    count_lsf = 0
    for fname in files:
        if (fname.find('undetermined') > -1):
            # skip undeterminded for now
            cmd = 'echo skipping undetermined files'
        elif not args.bsub:
            cmd = get_cmd(fname, args.cmd, p)
        else:
            cmd = bsub_cmd(args.bsub_user, args.log, args.bsub_mod) + get_cmd(fname, args.cmd, p)
            # Keep track of lsf job for listener
            count_lsf += expected_output  # stays 0 if the process will not wait until completion
        if args.verbose:
            print 'batch command helper running command:\n' + cmd
        os.system(cmd)

    if count_lsf > 0:
        check_done(count_lsf, p['out'])

    if args.verbose:
        print 'batch_process done'
