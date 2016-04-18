import sys
import os
import glob
import subprocess


def file_by_file_diff(dir1, dir2, ext):
    """
    Compare files with the same extensions to each other from different directories
    """
    if ext[0] == '.':
        ext = ext[1:]
    files1 = glob.glob(os.path.join(dir1, '*.{}'.format(ext)))
    all_same = True
    for fname in files1:
        lfile = os.path.join(dir1, fname)
        rfile = os.path.join(dir2, fname)
        cmd = 'diff {} {}'.format(lfile, rfile)
        print cmd
        response = subprocess.Popen(cmd.split(' '), stdout=subprocess.PIPE, stderr=subprocess.STDOUT).stdout.readlines()
        if len(response) == 0:
            print "files are identical"
        else:
            print "\nThese contain differences:\n{}\n{}\n"
            all_same = False
    return all_same

if __name__ == "__main__":
    usage = 'python dir_diff.py directory1 directory2 extension\n\n' \
            'Compare files with the same extensions to each other from different directories'
    if len(sys.argv) < 4:
        print usage
    dir1 = sys.argv[1]
    dir2 = sys.argv[2]
    ext = sys.argv[3]
    if not file_by_file_diff(dir1, dir2, ext):
        print "\n\nFiles compared did not all match!\n"
    else:
        print "\n\nFiles all matched.\n"
