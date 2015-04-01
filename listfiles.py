#!/usr/bin/python

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

__author__= 'Allison MacLeay'

import sys
import os


def get_names(dir):
	
	files=next(os.walk(dir))[2]
	files=map(lambda x: x.split('.')[0], files)
	dfiles=set(files)
	return dfiles

if __name__ == '__main__':
	dir='.'
	if len(sys.argv)>1:
		dir=sys.argv[1]
	print 'dir is '+dir
	f=get_names(dir)
	cmd="echo "+'\t'.join(f)
	os.system(cmd)
