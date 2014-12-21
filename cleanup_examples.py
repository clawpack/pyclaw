#!/usr/bin/env python

r"""
Cleans up by deleting object files, executable, and _output directory.
"""

import os,sys,glob
import shutil

if __name__ == "__main__":

    examples_dir = os.path.abspath('./examples')
    print "Will remove all _output, _plots, and build directories from ",examples_dir
    ans = raw_input("Ok? ")
    if ans.lower() not in ['y','yes']:
        print "Aborting."
        sys.exit()

    os.chdir(examples_dir)
    exdirlist = []
    for (dirpath, subdirs, files) in os.walk('.'):
        currentdir = os.path.abspath(os.getcwd())
        os.chdir(os.path.abspath(dirpath))
        
        print 'In directory ',dirpath
        
        if os.path.isdir('_output'):
            shutil.rmtree('./_output')
        if os.path.isdir('_plots'):
            shutil.rmtree('./_plots')
        if os.path.isdir('build'):
            shutil.rmtree('./build')
            
        os.chdir(currentdir)

