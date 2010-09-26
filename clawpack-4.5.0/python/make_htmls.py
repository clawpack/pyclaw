"""
Performs 'make .htmls' in each directory to create html files from 
source files using $CLAW/doc/clawcode2html.py.

This is useful so that links will work, e.g. if you build the Sphinx
documentation locally.

Use make_clean.py to first remove old versions if desired, or
set remove_first to True below.
"""

import os,sys,glob

remove_first = False   # to force remaking of htmls
                       # only set to True if you're sure

def make_htmls(rootdir):

    if rootdir==[]:   
        # if called from command line with no argument
        clawdir = os.path.expandvars('$CLAW')
        rootdir = clawdir
    else:
        # called with an argument, try to use this for rootdir:
        rootdir = rootdir[0]
        rootdir = os.path.abspath(rootdir)
    
    print "Will make htmls in all of ",rootdir
    ans = raw_input("Ok? ")
    if ans.lower() not in ['y','yes']:
        print "Aborting."
        sys.exit()
    
    os.chdir(rootdir)
    goodlist = []
    badlist = []
    
    
    for (dirpath, subdirs, files) in os.walk('.'):
        currentdir = os.path.abspath(os.getcwd())
        os.chdir(os.path.abspath(dirpath))
        
        if os.path.isfile('Makefile'):
            print 'In directory ',dirpath
            if remove_first:
                os.system('rm -f *.html')
            try:
                os.system('make .htmls')
                goodlist.append(dirpath)
            except:
                print "*** make .htmls Failed"
                badlist.append(dirpath)
            
        os.chdir(currentdir)
        
    print ' '
    print 'Created htmls in directories:'
    for d in goodlist:
        print '   ',d
    print ' '
    
    print 'Failed in the following directories:'
    for d in badlist:
        print '   ',d
    print ' '
    

if __name__=='__main__':
    make_htmls(sys.argv[1:])
