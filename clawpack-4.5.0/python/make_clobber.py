"""
Performs 'make clobber' in each subdirectory of the specified directory.
"""

import os,sys,glob

def make_clobber(rootdir):

    if rootdir==[]:   
        # if called from command line with no argument
        clawdir = os.path.expandvars('$CLAW')
        rootdir = clawdir
    else:
        # called with an argument, try to use this for rootdir:
        rootdir = rootdir[0]
        rootdir = os.path.abspath(rootdir)
    
    print "Will execute 'make clobber' in every subdirectory of "
    print "    ", rootdir
    ans = raw_input("Ok? ")
    if ans.lower() not in ['y','yes']:
        print "Aborting."
        sys.exit()
    
    os.chdir(rootdir)
    
    for (dirpath, subdirs, files) in os.walk('.'):
        currentdir = os.path.abspath(os.getcwd())
        os.chdir(os.path.abspath(dirpath))
        
        if os.path.isfile('Makefile'):
            print 'In directory ',dirpath
            try:
                os.system('make clobber')
            except:
                pass
            
        os.chdir(currentdir)
        

if __name__=='__main__':
    make_clobber(sys.argv[1:])
