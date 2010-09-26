"""
Performs 'make all' in each directory with a setrun.py file.
Usually this makes .plots and .htmls but might be overwritten by
a different set of things in the Makefile,  for example if a reference
solution must be created.

Sends output and errors to separate files to simplify looking for errors.
"""

import os,sys,glob
try:
    import subprocess
except:
    print '*** Error: require subprocess module from Python 2.4 or greater'
    raise ImportError()


def make_all(rootdir):

    if rootdir==[]:   
        # if called from command line with no argument
        clawdir = os.path.expandvars('$CLAW')
        rootdir = clawdir
    else:
        # called with an argument, try to use this for rootdir:
        rootdir = rootdir[0]
        rootdir = os.path.abspath(rootdir)

    print "Will 'make all' in every subdirectory of "
    print "    ", rootdir
    ans = raw_input("Ok? ")
    if ans.lower() not in ['y','yes']:
        print "Aborting."
        sys.exit()
    
    fname_output = 'make_all_output.txt'
    fout = open(fname_output, 'w')
    fout.write("ALL OUTPUT FROM RUNNING EXAMPLES\n\n")

    fname_errors = 'make_all_errors.txt'
    ferr = open(fname_errors, 'w')
    ferr.write("ALL ERRORS FROM RUNNING EXAMPLES\n\n")

    # Set environment variable to allow downloading topography for
    # GeoClaw examples:
    os.environ['CLAW_TOPO_DOWNLOAD'] = 'True'

    #os.chdir(rootdir)
    goodlist = []
    badlist = []
    
    # Traverse directories depth-first (topdown=False) to insure e.g. that code in
    # book/chap21/radialdam/1drad is run before code in book/chap21/radialdam
    
    for (dirpath, subdirs, files) in os.walk(rootdir,topdown=False):
        currentdir = os.path.abspath(os.getcwd())
        os.chdir(os.path.abspath(dirpath))
        rootdirpath = os.path.join(os.path.split(rootdir)[1],dirpath)
        

        if os.path.isfile('setrun.py'):

            fout.write("\n=============================================\n")
            fout.write(rootdirpath)
            fout.write("\n=============================================\n")
            ferr.write("\n=============================================\n")
            ferr.write(rootdirpath)
            ferr.write("\n=============================================\n")

            # flush I/O buffers:
            fout.flush()
            ferr.flush()

            print "Running 'make all' in ",rootdirpath
            job = subprocess.Popen(['make','all'], \
                             stdout=fout, stderr=ferr)
            return_code = job.wait()
            if return_code == 0:
                print "   Successful completion"
                goodlist.append(dirpath)
            else:
                print "   *** Errors encountered: see ", fname_errors
                badlist.append(dirpath)
            
        os.chdir(currentdir)
        
    
    print ' '
    print 'Ran Clawpack and created output in directories:'
    for d in goodlist:
        print '   ',d
    print ' '
    
    print 'Errors encountered in the following directories:'
    for d in badlist:
        print '   ',d
    print ' '
    
    fout.close()
    ferr.close()
    print 'For all output see ', fname_output
    print 'For all errors see ', fname_errors

if __name__=='__main__':
    make_all(sys.argv[1:])
