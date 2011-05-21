"""
Performs 'make .plots'
in each directory to create sample results for webpage.

First try 'make all', which might do other things as well.

Sends output and errors to separate files to simplify looking for errors.
"""

import os,sys,glob
try:
    import subprocess
except:
    print '*** Error: require subprocess module from Python 2.4 or greater'
    raise ImportError()


def make_plots(rootdir):

    if rootdir==[]:   
        # if called from command line with no argument
        clawdir = os.path.expandvars('$PYCLAW/apps')
        rootdir = clawdir
    else:
        # called with an argument, try to use this for rootdir:
        rootdir = rootdir[0]
        rootdir = os.path.abspath(rootdir)

    print "Will run code and make plots in every subdirectory of "
    print "    ", rootdir
    ans = raw_input("Ok? ")
    if ans.lower() not in ['y','yes']:
        print "Aborting."
        sys.exit()
    
    fname_output = 'make_plots_output.txt'
    fout = open(fname_output, 'w')
    fout.write("ALL OUTPUT FROM RUNNING EXAMPLES\n\n")

    fname_errors = 'make_plots_errors.txt'
    ferr = open(fname_errors, 'w')
    ferr.write("ALL ERRORS FROM RUNNING EXAMPLES\n\n")

    os.chdir(rootdir)
    goodlist_build = []
    badlist_build = []

    goodlist_run = []
    badlist_run = []
    
    # Traverse directories depth-first (topdown=False) to insure e.g. that code in
    # book/chap21/radialdam/1drad is run before code in book/chap21/radialdam
    
    for (dirpath, subdirs, files) in os.walk('.',topdown=False):
        currentdir = os.path.abspath(os.getcwd())
        os.chdir(os.path.abspath(dirpath))
        rootdirpath = os.path.join(os.path.split(rootdir)[1],dirpath)
        

        #By convention we assume that all python scripts are applications
        #unless they are named 'setup.py' or 'setplot.py'.
        files = os.listdir('.')
        pyfiles=[file for file in files if file.split('.')[-1]=='py']
        appfiles=[file for file in pyfiles if file.split('.')[0] not in ('setup','setplot')]
        if appfiles!=[]:

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
            job = subprocess.Popen(['make','all'], stdout=fout, stderr=ferr)
            return_code = job.wait()
            if return_code == 0:
                print "   Successful build"
                goodlist_build.append(dirpath)
            else:
                print "   *** Build errors encountered: see ", fname_errors
                badlist_build.append(dirpath)

            for appname in appfiles:
                if not os.access(appname, os.X_OK): 
                    ferr.write("Error: File is not executable")
                    continue
                appname_short=appname.split('.')[0]
                print appname_short
                job = subprocess.Popen(['./'+appname,'htmlplot=True'],stdout=fout,stderr=ferr)
                return_code = job.wait()
                if return_code == 0:
                    print "   Successful run"
                    goodlist_run.append(dirpath)
                else:
                    print "   *** Run errors encountered: see ", fname_errors
                    badlist_run.append(dirpath)

           
        os.chdir(currentdir)
        
    
    print ' '
    print 'Built successfully in directories:'
    for d in goodlist_build:
        print '   ',d
    print ' '
    
    print 'Build errors encountered in the following directories:'
    for d in badlist_build:
        print '   ',d
    print ' '

    print ' '
    print 'Ran PyClaw and created output and plots in directories:'
    for d in goodlist_run:
        print '   ',d
    print ' '
    
    print 'Run or plot errors encountered in the following directories:'
    for d in badlist_run:
        print '   ',d
    print ' '
    
    fout.close()
    ferr.close()
    print 'For all output see ', fname_output
    print 'For all errors see ', fname_errors

if __name__=='__main__':
    make_plots(sys.argv[1:])
