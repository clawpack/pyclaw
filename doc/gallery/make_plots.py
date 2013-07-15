"""
Performs 'make .plots'
in each directory to create sample results for webpage.

First try 'make all', which might do other things as well.

Sends output and errors to separate files to simplify looking for errors.
"""


def list_apps(apps_dir=None):
    """
    Searches all subdirectories of apps_dir for apps and prints out a list.
    """
    import os

    if apps_dir is None:
        from clawpack import pyclaw
        apps_dir = '/'.join(pyclaw.__path__[0].split('/')[:-2])+'/pyclaw/apps/'

    apps_dir = os.path.abspath(apps_dir)
    current_dir = os.getcwd()
    os.chdir(apps_dir)
    
    dirlist = []
    applist = []

    # Traverse directories depth-first (topdown=False) to insure e.g. that code in
    # book/chap21/radialdam/1drad is run before code in book/chap21/radialdam
    for (dirpath, subdirs, files) in os.walk('.',topdown=False):
        #By convention we assume that all python scripts are applications
        #unless they are named 'setup.py' or 'setplot.py' or have 'test' in their name.
        files = os.listdir(os.path.abspath(dirpath))
        pyfiles=[f for f in files if f.split('.')[-1]=='py']
        appfiles=[f for f in pyfiles if f.split('.')[0] not in ('setup','setplot','__init__')]
        appfiles=[f for f in appfiles if 'test' not in f]
        #Skip 3d for now because visclaw plotting is not set up for it
        appfiles=[file for file in appfiles if '3d' not in os.path.abspath(dirpath)]

        for filename in appfiles:
            dirlist.append(os.path.abspath(dirpath))
            applist.append(filename)

    os.chdir(current_dir)

    return applist, dirlist
        
def run_apps(apps_dir = None):
    """
    Runs all apps in subdirectories of apps_dir.
    """
    import os
    import subprocess

    current_dir = os.getcwd()

    app_list, dir_list = list_apps(apps_dir)
    for app, directory in zip(app_list,dir_list):
        print directory, app
        os.chdir(directory)
        process = subprocess.Popen(['python',app])#, stdout = subprocess.PIPE)
        stdout, stderr = process.communicate()

    os.chdir(current_dir)


def make_plots(apps_dir = None):
    import os,sys

    if apps_dir is None:
        print 'hello'
        from clawpack import pyclaw
        apps_dir = '/'.join(pyclaw.__path__[0].split('/')[:-2])+'/pyclaw/apps/'
        print pyclaw.__path__

    print apps_dir
    apps_dir = os.path.abspath(apps_dir)
    current_dir = os.getcwd()
 
    print "Will run code and make plots in every subdirectory of "
    print "    ", apps_dir
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

    os.chdir(apps_dir)

    goodlist_run = []
    badlist_run = []
    
    app_list, dir_list = list_apps(apps_dir)
    import subprocess
    for appname, directory in zip(app_list,dir_list):
        print 'Running ', directory, appname

        fout.write("\n=============================================\n")
        fout.write(directory)
        fout.write("\n=============================================\n")
        ferr.write("\n=============================================\n")
        ferr.write(directory)
        ferr.write("\n=============================================\n")

        os.chdir(directory)
        if not os.access(appname, os.X_OK): 
            ferr.write("Error: File is not executable")
            continue
        job = subprocess.Popen(['python','./'+appname,'htmlplot=True'], stdout = fout, stderr = ferr)
        return_code = job.wait()
        if return_code == 0:
            print "   Successful run"
            goodlist_run.append(directory)
        else:
            print "   *** Run errors encountered: see ", fname_errors
            badlist_run.append(directory)

        # flush I/O buffers:
        fout.flush()
        ferr.flush()
    
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

    os.chdir(current_dir)

if __name__=='__main__':
    make_plots()
