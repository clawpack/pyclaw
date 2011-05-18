"""
Searches all subdirectories for apps and prints out a list.
"""

import os,sys,glob
try:
    import subprocess
except:
    print '*** Error: require subprocess module from Python 2.4 or greater'
    raise ImportError()


def list_apps(rootdir):

    if rootdir==[]:   
        # if called from command line with no argument
        clawdir = os.path.expandvars('$PYCLAW/apps')
        rootdir = clawdir
    else:
        # called with an argument, try to use this for rootdir:
        rootdir = rootdir[0]
        rootdir = os.path.abspath(rootdir)

    os.chdir(rootdir)
    applist = []
    
    # Traverse directories depth-first (topdown=False) to insure e.g. that code in
    # book/chap21/radialdam/1drad is run before code in book/chap21/radialdam
    
    for (dirpath, subdirs, files) in os.walk('.',topdown=False):
        currentdir = os.path.abspath(os.getcwd())
        os.chdir(os.path.abspath(dirpath))
        rootdirpath = os.path.join(os.path.split(rootdir)[1],dirpath[2:])
        

        #By convention we assume that all python scripts are applications
        #unless they are named 'setup.py' or 'setplot.py'.
        files = os.listdir('.')
        pyfiles=[file for file in files if file.split('.')[-1]=='py']
        appfiles=[file for file in pyfiles if file.split('.')[0] not in ('setup','setplot')]

        if appfiles!=[]:
            print rootdirpath
            for appname in appfiles: print '     ',appname
           
        os.chdir(currentdir)
        
if __name__=='__main__':
    list_apps(sys.argv[1:])
