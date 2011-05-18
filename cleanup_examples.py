
#
# Cleans up by deleting object files, executable, and _output directory.
# Use after run_examples to clean up stuff not needed on webpages.

import os,sys,glob

clawdir = os.path.expandvars('$CLAW')
print "Will remove all .o, .so, _output from ",clawdir
ans = raw_input("Ok? ")
if ans.lower() not in ['y','yes']:
    print "Aborting."
    sys.exit()

os.chdir(clawdir)
exdirlist = []
for (dirpath, subdirs, files) in os.walk('.'):
    currentdir = os.path.abspath(os.getcwd())
    os.chdir(os.path.abspath(dirpath))
    
    if os.path.isfile('setrun.py'):
        print 'In directory ',dirpath
        
        for f in glob.glob('*.o'):
            os.remove(f)
        for f in glob.glob('*.so'):
            os.remove(f)
        if os.path.isdir('_output'):
            for f in glob.glob('_output/*'):
                os.remove(f)
            os.rmdir('_output')
        
    os.chdir(currentdir)

