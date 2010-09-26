#
# setenv.py 
#
# Customize your installation of claw.
# Set parameters in this file and execute it before running setup.py.
#
# CLAW = path to this directory where claw lives
# PYTHONPATH = any path needed to access modules required by python.
#    The $CLAW/python directory will be appended to this path.
# IPYTHONDIR = path to ipython profiles, icluding ipythonrc-claw.
# MATLABPATH = any path needed to use Matlab.
#    The $CLAW/matlab directory will be appended to this path.
# LD_LIBRARY_PATH or DYLD_LIBRARY_PATH = Path to shared libraries
# FC = fortran command, e.g. 'gfortran'
# Note that gfortran or some other flavor of f90/95 is now required
# due to dynamic memory allocation.

import os
import sys

clawdir = os.path.abspath('.')   
CLAW = clawdir

print " "
print "------------------------------------------------------------"

try:
    # check if the Fortran Compiler is already set:
    FC = os.environ['FC']
except:
    FC = 'gfortran'

if FC in ['f77','g77']:
    print '*** FC = ',FC,' will not work with this version.'  
    print '    gfortran or other flavor of f90/95 required'
    print '*** resetting FC to gfortran\n'
    FC = 'gfortran'

clawpythondir = os.path.join(clawdir,'python')
try:
    PYTHONPATH = os.environ['PYTHONPATH']
except:
    PYTHONPATH = clawpythondir

if clawpythondir not in PYTHONPATH:
    PYTHONPATH = clawpythondir +":"+ PYTHONPATH


try:
    IPYTHONDIR = os.environ['IPYTHONDIR']
    print 'IPYTHONDIR already set: you may want to move %s ' \
            % os.path.join(ipythondir,'ipythonrc-claw')
    print '   to this directory: ',IPYTHONDIR
except:
    IPYTHONDIR = os.path.join(clawpythondir,'ipythondir')


clawmatlabdir = os.path.join(clawdir,'matlab')
try:
    MATLABPATH = os.environ['MATLABPATH']
except:
    MATLABPATH = clawmatlabdir

if clawmatlabdir not in MATLABPATH:
    MATLABPATH = MATLABPATH +":"+ clawmatlabdir

# Possible platforms, 
#  win32 - Windows
#  cygwin - Windows with cygwin
#  linux - linux2 <- supported
#  sunos5 - Sun OS
#  darwin - Mac OS X <- supported
#  os2 - OS/2
#  os2emx - OS/2 EMX
#  RiscOS - riscos
#  atheos - AtheOS
if sys.platform.lower() == 'linux2':
    dylib_string = "LD_LIBRARY_PATH"
elif sys.platform.lower() == 'darwin':
    dylib_string = "DYLD_LIBRARY_PATH"
else:
    raise Exception("Unsupported system type %s" % sys.platform.lower())
dylib_path = os.path.join(clawdir,'lib')

print "Full path to claw directory should be:"
print "      $CLAW = ",clawdir
    

setenvcsh = open("setenv.csh","w")
setenvcsh.write("setenv CLAW '%s'\n" % CLAW)
setenvcsh.write("setenv FC '%s'\n\n" % FC)
setenvcsh.write("setenv MATLABPATH '%s'\n\n" % MATLABPATH)
setenvcsh.write("setenv PYTHONPATH '%s'\n" % PYTHONPATH)
setenvcsh.write("setenv IPYTHONDIR '%s'\n" % IPYTHONDIR)
setenvcsh.write('if ($?%s == 0) then\n' % dylib_string)
setenvcsh.write('    setenv %s "%s"\n' % (dylib_string,dylib_path))
setenvcsh.write('else\n')
setenvcsh.write('    setenv %s "%s:$%s"\n' % (dylib_string,dylib_path,dylib_string))
setenvcsh.write('endif\n')
setenvcsh.write("alias ipyclaw 'ipython -profile claw' \n")
setenvcsh.write("alias clawserver 'xterm -e python python/startserver.py &' \n")
setenvcsh.close()

setenvbash = open("setenv.bash","w")
setenvbash.write("export CLAW='%s'\n" % CLAW)
setenvbash.write("export FC='%s'\n\n" % FC)
setenvbash.write("export MATLABPATH='%s'\n\n" % MATLABPATH)
setenvbash.write("export PYTHONPATH='%s'\n" % PYTHONPATH)
setenvbash.write("export IPYTHONDIR='%s'\n" % IPYTHONDIR)
setenvbash.write('if [ -z "${%s}" ]; then\n' % (dylib_string))
setenvbash.write('    %s="%s"\n' % (dylib_string,dylib_path))
setenvbash.write('else\n')
setenvbash.write('    %s="%s:${%s}"\n' % (dylib_string,dylib_path,dylib_string))
setenvbash.write('fi\n')
setenvbash.write("alias ipyclaw='ipython -profile claw' \n")
setenvbash.write("alias clawserver='xterm -e python python/startserver.py &' \n")
setenvbash.close()

print "------------------------------------------------------------"
print "The files setenv.csh and setenv.bash contain the appropriate"
print "commands to set environment variables for csh or bash shells"
print "  and also some aliases you may find convenient             "
print "------------------------------------------------------------"
print " "


