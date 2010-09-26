
#
# Fix a set of target files in directory tree rootdir by replacing
# oldpat with newpat. 
#
# Now supports wildcards in list of targetfiles.
#

import os,sys,glob

rootdir = '..'
targetfiles = ['README.txt']

oldpat = "www.clawpack.org/doc.html"
newpat = "http://kingkong.amath.washington.edu/clawpack/users"

for (dirpath, subdirs, files) in os.walk(rootdir):
    currentdir = os.path.abspath(os.getcwd())
    os.chdir(os.path.abspath(dirpath))
    tfiles = []
    for fpat in targetfiles:
        for f in glob.glob(fpat):
            tfiles.append(f)
    for file in tfiles:

        infile = open(file,'r')
        lines = infile.read()
        infile.close()

        if lines.find(oldpat) > -1:
            lines = lines.replace(oldpat, newpat)
            print "Fixed file   ",dirpath + '/' + file
        else:
            print "No change to ",dirpath + '/' + file

        outfile = open(file,'w')
        outfile.write(lines)
        outfile.close()

    os.chdir(currentdir)

