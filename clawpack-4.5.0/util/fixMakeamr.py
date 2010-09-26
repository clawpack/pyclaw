#!/usr/bin/env python
#
# Fix amrclaw Makefile from 4.2 to work with Version 4.3 of CLAWPACK.
# Deletes the references to errsp.f and errsp.o and adds in new routines
# flag2refine.f, allowflag.f, bufnst.f, spest.f. 
# Also adds a "make clean" option if this doesn't exist.
#
# Works in 2d or 3d amrclaw applications.
#
# Writes output to Makefile, and renames original Makefile.original.

import sys,string,os

os.rename('Makefile', 'Makefile.original')
ifile = open('Makefile.original','r')
ofile = open('Makefile','w')

nclean = 0
for line in ifile:
    if string.find(line,'clean:') > -1:
        nclean = 1
    if string.find(line,'DO NOT remove') > -1:
        if nclean==0:
            # add a make clean option
            ofile.write('clean:\n')
            ofile.write('\t-rm -f $(OBJECTS) xamr xamrhdf\n')
            ofile.write('\n')

    if string.find(line,'$(CLAW)/amrclaw/2d/lib/errsp.o') > -1:
        print 'replacing amrclaw/2d/lib/errsp.o'
        ofile.write('  $(CLAW)/amrclaw/2d/lib/flag2refine.o \\\n')
        ofile.write('  $(CLAW)/amrclaw/2d/lib/allowflag.o \\\n')
        ofile.write('  $(CLAW)/amrclaw/2d/lib/bufnst.o \\\n')
        ofile.write('  $(CLAW)/amrclaw/2d/lib/spest.o \\\n')
    elif string.find(line,'$(CLAW)/amrclaw/2d/lib/errsp.f') > -1:
        print 'replacing amrclaw/2d/lib/errsp.f'
        ofile.write('  $(CLAW)/amrclaw/2d/lib/flag2refine.f \\\n')
        ofile.write('  $(CLAW)/amrclaw/2d/lib/allowflag.f \\\n')
        ofile.write('  $(CLAW)/amrclaw/2d/lib/bufnst.f \\\n')
        ofile.write('  $(CLAW)/amrclaw/2d/lib/spest.f \\\n')
    elif string.find(line,'$(CLAW)/amrclaw/3d/lib/errsp.o') > -1:
        print 'replacing amrclaw/3d/lib/errsp.o'
        ofile.write('  $(CLAW)/amrclaw/3d/lib/flag2refine.o \\\n')
        ofile.write('  $(CLAW)/amrclaw/3d/lib/allowflag.o \\\n')
        ofile.write('  $(CLAW)/amrclaw/3d/lib/bufnst.o \\\n')
        ofile.write('  $(CLAW)/amrclaw/3d/lib/spest.o \\\n')
    elif string.find(line,'$(CLAW)/amrclaw/3d/lib/errsp.f') > -1:
        print 'replacing amrclaw/3d/lib/errsp.f'
        ofile.write('  $(CLAW)/amrclaw/3d/lib/flag2refine.f \\\n')
        ofile.write('  $(CLAW)/amrclaw/3d/lib/allowflag.f \\\n')
        ofile.write('  $(CLAW)/amrclaw/3d/lib/bufnst.f \\\n')
        ofile.write('  $(CLAW)/amrclaw/3d/lib/spest.f \\\n')
    else:
        ofile.write('%s' % line)

print 'new version is now in Makefile (old version in Makefile.original)'

ifile.close()
ofile.close()

