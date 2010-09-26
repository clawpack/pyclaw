
"""
Create plots corresponding to each sample setplot function.  Search for all
files of the form setplot_*.py and loop over them.

Also create .rst files for each example.  The doc string for each setplot
file should start with at title underlined with ===, followed by a brief
description.  These are used in the rst file, which also includes the
setplot function itself and a pointer to the plots directory.
"""

def makeplots(filenames=[]):
    import os, glob, re
    from pyclaw.plotters.plotclaw import plotclaw

    thisdir = os.path.split(os.getcwd())[1]

    #os.system('make .plots')   # ensure output files and sample plots exist
    #os.system('make .htmls')   # ensure html files exist

    if filenames==[]:
        filenames = glob.glob('setplot_*.py')

    spnames = []
    dirname = os.path.split(os.getcwd())[1]
    rstname = "plot%s-" % dirname

    for setplotfile in filenames:
        print '=== Making plots using ',setplotfile 

        regexp = re.compile(r'setplot_(?P<spname>.*).py')
        result = regexp.search(setplotfile)
        spname = result.group('spname')
        spnames.append(spname)

        plotdir = 'plots_%s' % spname
        plotclaw(outdir="_output", plotdir=plotdir, setplot=setplotfile)

    for spname in spnames:
        setplotfile = 'setplot_%s.py' % spname
        rstfile_name = rstname + spname
        print '=== Making rst file %s.rst' % rstfile_name
        rstfile = open('../%s.rst' % rstfile_name, 'w')
        setplot_lines = open(setplotfile,'r').read()

        regexp = re.compile(r'"""(?P<descr>.*?)""" (?P<rest>.*)', \
                 re.DOTALL)
        result = regexp.search(setplot_lines)
        setplot_descr = result.group('descr')
        setplot_rest = result.group('rest')
        setplot_rest = setplot_rest.replace('\n','\n    ',1000)
        rstfile.write(""".. _%s: \n%s \n\n""" % (rstfile_name, setplot_descr))
	rstfile.write("Example generating data: `$CLAW/doc/sphinx/%s/README.html <../%s/README.html>`_\n\n" \
			% (thisdir, thisdir))
	rstfile.write("Resulting plots: `$CLAW/doc/sphinx/%s/plots_%s/_PlotIndex.html <../%s/plots_%s/_PlotIndex.html>`_\n\n::\n" \
			% (thisdir, spname, thisdir, spname))
        rstfile.write(setplot_rest)
        rstfile.close()
        os.system('chmod og+r ../%s.rst' % rstfile_name)

if __name__=='__main__':
    import sys
    makeplots(sys.argv[1:])
    
