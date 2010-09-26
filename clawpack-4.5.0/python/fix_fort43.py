
def fix_fortt(ndim=None):
    import glob, os
    if ndim==None:
        if os.path.exists('claw1ez.data'): 
            ndim = 1
        elif os.path.exists('claw2ez.data'): 
            ndim = 2
        elif os.path.exists('claw3ez.data'): 
            ndim = 3
        else:
            print 'No clawNez.data file, cannot deduce ndim'
            return
    fortfiles = glob.glob('fort.t*')
    for file in fortfiles:
        infile = open(file,'r')
        lines = infile.read()
        infile.close()

        if 'ndim' not in lines:
            lines = lines.replace('ngrids\n', \
                'ngrids\n    %s                 ndim\n' % ndim)

        outfile = open(file,'w')
        outfile.write(lines)
        outfile.close()


if __name__=='__main__':
    fix_fortt()
