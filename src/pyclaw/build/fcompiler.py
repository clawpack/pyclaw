def get_fcompiler():
    import numpy.distutils.fcompiler
    numpy.distutils.log.set_verbosity(-1)
    fc = numpy.distutils.fcompiler.new_fcompiler()
    fc.customize()
    return fc

if __name__ == "__main__":
    fc = get_fcompiler()
    import sys
    if sys.argv[1] == 'get_compiler':
        print(fc.compiler_f77[0])
    elif sys.argv[1] == 'get_flags':
        print(' '.join(fc.compiler_f77[1:]))
