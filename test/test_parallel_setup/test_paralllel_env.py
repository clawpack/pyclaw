from nose.plugins.attrib import attr

@attr(cat='sanity')
def test_sane():
    pass

def test_parallel_env_sane():
    """
    Test whether PETSc/MPI are interacting with each other correctly
    """
    import subprocess
    import json
    import os
    # run subprocess with the required number of processes
    try:
        from petsc4py import PETSc
    except ImportError:
        print "petsc4py.PETSc module is not found."
        PETSc = None
    
    if PETSc:
        np = 4
        file_name = os.path.join(os.path.dirname(os.path.abspath( __file__ )),
                                 'parallel_env_sanity.py')
        mpi_cmd = "mpiexec"
        run_command =str(mpi_cmd) , "-n", str(np) ,"python",file_name, "np="+str(np)
        
        p = subprocess.Popen(run_command, stdout=subprocess.PIPE ,stderr=subprocess.STDOUT)
        size = p.__dict__['stdout'].read()
    
        if size.rstrip().lstrip() == '':
            size=1
        else:
            try:
                size=int(size)
            except:
                raise Exception("size="+str(size)+" can not be parsed as integer.")   
        if size!=np:
            import procutils
            raise Exception('UNEXPECTED SIZE OF PETSc.COMM_WORLD: size = '+str(size)+
            ' while number of processes = '+str(np)+'. Your mpi library associated with the default '+
            mpi_cmd+' in ' +os.path.dirname(procutils.which(mpi_cmd)[0])+
            ' is probably different from the mpi library used for configuring your PETSc install.'+
            ' Fixing the path of the '+str(mpi_cmd)+
            ' command is necessary for using the feature of running pyclaw in parallel.')
        else:
            pass
