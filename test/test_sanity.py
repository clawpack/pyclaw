from nose.plugins.attrib import attr

@attr(cat='sanity')
def test_sane():
    pass

@attr(petsc=True)
def test_parallel_env_sane():
    """
    Test whether PETSc/MPI are interacting with each other correctly
    """
    import subprocess
    
    # run subprocess with the required number of processes
    np = 4
    file_name = "./test/parallel_env_sanity.py" 
    run_command = "mpiexec", "-n", str(np) ,"python",file_name, "np="+str(np)
    
    p = subprocess.Popen(run_command, stdout=subprocess.PIPE ,stderr=subprocess.STDOUT)
    (stdout_data, ignore) = p.communicate()

    s = sum([int(i) for i in stdout_data.split()])    

    if s!=np:
        raise Exception('Unexpected size of PETSc.COMM_WORLD, size = '+stdout_data.strip()+' while number of processes = '+str(np)+'. Check your PETSc/MPI installation')
    else:
        pass
