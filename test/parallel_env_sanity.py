

def parallel_env_sanity(np =4):
    from petsc4py import PETSc
    size = PETSc.Comm.getSize(PETSc.COMM_WORLD)
    rank = PETSc.Comm.getRank(PETSc.COMM_WORLD)
    if rank == 0:
        print size


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from pyclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        parallel_env_sanity(*args,**kwargs)
    else: parallel_env_sanity()

