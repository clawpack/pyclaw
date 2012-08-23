from nose.plugins.attrib import attr

if __name__=="__main__":
    import nose
    nose.main()
    
def qinit(state):
    pass
    
def test_initialization():
    from clawpack.pyclaw.classic.clawpack import ClawSolver2D
    
    solver = ClawSolver2D()
    
    from clawpack import peanoclaw
    peano_solver = peanoclaw.Solver(solver, 1.0, qinit)
    
    import inspect
    for member in inspect.getmembers(peano_solver):
        if(not member[0].startswith("_") and not inspect.ismethod(member[1])):
            print(member[0])
    for member in inspect.getmembers(solver):
        if(not member[0].startswith("_") and not inspect.ismethod(member[1])):
            print(member[0])
    
