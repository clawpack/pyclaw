import pytest

@pytest.mark.skip(reason="Test not working yet PeanoClaw.")
def test_initialization():
    from pyclaw.clawpack.clawpack import ClawSolver2D
    
    solver = ClawSolver2D()
    
    import peanoclaw
    peano_solver = peanoclaw.Solver(solver, 1.0)
    
    import inspect
    for member in inspect.getmembers(peano_solver):
        if(not member[0].startswith("_") and not inspect.ismethod(member[1])):
            print((member[0]))
    for member in inspect.getmembers(solver):
        if(not member[0].startswith("_") and not inspect.ismethod(member[1])):
            print((member[0]))
    

# if __name__=="__main__":
#     import nose
#     nose.main()
