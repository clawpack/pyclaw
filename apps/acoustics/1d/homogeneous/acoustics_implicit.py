#!/usr/bin/env python
# encoding: utf-8
import petsc4py, sys
petsc4py.init(sys.argv)
from petsc4py import PETSc

class AcousticEquation:
    def __init__(self, grid, solver):
        from petclaw.evolve.sharpclaw import RKStageState
        self.grid = grid
        self.solver = solver
        self.wrapper = RKStageState(grid)
        return

    def rhsfunction(self, ts, t, x, f):
        self.wrapper.gqVec[:] = x[:]
        f[:] = self.solver.dqdt(self.grid,self.wrapper)

    def ifunction(self, ts, t, x, xdot, f):
        self.rhsfunction(ts,t,x,f)
        f.aypx(-1.0,xdot)

    def ijacobian(self, ts, t, X, Xdot, a, J, P):
        from assembly import formIJacobian
        cc = self.grid.aux_global['cc']
        zz = self.grid.aux_global['zz']
        xmin = self.grid.x.lower
        xmax = self.grid.x.upper
        return formIJacobian(ts, t, X, Xdot, a, J, P, cc, zz, xmin, xmax)
    
def acoustics(kernel_language='Python',petscPlot=False,iplot=False,htmlplot=False,myplot=False,outdir='./_output',sclaw=False,**kwargs):
    import numpy as np
    """
    1D acoustics example.
    """

    from petclaw.grid import Grid
    from petclaw.grid import Dimension
    from pyclaw.solution import Solution
    from pyclaw.controller import Controller
    from petclaw import plot

    petscts = kwargs.get('petscts', False)

    # Initialize grids and solutions
    xmin = 0.0
    xmax = 1.0
    ncells = kwargs.get('ncells',100)
    x = Dimension('x',xmin,xmax,ncells,mthbc_lower=2,mthbc_upper=2)
    grid = Grid(x)
    rho = 1.0
    bulk = 1.0
    grid.aux_global['rho']=rho
    grid.aux_global['bulk']=bulk
    grid.aux_global['zz']=np.sqrt(rho*bulk)
    grid.aux_global['cc']=np.sqrt(rho/bulk)
    from classic1 import cparam 
    for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)
    grid.meqn=2
    if sclaw:
        grid.mbc=3
    grid.t = 0.0

    # init_q_petsc_structures must be called 
    # before grid.x.center and such can be accessed.
    grid.init_q_petsc_structures()
    xc=grid.x.center
    q=np.zeros([grid.meqn,len(xc)], order = 'F')
    beta=100; gamma=0; x0=0.75
    q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))
    q[1,:]=0.
    grid.q=q
    
    init_solution = Solution(grid)

    if sclaw:
        from petclaw.evolve.sharpclaw import SharpClawSolver1D
        solver = SharpClawSolver1D()
        solver.lim_type = kwargs.get('lim_type',2)
        solver.time_integrator = 'SSP33'
        solver.mwaves=2
        solver.char_decomp=0
    else:
        from petclaw.evolve.clawpack import ClawSolver1D
        solver = ClawSolver1D()
        solver.mwaves=2
        from pyclaw.evolve import limiters
        solver.mthlim = limiters.MC

    if kernel_language=='Python': solver.set_riemann_solver('acoustics')
    solver.dt=grid.d[0]/grid.aux_global['cc']*0.1
    solver.kernel_language=kernel_language

    claw = Controller()
    #claw.outstyle = 3
    claw.keep_copy = True
    claw.nout = 5
    # The output format MUST be set to petsc!
    claw.output_format = 'ascii'
    claw.outdir = outdir
    claw.tfinal = 1.0
    claw.solutions['n'] = init_solution
    claw.solver = solver

    # Solve
    if petscts:
        # Uses PETSc time integrator. Needs sclaw=1 to get semidiscrete form. Examples:
        # ./acoustics.py sclaw=1 petscts=1 -ts_monitor_solution -ts_type theta -snes_mf_operator -ts_dt 0.02 -ts_max_steps 50 -ts_max_time 10 -draw_pause 0.05
        # ./acoustics.py sclaw=1 petscts=1 -ts_monitor_solution -ts_type ssp -ts_ssp_type rk104
        x = grid.q_da.createGlobalVector()
        f = x.duplicate()

        ts = PETSc.TS().create(PETSc.COMM_WORLD)
        ts.setDM(grid.q_da)
        ts.setProblemType(PETSc.TS.ProblemType.NONLINEAR)
        ts.setType(ts.Type.THETA)

        eqn = AcousticEquation(grid, solver)
        x[:] = grid.q
        ts.setRHSFunction(eqn.rhsfunction, f)
        ts.setIFunction(eqn.ifunction, f)
        ts.setIJacobian(eqn.ijacobian, grid.q_da.createMat())

        ts.setTimeStep(solver.dt)
        ts.setDuration(claw.tfinal, max_steps=1000)
        ts.setFromOptions()
        q0 = x[:].copy()
        ts.solve(x)
        qfinal = x[:].copy()
        dx = grid.d[0]
    else:
        status = claw.run()
        #This test is set up so that the waves pass through the domain
        #exactly once, and the final solution should be equal to the
        #initial condition.  Here we output the 1-norm of their difference.
        q0=claw.frames[0].grid.gqVec.getArray().reshape([-1])
        qfinal=claw.frames[claw.nout].grid.gqVec.getArray().reshape([-1])
        dx=claw.frames[0].grid.d[0]

    if htmlplot:  plot.plotHTML(outdir=outdir,format=output_format)
    if iplot:     plot.plotInteractive(outdir=outdir,format=output_format)
    if petscPlot: plot.plotPetsc(output_object)

    print 'Max error:', np.max(qfinal - q0)

    return dx*np.sum(np.abs(qfinal-q0))


if __name__=="__main__":
    import sys
    from petclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    error=acoustics(**kwargs)
    print 'Error: ',error
