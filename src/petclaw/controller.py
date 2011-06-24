"""
Module for PetClaw controller class.  The PetClaw controller is identical to the
PyClaw controller except for the default value of output_format.
"""

from pyclaw.controller import Controller as pyclawController

class Controller(pyclawController):
    def __init__(self):
        super(Controller,self).__init__()

        self.output_format = 'petsc'


    def set_gauges(self):
        import os
        from numpy import floor
        if not os.path.exists(self.gauge_path):
            try:
                os.makedirs(self.gauge_path)
            except OSError:
                print "gauge directory already exists, ignoring"
        mygauges=[] #List of gauges: for each processor that owns a gauge
        gauge_files=[]  #List of files: for each processor that owns a gauge
        
        grid = self.solution.state.grid
        gauge_ind = [0]*self.solution.ndim

        for gauge in self.gauges: 
            is_mine = True
            for n in xrange(self.solution.ndim): 
                # Determine gauge locations in units of grid spacing
                gauge_ind[n] = int(floor(gauge[n]/grid.d[n]))
                istart = grid.dimensions[n].nstart
                iend   = grid.dimensions[n].nend
                if not (istart <= gauge_ind[n] <= iend):
                    is_mine = False
                else:
                    gauge_ind[n] = gauge_ind[n] - istart
            if is_mine:
                # append the gauge to mygauges and create a file for it
                # mygauges and gauge_files  are Lists because 
                #      one processor may have several gauges
                mygauges.append(list(gauge_ind))
                gauge_path = self.gauge_path+'gauge'+'_'.join(str(coord) for coord in gauge)+'.txt'
                # Is this a good idea???  It should depend on self.overwrite.
                if os.path.isfile(gauge_path): os.remove(gauge_path)
                gauge_files.append(open(gauge_path,'a'))

        self.solver.gauges = mygauges #pass to the solver
        self.solver.gauge_files = gauge_files
        self.solver.gauge_pfunction = self.gauge_pfunction
        self.solver.write_gauges(self.solution)

    def write_F(self):
        if self.compute_F is not None:
            self.compute_F(self.solution.state)
            F = [0]*self.solution.state.mF
            for i in xrange(self.solution.state.mF):
                #Sum across all processors
                #Warning: this actually computes the absolute sum!
                #Can this line be moved to after the if?
                F[i] = self.solution.state.gFVec.strideNorm(i,0)

            from petsc4py import PETSc
            rank = PETSc.Comm.getRank(PETSc.COMM_WORLD)
            if rank == 0:
                F_file = open(self.F_path,'a')
                F_file.write(' '.join(str(j) for j in F) + '\n')
                F_file.close()
