#!/usr/bin/env python
# encoding: utf-8
r"""
Controller for basic computation and plotting setup

This module defines the Pyclaw controller class.  It can be used to perform
simulations similar to previous versions of clawpack, i.e. with outstyle and
output time specification.  It also can be used to setup easy plotting and 
running of compiled fortran binaries.

:Authors:
    Kyle T. Mandli (2008-02-15) Initial version
    
    Randall J. LeVeque and Kyle T Mandli (2009) Plotting and run updates

    David I. Ketcheson (2011) Minor additional functionality
"""

import logging
import sys
import os
import copy
import shutil
import time

from data import Data
from solution import Solution
from solver import Solver
from util import FrameCounter

class Controller(object):
    r"""Controller for pyclaw simulation runs and plotting
            
    :Initialization:
    
        Input: None
    
    :Version: 1.0 (2009-06-01)
    """
    #  ======================================================================
    #   Initialization routines
    #  ======================================================================
    def __init__(self):
        r"""
        Initialization routine for a Controller object.
        
        See :class:`Controller` for full documentation.
        """
        
        import numpy as np

        self.viewable_attributes = ['xdir','rundir','outdir','overwrite',
                        'xclawcmd','xclawout','xclawerr','runmake','savecode',
                        'solver','keep_copy','write_aux_init',
                        'write_aux_always','output_format',
                        'output_file_prefix','output_options','nout',
                        'outstyle','verbosity']
        r"""(list) - Viewable attributes of the `:class:`~pyclaw.controller.Controller`"""

        # Global information for running and/or plotting
        self.xdir = os.getcwd()
        r"""(string) - Executable path, executes xclawcmd in xdir"""
        self.rundir = os.getcwd()
        r"""(string) - Directory to run from (containing \*.data files), uses 
        \*.data from rundir"""
        self.outdir = os.getcwd()+'/_output'
        r"""(string) - Output directory, directs output files to outdir"""
        self.overwrite = True
        r"""(bool) - Ok to overwrite old result in outdir, ``default = True``"""

        self.xclawcmd = 'xclaw'
        r"""(string) - Command to execute (if using fortran), defaults to xclaw or
        xclaw.exe if cygwin is being used (which it checks vis sys.platform)"""
        if sys.platform == 'cygwin':
             self.xclawcmd = 'xclaw.exe'

        self.start_frame = 0
        self.xclawout = None
        r"""(string) - Where to write timestep messages"""
        self.xclawerr = None
        r"""(string) - Where to write error messages"""
        self.runmake = False
        r"""(bool) - Run make in xdir before xclawcmd"""
        self.savecode = False
        r"""(bool) - Save a copy of \*.f files in outdir"""
        
        # Solver information
        self.solutions = {}              # Solutions dictionary
        self.solver = None
        r"""(:class:`~pyclaw.solver.Solver`) - Solver object"""
        
        # Output parameters for run convenience method
        self.keep_copy = False 
        r"""(bool) - Keep a copy in memory of every output time, 
        ``default = False``"""
        self.frames = []
        r"""(list) - List of saved frames if ``keep_copy`` is set to ``True``"""
        self.write_aux_init = False
        r"""(bool) - Write out initial auxiliary array, ``default = False``"""
        self.write_aux_always = False
        r"""(bool) - Write out auxiliary array at every time step, 
        ``default = False``"""
        self.output_format = 'ascii'
        r"""(list of strings) - Format or list of formats to output the data, 
        if this is None, no output is performed.  See _pyclaw_io for more info
        on available formats.  ``default = 'ascii'``"""
        self.output_file_prefix = None
        r"""(string) - File prefix to be appended to output files, 
        ``default = None``"""
        self.output_options = {}
        r"""(dict) - Output options passed to function writing and reading 
        data in output_format's format.  ``default = {}``"""
        
        # Classic output parameters, used in run convenience method
        self.tfinal = 1.0
        r"""(float) - Final time output, ``default = 1.0``"""
        self.outstyle = 1
        r"""(int) - Time output style, ``default = 1``"""
        self.verbosity = 0 
        r"""(int) - Level of output, ``default = 0``"""
        self.nout = 10                  # Outstyle 1 defaults
        r"""(int) - Number of output times, only used with ``outstyle = 1``,
        ``default = 10``"""
        self.out_times = np.linspace(0.0,self.tfinal,self.nout) # Outstyle 2
        r"""(int) - Output time list, only used with ``outstyle = 2``,
        ``default = numpy.linspace(0.0,tfinal,nout)``"""
        
        self.nstepout = 1               # Outstyle 3 defaults
        r"""(int) - Number of steps between output, only used with 
        ``outstyle = 3``, ``default = 1``"""
        
        # Data objects
        self.plotdata = None
        r"""(:class:`~pyclaw.plotters.data.ClawPlotData`) - An instance of a 
        :class:`~pyclaw.plotters.data.ClawPlotData` object defining the 
        objects plot parameters."""
        
        # Derived quantity p
        self.file_prefix_p = 'claw_p'
        r"""(string) - File prefix to be prepended to derived quantity output files"""
        self.compute_p = None
        r"""(function) - function that computes derived quantities"""
        self.outdir_p = self.outdir+'/_p'
        r"""(string) - Directory to use for writing derived quantity files"""

        # functionals
        self.compute_F = None
        r"""(function) - Function that computes F"""
        self.F_file_name = 'F'
        r"""(string) - Name of text file containing functionals"""
        self.F_path = './_output/'+self.F_file_name+'.txt'
        r"""(string) - Full path to output file for functionals"""

        # gauges
        self.gauges = []
        r"""(list) - List of gauges"""
        self.gauge_path = './_output/_gauges/'
        r"""(string) - Full path to output directory for gauges"""
        self.gauge_pfunction = None
        r"""(function) - Function that computes quantities to be recorded at gaugues"""

    # ========== Access methods ===============================================
    def __str__(self):        
        output = "Controller attributes:\n"
        for attr in self.viewable_attributes:
            value = getattr(self,attr)
            output = output + "  %s = %s \n" % (attr,value)
        output = output + '\n'
        if self.plotdata is not None:
            output = output + "  Data "+str(self.plotdata)+"\n"
        if self.solver is not None:
            output = output + "  Solver "+str(self.solver)+"\n"
        if len(self.frames) > 0:
            output = output + "  Frames \n"
            for frame in self.frames:
                output = output + "    " + str(frame) + "\n"
        return output
        
    # ========== Properties ==================================================
    def solution():
        def fget(self): return self.solutions['n']
        def fset(self, value): self.solutions['n'] = value
        return locals()
    solution = property(**solution())
    
    # ========== Plotting methods ============================================        
    # ========== Solver convenience methods ==================================
    def run(self):
        r"""
        Convenience routine that will evolve solutions['n'] based on the 
        traditional clawpack output and run parameters.
        
        This function uses the run parameters and solver parameters to evolve
        the solution to the end time specified in run_data, outputting at the
        appropriate times.
        
        :Input:
            None
            
        :Ouput:
            (dict) - Return a dictionary of the status of the solver.
            
        :Version: 1.0 (2009-05-01)
        """
        
        import numpy as np

        frame = FrameCounter()
        frame.set_counter(self.start_frame)
        if self.keep_copy:
            self.frames = []
                    
        # Check to make sure we have a valid solver to use
        if self.solver is None:
            raise Exception("No solver set in controller.")
        if not isinstance(self.solver,Solver):
            raise Exception("Solver is not of correct type.")
        if not self.solver.is_valid():
            raise Exception("The solver failed to initialize properly.") 
            
        # Call solver's setup routine
        self.solver.setup(self.solutions)
            
        # Check to make sure the initial solutions are valid
        if not all([sol.is_valid() for sol in self.solutions.itervalues()]):
            raise Exception("Initial solutions are not valid.")
        
        # Assign gauges and determine their grid indices
        if self.gauges is not None:
            self.set_gauges()

        # Output styles
        # We ought to have a check confirming that all solutions
        # are at the same time t
        if self.outstyle == 1:
            output_times = np.linspace(self.solution.t,
                    self.tfinal,self.nout+1)
        elif self.outstyle == 2:
            output_times = self.out_times
        elif self.outstyle == 3:
            output_times = np.ones((self.nout+1))
        else:
            raise Exception("Invalid output style %s" % self.outstyle)  
         
        # Output and save initial time
        if self.keep_copy:
            self.frames.append(copy.deepcopy(self.solutions['n']))
        if self.output_format is not None:
            if self.compute_p is not None:
                self.compute_p(self.solution.state)
                self.solutions['n'].write(0,self.outdir_p,
                                        self.output_format,
                                        self.file_prefix_p,
                                        write_aux = False,
                                        options = self.output_options,
                                        write_p = True) 

            self.solutions['n'].write(0,self.outdir,
                                        self.output_format,
                                        self.output_file_prefix,
                                        self.write_aux_init,
                                        self.output_options)

        if self.compute_F is not None:
            from petsc4py import PETSc
            rank = PETSc.Comm.getRank(PETSc.COMM_WORLD)
            self.compute_F(self.solution.state)
            F = np.zeros(self.solution.state.mF)
            for i in xrange(self.solution.state.mF):
                #Sum across all processors
                #Warning: this actually computes the absolute sum!
                #Can this line be moved to after the if?
                F[i] = self.solution.state.gFVec.strideNorm(i,0)
            if rank == 0:
                F_file = open(self.F_path,'w')
                F_file.write(' '.join(str(j) for j in F) + '\n')
                F_file.close()


        logging.info("Solution %s computed for time t=%f" % 
                        (0,self.solutions['n'].t) )

        for t in output_times[1:]:                
            if self.outstyle < 3:
                status = self.solver.evolve_to_time(self.solutions,t)
            else:
                # Take nstepout steps and output
                for n in xrange(self.nstepout):
                    status = self.solver.evolve_to_time(self.solutions)
            frame.increment()
            if self.keep_copy:
                # Save current solution to dictionary with frame as key
                self.frames.append(copy.deepcopy(self.solutions['n']))
            if self.output_format is not None:
                if self.compute_p is not None:
                    self.compute_p(self.solution.state)
                    self.solutions['n'].write(frame,self.outdir_p,
                                            self.output_format,
                                            self.file_prefix_p,
                                            write_aux = False, 
                                            options = self.output_options,
                                            write_p = True) 
                
                self.solutions['n'].write(frame,self.outdir,
                                            self.output_format,
                                            self.output_file_prefix,
                                            self.write_aux_always,
                                            self.output_options)
            if self.compute_F is not None:
                rank = PETSc.Comm.getRank(PETSc.COMM_WORLD) 
                self.compute_F(self.solution.state) 
                F = np.zeros(self.solution.state.mF)
                for i in xrange(self.solution.state.mF):
                    F[i] = self.solution.state.gFVec.strideNorm(i,0)
                if rank == 0:
                    F_file = open(self.F_path,'a')
                    F_file.write(' '.join(str(j) for j in F) + '\n')
                    F_file.close()

            logging.info("Solution %s computed for time t=%f"
                % (frame,self.solutions['n'].t))
            
        self.solver.teardown()
        for file in self.solver.gauge_files: file.close()

        # Return the current status of the solver
        return status
    
    def set_gauges(self):
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
            for n in xrange(self.solution.ndim): 
                # Determine gauge locations in units of grid spacing
                gauge_ind[n] = int(floor(gauge[n]/grid.d[n]))
                istart = grid.dimensions[n].nstart
                iend   = grid.dimensions[n].nend
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


    # ========== Output Data object based on solver and solutions['n'] =======
    def get_data(self,claw_path=None):
        r"""
        Create a data object from this controller's solver and solutions
        
        This function will take the current solver and solutions['n'] and
        create a data object that can be read in via classic clawpack.
        
        If claw_path is provided, then the data that should be written to the
        claw.data file will be written to that path.
        
        :Input:
            - *claw_path* - (string) Path to write data file to
            
        :Output:
            - (:class:`~pyclaw.data.Data`) - Data object claw_data containing 
              the appropriate data for a claw.data file.
        """
        
        # Check to make sure we have a valid solver and solution
        if not self.solver.is_valid() or not self.solution.is_valid():
            raise Exception("Invalid solver or solution.")
        
        claw_data = Data()
        
        claw_data.add_attribute('ndim',value=self.solution.ndim)
        claw_data.add_attribute('mx',value=self.solution.dimensions[0].n)
        if claw_data.ndim > 1:
            claw_data.add_attribute('my',value=self.solution.dimensions[1].n)
        if claw_data.ndim > 2:
            claw_data.add_attribute('mz',value=self.solution.dimensions[2].n)
            
        claw_data.add_attribute('nout',value=self.nout)
        claw_data.add_attribute('outstyle',value=self.outstyle)
        if claw_data.outstyle == 2:
            claw_data.add_attribute('out_times',value=self.out_times)
        elif claw_data.outstyle == 3:
            claw_data.add_attribute('nstepout',value=self.nstepout)
            
        claw_data.add_attribute('dt_initial',value=self.solver.dt)
        claw_data.add_attribute('dt_max',value=self.solver.dt_max)
        claw_data.add_attribute('cfl_max',value=self.solver.cfl_max)
        claw_data.add_attribute('cfl_desired',value=self.solver.cfl_desired)
        claw_data.add_attribute('max_steps',value=self.solver.max_steps)
        
        if self.solver.dt_variable:
            claw_data.add_attribute('dt_variable',value=1)
        else:
            claw_data.add_attribute('dt_variable',value=0)
        claw_data.add_attribute('order',value=self.solver.order)
        if claw_data.ndim == 1:
            claw_data.add_attribute('order_trans',value=0)
        else:
            claw_data.add_attribute('order_trans',value=self.solver.order_trans)
        claw_data.add_attribute('verbosity',value=self.verbosity)
        claw_data.add_attribute('src_split',value=self.solver.src_split)
        if self.solution.capa is not None:
            raise Exception("Not sure what to do here, they're are different!")
        else:
            claw_data.add_attribute('mcapa',value=0)
        claw_data.add_attribute('maux',value=self.solution.maux)
        
        claw_data.add_attribute('meqn',value=self.solution.meqn)
        # claw_data.add_attribute('mwaves',value=self.solution.mwaves)
        claw_data.add_attribute('mthlim',value=self.solver.mthlim)
        
        claw_data.add_attribute('t0',value=self.solution.t)
        claw_data.add_attribute('xlower',
                                    value=self.solution.dimensions[0].lower)
        claw_data.add_attribute('xupper',
                                    value=self.solution.dimensions[0].upper)
        if claw_data.ndim > 1:
            claw_data.add_attribute('ylower',
                                    value=self.solution.dimensions[1].lower)
            claw_data.add_attribute('yupper',
                                    value=self.solution.dimensions[1].upper)
        if claw_data.ndim > 2:
            claw_data.add_attribute('zlower',
                                    value=self.solution.dimensions[2].lower)
            claw_data.add_attribute('zupper',
                                    value=self.solution.dimensions[2].upper)
        
        claw_data.add_attribute('mbc',value=self.solution.mbc)
        claw_data.add_attribute('mthbc_xlower',
                                value=self.solution.dimensions[0].mthbc_lower)
        claw_data.add_attribute('mthbc_xupper',
                                value=self.solution.dimensions[0].mthbc_upper)
        if claw_data.ndim > 1:
            claw_data.add_attribute('mthbc_ylower',
                                value=self.solution.dimensions[1].mthbc_lower)
            claw_data.add_attribute('mthbc_yupper',
                                value=self.solution.dimensions[1].mthbc_upper)
        if claw_data.ndim > 2:
            claw_data.add_attribute('mthbc_zlower',
                                value=self.solution.dimensions[2].mthbc_lower)
            claw_data.add_attribute('mthbc_zupper',
                                value=self.solution.dimensions[2].mthbc_upper) 
            
        if claw_path is not None:
            # Write out this data object
            pass

        return claw_data
        
    def read_data(self,path):
        r"""Read in a claw.data file and initialize accordingly
        
        .. warning::
            
            Not implemented!
        """
        raise NotImplementedException()
