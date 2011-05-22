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
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#      Copyright (C) 2009 Randall J. LeVeque <rjl@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import logging
import sys
import os
import copy
import shutil
import time

from data import Data
from solution import Solution
from evolve.solver import Solver
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
        r"""(:class:`~pyclaw.evolve.solver.Solver`) - Solver object"""
        
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
        self.t0 = 0.0
        r"""(float) - Starting time, ``default = 0.0``"""
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
    def plotclaw(self, datadir='.'):
        pydir = '/home/rjl/claw/trunk/claw/python'
        if sys.platform in ['cygwin', 'win32']:
            syscmd = " C:/Python25/python.exe C:/cygwin%s/pyclaw/plotclaw.py  %s" \
                   % (pydir, datadir) 
        else:
            syscmd = " python %s/pyclaw/plotclaw.py  %s" \
                   % (pydir, datadir) 
        os.system(syscmd)
    
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
        if not reduce(lambda x,y: x*y,[sol.is_valid() for sol in 
                        self.solutions.itervalues()]):
            raise Exception("Initial solutions are not valid.")
        
        # Output styles
        if self.outstyle == 1:
            output_times = np.linspace(self.t0,
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
            self.solutions['n'].write(0,self.outdir,
                                        self.output_format,
                                        self.output_file_prefix,
                                        self.write_aux_init,
                                        self.output_options)                            
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
            # Save current solution to dictionary with frame as key
            if self.keep_copy:
                self.frames.append(copy.deepcopy(self.solutions['n']))
            if self.output_format is not None:
                self.solutions['n'].write(frame,self.outdir,
                                            self.output_format,
                                            self.output_file_prefix,
                                            self.write_aux_always,
                                            self.output_options)

            logging.info("Solution %s computed for time t=%f"
                % (frame,self.solutions['n'].t))
            
        self.solver.teardown()
        # Return the current status of the solver
        return status
    
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
        
    #------------------------------------------------------------
    def runxclaw(self, verbose=True):
    #------------------------------------------------------------
        r'''
        Run the command xdir/xclawcmd, directing the output fort.*
        files to outdir, writing unit 6 timestepping info to file xclawout.
        Runtime error messages are written to file xclawerr.
        If xclawout(xclawerr) is None, then output to stdout(stderr).
    
        If savecode==True, archive a copy of the code into directory outdir.
    
        This function returns the returncode from the process running xclawcmd,
        which will be nonzero if a runtime error occurs.
        '''
        import os, shutil, glob
    
        try:
            import subprocess
        except:
            print 'Must use a more recent version of Python with module subprocess'
            print 'Python 2.5 is recommended'
            raise Exception('subprocess module missing')
            return
    
        debug = False
        
        try:
            # the following attributes are needed:
            xdir = self.xdir
            rundir = self.rundir
            outdir = self.outdir
            overwrite = self.overwrite
            xclawcmd = self.xclawcmd
            xclawout = self.xclawout
            xclawerr = self.xclawerr
            runmake = self.runmake
            savecode = self.savecode
        except:
            raise Exception('missing attributes in runxclaw')
            return

        startdir = os.getcwd()
        xdir = os.path.abspath(xdir)
        outdir = os.path.abspath(outdir)
        rundir = os.path.abspath(rundir)
        xclawcmd = os.path.join(xdir,xclawcmd)
    
        #import pdb; pdb.set_trace()
    
        try:
            os.chdir(xdir)
        except:
            raise Exception( "Cannot change to directory xdir = %s" %xdir)
            return 
    
        if runmake:
            try:
                pass
                #os.system('make')
            except:
                print 'Warning: no make file in directory xdir = ',xdir
    

    
        if os.path.isfile(outdir):
            print "Error: outdir specified is a file"
            return
        
        if (os.path.isdir(outdir) & (not overwrite)):
            # move the old outdir instead of just clobbering it:
            tm = time.localtime(os.path.getmtime(outdir))
            year = str(tm[0]).zfill(4)
            month = str(tm[1]).zfill(2)
            day = str(tm[2]).zfill(2)
            hour = str(tm[3]).zfill(2)
            minute = str(tm[4]).zfill(2)
            second = str(tm[5]).zfill(2)
            outdir_backup = outdir + '_%s%s%s-%s%s%s' \
                  % (year,month,day,hour,minute,second)
            if verbose:
                print "*** Directory already exists: ",os.path.split(outdir)[1]
                print "*** Moving directory to:      ",os.path.split(outdir_backup)[1]
                time.sleep(1)
            
            shutil.move(outdir,outdir_backup)
            
            
        if (not os.path.isdir(outdir)):
            try:
                os.mkdir(outdir)
            except:
                print "Cannot make directory ",outdir
                return
    
        try:
            os.chdir(outdir)
        except:
            print '*** Error in runxclaw: cannot move to outdir = ',\
                  outdir
            raise
            return
    
        fortfiles = glob.glob(os.path.join(outdir,'fort.*'))
        if overwrite:
            # remove any old versions:
            for file in fortfiles:
                os.remove(file)
        else:
            if len(fortfiles) > 1:
                print "*** Remove fort.* and try again,"
                print "  or use overwrite=True in call to runxclaw"
                return
            
        
        try:
            os.chdir(rundir)
        except:
            raise Exception("Cannot change to directory %s" % rundir)
            return 
    
        datafiles = glob.glob('*.data')
        if datafiles == ():
            print "Warning: no data files found in directory ",rundir
        else:
            if rundir != outdir:
                for file in datafiles:
                    shutil.copy(file,os.path.join(outdir,file))
    
        if xclawout:
            xclawout = open(xclawout,'wb')
        if xclawerr:
            xclawerr = open(xclawerr,'wb')
    
        os.chdir(outdir)
    
        #print "\nIn directory outdir = ",outdir,"\n"
    
        # execute command to run fortran program:
    
        try:
            #print "\nExecuting ",xclawcmd, "  ...  "
            #pclaw = subprocess.Popen(xclawcmd,stdout=xclawout,stderr=xclawerr)
            #print '+++ pclaw started'
                #pclaw.wait()   # wait for code to run
            #returncode = pclaw.returncode
            #print '+++ pclaw done'
            
            returncode = os.system(xclawcmd)
    
            if returncode == 0:
                print "\nFinished executing\n"
            else:
                print "\n *** Runtime error: return code = %s\n " % returncode
        except:
            raise Exception("Could not execute command %s" % xclawcmd)
            return
    
        os.chdir(startdir)
        return returncode
