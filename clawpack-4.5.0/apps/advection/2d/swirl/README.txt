
begin_html [use: jsMath]  [use: doc/doc.css]

<!----------------------------------------------------------------      
     For a more readable version of this file, execute
                  unix>  make htmls
     in this directory and then point your browser to README.html
--------------------------------------------------------------------->

<h2>
CLAWPACK Sample Code
</h2>

<h4>Description</h4>
<p>


Swirling flow in a box

Edge velocities are stored in aux array (see [code:setaux.f]) with
velocity specified by differencing the streamfunction psi.f

The velocities are time-dependent giving reversing flow.  These
velocities are computed in [code:b4step2.f].   The period is specified
by the parameter tperiod in setprob.data.   

In theory the solution should agree with the initial data at times 
t = N*tperiod/2 for all integers N, but because of the numerical 
diffusion this won't happen.

As a special case, if tperiod = 0, then 
the velocities are constant in time and b4step2 does nothing.
The velocities specified in setaux.f are then used at all times.

An example using AMRCLAW is also included in directory 
[link: amr/README.html ./amr]

This test problem was used in:

<pre>
@article{rjl:advect,
  author="R. J. LeVeque",
  title="High-resolution conservative algorithms for advection in
  incompressible flow",
  journal="SIAM J. Numer. Anal.",
  volume="33",
  year="1996",
  pages="627-665"
}
</pre>

See [http://www.amath.washington.edu/~rjl/pubs/hiresadv]

See [link: #instructions Instructions]

<h4>
Plots of results
</h4>
After running this code and creating plots via "make .plots", you should be
able to view the plots in [link: _plots/_PlotIndex.html].



<h4>
Fortran files
</h4>


<p>
<dl>
<dt>[code: Makefile]
<dd> Determines which version of fortran files
are used when compiling the code with make.


<p>
<dt>[code: qinit.f95]
<dd>
Sets initial conditions.

<p>
<dt>[code: setprob.f]
<dd>
A routine by this name is called by the library routine
[clawcode: clawpack/1d/lib/claw1ez.f]
and is generally used to set any values needed for the specific problem
being solved.  In this example, the value of the advection velocity $u$  and
a parameter $\beta$ used in setting the initial conditions are
read from the file [code: setprob.data].


<dt>[code: setaux.f]
<dd> Set values in the aux arrays.

<dt>[code: psi.f]
<dd> Stream function for advection velocities.

<dt>[code: b4step2.f]
<dd> Called before each time step.  Used here to make advection velocity time
dependent.

<p>
<dt><b>Riemann solvers:</b>

<p>
<dt>[code: rpn2ad1.f]
<dd> Normal Riemann solver (normal to cell interface).

<dt>[code: rpt2ad1.f]
<dd> Transverse Riemann solver.

</dl>

<p>
<h4>
Data files
</h4>
<p>
<dt>[code: claw2ez.data]
<dd> 
Standard parameter values that are read by library routine 
[clawcode:clawpack/2d/lib/claw2ez.f]
Each line contains one or more values to be read
in, followed by comments that are ignored by the Fortran code but
used by the Python read or write methods of class clawtools.ClawData.

Some parameters that you might want to modify are described in the
[http://kingkong.amath.washington.edu/clawpack/users documentation].


<p>
<dt>  [code: setprob.data]
<dd> 
Problem-specific parameters.

<p>
    
<dt> [code: setplot.data]
<dd> This file contains plotting parameters.

</dl>
<h4>
Python files
</h4>
<dl>
<dt>[code:mapc2p.py]
<dd> Maps the computational rectangular domain to the physical grid.

<dt> [code: setplot.py]
<dd> This file contains commands that are executed 
on start-up, and sets various plotting parameters. 

</dl>


<p>
<h3>Library routines</h3>
<p>

In addition to the Fortran routines in this library, several library
routines from [claw:clawpack/2d/lib] are used.  See the [code: Makefile]
to determine which ones are used.

[name:  instructions]
<h4>Instructions</h4> 


To run code:

{{{
    $ make .output
}}}

View plots interactively with Iplotclaw or use "make .plots" to create html
files.

end_html
