

begin_html [use: jsMath] [use: doc/doc.css]

<!----------------------------------------------------------------      
     For a more readable version of this file, execute
                  unix>  make htmls
     in this directory and then point your browser to README.html
--------------------------------------------------------------------->

<h2>
AMRCLAW Sample Code
</h2>


<h2>Description</h2>
<p>

Advection at constant velocity with periodic boundary conditions.
The velocity is specified in [code:setrun.py].

This example illustrates how to set gauge locations (see [code:setrun.py])
and plot the results as a function of time at each gauge 
(see [code: setplot.py]).

For more details, see the parent directory for the code without AMR:
[link: ../README.html]


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

<dt>Other Fortran files are either up one level (see [link: ../README.html])
or in the  AMRCLAW library [claw:amrclaw/2d/lib].

</dl>

<p>
<h4>
Data files
</h4>

<p>
<dl>
<dt>[code: amr2ez.data]
<dd> 
Standard parameter values that are read by library routine 
[clawcode:amrclaw/2d/lib/amr2ez.f]
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
    

</dl>
<h4>
Python files
</h4>
<dl>

<dt>[code: setrun.py]
<dd> This file contains a function that
specifies what run-time parameters will be used.

<dt>[code: setplot.py]
<dd> This file contains a function that
specifies what plots will be done and
sets various plotting parameters.


</dl>

<p>
<h3>Library routines</h3>
<p>

In addition to the Fortran routines in this library, several library
routines from [claw:amrclaw/2d/lib] are used.  See the [code: Makefile]
to determine which ones are used.



end_html
