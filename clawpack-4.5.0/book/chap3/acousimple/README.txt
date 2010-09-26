begin_html  [use: doc/doc.css] [use: jsMath]
<!--   For a more readable version of this file, execute
                  unix>  make htmls
       in this directory and then point your browser to README.html -->

<h2>
Clawpack Sample Code
</h2>

Example [book/chap3/acousimple]
to accompany Figure 3.1 of the book <br> &nbsp;&nbsp;
  [www.clawpack.org/book.html Finite Volume Methods for Hyperbolic Problems] 
  by R. J. LeVeque.

Converted to [www.clawpack.org Clawpack 4.4] form by RJL in 2009.


<h4>Description:</h4>
<p>

1d acoustics in a constant medium.
          \[   q_t + A q_x = 0  \]
where
          \[   q(x,t) = \vector{ p(x,t)\\ u(x,t)}           \]
and the coefficient matrix is
          \[   A = \begin{matrix}
                        0         & K\\
                        1/\rho & 0
                        \end{matrix}.
          \]



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

Also provides options for running the code and plotting results.
See the documentation at [claw: doc/?]


<p>
<dt>[code: setprob.f95]
<dd>
A routine by this name is called by the library routine 
[clawcode: clawpack/1d/lib/main.f95]
and is generally used to set any values needed for the specific problem
being solved.  In this example, values specifying the material and the
intitial conditions are read from the file [code: setprob.data]. 

<p>
<dt>[code: rp1.f]
<dd>
This is the Riemann solver, which takes the $q$ values stored in the
arrays <tt>ql</tt> and <tt>qr</tt> and returns the waves in the array
<tt>wave</tt> and speeds in the array <tt>s</tt> that result in solving the
Riemann problem at each cell interface, and the fluctuations <tt>amdq</tt>
and <tt>apdq</tt>.  See [claw:doc/rp1.html] for more information about 1d
Riemann solvers.
<p>

<p>
<dt>[code: qinit.f]
<dd>
This subroutine sets the initial data $q(x,0)$ at time $t=0$.

</dl>

<h4>
Python files
</h4>
<dl>

<dt>[code: setrun.py]
<dd> This file contains a function that specifies the run-time parameters to
be used.  "python setrun.py" creates data files *.data needed by the
Fortran codes.
<p>
Some parameters that you might want to modify are described in the
[http://kingkong.amath.washington.edu/clawpack/users documentation].

<dt>[code: setplot.py]
<dd> This file contains a function that
specifies what plots will be done and
sets various plotting parameters.



</dl>


<h4>
Data files
</h4>

<p>
<font color='red'>Note:</font> The files *.data are generated
when [code: setrun.py] is executed.

<p>
<dl>
<dt>[code: claw.data]
<dd> This file contains a number of
parameter values that are used by CLAWPACK.
The values in this file are read by the library routine
[clawcode: clawpack/1d/lib/main.f].


<p>
<dt> [code: setprob.data]
<dd> This file contains the advection velocity $u$ and various other
parameters used in setting the initial conditions.
Values in this file are read  in by the
subroutine [code: setprob.f95].



</dl>

<p>
<h4>Library routines</h4>
<p>

In addition to the Fortran routines in this library, several library
routines from [claw:clawpack/1d/lib] are used.  See the [code: Makefile]
to determine which ones are used.

end_html
