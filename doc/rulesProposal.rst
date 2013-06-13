=============================================
Conventions for coding and documenting PyClaw
=============================================
In order to improve and maintain PyClaw and allow more developers or simply 
users to acces and use it all Python/Fortran code should follow some 
conventions. As stressed in many books, reports, notes, etc. (see for instance
`Wikipedia <http://en.wikipedia.org/wiki/Coding_conventions>`_) code conventions
are important to programmers for a number of reasons:

    * A big part of the lifetime cost of a piece of software goes to maintenance
    * Rarely is software maintained for its whole life by the original authors
    * Code conventions improve the readability of the software, allowing new developers and users to understand new code more quickly and thoroughly (efficiency increases!)
    * Standards makes the source code look consistent, well-organized, and professional
    * etc.

This page provides some suggested coding conventions for PyClaw developers.
Most of the ideas/proposed rules listed here are adapted
from some conventions already used by most of the Python, Fortran90/95, C, C++
communities. 

Most of the proposed rules listed here have been extracted from the
documentation available at following links:
    
    * `pep-0008 <http://www.python.org/dev/peps/pep-0008/>`_
    * `pep-0257 <http://www.python.org/dev/peps/pep-0257/>`_
    * `numpy <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
    * `sage <http://www.sagemath.org/doc/developer/conventions.html>`_


Code layout
===========
Indentation:
    * Use 4 spaces per indentation level
    * Don't use tabs
    * Confine your line width to 80 characters, so that your code can be displayed on screen for side-by-side comparison and printed easily on A4 paper.
    * For flowing long blocks of text (docstrings or comments), limiting the length to 72 characters is recommended.
    * Separate top-level function and class definitions with two blank lines.
    * Method definitions inside a class are separated by a single blank line.


**NOTE:** For Python 3.0 and beyond, the following policy is prescribed for the
standard library (see PEP 3131): All identifiers in the Python standard library
MUST use ASCII-only identifiers, and SHOULD use English words wherever feasible
(in many cases, abbreviations and technical terms are used which aren't
English). In addition, string literals and comments must also be in ASCII. The
only exceptions are (a) test cases testing the non-ASCII features, and (b)
names of authors. Authors whose names are not based on the latin alphabet MUST
provide a latin transliteration of their names. **Open source projects with a
global audience are encouraged to adopt a similar policy.**

Docstring conventions
=====================
A docstring is a string literal that occurs as the first statement in a module, function, class, or method definition. Such a docstring becomes the ``__doc__ special`` attribute of that object.

All modules should normally have docstrings, and all functions and classes exported by a module should also have docstrings. Public methods (including the __init__ constructor) should also have docstrings. A package may be documented in the module docstring of the ``__init__.py`` file in the package directory.

**One-line Docstrings:** One-liners are for really obvious cases. They should really fit on one line. For example:

.. code-block:: python

    def limiter_type():
        """Return the type of limiter."""
        global lim_type
        if lim_type: return lim_type
        ...

.. code-block:: fortran

    subroutine limiter_type(method_data)
    ! Return the type of limiter.
    ...
    ...
    return lim_type
    end subroutine limiter_type


**NOTES:**
    * Triple quotes are used even though the string fits on one line. This makes it easy to later expand it
    * The closing quotes are on the same line as the opening quotes. This looks better for one-liners
    * There's no blank line either before or after the docstring
    * The docstring is a phrase ending in a period. It prescribes the function or method's effect as a command ("Do this", "Return that"), not as a description
    * The one-line docstring should NOT be a "signature" reiterating the function/method parameters (which can be obtained by introspection)


**Multi-line docstrings:**

Multi-line docstrings consist of a summary line just like a one-line docstring, followed by a blank line, followed by a more elaborate description. The summary line may be used by automatic indexing tools; it is important that it fits on one line and is separated from the rest of the docstring by a blank line. The summary line may be on the same line as the opening quotes or on the next line. The entire docstring is indented the same as the quotes at its first line.

The docstring of a script (a stand-alone program) should be usable as its "usage" message, printed when the script is invoked with incorrect or missing arguments (or perhaps with a "-h" option, for "help"). Such a docstring should document the script's function and command line syntax, environment variables, and files. Usage messages can be fairly elaborate (several screens full) and should be sufficient for a new user to use the command properly, as well as a complete quick reference to all options and arguments for the sophisticated user.

The docstring for a module should generally list the classes, exceptions and functions (and any other objects) that are exported by the module, with a one-line summary of each. (These summaries generally give less detail than the summary line in the object's docstring.) The docstring for a package (i.e., the docstring of the package's ``__init__.py`` module) should also list the modules and subpackages exported by the package.

The docstring for a function or method should summarize its behavior and document its arguments, return value(s), side effects, exceptions raised, and restrictions on when it can be called (all if applicable). Optional arguments should be indicated. It should be documented whether keyword arguments are part of the interface.

The docstring for a class should summarize its behavior and list the public methods and instance variables. If the class is intended to be subclassed, and has an additional interface for subclasses, this interface should be listed separately (in the docstring). The class constructor should be documented in the docstring for its ``__init__`` method. Individual methods should be documented by their own docstring.

The extended summary should be used to clarify functionality, not to discuss implementation detail or background theory, which should rather be explored in the notes section below. You may refer to the parameters and the function name, but parameter descriptions still belong in the parameters section.

.. code-block:: python

    def complex(real=0.0, imag=0.0):
        """Form a complex number.

        Keyword arguments:
        real -- the real part (default 0.0)
        imag -- the imaginary part (default 0.0)

        """
        if imag == 0.0 and real == 0.0: return complex_zero
        ...

.. code-block:: fortran
    
    subroutine tfluct(ixy,maxmx,num_eqn,num_waves,num_ghost,mx,ql,qr,auxl,auxr,s,adq)

    ! Solve Riemann problems for the 2D shallow water equations
    ! using f-wave algorithm and Roe's approximate Riemann solver.  
    ! 
    ! Input arguments:
    ! ql -- left state vector at the left edge of each cell
    ! qr -- right state vector at the right edge of each cell
    !
    ! Output arguments:
    ! wave -- Riemann problem waves, 
    ! s    -- Waves speed, 
    ! amdq -- left-going flux difference  A^- \Delta q
    ! apdq -- right-going flux difference  A^+ \Delta q
    !
    !
    ! Note that the i'th Riemann problem has left state qr(i-1,:)
    !                                     and right state ql(i,:)
    ! From the basic clawpack routine step1, rp is called with ql = qr = q.


If it is not necessary to specify a keyword argument, use optional:

.. code-block:: python

    x : int, optional

An optional section detailing which errors get raised and under what conditions:

.. code-block:: python

    Raises
    ------
        LinAlgException
            If the matrix is not numerically invertible.

References cited in the notes section may be listed here, e.g. if you cited the article below using the text [Ref]_, include it as in the list as follows:

.. [Ref] Amal Alghamdi, Aron Ahmadia, David I. Ketcheson, Matthew G. Knepley, Kyle T. Mandli, and Lisandro Dalcin, PetClaw: A Scalable Parallel Nonlinear Wave Propagation Solver for Python.

Naming conventions
==================
Reasons for using a naming convention (as opposed to allowing programmers to choose any character sequence) include the following:
    * to reduce the effort needed to read and understand source code
    * to enhance source code appearance (for example, by disallowing overly long names or abbreviations)

Some of the potential benefits that can be obtained by adopting a naming convention include the following:
    * to provide additional information (i.e., metadata) about the use to which an identifier is put
    * to help formalize expectations and promote consistency within a development team
    * to enable the use of automated refactoring or search and replace tools with minimal potential for error
    * to help avoid "naming collisions"
    * to provide better understanding in case of code reuse after a long interval of time

Thus, I would use "self-explaining names for variables, procedures etc."


**Multiple-word identifiers:**
    * Delimiter-separated words: hyphen ('-') and the underscore ('_')
    * Letter-case separated words: indicate word boundaries using medial capitalization (here we can have the first word in lower case), e.g. limiterType


Order of the test cases instruction
===================================
It would be useful to follow also some rules when preparing the Python script
of a new test case. Listing  phases and instructions in a logical order could
improve the readability of the set-up. One idea could be:

    * Import libraries needed by all the functions
    * Define the functions use by the main program, e.g. qinit, setaux, etc.
      Here the conventions introduce previously for the docstrings should be used
    * Main function
        * Import libraries 
        * Initialize grid, solution and aux array
        * Setup the solver and solver parameters
        * Setup controller and controller parameters
        * Solve problem
        * Plot results

Add regression test to check new piece of code
==============================================
    * Add one or more regression test to check the functionality of the new code
    * Check with nose if all the tests pass before to commit

Add documentation when new code is added
========================================
    * What the new code does
    * How to use it
    * Document inputs outputs and  default parameters

Write Comments as You Code
==========================
You won't ever go back later and document your code. You just won't. So when
you do something document it right then and there. When you create a class-
document it. When you create a method- document it. And so on. That way when
you finish coding you will also be finished documenting.
Won't this break the flow? No, I think it improves flow because it keeps you
mindful of what you are doing, why you are doing, and how it fits in the big
picture.

Some fortran tips
======================
    * Always use **IMPLICIT NONE**
    * Always allocate and deallocate memory
