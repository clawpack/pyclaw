c
c
c
c     ==============================================================
      subroutine claw2(maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &           q,aux,xlower,ylower,dx,dy,tstart,tend,dtv,
     &           cflv,nv,method,mthlim,mthbc,
     &           work,mwork,info,bc2,rpn2,rpt2,src2,b4step2)
c     ==============================================================
c
c
c
c  Solves a hyperbolic system of conservation laws in two space dimensions
c  of the general form
c  
c     capa * q_t + A q_x + B q_y = psi
c
c  The "capacity function" capa(x,y) and source term psi are optional 
c  (see below).
c
c  For a more complete description see the documentation at
c      http://www.amath.washington.edu/~claw
c
c  Sample driver programs and user-supplied subroutines are available.
c  See the the directories claw/clawpack/2d/example* for some examples, and
c  codes in claw/applications for more extensive examples.
c
c  --------------------------------------------------------
c
c  The user must supply the following subroutines:
c
c    bc2, rpn2, rpt2,        subroutines specifying the boundary conditions
c                            and Riemann solvers.
c
c    b4step2            The routine b4step2 is called each time step and
c                       can be supplied by the user in order to perform
c                       other operations that are necessary every time
c                       step.  For example, if the variables stored in
c                       the aux arrays are time-dependent then these
c                       values can be set.   
c
c  In addition, if the equation contains source terms psi, then the user
c  must provide:
c
c    src2               subroutine that solves capa * q_t = psi
c                       over a single time step.
c
c  These routines must be declared EXTERNAL in the main program.
c  For description of the calling sequences, see below.
c
c  Dummy routines b4step1.f and src1.f are available in
c       claw/clawpack/1d/lib
c
c  A subroutine implementing many standard boundary conditions is
c  available in claw/clawpack/2d/lib/bc2.f.
c
c
c
c  Description of parameters...
c  ----------------------------
c
c    maxmx is the maximum number of interior grid points in x, 
c          and is used in declaration of the array q
c
c    maxmy is the maximum number of interior grid points in y, 
c          and is used in declaration of the array q
c
c    meqn is the number of equations in the system of
c         conservation laws.
c
c    mwaves is the number of waves that result from the
c           solution of each Riemann problem.  Often mwaves = meqn but
c           for some problems these may be different, e.g. for the Euler
c           equations meqn = 4 but mwaves = 3 since there are only 3
c           distinct wave speeds.
c
c    mbc is the number of "ghost cells" that must be added on to each
c       side of the domain to handle boundary conditions.  The cells
c       actually in the physical domain are labelled from 1 to mx in x and
c       from 1 to my in y.  The arrays are dimensioned actually indexed
c       from 1-mbc to mx+mbc and from 1-mbc to my+mbc.
c       For the methods currently implemented, mbc = 2 should be used.
c       If the user implements another method that has a larger stencil and
c       hence requires more ghost cells, a larger value of mbc could be used.
c       q is extended from the physical domain to the ghost cells by the
c       user-supplied routine bc2.
c
c    mx is the number of grid cells in the x-direction, in the
c       physical domain.  In addition there are mbc grid cells
c       along each edge of the grid that are used for boundary
c       conditions.
c       Must have mx .le. maxmx
c 
c    my is the number of grid cells in the y-direction, in the
c       physical domain.  In addition there are mbc grid cells
c       along each edge of the grid that are used for boundary
c       conditions.
c       Must have my .le. maxmy
c 
c    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn) 
c        On input:  initial data at time tstart.
c        On output: final solution at time tend.
c        q(i,j,m) = value of mth component in the (i,j) cell.
c        Values within the physical domain are in q(i,j,m) 
c                for i = 1,2,...,mx   and j = 1,2,...,my.
c        mbc extra cells on each end are needed for boundary conditions
c        as specified in the routine bc2.
c
c    aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)
c        Array of auxiliary variables that are used in specifying the problem.
c        If method(7) = 0 then there are no auxiliary variables and aux
c                         can be a dummy variable.
c        If method(7) = maux > 0 then there are maux auxiliary variables
c                         and aux must be dimensioned as above.
c
c        Capacity functions are one particular form of auxiliary variable.
c        These arise in some applications, e.g. the
c        determinant of the Jacobian if a mapped grid is used, or a density
c        or porosity function in some advection problems.  
c        See Clawpack Note # 5 for examples.
c
c        If method(6) = 0 then there is no capacity function.
c        If method(6) = mcapa > 0  then there is a capacity function and 
c            capa(i,j), the "capacity" of the (i,j) cell, is assumed to be 
c            stored in aux(i,j,mcapa).
c            In this case we require method(7).ge.mcapa.
c
c    dx = grid spacing in x.  
c         (for a computation in ax <= x <= bx,  set dx = (bx-ax)/mx.)
c
c    dy = grid spacing in y.  
c         (for a computation in ay <= y <= by,  set dy = (by-ay)/my.)
c
c    tstart = initial time.
c
c    tend = Desired final time (on input).
c              If tend<tstart, then claw2 returns after a single successful
c                 time step has been taken (single-step mode).  
c              Otherwise, as many steps are taken as needed to reach tend, 
c                 up to a maximum of nv(1).
c         = Actual time reached (on output).
c
c    dtv(1:5) = array of values related to the time step:
c               (Note: method(1)=1 indicates variable size time steps)
c         dtv(1) = value of dt to be used in all steps if method(1) = 0
c                = value of dt to use in first step if method(1) = 1
c         dtv(2) = unused if method(1) = 0.
c                = maximum dt allowed if method(1) = 1.
c         dtv(3) = smallest dt used (on output)
c         dtv(4) = largest dt used (on output)
c         dtv(5) = dt used in last step (on output)
c
c    cflv(1:4) = array of values related to Courant number:
c         cflv(1) = maximum Courant number to be allowed.  
c                   With variable time steps the step is retracted and a 
c                   smaller step taken if the Courant
c                   number is larger than this value.  
c                   With fixed time steps the routine aborts.
c                   Usually cflv(1)=1.0 should work 
c                   (or cflv(1)=0.5 if method(3)=0).
c         cflv(2) = unused if method(1) = 0.
c                 = desired Courant number if method(1) = 1.
c                   Should be somewhat less than cflv(1), e.g. 0.9
c         cflv(3) = largest Courant number observed (on output).
c         cflv(4) = Courant number in last step (on output).
c
c    nv(1:2) = array of values related to the number of time steps:
c         nv(1) = unused if method(1) = 0
c               = maximum number of time steps allowed if method(1) = 1
c         nv(2) = number of time steps taken (on output).
c
c    method(1:7) = array of values specifying the numerical method to use
c                  and also indicating whether source terms, capacity
c                  function, auxiliary variables are present in the equation.
c
c         method(1) = 0 if fixed size time steps are to be taken.
c                       In this case, dt = dtv(1) in all steps.
c                   = 1 if variable time steps are to be used.
c                       In this case, dt = dtv(1) in the first step and
c                       thereafter the value cflv(2) is used to choose the
c                       next time step based on the maximum wave speed seen
c                       in the previous step.  Note that since this value
c                       comes from the previous step, the Courant number will
c                       not in general be exactly equal to the desired value
c                       If the actual Courant number in the next step is
c                       greater than cflv(1), then this step is redone with a 
c                       smaller dt.
c
c         method(2) = 1 if only first order increment waves are to be used.
c                   = 2 if second order correction terms are to be added, with
c                       a flux limiter as specified by mthlim.  
c
c
c         method(3) = 0 if no transverse propagation is to be applied.
c                       Increment and perhaps correction waves are propagated
c                       normal to the interface.
c                   = 1 if transverse propagation of increment waves 
c                       (but not correction waves, if any) is to be applied.
c                   = 2 if transverse propagation of correction waves is also
c                       to be included.  
c
c                   = -1 if dimensional splitting is to be used instead
c                        of the multi-dimensional wave-propagation.  The
c                        Godunov splitting is used which consists of
c                        sweeping first in x and then in y, with a step of
c                        length dt in each.  The routine bc2 is called
c                        before either sweep to set boundary data, and in
c                        the x-sweep goes over the rows of ghost cells too
c                        so that proper boundary conditions should be set
c                        for the y-sweeps by this process.  Dimensional
c                        splitting is somewhat faster than the unsplit
c                        method and works as well for many (though not all)
c                        problems.
c
c                   = -2 if dimensional splitting is to be used with the
c                        Strang splitting, consisting of 
c                           sweep in x over time dt/2
c                           sweep in y over time dt
c                           sweep in x over time dt/2
c                        This is not recommended because it is slower than
c                        the Godunov splitting and does not appear to be
c                        appreciably better.  Moreover, the boundary
c                        conditions will not be properly set for the final
c                        x-sweep.  (The code could be modified to achieve
c                        this by sweeping over more ghost cells.)
c
c         method(4) = 0 to suppress printing
c                   = 1 to print dt and Courant number every time step
c
c         method(5) = 0 if there is no source term psi.  In this case
c                       the subroutine src2 is never called so a dummy
c                       parameter can be given.
c                   = 1 if there is a source term.  In this case 
c                       the subroutine src2 must be provided and a 
c                       fractional step method is used.
c                       In each time step the following sequence is followed:
c                            call bc to extend data to ghost cells
c                            call step2 to advance hyperbolic eqn by dt
c                            call src2 to advance source terms by dt
c                   = 2 if there is a source term and Strang splitting is to
c                       be used instead of the Godunov splitting above.
c                       In each time step the following sequence is followed:
c                            call bc to extend data to ghost cells
c                            call src2 to advance source terms by dt/2
c                            call step2 to advance hyperbolic equation by dt
c                            call src2 to advance source terms by dt/2
c                       For most problems 1 is recommended rather than 2
c                       since it is less expensive and works essentially as 
c                       well on most problems.  
c                           
c
c         method(6) = 0 if there is no capacity function capa.  
c                   = mcapa > 0 if there is a capacity function.  In this case 
c                       aux(i,j,mcapa) is the capacity of cell (i,j) and you
c                       must also specify method(7) .ge. mcapa and set aux.
c
c         method(7) = 0 if there is no aux array used.
c                   = maux > 0  if there are maux auxiliary variables.
c
c         The recommended choice of methods for most problems is 
c            method(1) = 1,  method(2) = 2,  method(3) = 2.
c
c
c    mthlim(1:mwaves) = array of values specifying the flux limiter to be used
c                     in each wave family mw.  Often the same value will be used
c                     for each value of mw, but in some cases it may be
c                     desirable to use different limiters.  For example,
c                     for the Euler equations the superbee limiter might be
c                     used for the contact discontinuity (mw=2) while another
c                     limiter is used for the nonlinear waves.  Several limiters
c                     are built in and others can be added by modifying the
c                     subroutine philim.
c
c        mthlim(mw) = 0 for no limiter
c                   = 1 for minmod
c                   = 2 for superbee
c                   = 3 for van Leer
c                   = 4 for monotonized centered
c
c    mthbc(1:4) = array of values specifying what boundary conditions should
c                 be used at each edge of the domain, if the standard
c                 bc2.f routine is used.  Passed to bc2.
c
c    work(mwork) = double precision work array of length at least mwork
c
c    mwork = length of work array.  Must be at least
c               N * (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn   
c               + (max(mx,my) + 2*mbc) * (10*meqn + mwaves + meqn*mwaves 
c                                          + 3*maux + 2) 
c            where N = 1 if method(5)<2  (no source term or Godunov splitting)
c                  N = 2 if method(5)=2  (source term with Strang splitting)
c            If mwork is too small then the program returns with info = 4
c            and also prints the required value of mwork to unit 6.
c
c            
c    info = output value yielding error information:
c         = 0 if normal return.
c         = 1 if mx.gt.maxmx  or   my.gt.maxmy  or  mbc.lt.2
c         = 2 if method(1)=0 and dt doesn't divide (tend - tstart).
c         = 3 if method(1)=1 and cflv(2) > cflv(1).
c         = 4 if mwork is too small.
c         = 5 if method(6) > method(7)
c         = 6 if method(3) > method(2)
c         = 11 if the code attempted to take too many time steps, n > nv(1).
c              This could only happen if method(1) = 1 (variable time steps).
c         = 12 if the method(1)=0 and the Courant number is greater than 1
c              in some time step.
c
c           Note: if info.ne.0, then tend is reset to the value of t actually
c           reached and q contains the value of the solution at this time.
c
c    User-supplied subroutines
c    -------------------------
c
c    bc2 = subroutine that specifies the boundary conditions.  
c         This subroutine should extend the values of q from cells
c         (1:mx, 1:my) to the mbc ghost cells along each edge of the domain.
c
c
c    rpn2 = user-supplied subroutine that implements the Riemann solver
c           along a one-dimensional slice of data.
c
c          The form of this subroutine is
c  -------------------------------------------------
c     subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
c    &                  auxl,auxr,wave,s,amdq,apdq)
c
c     implicit double precision (a-h,o-z)
c     dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
c     dimension    s(1-mbc:maxm+mbc, mwaves)
c     dimension   ql(1-mbc:maxm+mbc, meqn)
c     dimension   qr(1-mbc:maxm+mbc, meqn)
c     dimension auxl(1-mbc:maxm+mbc, *)
c     dimension auxr(1-mbc:maxm+mbc, *)
c     dimension amdq(1-mbc:maxm+mbc, meqn)
c     dimension apdq(1-mbc:maxm+mbc, meqn)
c  -------------------------------------------------
c
c         On input, ql contains the state vector at the left edge of each cell
c                   qr contains the state vector at the right edge of each cell
c                 auxl contains auxiliary values at the left edge of each cell
c                 auxr contains auxiliary values at the right edge of each cell
c
c         This data is along a slice in the x-direction if ixy=1
c                                    or the y-direction if ixy=2.
c
c         Note that the i'th Riemann problem has left state qr(i-1,:)
c                                            and right state ql(i,:)
c         In the standard clawpack routines, this Riemann solver is 
c         called with ql=qr=q along this slice.  More flexibility is allowed
c         in case the user wishes to implement another solution method
c         that requires left and rate states at each interface.

c         If method(7)=maux > 0 then the auxiliary variables along this slice
c         are passed in using auxl and auxr.  Again, in the standard routines
c         auxl=auxr is just the values of aux along this slice.

c          On output, 
c             wave(i,m,mw) is the mth component of the jump across
c                              wave number mw in the ith Riemann problem.
c             s(i,mw) is the wave speed of wave number mw in the
c                              ith Riemann problem.
c             amdq(i,m) is the m'th component of the left-going flux difference.
c             apdq(i,m) is the m'th component of the right-going flux difference.
c           It is assumed that each wave consists of a jump discontinuity
c           propagating at a single speed, as results, for example, from a
c           Roe approximate Riemann solver.  An entropy fix can be included
c           into the specification of amdq and apdq.
c
c
c    rpt2 = user-supplied subroutine that implements the splitting of
c           a flux difference asdq into waves in the transverse direction.
c           The form of this subroutine is
c  -------------------------------------------------
c     subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,aux1,aux2,aux3,
c                     imp,asdq,bmasdq,bpasdq)
c
c     implicit double precision (a-h,o-z)
c     dimension     ql(1-mbc:maxm+mbc, meqn)
c     dimension     qr(1-mbc:maxm+mbc, meqn)
c     dimension   aux1(1-mbc:maxm+mbc, maux)
c     dimension   aux2(1-mbc:maxm+mbc, maux)
c     dimension   aux3(1-mbc:maxm+mbc, maux)
c     dimension   asdq(1-mbc:maxm+mbc, meqn)
c     dimension bmasdq(1-mbc:maxm+mbc, meqn)
c     dimension bpasdq(1-mbc:maxm+mbc, meqn)
c  -------------------------------------------------
c          On input, 
c              ql,qr is the data along some one-dimensional slice, as in rpn2
c                   This slice is in the x-direction
c                   if ixy=1, or in the y-direction if ixy=2.  
c              aux2 is the auxiliary array (if method(6)=maux>0) along
c                   this slice, say at j=J if ixy=1.
c              aux1 is the auxiliary array along the adjacent slice J-1
c              aux3 is the auxiliary array along the adjacent slice J+1
c          
c              asdq is an array of flux differences (A^* \Delta q).  
c                   asdq(i,:) is the flux difference propagating away from
c                   the interface between cells i-1 and i.
c              imp = 1 if asdq = A^- \Delta q,  the left-going flux difference
c                    2 if asdq = A^+ \Delta q, the right-going flux difference
c          On output, 
c              bmasdq is the down-going portion of the flux difference
c                   determined by solving a Riemann problem in the transverse
c                   direction using asdq as data.  
c              bpasdq is the up-going portion of the flux difference.
c           For example, for a linear system q_t + Aq_x + Bq_y = 0,  
c                   asdq = A^+ dq  or  A^- dq
c                   and this is then split into
c                       bmasdq = B^- asdq   and   bpasdq = B^+ asdq
c
c
c    src2 = user-supplied subroutine that takes one time step on the 
c           source terms alone, solving
c               capa * q_t = psi
c           over time dt.
c
c           If method(5)=0 then the equation does not contain a source
c           term and this routine is never called.  A dummy argument can
c           be used with many compilers, or provide a dummy subroutine that
c           does nothing (such a subroutine can be found in 
c           clawpack/2d/lib/src2.f)
c
c           The form of this subroutine is
c  -------------------------------------------------
c      subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
c    &                 dx,dy,q,maux,aux,told,dt2)
c      implicit double precision (a-h,o-z)
c      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
c      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
c  -------------------------------------------------
c      If method(7)=0  or the auxiliary variables are not needed in this solver,
c      then the latter dimension statement can be omitted, but aux should
c      still appear in the argument list.
c
c      On input, q(i,j,m) contains the data for solving the 
c                source term equation.
c      On output, q(i,j,m) should have been replaced by the solution to
c                 the source term equation after a step of length dt.
c
c
c
c      b4step2 = subroutine that is called from claw2 before each call to
c                step2.  Use to set time-dependent aux arrays or perform
c                other tasks which must be done every time step.
c
c          The form of this subroutine is
c
c  -------------------------------------------------
c      subroutine b4step2(maxmx,maxmy,mbc,mx,my,meqn,q,
c    &            xlower,ylower,dx,dy,time,dt,maux,aux)
c      implicit double precision (a-h,o-z)
c      dimension   q(1-mbc:maxmx+mbc, meqn)
c      dimension aux(1-mbc:maxmx+mbc, *)
c  -------------------------------------------------
c
c  
c
c =========================================================================
c
c  Copyright 1994 -- 1999 R. J. LeVeque
c
c  This software is made available for research and instructional use only. 
c  You may copy and use this software without charge for these non-commercial
c  purposes, provided that the copyright notice and associated text is
c  reproduced on all copies.  For all other uses (including distribution of
c  modified versions), please contact the author at the address given below. 
c  
c  *** This software is made available "as is" without any assurance that it
c  *** will work for your purposes.  The software may in fact have defects, so
c  *** use the software at your own risk.
c
c  --------------------------------------
c    CLAWPACK Version 4.1,  August, 2002
c    Webpage: http://www.amath.washington.edu/~claw
c  --------------------------------------
c
c    Author:  Randall J. LeVeque
c             Applied Mathematics
c             Box 352420
c             University of Washington, 
c             Seattle, WA 98195-2420
c             rjl@amath.washington.edu
c =========================================================================
c
c
c
c            
c    ======================================================================
c    Beginning of claw2 code
c    ======================================================================
c 
      implicit double precision (a-h,o-z)
      external bc2,rpn2,rpt2,src2,b4step2
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
      dimension work(mwork)
      dimension mthlim(mwaves),method(7),dtv(5),cflv(4),nv(2),mthbc(4)
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom
c
      maxm = max0(maxmx, maxmy)
      info = 0
      t = tstart
      maxn = nv(1)
      dt = dtv(1)   !# initial dt
      cflmax = 0.d0
      dtmin = dt
      dtmax = dt
      nv(2) = 0
      maux = method(7)
c
c     # check for errors in data:
c     ---------------------------
c
      if (mx.gt.maxmx .or. my.gt.maxmy .or. mbc.lt.2) then
         info = 1
         write(6,*) 'CLAW2 ERROR...  check mx,maxmx,my,maxmy,mbc'
         go to 900
         endif
c
      if (method(1) .eq. 0) then
c        # fixed size time steps.  Compute the number of steps:
         if (tend .lt. tstart) then
c             # single step mode 
              maxn = 1
           else
              maxn = (tend - tstart + 1d-10) / dt
              if (dabs(maxn*dt - (tend-tstart)) .gt.
     &                          1d-5*(tend-tstart)) then
c                # dt doesn't divide time interval integer number of times
                 info = 2
                 write(6,*) 
     &               'CLAW2 ERROR... dt does not divide (tend-tstart)'
                 go to 900
                 endif
           endif
         endif

c
      if (method(1).eq.1 .and. cflv(2).gt.cflv(1)) then
         info = 3
         write(6,*) 'CLAW2 ERROR...  cflv(2) > cflv(1)'
         go to 900
         endif
c
      if (method(6).gt.method(7)) then
         info = 5
         write(6,*) 'CLAW2 ERROR...  method(6) > method(7)'
         go to 900
         endif
c
      if (method(2) .lt. method(3)) then
         info = 6
         write(6,*) 'CLAW2 ERROR...  method(3) > method(2)'
         go to 900
         endif
c
      if (method(5).lt.2) then
          narray = 1   !# only need one qwork array
        else
          narray = 2   !# need two qwork arrays for Strang splitting
        endif
c
      mwork0 = (maxm+2*mbc)*(10*meqn + mwaves + meqn*mwaves 
     &                      + 3*maux + 2) 
     &          + narray * (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn   
c
      if (mwork .lt. mwork0) then
         info = 4
         write(6,*) 'CLAW2 ERROR... mwork should be increased to ',
     &               mwork0
         go to 900
         endif
c
c     # partition work array into pieces needed for local storage in 
c     # step2 routine. Find starting index of each piece:
c
      i0qadd = 1
      i0fadd = i0qadd + (maxm+2*mbc)*meqn
      i0gadd = i0fadd + (maxm+2*mbc)*meqn
      i0q1d = i0gadd + 2*(maxm+2*mbc)*meqn 
      i0dtdx1 = i0q1d + (maxm+2*mbc)*meqn  
      i0dtdy1 = i0dtdx1 + (maxm+2*mbc)
      i0qwrk1 = i0dtdy1 + (maxm+2*mbc)
c
      nqwork = (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn  !# size of q array
      if (method(5).lt.2) then
          i0qwrk2 = i0qwrk1  !# qwrk2 points to same storage as qwrk1
        else
          i0qwrk2 = i0qwrk1 + nqwork  !# second qwork array is needed for
                                      !# Strang spliting
        endif
c
      i0aux1 = i0qwrk2 + nqwork
      i0aux2 = i0aux1 + (maxm+2*mbc)*maux
      i0aux3 = i0aux2 + (maxm+2*mbc)*maux
c
      i0next = i0aux3 + (maxm+2*mbc)*maux  !# next free space
      mused = i0next - 1                  !# space already used
      mwork1 = mwork - mused              !# remaining space (passed to step2)
c
c
c
c     -----------
c     # main loop
c     -----------
c
      if (maxn.eq.0) go to 900
      do 100 n=1,maxn
         told = t   !# time at beginning of time step.

c        # adjust dt to hit tend exactly if we're near end of computation
c        #  (unless tend < tstart, which is a flag to take only a single step)
         if (told+dt.gt.tend .and. tstart.lt.tend) dt = tend - told

c
   40    continue
c
c        # store dt and t in the common block comxyt in case they are needed
c        # in the Riemann solvers (for variable coefficients)
         tcom = told
         dtcom = dt
         dxcom = dx
         dycom = dy
c
c
c
c        ================================================================
c
c        -------------------------
c        # main steps in algorithm
c        -------------------------
c
c        # extend data from grid to bordering boundary cells:
         call bc2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &               dx,dy,q,maux,aux,told,dt,mthbc)
c
c
c
c        # call user-supplied routine which might set aux arrays
c        # for this time step, for example.

         call b4step2(maxmx,maxmy,mbc,mx,my,meqn,q,
     &                xlower,ylower,dx,dy,told,dt,maux,aux)
c
c
c
         if (method(5).eq.2) then
c            # with Strang splitting for source term:
c            # First need to store solution before taking
c            # step on source terms in case we need to redo everything with
c            # a smaller time step if the Courant number is too large in 
c            # subroutine step2.
             call copyq2(maxmx,maxmy,meqn,mbc,mx,my,q,work(i0qwrk2))
c
c            # source terms over a half time step:
             dt2 = dt / 2.d0
             call src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,maux,aux,told,dt2)
             endif
c
c        # copy q into qwork1.  q is updated in step2 and qwork1 is
c        # preserved to provide data for Riemann problems.
c        # qwork1 can also be used to restart if the Courant number is 
c        # too large, unless Strang splitting is used in which case we
c        # must used the values already stored above before
c        # taking the source term step.
         call copyq2(maxmx,maxmy,meqn,mbc,mx,my,q,work(i0qwrk1))
c
c        # take one step on the conservation law:
c
         if( method(3) .ge. 0 )then
c            # unsplit version
c             
             call step2(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &                  work(i0qwrk1),q,aux,
     &                  dx,dy,dt,method,mthlim,cfl,
     &                  work(i0qadd),work(i0fadd),work(i0gadd),
     &                  work(i0q1d),work(i0dtdx1),work(i0dtdy1),
     &                  work(i0aux1),work(i0aux2),work(i0aux3),
     &                  work(i0next),mwork1,rpn2,rpt2)
c
         else
c           # dimensional splitting (fractional steps)
c
            call dimsp2(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &                  work(i0qwrk1),q,aux,
     &                  dx,dy,dt,method,mthlim,cfl,cflv,
     &                  work(i0qadd),work(i0fadd),work(i0gadd),
     &                  work(i0q1d),work(i0dtdx1),work(i0dtdy1),
     &                  work(i0aux1),work(i0aux2),work(i0aux3),
     &                  work(i0next),mwork1,rpn2,rpt2)
c
         endif
c
         t = told + dt
c
         if (method(4).eq.1) then
c            # verbose mode
             write(6,601) n,cfl,dt,t
  601        format('CLAW2... Step',i6,
     &                   '   Courant number =',f6.3,'  dt =',d12.4,
     &                   '  t =',d12.4)
             endif
c
c
c        # check to see if the Courant number was too large:
         if (cfl .le. cflv(1)) then
c             # accept this step
              cflmax = dmax1(cfl,cflmax)
            else
c             # Reject this step.  Reset q to qwork from previous time:
c             # Note that i0qwrk2 points to work space where previous
c             # solution is stored in all cases method(5) = 0,1, or 2. 
              t = told
              call copyq2(maxmx,maxmy,meqn,mbc,mx,my,work(i0qwrk2),q)
c
              if (method(4).eq.1) then
c                # verbose mode
                 write(6,602)
                 endif
  602         format('CLAW2 rejecting step... Courant number too large')
c
              if (method(1).eq.1) then
c                 # if variable dt, go back and take a smaller step.
                  dt = dmin1(dtv(2), dt * cflv(2)/cfl)
                  go to 40
                else
c                 # if fixed dt, give up and return
                  cflmax = dmax1(cfl,cflmax)
                  go to 900
                endif
            endif
c
c
c        # claw2 step is accepted
c        # now apply source terms:
c
         if (method(5).eq.2) then
c            # source terms over a second half time step for Strang splitting:
c            # Note it is not so clear what time t should be used here if
c            # the source terms are time-dependent!
             call src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,maux,aux,t,dt2)
             endif
c
         if (method(5).eq.1) then
c            # source terms over a full time step:
             call src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,maux,aux,t,dt)
             endif
c
c        ================================================================
c
c
c
         if (method(1) .eq. 1) then
c           # choose new time step if variable time step
            if (cfl.eq.0.d0) then 
                dt = dtv(2)
              else
                dt = dmin1(dtv(2), dt * cflv(2)/cfl)
              endif
            dtmin = dmin1(dt,dtmin)
            dtmax = dmax1(dt,dtmax)
            endif
c
c
c        # see if we are done:
c
         nv(2) = nv(2) + 1
         if (t .ge. tend) go to 900
  100    continue
c
  900  continue
c 
c      # return information
c
       if (method(1).eq.1 .and. t.lt.tend .and. nv(2) .eq. maxn) then
c         # too many timesteps
          write(6,*) 'CLAW2 ERROR...  too many timesteps'
          info = 11
          endif
       if (method(1).eq.0 .and. cflmax .gt. cflv(1)) then
c         # Courant number too large with fixed dt
          write(6,*) 'CLAW2 ERROR...  Courant number too large'
          info = 12
          endif
       tend = t
       cflv(3) = cflmax
       cflv(4) = cfl
       dtv(3) = dtmin
       dtv(4) = dtmax
       dtv(5) = dt
c
       return 
       end
