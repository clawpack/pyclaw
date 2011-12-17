c
c
c ===================================================================
      subroutine step1(meqn,mwaves,mbc,maux,mx,q,aux,dx,dt,
     &              method,mthlim,cfl,f,fwave,s,amdq,apdq,dtdx)
c ===================================================================
c     # This version uses interleaved arrays.
c
c     # Take one time step, updating q.
c
c     ----------------------------------------------------------------------
c     # step1fw is a modified version of step1 to use fwave instead of wave.
c     # A modified Riemann solver rp1 must be used in conjunction with this
c     # routine, which returns fwave's instead of wave's.
c     # See http://amath.washington.edu/~claw/fwave.html
c
c     # Limiters are applied to the fwave's, and the only significant
c     # modification of this code is in the "do 110" loop, for the 
c     # second order corrections.
c     ----------------------------------------------------------------------
c
c     method(1) = 1   ==>  Godunov method
c     method(1) = 2   ==>  Slope limiter method
c     mthlim(p)  controls what limiter is used in the pth family
c
c
c     amdq, apdq, fwave, s, and f are used locally:
c
c     amdq(meqn,1-mbc:mx+mbc) = left-going flux-differences
c     apdq(meqn,1-mbc:mx+mbc) = right-going flux-differences
c        e.g. amdq(m,i) = m'th component of A^- \Delta q from i'th Riemann
c                         problem (between cells i-1 and i).
c
c     fwave( meqn, mwaves,1-mbc:mx+mbc) = waves from solution of
c                                           Riemann problems,
c            fwave(m,mw,i) = mth component of jump in f across
c                           wave in family mw in Riemann problem between
c                           states i-1 and i.
c
c     s(1-mbc:mx+mbc, mwaves) = wave speeds,
c            s(mw,i) = speed of wave in family mw in Riemann problem between
c                      states i-1 and i.
c
c     f(meqn,1-mbc:mx+mbc) = correction fluxes for second order method
c            f(m,i) = mth component of flux at left edge of ith cell 
c     --------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      dimension    q(meqn,1-mbc:mx+mbc)
      dimension  aux(maux,1-mbc:mx+mbc)
      dimension    f(meqn,1-mbc:mx+mbc)
      dimension    s(mwaves,1-mbc:mx+mbc)
      dimension fwave(meqn, mwaves,1-mbc:mx+mbc)
      dimension amdq(meqn,1-mbc:mx+mbc)
      dimension apdq(meqn,1-mbc:mx+mbc)
      dimension dtdx(1-mbc:mx+mbc)
      dimension method(7),mthlim(mwaves)
      logical limit


cf2py intent(in,out) q  
cf2py intent(out) cfl  
cf2py intent(in) meqn  
cf2py intent(in) mbc  
cf2py intent(in) mx  
cf2py optional f, amdq, apdq, dtdx, s, fwave

c
c     # check if any limiters are used:
      limit = .false.
      do 5 mw=1,mwaves
	 if (mthlim(mw) .gt. 0) limit = .true.
   5     continue
c
      mcapa = method(6)
      do 10 i=1-mbc,mx+mbc
	 if (mcapa.gt.0) then
	     if (aux(mcapa,i) .le. 0.d0) then
		write(6,*) 'Error -- capa must be positive'
		stop
		endif
             dtdx(i) = dt / (dx*aux(mcapa,i))
	    else
             dtdx(i) = dt/dx
	    endif
   10	 continue
c
c
c
c     # solve Riemann problem at each interface 
c     -----------------------------------------
c
      call rp1(mx,meqn,mwaves,mbc,mx,q,q,aux,aux,fwave,s,amdq,apdq)
c
c     # Modify q for Godunov update:
c     # Note this may not correspond to a conservative flux-differencing
c     # for equations not in conservation form.  It is conservative if
c     # amdq + apdq = f(q(i)) - f(q(i-1)).
c
      do 40 i=1,mx+1
         do 39 m=1,meqn
            q(m,i) = q(m,i) - dtdx(i)*apdq(m,i)
            q(m,i-1) = q(m,i-1) - dtdx(i-1)*amdq(m,i)
   39       continue
   40    continue

c
c     # compute maximum wave speed:
      cfl = 0.d0
      do 50 mw=1,mwaves
	 do 45 i=1,mx+1
c          # if s>0 use dtdx(i) to compute CFL,
c          # if s<0 use dtdx(i-1) to compute CFL:
	   cfl = dmax1(cfl, dtdx(i)*s(mw,i), -dtdx(i-1)*s(mw,i))
   45      continue
   50    continue
c
      if (method(2) .eq. 1) go to 900
c
c     # compute correction fluxes for second order q_{xx} terms:
c     ----------------------------------------------------------
c
      do 100 m = 1, meqn
            do 100 i = 1-mbc, mx+mbc
               f(m,i) = 0.d0
  100          continue
c
c     # apply limiter to waves:
      if (limit) then
         call limiter(mx,meqn,mwaves,mbc,mx,fwave,s,mthlim)
         endif

c
      do 120 i=1,mx+1
	 do 120 m=1,meqn
	    do 110 mw=1,mwaves
	       dtdxave = 0.5d0 * (dtdx(i-1) + dtdx(i))
	       f(m,i) = f(m,i) + 0.5d0 * dsign(1.d0,s(mw,i))
     &           * (1.d0 - dabs(s(mw,i))*dtdxave) * fwave(m,mw,i)
c

  110          continue
  120       continue
c
c
c     # update q by differencing correction fluxes 
c     ============================================
c
c     # (Note:  Godunov update has already been performed above)
c
      do 150 m=1,meqn
	 do 150 i=1,mx
	    q(m,i) = q(m,i) - dtdx(i) * (f(m,i+1) - f(m,i))
  150       continue
c
  900 continue
      return
      end
