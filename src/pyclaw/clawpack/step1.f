c
c
c ===================================================================
      subroutine step1(num_eqn,num_waves,num_ghost,num_aux,mx,q,aux,dx,
     &           dt,method,mthlim,cfl,f,wave,s,amdq,apdq,dtdx,use_fwave)
c ===================================================================
c
c     # Take one time step, updating q.
c
c     method(1) = 1   ==>  Godunov method
c     method(1) = 2   ==>  Slope limiter method
c     mthlim(p)  controls what limiter is used in the pth family
c
c
c     amdq, apdq, wave, s, and f are used locally:
c
c     amdq(1-num_ghost:mx+num_ghost, num_eqn) = left-going flux-differences
c     apdq(1-num_ghost:mx+num_ghost, num_eqn) = right-going flux-differences
c        e.g. amdq(m,i) = m'th component of A^- \Delta q from i'th Riemann
c                         problem (between cells i-1 and i).
c
c     wave(1-num_ghost:mx+num_ghost, num_eqn, num_waves) = waves from solution of
c                                           Riemann problems,
c            wave(m,mw,i) = mth component of jump in q across
c                           wave in family mw in Riemann problem between
c                           states i-1 and i.
c
c     s(1-num_ghost:mx+num_ghost, num_waves) = wave speeds,
c            s(m,iw) = speed of wave in family mw in Riemann problem between
c                      states i-1 and i.
c
c     f(1-num_ghost:mx+num_ghost, num_eqn) = correction fluxes for second order method
c            f(m,i) = mth component of flux at left edge of ith cell 
c     --------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)    
      double precision q(num_eqn,1-num_ghost:mx+num_ghost)
      double precision aux(num_aux,1-num_ghost:mx+num_ghost)
      double precision f(num_eqn,1-num_ghost:mx+num_ghost)
      double precision s(num_waves,1-num_ghost:mx+num_ghost)
      double precision wave(num_eqn, num_waves,1-num_ghost:mx+num_ghost)
      double precision amdq(num_eqn,1-num_ghost:mx+num_ghost)
      double precision apdq(num_eqn,1-num_ghost:mx+num_ghost)
      double precision dtdx(1-num_ghost:mx+num_ghost)
      integer          method(7),mthlim(num_waves)
      logical          use_fwave
      logical limit

cf2py intent(in,out) q  
cf2py intent(out) cfl  
cf2py intent(in) num_eqn  
cf2py intent(in) num_ghost  
cf2py intent(in) mx  
cf2py optional f, amdq, apdq, dtdx, s, wave

c
c     # check if any limiters are used:
      limit = .false.
      do 5 mw=1,num_waves
         if (mthlim(mw) .gt. 0) limit = .true.
   5     continue
c
      index_capa = method(6)
      do 10 i=1-num_ghost,mx+num_ghost
         if (index_capa.gt.0) then
             if (aux(index_capa,i) .le. 0.d0) then
                write(6,*) 'Error -- capa must be positive'
                stop
                endif
             dtdx(i) = dt / (dx*aux(index_capa,i))
            else
             dtdx(i) = dt/dx
            endif
   10    continue
c
c
c
c     # solve Riemann problem at each interface 
c     -----------------------------------------
c
      call rp1(mx,num_eqn,num_waves,num_ghost,mx,q,q,aux,aux,wave,
     &              s,amdq,apdq)
c
c     # Modify q for Godunov update:
c     # Note this may not correspond to a conservative flux-differencing
c     # for equations not in conservation form.  It is conservative if
c     # amdq + apdq = f(q(i)) - f(q(i-1)).
c



      
c     print *,"q before, from fortran",q(24,1)

      forall(i=1:mx+1, m=1:num_eqn)
         q(m,i) = q(m,i) - dtdx(i)*apdq(m,i)
         q(m,i-1) = q(m,i-1) - dtdx(i-1)*amdq(m,i)
      end forall

c
c     # compute maximum wave speed:
      cfl = 0.d0
      do 50 mw=1,num_waves
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
      forall(i=1-num_ghost:mx+num_ghost, m=1:num_eqn)
         f(m,i) = 0.d0
      end forall
c
c      # apply limiter to waves:
      if (limit) call limiter(mx,num_eqn,num_waves,num_ghost,mx,
     &                          wave,s,mthlim)
c
      if (use_fwave.eqv..false.) then
          do 120 i=1,mx+1
             do 120 m=1,num_eqn
                do 110 mw=1,num_waves
             dtdxave = 0.5d0 * (dtdx(i-1) + dtdx(i))
             f(m,i) = f(m,i) + 0.5d0 * dabs(s(mw,i))
     &                  * (1.d0 - dabs(s(mw,i))*dtdxave) * wave(m,mw,i)
  110           continue
  120       continue
      else
          do 121 i=1,mx+1
             do 121 m=1,num_eqn
                do 111 mw=1,num_waves
             dtdxave = 0.5d0 * (dtdx(i-1) + dtdx(i))
             f(m,i) = f(m,i) + 0.5d0 * dsign(1.d0,s(mw,i))
     &                  * (1.d0 - dabs(s(mw,i))*dtdxave) * wave(m,mw,i)
  111           continue
  121       continue
      endif

c
c
c     # update q by differencing correction fluxes 
c     ============================================
c
c     # (Note:  Godunov update has already been performed above)
c
      forall(i=1:mx+1, m=1:num_eqn)
            q(m,i) = q(m,i) - dtdx(i) * (f(m,i+1) - f(m,i))
      end forall
c
  900 continue
      return
      end
