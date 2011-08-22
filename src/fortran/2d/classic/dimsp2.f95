c     ==========================================================
      subroutine dimsp2(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,mx,my,
     &                  qold,qnew,aux,dx,dy,dt,method,mthlim,cfl,
     &                  cflv,qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                  aux1,aux2,aux3,work,mwork)
c     ==========================================================
c
c     # Take one time step, updating q, using dimensional
c     # splitting. Two choices are available:
c     #
c     # method(3) = -1   gives Godunov splitting:
c     #    time step dt in x-direction
c     #    time step dt in y-direction
c
c     # method(3) = -2   gives Strang splitting
c     #    time step dt/2 in x-direction
c     #    time step dt   in y-direction
c     #    time step dt/2 in x-direction
c
c     # Godunov splitting is recommended over Strang splitting normally
c     # since it typically works as well, is faster, and boundary
c     # conditions are handled properly.
c
      implicit double precision (a-h,o-z)
      external rpn2,rpt2
      double precision qold(meqn, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      double precision qnew(meqn, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      double precision  q1d(meqn, 1-mbc:maxm+mbc)
      double precision cflv(4)
      double precision qadd(meqn, 1-mbc:maxm+mbc)
      double precision fadd(meqn, 1-mbc:maxm+mbc)
      double precision gadd(meqn, 2, 1-mbc:maxm+mbc)
      double precision aux(maux,1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      double precision aux1(maux,1-mbc:maxm+mbc)
      double precision aux2(maux,1-mbc:maxm+mbc)
      double precision aux3(maux,1-mbc:maxm+mbc)

      double precision dtdx1d(1-mbc:maxm+mbc)
      double precision dtdy1d(1-mbc:maxm+mbc)
      integer method(7),mthlim(mwaves)
      double precision work(mwork)

cf2py intent(in,out) cfl
cf2py intent(in,out) qnew  
cf2py optional q1d, qadd, fadd, gadd, dtdx1d, dtdy1d

c
c     # If method(3) = -1, take a full time step in x.
c     # If method(3) = -2, take a half time step in x.
c
      dt2 = dt/2.d0
c
      if( method(3) .eq. -2 )then
          call step2ds(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,mx,my,
     &                 qold,qnew,aux,dx,dy,dt2,method,mthlim,cflx,
     &                 qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,rpn2,rpt2,1)
      else
          call step2ds(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,mx,my,
     &                 qold,qnew,aux,dx,dy,dt,method,mthlim,cflx,
     &                 qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,rpn2,rpt2,1)
      endif
c
      if (cflx .gt. cflv(1)) then
c        # Abort if the Courant number was too large in x-sweep
         cfl = cflx
         return
         endif
c
c     # Take full step in y-direction
c
      call step2ds(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,mx,my,
     &             qnew,qnew,aux,dx,dy,dt,method,mthlim,cfly,
     &             qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &             aux1,aux2,aux3,work,mwork,rpn2,rpt2,2)
c
      cfl = dmax1(cflx,cfly)
c
c     # Finally, take a half time step in the x-direction
c     # if Strang splitting is used.  NOTE: boundary conditions may
c     # not be set properly for this sweep.
c
      if( method(3) .eq. -2 )then
         if (cfly .gt. cflv(1)) then
c           # Abort if the Courant number was too large in y-sweep
            cfl = cfly
            return
            endif
          call step2ds(maxm,maxmx,maxmy,meqn,mwaves,maux,mbc,mx,my,
     &                 qnew,qnew,aux,dx,dy,dt2,method,mthlim,cflx,
     &                 qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,rpn2,rpt2,1)
          cfl = dmax1(cfl,cflx)
      endif
c
      return
      end
