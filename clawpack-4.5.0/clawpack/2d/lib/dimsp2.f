c
c
c     ==========================================================
      subroutine dimsp2(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &                  qold,qnew,aux,dx,dy,dt,method,mthlim,cfl,
     &                  cflv,qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                  aux1,aux2,aux3,work,mwork,rpn2,rpt2)
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
      dimension qold(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension qnew(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      dimension  q1d(1-mbc:maxm+mbc, meqn)
      dimension cflv(4)
      dimension qadd(1-mbc:maxm+mbc, meqn)
      dimension fadd(1-mbc:maxm+mbc, meqn)
      dimension gadd(1-mbc:maxm+mbc, meqn, 2)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
      dimension aux1(1-mbc:maxm+mbc, *)
      dimension aux2(1-mbc:maxm+mbc, *)
      dimension aux3(1-mbc:maxm+mbc, *)

      dimension dtdx1d(1-mbc:maxmx+mbc)
      dimension dtdy1d(1-mbc:maxmx+mbc)
      dimension method(7),mthlim(mwaves)
      dimension work(mwork)

c
c     # If method(3) = -1, take a full time step in x.
c     # If method(3) = -2, take a half time step in x.
c
      dt2 = dt/2.d0
c
      if( method(3) .eq. -2 )then
          call step2ds(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &                 qold,qnew,aux,dx,dy,dt2,method,mthlim,cflx,
     &                 qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,rpn2,rpt2,1)
      else
          call step2ds(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
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
      call step2ds(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
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
          call step2ds(maxm,maxmx,maxmy,meqn,mwaves,mbc,mx,my,
     &                 qnew,qnew,aux,dx,dy,dt2,method,mthlim,cflx,
     &                 qadd,fadd,gadd,q1d,dtdx1d,dtdy1d,
     &                 aux1,aux2,aux3,work,mwork,rpn2,rpt2,1)
          cfl = dmax1(cfl,cflx)
      endif
c
      return
      end
