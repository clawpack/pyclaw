c
c
c =========================================================
      subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,maux,aux,t,dt)
c =========================================================
      implicit double precision(a-h,o-z)
      dimension    q(meqn, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
      dimension  aux(maux, 1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc)
cf2py intent(in,out) q 
c
c     # source terms for cylindrical symmetry in 2d Euler equations
c     # about y = 0, so radius = y
c     # y value at cell center is stored in aux(i,j,1)
c
      dimension qstar(4)
      common /param/  gamma,gamma1
c
c     # 2-stage Runge-Kutta method
c
      dt2    = dt/2.d0
      press  = 0.d0
      ndim   = 2

c      forall (i=1:mx, j=1:my)
      do 10 i=1,mx
        do 10 j=1,my
         rad      = aux(1,i,j)
         rho      = q(1,i,j)
         u        = q(2,i,j)/q(1,i,j)
         v        = q(3,i,j)/q(1,i,j)
         press    = gamma1*(q(4,i,j) - 0.5d0*rho*(u**2 + v**2))

         if (rad.eq.0.d0) write(6,*) 'rad = 0 in src2'
         if (rho.eq.0.d0) write(6,*) 'rho = 0 in q'

         qstar(1) = q(1,i,j) - dt2*(ndim-1)/rad * q(3,i,j)
         qstar(2) = q(2,i,j) - dt2*(ndim-1)/rad * 
     &                          (rho*u*v)
         qstar(3) = q(3,i,j) - dt2*(ndim-1)/rad * 
     &                          (rho*v*v)
         qstar(4) = q(4,i,j) - dt2*(ndim-1)/rad * 
     &                          v*(q(4,i,j) + press)
c
c        # second stage
c
         rho      = qstar(1)
         u        = qstar(2)/qstar(1)
         v        = qstar(3)/qstar(1)
         press    = gamma1*(qstar(4) - 0.5d0*rho*(u**2 + v**2))
         if (rho.eq.0.d0) write(6,*) 'rho = 0 in qstar'

         q(1,i,j) = q(1,i,j) - dt*(ndim-1)/rad * qstar(3)
         q(2,i,j) = q(2,i,j) - dt*(ndim-1)/rad * 
     &                          (rho*u*v)
         q(3,i,j) = q(3,i,j) - dt*(ndim-1)/rad * 
     &                          (rho*v*v)
         q(4,i,j) = q(4,i,j) - dt*(ndim-1)/rad * 
     &                          v*(qstar(4) + press)
c      end forall
   10    continue
      return
      end
