c      =======================================================
       subroutine src2(maxmx,maxmy,meqn,mbc,mx,my,xlower,ylower,
     &                 dx,dy,q,maux,aux,t,dt,Rsphere)
c      =======================================================
c
c
c     # Compute source term for Rossby-Haurwitz wave.
c     # The source term models the Coriolis force using a 4-stage RK method
c     # and the projection of the velocity components to the tangent plane.
c

       implicit double precision (a-h,o-z)
       dimension    q(meqn, 1:mx, 1:my)
       dimension aux(maux,  1:mx, 1:my)
       double precision RK(4,3) 

cf2py integer intent(in) maxmx
cf2py integer intent(in) maxmy
cf2py integer optional, intent(in) meqn
cf2py integer intent(in) mbc
cf2py integer intent(in) mx
cf2py integer intent(in) my
cf2py double precision intent(in) xlower
cf2py double precision intent(in) ylower
cf2py double precision intent(in) dx
cf2py double precision intent(in) dy
cf2py intent(in,out) q
cf2py integer optional, intent(in)  maux
cf2py intent(in) aux
cf2py double precision intent(in) t
cf2py double precision intent(in) dt
cf2py double precision intent(in) Rsphere


c     # Parameters
      a = 6.37122d6    
      Omega = 7.292d-5  
      df=12.600576e0 

c     # project momentum components of q onto tangent plane:

      do i=1,mx
        do j=1, my
            erx = aux(14,i,j)
            ery = aux(15,i,j)
            erz = aux(16,i,j)
            qn = erx*q(2,i,j) + ery*q(3,i,j) + erz*q(4,i,j)

            q(2,i,j) = q(2,i,j) - qn*erx
            q(3,i,j) = q(3,i,j) - qn*ery
            q(4,i,j) = q(4,i,j) - qn*erz

        enddo
      enddo

c     # calculate Coriolis term
      do i=1, mx
        xc = xlower + (i-0.5d0)*dx
        do j=1, my
            yc = ylower + (j-0.5d0)*dy
c
            call mapc2p(xc,yc,xp,yp,zp,Rsphere)
            erx = xp 
            ery = yp 
            erz = zp 
c
            rad = dsqrt(xp**2. + yp**2.)
            if (rad.gt.1.d-4) then
              theta = dacos(xp/rad)
            elseif (yp.gt.0.d0) then
              theta = 0.d0
            endif
            if (yp.lt.0.d0) then
              theta = -theta
            endif

c           # compute phi, at north pole: pi/2 at south pool: -pi/2
            if (zp.gt. 0.d0) then
              phi =   dacos(rad/Rsphere)
            else
              phi = - dacos(rad/Rsphere)
            endif    
c
            fcor = df*erz 

c           stage 1
            hu = q(2,i,j)
            hv = q(3,i,j)
            hw = q(4,i,j)
                              
            RK(1,1) = fcor*dt*(erz*hv-ery*hw)
            RK(1,2) = dt*fcor*(erx*hw-erz*hu)
            RK(1,3) = dt*fcor*(ery*hu-erx*hv)

c           stage 2
            hu = q(2,i,j) + 0.5d0*RK(1,1)
            hv = q(3,i,j) + 0.5d0*RK(1,2)
            hw = q(4,i,j) + 0.5d0*RK(1,3)

            RK(2,1) = fcor*dt*(erz*hv-ery*hw)
            RK(2,2) = dt*fcor*(erx*hw-erz*hu)
            RK(2,3) = dt*fcor*(ery*hu-erx*hv)

c           stage 3
            hu = q(2,i,j) + 0.5d0*RK(2,1)
            hv = q(3,i,j) + 0.5d0*RK(2,2)
            hw = q(4,i,j) + 0.5d0*RK(2,3)

            RK(3,1) = fcor*dt*(erz*hv-ery*hw)
            RK(3,2) = dt*fcor*(erx*hw-erz*hu)
            RK(3,3) = dt*fcor*(ery*hu-erx*hv)

c           stage 4
            hu = q(2,i,j) + 0.5d0*RK(3,1)
            hv = q(3,i,j) + 0.5d0*RK(3,2)
            hw = q(4,i,j) + 0.5d0*RK(3,3)

            RK(4,1) = fcor*dt*(erz*hv-ery*hw)
            RK(4,2) = dt*fcor*(erx*hw-erz*hu)
            RK(4,3) = dt*fcor*(ery*hu-erx*hv)
            
            do m=2,meqn
               q(m,i,j) = q(m,i,j) 
     &                 + (RK(1,m-1) + 2.d0*RK(2,m-1)+
     &                 2.d0*RK(3,m-1) + RK(4,m-1))/6.d0
            enddo

        enddo
      enddo
 
c     # project momentum components of q onto tangent plane:

      do i=1,mx
        do j=1, my
            erx = aux(14,i,j)
            ery = aux(15,i,j)
            erz = aux(16,i,j)
            qn = erx*q(2,i,j) + ery*q(3,i,j) + erz*q(4,i,j)
            q(2,i,j) = q(2,i,j) - qn*erx
            q(3,i,j) = q(3,i,j) - qn*ery
            q(4,i,j) = q(4,i,j) - qn*erz

        enddo
      enddo

      return
      end
