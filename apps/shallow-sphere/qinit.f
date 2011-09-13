c
c
c
c     =====================================================
       subroutine qinit(maxmx,maxmy,meqn,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux,Rsphere)
c     =====================================================
c
c      # Set initial conditions for q 

c      # -------4-Rossby-Haurwitz wave-----------------------
c
       implicit double precision (a-h,o-z)
       dimension q(meqn,1:maxmx, 1:maxmy)
       dimension aux(maux, 1:maxmx, 1:maxmy)
       double precision Uin(3),Uout(3)
       double precision K, Rsphere
cf2py intent(in,out) q
cf2py integer optional,intent(in) maxmx
cf2py integer optional,intent(in) maxmy
cf2py integer optional, intent(in) meqn
cf2py integer intent(in) mx
cf2py integer intent(in) my
cf2py double precision intent(in) xlower
cf2py double precision intent(in) ylower
cf2py double precision intent(in) dx
cf2py double precision intent(in) dy
cf2py integer optional, intent(in)  maux
cf2py intent(in) aux
cf2py double precision intent(in) Rsphere



c
       pi = 4.d0*datan(1.d0)

       a = 6.37122d6    
       K = 7.848d-6     
       Omega = 7.292d-5 
       G = 9.80616d0    
       t0 = 86400.d0    
       h0 = 8.d3        
       R = 4.d0

      
       do 20 i=1,mx
        xc = xlower + (i-0.5d0)*dx
          do 20 j=1,my
        yc = ylower + (j-0.5d0)*dy
        call mapc2m(xc,yc,xp,yp,zp)
c            # compute longitude theta from positive x axis:
             rad = dmax1(dsqrt(xp**2 + yp**2),1.d-6)
             theta = 0.d0

              if(xp.gt.0.d0.and.yp.gt.0.d0) then
                 theta = dasin(yp/rad) 
              elseif(xp.lt.0.d0.and.yp.gt.0.d0) then
                 theta = pi - dasin(yp/rad)
              elseif(xp.lt.0.d0.and.yp.lt.0.d0) then
                 theta = -pi+dasin(-yp/rad)
              elseif(xp.gt.0.d0.and.yp.lt.0.d0) then
                 theta = -dasin(-yp/rad)
              endif 

c            # compute phi, at north pole: pi/2 at south pool: -pi/2
             if (zp.gt. 0.d0) then
               phi =  dasin(zp/Rsphere) 
             else
               phi = -dasin(-zp/Rsphere)  
             endif 
c
             xp = theta 
             yp = phi 
c
             bigA = 0.5d0*K*(2.d0*Omega+K)*dcos(yp)**2.d0 + 
     &               0.25d0*K*K*dcos(yp)**(2.d0*R)*( 
     &               (1.d0*R+1.d0)*dcos(yp)**2.d0
     &               +(2.d0*R*R - 1.d0*R - 2.d0) 
     &               - 2.d0*R*R*(dcos(yp))**(-2.d0))
             bigB = (2.d0*(Omega+K)*K)/((1.d0*R+1.d0)*(1.d0*R+2.d0))
     &               *dcos(yp)**R*( (1.d0*R*R + 2.d0*R + 2.d0) 
     &               - (1.d0*R+1.d0)**(2)*dcos(yp)**2 )
             bigC = 0.25d0*K*K*dcos(yp)**(2*R)*( (1.d0*R + 1.d0)*
     &               dcos(yp)**2 - (1.d0*R + 2.d0))
             
c            # longitude (angular) velocity component
             Uin(1) = (K*dcos(yp)+K*dcos(yp)**(R-1.)*( R*dsin(yp)**2. 
     &              - dcos(yp)**2. )*dcos(R*xp))*t0

c            # latitude (angular) velocity component
             Uin(2) = (-K*R*dcos(yp)**(R-1.)*dsin(yp)*dsin(R*xp))*t0

c            # radial velocity component
             Uin(3) = 0.d0
c

c            # calculate velocity vetor in cartesian coordinates
             Uout(1) = (-dsin(xp)*Uin(1)-dsin(yp)*dcos(xp)*Uin(2))
             Uout(2) = (dcos(xp)*Uin(1)-dsin(yp)*dsin(xp)*Uin(2))
             Uout(3) = dcos(yp)*Uin(2)
c
c            # set the clawpack initial values:             
             q(1,i,j) =  h0/a + (a/G)*( bigA + bigB*dcos(R*xp) 
     &              + bigC*dcos(2.d0*R*xp))
             q(2,i,j) = q(1,i,j)*Uout(1) 
             q(3,i,j) = q(1,i,j)*Uout(2) 
             q(4,i,j) = q(1,i,j)*Uout(3) 

            
  20         continue


       return
       end
