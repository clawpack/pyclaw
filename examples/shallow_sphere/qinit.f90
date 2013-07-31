!     =====================================================
    subroutine qinit(maxmx,maxmy,num_eqn,num_ghost,mx,my,xlower,ylower, &
    dx,dy,q,num_aux,aux,Rsphere)
!     =====================================================

!      # Set initial conditions for q.

!      # -------4-Rossby-Haurwitz wave-----------------------

    implicit double precision (a-h,o-z)
    dimension q(num_eqn, 1-num_ghost:maxmx+num_ghost, 1-num_ghost:maxmy+num_ghost)
    dimension aux(num_aux, 1-num_ghost:maxmx+num_ghost, 1-num_ghost:maxmy+num_ghost)
    double precision :: Uin(3),Uout(3)
    double precision :: K
!f2py integer intent(in) maxmx
!f2py integer intent(in) maxmy
!f2py integer optional,intent(in) num_eqn
!f2py integer intent(in) num_ghost
!f2py integer intent(in) mx
!f2py integer intent(in) my
!f2py double precision intent(in) xlower
!f2py double precision intent(in) ylower
!f2py double precision intent(in) dx
!f2py double precision intent(in) dy
!f2py intent(in,out) q
!f2py integer optional, intent(in)  num_aux
!f2py intent(in) aux
!f2py double precision intent(in) Rsphere

    pi = 4.d0*datan(1.d0)

!      # Parameters
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
            call mapc2p(xc,yc,xp,yp,zp,Rsphere)
        !            # compute longitude theta from positive x axis:
            rad = dmax1(dsqrt(xp**2 + yp**2),1.d-6)

            if(xp > 0.d0 .AND. yp > 0.d0) then
                theta = dasin(yp/rad)
            elseif(xp < 0.d0 .AND. yp > 0.d0) then
                theta = pi - dasin(yp/rad)
            elseif(xp < 0.d0 .AND. yp < 0.d0) then
                theta = -pi+dasin(-yp/rad)
            elseif(xp > 0.d0 .AND. yp < 0.d0) then
                theta = -dasin(-yp/rad)
            endif

        !            # compute phi, at north pole: pi/2 at south pool: -pi/2
            if (zp > 0.d0) then
                phi =  dasin(zp/Rsphere)
            else
                phi = -dasin(-zp/Rsphere)
            endif
        
            xp = theta
            yp = phi
        
            bigA = 0.5d0*K*(2.d0*Omega+K)*dcos(yp)**2.d0 + &
            0.25d0*K*K*dcos(yp)**(2.d0*R)*( &
            (1.d0*R+1.d0)*dcos(yp)**2.d0 &
            +(2.d0*R*R - 1.d0*R - 2.d0) &
            - 2.d0*R*R*(dcos(yp))**(-2.d0))
            bigB = (2.d0*(Omega+K)*K)/((1.d0*R+1.d0)*(1.d0*R+2.d0)) &
            *dcos(yp)**R*( (1.d0*R*R + 2.d0*R + 2.d0) &
            - (1.d0*R+1.d0)**(2)*dcos(yp)**2 )
            bigC = 0.25d0*K*K*dcos(yp)**(2*R)*( (1.d0*R + 1.d0)* &
            dcos(yp)**2 - (1.d0*R + 2.d0))
                         
        !            # longitude (angular) velocity component
            Uin(1) = (K*dcos(yp)+K*dcos(yp)**(R-1.)*( R*dsin(yp)**2. &
            - dcos(yp)**2. )*dcos(R*xp))*t0

        !            # latitude (angular) velocity component
            Uin(2) = (-K*R*dcos(yp)**(R-1.)*dsin(yp)*dsin(R*xp))*t0

        !            # radial velocity component
            Uin(3) = 0.d0
        

        !            # calculate velocity vetor in cartesian coordinates
            Uout(1) = (-dsin(xp)*Uin(1)-dsin(yp)*dcos(xp)*Uin(2))
            Uout(2) = (dcos(xp)*Uin(1)-dsin(yp)*dsin(xp)*Uin(2))
            Uout(3) = dcos(yp)*Uin(2)
        
        !            # set the clawpack initial values:
            q(1,i,j) =  h0/a + (a/G)*( bigA + bigB*dcos(R*xp) &
            + bigC*dcos(2.d0*R*xp))
            q(2,i,j) = q(1,i,j)*Uout(1)
            q(3,i,j) = q(1,i,j)*Uout(2)
            q(4,i,j) = q(1,i,j)*Uout(3)

                        
    20 END DO

    return
    end subroutine qinit
