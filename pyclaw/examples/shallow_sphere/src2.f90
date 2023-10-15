!      =======================================================
    subroutine src2(maxmx,maxmy,num_eqn,num_ghost,mx,my,xlower,ylower, &
    dx,dy,q,num_aux,aux,t,dt,Rsphere)
!      =======================================================


!     # Compute source term for Rossby-Haurwitz wave.
!     # The source term models the Coriolis force using a 4-stage RK method
!     # and the projection of the velocity components to the tangent plane.


    implicit double precision (a-h,o-z)
    dimension    q(num_eqn, 1:mx, 1:my)
    dimension aux(num_aux,  1:mx, 1:my)
    double precision :: RK(4,3)

!f2py integer intent(in) maxmx
!f2py integer intent(in) maxmy
!f2py integer optional, intent(in) num_eqn
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
!f2py double precision intent(in) t
!f2py double precision intent(in) dt
!f2py double precision intent(in) Rsphere


!     # Parameters
    a = 6.37122d6
    Omega = 7.292d-5
    df=12.600576e0

!     # project momentum components of q onto tangent plane:

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

!     # calculate Coriolis term
    do i=1, mx
        xc = xlower + (i-0.5d0)*dx
        do j=1, my
            yc = ylower + (j-0.5d0)*dy
        
            call mapc2p(xc,yc,xp,yp,zp,Rsphere)
            erx = xp
            ery = yp
            erz = zp
        
            rad = dsqrt(xp**2. + yp**2.)
            if (rad > 1.d-4) then
                theta = dacos(xp/rad)
            elseif (yp > 0.d0) then
                theta = 0.d0
            endif
            if (yp < 0.d0) then
                theta = -theta
            endif

        !           # compute phi, at north pole: pi/2 at south pool: -pi/2
            if (zp > 0.d0) then
                phi =   dacos(rad/Rsphere)
            else
                phi = - dacos(rad/Rsphere)
            endif
        
            fcor = df*erz

        !           stage 1
            hu = q(2,i,j)
            hv = q(3,i,j)
            hw = q(4,i,j)
                                          
            RK(1,1) = fcor*dt*(erz*hv-ery*hw)
            RK(1,2) = dt*fcor*(erx*hw-erz*hu)
            RK(1,3) = dt*fcor*(ery*hu-erx*hv)

        !           stage 2
            hu = q(2,i,j) + 0.5d0*RK(1,1)
            hv = q(3,i,j) + 0.5d0*RK(1,2)
            hw = q(4,i,j) + 0.5d0*RK(1,3)

            RK(2,1) = fcor*dt*(erz*hv-ery*hw)
            RK(2,2) = dt*fcor*(erx*hw-erz*hu)
            RK(2,3) = dt*fcor*(ery*hu-erx*hv)

        !           stage 3
            hu = q(2,i,j) + 0.5d0*RK(2,1)
            hv = q(3,i,j) + 0.5d0*RK(2,2)
            hw = q(4,i,j) + 0.5d0*RK(2,3)

            RK(3,1) = fcor*dt*(erz*hv-ery*hw)
            RK(3,2) = dt*fcor*(erx*hw-erz*hu)
            RK(3,3) = dt*fcor*(ery*hu-erx*hv)

        !           stage 4
            hu = q(2,i,j) + 0.5d0*RK(3,1)
            hv = q(3,i,j) + 0.5d0*RK(3,2)
            hw = q(4,i,j) + 0.5d0*RK(3,3)

            RK(4,1) = fcor*dt*(erz*hv-ery*hw)
            RK(4,2) = dt*fcor*(erx*hw-erz*hu)
            RK(4,3) = dt*fcor*(ery*hu-erx*hv)
                        
            do m=2,num_eqn
                q(m,i,j) = q(m,i,j) &
                + (RK(1,m-1) + 2.d0*RK(2,m-1)+ &
                2.d0*RK(3,m-1) + RK(4,m-1))/6.d0
            enddo

        enddo
    enddo
     
!     # project momentum components of q onto tangent plane:

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
    end subroutine src2
