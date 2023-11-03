! ============================================
subroutine setaux(maxmx,maxmy,num_ghost,mx,my,xlower,ylower,dxc,dyc, &
                  num_aux,aux,Rsphere)
! ============================================

! Set auxiliary arrays for shallow water equations on the sphere.

! On input: (xc(i),yc(j)) gives uniformly spaced computational grid.
! On output: he aux array has the following elements:
!         1  kappa = ratio of cell area to dxc*dyc
!         2  enx = x-component of normal vector to left edges in tangent plane
!         3  eny = y-component of normal vector to left edges in tangent plane
!         4  enz = z-component of normal vector to left edges in tangent plane
!         5  etx = x-component of tangent vector to left edges in tangent plane
!         6  ety = y-component of tangent vector to left edges in tangent plane
!         7  etz = z-component of tangent vector to left edges in tangent plane
!         8  enx = x-component of normal vector to bottom edges in tangent plane
!         9  eny = y-component of normal vector to bottom edges in tangent plane
!        10  enz = z-component of normal vector to bottom edges in tangent plane
!        11  etx = x-component of tangent vector to bottom edges in tangent plane
!        12  ety = y-component of tangent vector to bottom edges in tangent plane
!        13  etz = z-component of tangent vector to bottom edges in tangent plane
!        14  erx = x-component of unit vector in radial direction at cell ctr
!        15  ery = y-component of unit vector in radial direction at cell ctr
!        16  erz = z-component of unit vector in radial direction at cell ctr


    implicit double precision (a-h,o-z)
    !parameter (maxm3 = 1005)
    dimension xc(-num_ghost:maxmx+num_ghost), yc(-num_ghost:maxmx+num_ghost)
    dimension xp(-num_ghost:maxmx+num_ghost,-num_ghost:maxmx+num_ghost)
    dimension yp(-num_ghost:maxmx+num_ghost,-num_ghost:maxmx+num_ghost)
    dimension zp(-num_ghost:maxmx+num_ghost,-num_ghost:maxmx+num_ghost)
    dimension theta(-num_ghost:maxmx+num_ghost,-num_ghost:maxmx+num_ghost)
    dimension phi(-num_ghost:maxmx+num_ghost,-num_ghost:maxmx+num_ghost)
    dimension aux(num_aux,1-num_ghost:maxmx+num_ghost,1-num_ghost:maxmy+num_ghost)

!f2py integer intent(in) maxmx
!f2py integer intent(in) maxmy
!f2py integer intent(in) num_ghost
!f2py integer intent(in) mx
!f2py integer intent(in) my
!f2py double precision intent(in) xlower
!f2py double precision intent(in) ylower
!f2py double precision intent(in) dx
!f2py double precision intent(in) dy
!f2py integer optional, intent(in)  num_aux
!f2py intent(in,out) aux
!f2py double precision intent(in) Rsphere

    pi = 4.d0*datan(1.d0)

    ! Set xc and yc so that (xc(i),yc(j)) is the
    ! lower left corner of (i,j) cell in computational space:

    do i=1-num_ghost,mx+num_ghost+1
        xc(i) = xlower + (i-1.d0) * dxc
    end do

    do j=1-num_ghost,my+num_ghost+1
        yc(j) = ylower + (j-1.d0) * dyc
    end do

    ! compute cell corners on sphere and angles phi, theta
    ! related to latitude and longitude

    do j=1-num_ghost,my+num_ghost+1
        do i=1-num_ghost,mx+num_ghost+1

            ! map computational point to (xp,yp,zp) on sphere:
            call mapc2p(xc(i),yc(j),xp(i,j),yp(i,j),zp(i,j),Rsphere)

            ! compute longitude theta from positive x axis:
            r = dsqrt(xp(i,j)**2 + yp(i,j)**2)

            ! theta(i,j) = atan2(xp(i,j),yp(i,j))

            if (r > 1.d-4) then
                theta(i,j) = dacos(xp(i,j)/r)
            elseif (yp(i,j) > 0.d0) then
                theta(i,j) = 0.d0
            endif
            if (yp(i,j) < 0.d0) then
                theta(i,j) = -theta(i,j)
            endif
            ! compute phi, angle down from north pole:
            if (zp(i,j) > 0.d0) then
                phi(i,j) = pi/2.d0 - dacos(r/Rsphere)
            else
                phi(i,j) = pi/2.d0 + dacos(r/Rsphere)
            endif
        end do
    end do


    do j=1-num_ghost,my+num_ghost
        do i=1-num_ghost,mx+num_ghost
        
            ! compute normal and tangent vectors to left edges (in tangent plane)
        
            ! tangent vector is edges vector:
            etx = xp(i,j+1) - xp(i,j)
            ety = yp(i,j+1) - yp(i,j)
            etz = zp(i,j+1) - zp(i,j)
            aux(5,i,j) = etx
            aux(6,i,j) = ety
            aux(7,i,j) = etz

            ! normal to sphere in radial direction at midpoint of edges:
            erx = 0.5d0*(xp(i,j) + xp(i,j+1))
            ery = 0.5d0*(yp(i,j) + yp(i,j+1))
            erz = 0.5d0*(zp(i,j) + zp(i,j+1))
        
            ! normal to edges in tangent plane is cross product of et and er:
            enx = ety*erz - etz*ery
            eny = etz*erx - etx*erz
            enz = etx*ery - ety*erx
            ennorm = dsqrt(enx**2 + eny**2 + enz**2)
            aux(2,i,j) = enx/ennorm
            aux(3,i,j) = eny/ennorm
            aux(4,i,j) = enz/ennorm


            ! compute normal and tangent vectors to bottom edges (in tang. pl.)
        
            ! tangent vector is edges vector:
            etx = xp(i+1,j) - xp(i,j)
            ety = yp(i+1,j) - yp(i,j)
            etz = zp(i+1,j) - zp(i,j)
            aux(11,i,j) = etx
            aux(12,i,j) = ety
            aux(13,i,j) = etz

            ! normal to sphere in radial direction at midpoint of edges:
            erx = 0.5d0*(xp(i,j) + xp(i+1,j))
            ery = 0.5d0*(yp(i,j) + yp(i+1,j))
            erz = 0.5d0*(zp(i,j) + zp(i+1,j))
        
            ! normal to edges in tangent plane is cross product of er and et:
            enx = ery*etz - erz*ety
            eny = erz*etx - erx*etz
            enz = erx*ety - ery*etx
            ennorm = dsqrt(enx**2 + eny**2 + enz**2)
            aux(8,i,j) = enx/ennorm
            aux(9,i,j) = eny/ennorm
            aux(10,i,j) = enz/ennorm

            ! normal to sphere in radial direction at cell centers:
            xcm = xlower+(i-0.5)*dxc
            ycm = ylower+(j-0.5)*dyc
            call mapc2p(xcm,ycm,xpm,ypm,zpm,Rsphere)

            aux(14,i,j) = xpm
            aux(15,i,j) = ypm
            aux(16,i,j) = zpm

            ! compute area of physical cell from four corners:
            ! find area on the sphere of two spherical triangles obtained
            ! by subdividing rectangle.  See
            !    http://mathforum.org/library/drmath/view/65316.html

            ! corners are labeled
            !      1: (i,j)      2: (i+1,j)    3: (i,j+1)    4: (i+1,j+1)

            beta12 = dsin(phi(i,j))*dsin(phi(i+1,j))* &
                dcos(theta(i,j) - theta(i+1,j)) &
                + dcos(phi(i,j))*dcos(phi(i+1,j))
            beta23 = dsin(phi(i,j+1))*dsin(phi(i+1,j))* &
                dcos(theta(i,j+1) - theta(i+1,j)) &
                + dcos(phi(i,j+1))*dcos(phi(i+1,j))
            beta13 = dsin(phi(i,j+1))*dsin(phi(i,j))* &
                dcos(theta(i,j+1) - theta(i,j)) &
                + dcos(phi(i,j+1))*dcos(phi(i,j))
            beta24 = dsin(phi(i+1,j+1))*dsin(phi(i+1,j))* &
                dcos(theta(i+1,j+1) - theta(i+1,j)) &
                + dcos(phi(i+1,j+1))*dcos(phi(i+1,j))
            beta34 = dsin(phi(i+1,j+1))*dsin(phi(i,j+1))* &
                dcos(theta(i+1,j+1) - theta(i,j+1)) &
                + dcos(phi(i+1,j+1))*dcos(phi(i,j+1))

            ! great circles distances between corners:
            d12 = Rsphere * dacos(beta12)
            d23 = Rsphere * dacos(beta23)
            d13 = Rsphere * dacos(beta13)
            d24 = Rsphere * dacos(beta24)
            d34 = Rsphere * dacos(beta34)

            s123 = 0.5d0 * (d12 + d23 + d13)
            s234 = 0.5d0 * (d23 + d34 + d24)

            ! spherical excess for each triangle:
            t123 = dtan(s123/2.d0)*dtan((s123-d12)/2.d0) &
                *dtan((s123-d23)/2.d0)*dtan((s123-d13)/2.d0)
            t123 = dmax1(t123,0.d0)
            E123 = 4.d0*datan(sqrt(t123))

            t234 = dtan(s234/2.d0)*dtan((s234-d23)/2.d0) &
                *dtan((s234-d34)/2.d0)*dtan((s234-d24)/2.d0)
            t234 = dmax1(t234,0.d0)
            E234 = 4.d0*datan(sqrt(t234))

            ! = Rsphere**2 * (E123 + E234)
            area = (E123 + E234)
        
            ! capacity kappa:
            aux(1,i,j) = area / (dxc*dyc)
                       
        end do
    end do

    return
end subroutine setaux
