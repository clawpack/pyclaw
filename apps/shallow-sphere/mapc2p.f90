!     ==========================================
    subroutine mapc2p(x1,y1,xp,yp,zp,Rsphere)
!     ==========================================

!     Map computational point (x1,y1) to the point (xp,yp,zp)
!     on the sphere.

!     Points in [-1,1]x[-1,1]  are mapped to the upper hemisphere
!     Points in [-3,-1]x[-1,1] are mapped to the lower hemisphere

    implicit real*8(a-h,o-z)
!f2py  double precision intent(in) x1
!f2py  double precision intent(in) y1
!f2py  double precision intent(in) xp
!f2py  double precision intent(in) yp
!f2py  double precision intent(in) zp
!f2py  double precision intent(in) Rsphere


    r1 = Rsphere
    xc = x1
    yc = y1
    pi = 4.d0*datan(1.d0)

!     # ghost cell values outside of [-3,1]x[-1,1] get mapped to other
!     # hemisphere:
    if (xc >= 1.d0) xc = xc-4.d0
    if (xc < -3.d0) xc = xc+4.d0

    if (yc >= 1.d0) then
        yc = 2.d0 - yc
        xc = -2.d0 - xc
    endif
          
    if (yc < -1.d0) then
        yc = -2.d0 - yc
        xc = -2.d0 - xc
    endif

    if (xc < -1.d0) then
    !         # points in [-3,-1] map to lower hemisphere - reflect about x=-1
    !         # to compute x,y mapping and set sgnz appropriately:
        xc = -2.d0 - xc
        sgnz = -1.d0
    else
        sgnz = 1.d0
    endif
          
    sgnxc = dsign(1.d0,xc)
    sgnyc = dsign(1.d0,yc)

    xc1 = dabs(xc)
    yc1 = dabs(yc)
    d = dmax1(dmax1(xc1,yc1), 1.d-10)

    DD = r1*d*(2.d0 - d) / dsqrt(2.d0)
    R = r1

    centers = DD - dsqrt(dmax1(R**2 - DD**2, 0.d0))
    xp = DD/d * xc1
    yp = DD/d * yc1

    if (yc1 > xc1) then
        yp = centers + dsqrt(dmax1(R**2 - xp**2, 0.d0))
    else
        xp = centers + dsqrt(dmax1(R**2 - yp**2, 0.d0))
    endif

    zp = dsqrt(dmax1(r1**2 - (xp**2 + yp**2), 0.d0))
    xp = xp*sgnxc
    yp = yp*sgnyc
    zp = zp*sgnz

    return
    end subroutine mapc2p

