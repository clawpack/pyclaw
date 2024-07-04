MODULE euler3d_mappedgrid

  REAL(kind=8), PARAMETER :: pi=3.1415926535897932385,twopi=2*pi
  REAL(kind=8), PARAMETER :: rEarth=637120000.0

CONTAINS
  !-------------------------------------------------------------------------
  ! Name:
  !   volumeOfHex(x,y,z):
  !
  ! Description:
  !   Compute Volume of an LD Hexahedron given the 8 vertices
  !    from paper: 
  !    @book{Grandy_1997, 
  !    title={Efficient computation of volume of hexahedral cells}, 
  !    url={http://www.osti.gov/scitech/servlets/purl/632793}, 
  !    DOI={10.2172/632793}, 
  !    abstractNote={This report describes an efficient method to compute 
  !                  the volume of hexahedral cells used in three-dimensional 
  !                  hydrodynamics simulation. Two common methods for 
  !                  creating the hexahedron using triangular boundaries are 
  !                  considered.}, 
  !    author={Grandy, J.}, year={1997}, month={Oct}}
  !    https://doi.org/10.2172/632793
  !
  !    6*V_LD = 
  ! |x6-x0  x3-x0  x2-x7|   |x6-x0  x4-x0  x7-x5|   |x6-x0  x1-x0  x5-x2|
  ! |y6-y0  y3-y0  y2-y7| + |y6-y0  y4-y0  y7-y5| + |y6-y0  y1-y0  y5-y2|
  ! |z6-z0  z3-z0  z2-z7|   |z6-z0  z4-z0  z7-z5|   |z6-z0  z1-z0  z5-z2|
  !
  ! Inputs:
  !   x[:],y[:],z[:]  : x,y, and z coordinates of eight hexahedron vertices
  ! 
  ! Output:
  !   volumeOfHex     : Volume of LD hexahedron
  !-------------------------------------------------------------------------
  REAL(kind=8) FUNCTION volumeOfHex(x,y,z)

    IMPLICIT NONE

    REAL(kind=8), DIMENSION(8), INTENT(IN) :: x,y,z
    
    ! Compute Volume of Hexahedron
    volumeOfHex = &
         (x(7)-x(1))*((y(4)-y(1))*(z(3)-z(8)) - (y(3)-y(8))*(z(4)-z(1)))+&
         (x(4)-x(1))*((y(3)-y(8))*(z(7)-z(1)) - (y(7)-y(1))*(z(3)-z(8)))+&
         (x(3)-x(8))*((y(7)-y(1))*(z(4)-z(1)) - (y(4)-y(1))*(z(7)-z(1)))&
         +&
         (x(7)-x(1))*((y(5)-y(1))*(z(8)-z(6)) - (y(8)-y(6))*(z(5)-z(1)))+&
         (x(5)-x(1))*((y(8)-y(6))*(z(7)-z(1)) - (y(7)-y(1))*(z(8)-z(6)))+&
         (x(8)-x(6))*((y(7)-y(1))*(z(5)-z(1)) - (y(5)-y(1))*(z(7)-z(1)))&
         +&
         (x(7)-x(1))*((y(2)-y(1))*(z(6)-z(3)) - (y(6)-y(3))*(z(2)-z(1)))+&
         (x(2)-x(1))*((y(6)-y(3))*(z(7)-z(1)) - (y(7)-y(1))*(z(6)-z(3)))+&
         (x(6)-x(3))*((y(7)-y(1))*(z(2)-z(1)) - (y(2)-y(1))*(z(7)-z(1)))
    volumeOfHex = ABS(volumeOfHex)/6.d0
    RETURN
  END FUNCTION volumeOfHex  
  
  !---------------------------------------------------------------------------
  ! Name: mapc2pWrapper(xc,yc,zc,mx,my,mz,xyz0,xyzN,mapType,xp,yp,zp)
  !
  ! Description: mapping from computational to physical space
  !
  ! Inputs: 
  !        xc,yc,zc <REAL>     : computational coordinate values
  !        mx,my,mz <INTEGER>  : grid sizes
  !        xyz0,xyzN <REAL>    : altitude at bottom and top of grid
  !        mapType <CHARACTER> : type of mapping to compute         
  !
  ! Outputs: 
  !        xp,yp,zp <REAL>     : physical coordinate values
  !---------------------------------------------------------------------------
  SUBROUTINE mapc2pWrapper(xc,yc,zc,mz,xyz0,xyzN,mapType,xp,yp,zp)

    IMPLICIT NONE
    
    ! Input
    REAL(kind=8), DIMENSION(:), INTENT(IN) :: xc,yc,zc
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: xyz0,xyzN
    INTEGER, INTENT(IN) :: mz
    CHARACTER(20) :: mapType
    
    ! Output
    REAL(kind=8), DIMENSION(SIZE(xc)), INTENT(OUT) :: xp
    REAL(kind=8), DIMENSION(SIZE(yc)), INTENT(OUT) :: yp
    REAL(kind=8), DIMENSION(SIZE(zc)), INTENT(OUT) :: zp
    
    ! Compute Physical Grid From Computational
    SELECT CASE(mapType)
    CASE("Identity")
       CALL mapc2pIdentity(xc,yc,zc,xp,yp,zp)
    CASE("ZeroToOne")
       CALL mapc2pZeroToOne(xc,yc,zc,xyz0,xyzN,xp,yp,zp)
    CASE("Rotate")
       CALL mapc2pRotate(xc,yc,zc,xyz0,xyzN,xp,yp,zp)
    CASE("Spherical")
       CALL mapc2pSpherical(xc,yc,zc,xyz0,xyzN,xp,yp,zp)
    CASE DEFAULT
       WRITE(*,*)"Invalid mapType: ",mapType
       STOP
    END SELECT
    
    RETURN
  END SUBROUTINE mapc2pWrapper

  !---------------------------------------------------------------------------
  ! Name: mapc2pIdentity(xc,yc,zc,xp,yp,zp)
  !
  ! Description: Identity mapping from computational to physical space
  !
  ! Inputs: 
  !        xc,yc,zc <REAL>     : computational coordinate values
  ! Outputs: 
  !        xp,yp,zp <REAL>     : physical coordinate values
  !---------------------------------------------------------------------------
  SUBROUTINE mapc2pIdentity(xc,yc,zc,xp,yp,zp)

    IMPLICIT NONE
    
    ! Input
    REAL(kind=8), DIMENSION(:), INTENT(IN) :: xc,yc,zc
    
    ! Output
    REAL(kind=8), DIMENSION(SIZE(xc)), INTENT(OUT) :: xp
    REAL(kind=8), DIMENSION(SIZE(yc)), INTENT(OUT) :: yp
    REAL(kind=8), DIMENSION(SIZE(zc)), INTENT(OUT) :: zp
    
    ! x,xi coordinate
    xp = xc
    
    ! y,eta coordinate
    yp = yc
    
    ! z,phi coordinate
    zp = zc
    
    RETURN
  END SUBROUTINE mapc2pIdentity
      
  !---------------------------------------------------------------------------
  ! Name: mapc2pZeroToOne(xc,yc,zc,xyz0,xyzN,xp,yp,zp)
  !
  ! Description: mapping from computational to physical space using a
  !                computational grid from 0 to 1
  !
  ! Inputs: 
  !        xc,yc,zc <REAL>     : computational coordinate values
  !        xyz0,xyzN <REAL>    : physical coordinates at bottom and top of grid
  !
  ! Outputs: 
  !        xp,yp,zp <REAL>     : physical coordinate values
  !---------------------------------------------------------------------------
  SUBROUTINE mapc2pZeroToOne(xc,yc,zc,xyz0,xyzN,xp,yp,zp)

    IMPLICIT NONE
    
    ! Input
    REAL(kind=8), DIMENSION(:), INTENT(IN) :: xc,yc,zc
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: xyz0,xyzN
    
    ! Output
    REAL(kind=8), DIMENSION(SIZE(xc)), INTENT(OUT) :: xp
    REAL(kind=8), DIMENSION(SIZE(yc)), INTENT(OUT) :: yp
    REAL(kind=8), DIMENSION(SIZE(zc)), INTENT(OUT) :: zp
    
    ! x,xi coordinate
    xp = xyz0(1) + xc*(xyzN(1)-xyz0(1))
    
    ! y,eta coordinate
    yp = xyz0(2) + yc*(xyzN(2)-xyz0(2))
    
    ! z,phi coordinate
    zp = xyz0(3) + zc*(xyzN(3)-xyz0(3))
    
    RETURN
  END SUBROUTINE mapc2pZeroToOne
    
  !---------------------------------------------------------------------------
  ! Name: mapc2pRotate(xc,yc,zc,xyz0,xyzN,xp,yp,zp)
  !
  ! Description: mapping from computational to physical space using a
  !                computational grid from 0 to 1 and rotated by angle theta
  !
  ! Inputs: 
  !        xc,yc,zc <REAL>     : computational coordinate values
  !        xyz0,xyzN <REAL>    : physical coordinates at bottom and top of grid
  !
  ! Outputs: 
  !        xp,yp,zp <REAL>     : physical coordinate values
  !---------------------------------------------------------------------------
  SUBROUTINE mapc2pRotate(xc,yc,zc,xyz0,xyzN,xp,yp,zp)

    IMPLICIT NONE
    
    ! Input
    REAL(kind=8), DIMENSION(:), INTENT(IN) :: xc,yc,zc
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: xyz0,xyzN
    
    ! Output
    REAL(kind=8), DIMENSION(SIZE(xc)), INTENT(OUT) :: xp
    REAL(kind=8), DIMENSION(SIZE(yc)), INTENT(OUT) :: yp
    REAL(kind=8), DIMENSION(SIZE(zc)), INTENT(OUT) :: zp
    
    ! Local
    REAL(kind=8) :: theta
    REAL(kind=8), DIMENSION(SIZE(xc)) :: xp0
    REAL(kind=8), DIMENSION(SIZE(yc)) :: yp0
    REAL(kind=8), DIMENSION(SIZE(zc)) :: zp0
    
    theta=0.d0*pi/180.d0

    ! x,xi coordinate
    xp = xyz0(1) + xc*(xyzN(1)-xyz0(1))

    ! y,eta coordinate
    yp = xyz0(2) + yc*(xyzN(2)-xyz0(2))

    ! z,phi coordinate
    zp = xyz0(3) + zc*(xyzN(3)-xyz0(3))

    xp0 = xp
    yp0 = yp
    zp0 = zp    
    ! Rotate about x-axis by angle theta
    xp = xp0
    yp = yp0*COS(theta) - zp0*SIN(theta)
    zp = yp0*SIN(theta) + zp0*COS(theta)

    xp0 = xp
    yp0 = yp
    zp0 = zp
    ! Rotate about y-axis by angle theta
    xp = xp0*COS(theta) + zp0*SIN(theta)
    yp = yp0
    zp =-xp0*SIN(theta) + zp0*COS(theta)

    xp0 = xp
    yp0 = yp
    zp0 = zp    
    ! Rotate about z-axis by angle theta
    xp = xp0*COS(theta) - yp0*SIN(theta)
    yp = xp0*SIN(theta) + yp0*COS(theta)
    zp = zp0

    RETURN
  END SUBROUTINE mapc2pRotate
  
  !---------------------------------------------------------------------------
  ! Name: mapc2pSpherical(xc,yc,zc,xyz0,xyzN,xp,yp,zp)
  !
  ! Description: mapping from computational to physical space using a
  !                computational grid from 0 to 1 and then using a
  !                spherical transformation
  !
  ! Inputs: 
  !        xc,yc,zc <REAL>     : computational coordinate values
  !        xyz0,xyzN <REAL>    : physical coordinates at bottom and top of grid
  !
  ! Outputs: 
  !        xp,yp,zp <REAL>     : physical coordinate values
  !---------------------------------------------------------------------------
  SUBROUTINE mapc2pSpherical(xc,yc,zc,xyz0,xyzN,xp,yp,zp)

    IMPLICIT NONE

    ! Input
    REAL(kind=8), DIMENSION(:), INTENT(IN) :: xc,yc,zc
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: xyz0,xyzN

    ! Output
    REAL(kind=8), DIMENSION(SIZE(xc)), INTENT(OUT) :: xp
    REAL(kind=8), DIMENSION(SIZE(yc)), INTENT(OUT) :: yp
    REAL(kind=8), DIMENSION(SIZE(zc)), INTENT(OUT) :: zp

    ! Local
    REAL(kind=8), DIMENSION(SIZE(xc)) :: theta
    REAL(kind=8), DIMENSION(SIZE(yc)) :: phi
    REAL(kind=8), DIMENSION(SIZE(zc)) :: r

    ! x,xi coordinate
    xp = xyz0(1) + xc*(xyzN(1)-xyz0(1))

    ! y,eta coordinate
    yp = xyz0(2) + yc*(xyzN(2)-xyz0(2))

    ! z,phi coordinate
    zp = xyz0(3) + zc*(xyzN(3)-xyz0(3))

    ! Spherical to Cartesian Coordinates
    theta = xp
    phi = yp
    r = zp    
    xp = r*SIN(theta)*COS(phi)
    yp = r*SIN(theta)*SIN(phi)
    zp = r*COS(theta)

    RETURN
  END SUBROUTINE mapc2pSpherical

  !---------------------------------------------------------------------------
  ! Name: determinant(A,detA)
  !
  ! Description: Compute the determinant of a 3x3 array
  !
  ! Inputs:  A <REAL>    : 3x3 Array
  !         
  ! Outputs: detA <REAL> : determinant of A
  !---------------------------------------------------------------------------
  SUBROUTINE determinant(A,detA)

    IMPLICIT NONE

    ! Input
    REAL(kind=8), DIMENSION(3,3), INTENT(IN) :: A

    ! Output
    REAL(kind=8), INTENT(OUT) :: detA

    detA = A(1,1)*A(2,2)*A(3,3) &
         + A(2,1)*A(3,2)*A(1,3) &
         + A(3,1)*A(1,2)*A(2,3) &
         - A(3,1)*A(2,2)*A(1,3) &
         - A(2,1)*A(1,2)*A(3,3) &
         - A(1,1)*A(3,2)*A(2,3)

    RETURN
  END SUBROUTINE determinant

  !---------------------------------------------------------------------------
  ! Name: unitNormal(a,b,c,n)
  !
  ! Description: Compute the unit normal of given 3 indices in 3D space
  !
  ! Inputs:  a,b,c <REAL> : indices with x,y, and z coordinates
  !         
  ! Outputs: n <REAL>     : unit normals of plane defined from points a,b,c
  !---------------------------------------------------------------------------
  SUBROUTINE unitNormal(a,b,c,n)

    IMPLICIT NONE
    
    ! Input
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: a,b,c
    
    ! Output
    REAL(kind=8), DIMENSION(3), INTENT(OUT) :: n

    ! Local
    REAL(kind=8), DIMENSION(3,3) :: D
    REAL(kind=8) :: magnitude

    D(:,1) = (/1.d0,a(2),a(3)/)
    D(:,2) = (/1.d0,b(2),b(3)/)
    D(:,3) = (/1.d0,c(2),c(3)/)
    CALL determinant(D,n(1))
    D(:,1) = (/a(1),1.d0,a(3)/)
    D(:,2) = (/b(1),1.d0,b(3)/)
    D(:,3) = (/c(1),1.d0,c(3)/)
    CALL determinant(D,n(2))
    D(:,1) = (/a(1),a(2),1.d0/)
    D(:,2) = (/b(1),b(2),1.d0/)
    D(:,3) = (/c(1),c(2),1.d0/)
    CALL determinant(D,n(3))

    magnitude = SQRT(n(1)**2 + n(2)**2 + n(3)**2)
    n(1) = n(1)/magnitude
    n(2) = n(2)/magnitude
    n(3) = n(3)/magnitude

    RETURN
  END SUBROUTINE unitNormal

  !---------------------------------------------------------------------------
  ! Name: dot(a,b,adotb)
  !
  ! Description: Compute the dot product of vectors a and b
  !
  ! Inputs:  a,b <REAL>   : 3D vectors
  !         
  ! Outputs: adotb <REAL> : dot product of a and b
  !---------------------------------------------------------------------------
  SUBROUTINE dot(a,b,adotb)

    IMPLICIT NONE
    
    ! Input
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: a,b
    
    ! Output
    REAL(kind=8), INTENT(OUT) :: adotb

    adotb = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

    RETURN
  END SUBROUTINE dot

  !---------------------------------------------------------------------------
  ! Name: cross(a,b,c)
  !
  ! Description: Compute the cross product of vectors a and b
  !
  ! Inputs:  a,b <REAL> : 3D vectors
  !         
  ! Outputs: c <REAL>   : cross product of a and b
  !---------------------------------------------------------------------------
  SUBROUTINE cross(a,b,c)

    IMPLICIT NONE
    
    ! Input
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: a,b
    
    ! Output
    REAL(kind=8), DIMENSION(3), INTENT(OUT) :: c

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

    RETURN
  END SUBROUTINE cross

  !---------------------------------------------------------------------------
  ! Name: areaPolygon(poly)
  !
  ! Description: Compute the area of a polygon given its vertices
  !
  ! Inputs:  poly <REAL> : vertices of polygon
  !         
  ! Outputs: area <REAL> : area of polygon
  !---------------------------------------------------------------------------
  SUBROUTINE areaPolygon(poly,polyArea)

    IMPLICIT NONE
    
    ! Input
    REAL(kind=8), DIMENSION(:,:), INTENT(IN) :: poly
    
    ! Output
    REAL(kind=8), INTENT(OUT) :: polyArea

    ! Local
    INTEGER :: i,np
    REAL(kind=8), DIMENSION(3) :: total,prod,vi1,vi2,n

    ! Compute Area of a Polygon Given its Vertices
    np = SIZE(poly,2)
    IF(np < 3)THEN
       WRITE(*,*)"This is not a polygon, no area"
       STOP
    END IF
    total=0
    DO i = 1,np
       vi1 = poly(:,i)
       IF(i == np)THEN
          vi2 = poly(:,1)
       ELSE
          vi2 = poly(:,i+1)
       END IF
       CALL cross(vi1,vi2,prod)
       total(1) = total(1) + prod(1)
       total(2) = total(2) + prod(2)
       total(3) = total(3) + prod(3)
    END DO
    CALL unitNormal(poly(:,1),poly(:,2),poly(:,3),n)
    CALL dot(total,n,polyArea)
    polyArea = ABS(polyArea*0.5d0)

    RETURN
  END SUBROUTINE areaPolygon

  !---------------------------------------------------------------------------
  ! Name: setAuxiliaryVariables(num_aux,mbc,mx,my,mz,xlower,ylower,zlower,&
  !                             dxc,dyc,dzc,xyz0,xyzN,aux)
  !
  ! Description: Compute auxiliary variables for non-uniform/mapped grids
  !
  ! Inputs: 
  !         xc,yc,zc <REAL>    : computational coordinate values
  !         mx,my,mz <INTEGER> : grid sizes
  !         xyz0,xyzN <REAL>   : altitude at bottom and top of grid
  !         
  ! Outputs: 
  !         xp,yp,zp <REAL>    : physical coordinate values
  !---------------------------------------------------------------------------
  SUBROUTINE setAuxiliaryVariables(num_aux,mbc,mx,my,mz,xlower,ylower,zlower,&
       dxc,dyc,dzc,xyz0,xyzN,mapType,aux)

    IMPLICIT NONE
    
    ! Input
    INTEGER, INTENT(IN) :: num_aux,mbc,mx,my,mz
    REAL(kind=8), INTENT(IN) :: xlower,ylower,zlower,dxc,dyc,dzc
    REAL(kind=8), DIMENSION(3), INTENT(IN) :: xyz0,xyzN
    CHARACTER(20) :: mapType

    ! Output
    REAL(kind=8), DIMENSION(num_aux,mx+2*mbc,my+2*mbc,mz+2*mbc), &
         INTENT(OUT) :: aux
    
    ! Local
    REAL(kind=8), DIMENSION(3) :: norm
    REAL(kind=8), DIMENSION(3,4) :: poly
    REAL(kind=8), DIMENSION(8) :: xcCorner,ycCorner,zcCorner,&
         xpCorner,ypCorner,zpCorner
    INTEGER :: i,j,k
    REAL(kind=8) :: vHex,area
    REAL(kind=8), DIMENSION(1) :: xc,yc,zc,xp,yp,zp

    ! Compute Auxiliary Variables
    DO i = 1,mx+2*mbc
       DO j = 1,my+2*mbc
          DO k = 1,mz+2*mbc
             
             ! Compute Physical Coordinates of Corners of Hexahedron

             ! i-1/2, j-1/2, k-1/2 Corner
             xcCorner(1) = xlower + (i-1)*dxc
             ycCorner(1) = ylower + (j-1)*dyc
             zcCorner(1) = zlower + (k-1)*dzc
             xc(1) = xcCorner(1)
             yc(1) = ycCorner(1)
             zc(1) = zcCorner(1)
             CALL mapc2pWrapper(xc,yc,zc,mz,xyz0,xyzN,mapType,xp,yp,zp)
             xpCorner(1) = xp(1)
             ypCorner(1) = yp(1)
             zpCorner(1) = zp(1)
             
             ! i+1/2, j-1/2, k-1/2 Corner
             xcCorner(2) = xcCorner(1) + dxc
             ycCorner(2) = ycCorner(1)
             zcCorner(2) = zcCorner(1)
             xc(1) = xcCorner(2)
             yc(1) = ycCorner(2)
             zc(1) = zcCorner(2)
             CALL mapc2pWrapper(xc,yc,zc,mz,xyz0,xyzN,mapType,xp,yp,zp)
             xpCorner(2) = xp(1)
             ypCorner(2) = yp(1)
             zpCorner(2) = zp(1)
             
             ! i+1/2, j+1/2, k-1/2 Corner
             xcCorner(3) = xcCorner(1) + dxc
             ycCorner(3) = ycCorner(1) + dyc
             zcCorner(3) = zcCorner(1)
             xc(1) = xcCorner(3)
             yc(1) = ycCorner(3)
             zc(1) = zcCorner(3)
             CALL mapc2pWrapper(xc,yc,zc,mz,xyz0,xyzN,mapType,xp,yp,zp)
             xpCorner(3) = xp(1)
             ypCorner(3) = yp(1)
             zpCorner(3) = zp(1)
             
             ! i-1/2, j+1/2, k-1/2 Corner
             xcCorner(4) = xcCorner(1)
             ycCorner(4) = ycCorner(1) + dyc
             zcCorner(4) = zcCorner(1)
             xc(1) = xcCorner(4)
             yc(1) = ycCorner(4)
             zc(1) = zcCorner(4)
             CALL mapc2pWrapper(xc,yc,zc,mz,xyz0,xyzN,mapType,xp,yp,zp)
             xpCorner(4) = xp(1)
             ypCorner(4) = yp(1)
             zpCorner(4) = zp(1)

             ! i-1/2, j-1/2, k+1/2 Corner
             xcCorner(5) = xcCorner(1)
             ycCorner(5) = ycCorner(1)
             zcCorner(5) = zcCorner(1) + dzc
             xc(1) = xcCorner(5)
             yc(1) = ycCorner(5)
             zc(1) = zcCorner(5)
             CALL mapc2pWrapper(xc,yc,zc,mz,xyz0,xyzN,mapType,xp,yp,zp)
             xpCorner(5) = xp(1)
             ypCorner(5) = yp(1)
             zpCorner(5) = zp(1)

             ! i+1/2, j-1/2, k+1/2 Corner
             xcCorner(6) = xcCorner(1) + dxc
             ycCorner(6) = ycCorner(1)
             zcCorner(6) = zcCorner(1) + dzc
             xc(1) = xcCorner(6)
             yc(1) = ycCorner(6)
             zc(1) = zcCorner(6)
             CALL mapc2pWrapper(xc,yc,zc,mz,xyz0,xyzN,mapType,xp,yp,zp)
             xpCorner(6) = xp(1)
             ypCorner(6) = yp(1)
             zpCorner(6) = zp(1)
             
             ! i+1/2, j+1/2, k+1/2 Corner
             xcCorner(7) = xcCorner(1) + dxc
             ycCorner(7) = ycCorner(1) + dyc
             zcCorner(7) = zcCorner(1) + dzc
             xc(1) = xcCorner(7)
             yc(1) = ycCorner(7)
             zc(1) = zcCorner(7)
             CALL mapc2pWrapper(xc,yc,zc,mz,xyz0,xyzN,mapType,xp,yp,zp)
             xpCorner(7) = xp(1)
             ypCorner(7) = yp(1)
             zpCorner(7) = zp(1)
             
             ! i-1/2, j+1/2, k+1/2 Corner
             xcCorner(8) = xcCorner(1)
             ycCorner(8) = ycCorner(1) + dyc
             zcCorner(8) = zcCorner(1) + dzc
             xc(1) = xcCorner(8)
             yc(1) = ycCorner(8)
             zc(1) = zcCorner(8)
             CALL mapc2pWrapper(xc,yc,zc,mz,xyz0,xyzN,mapType,xp,yp,zp)
             xpCorner(8) = xp(1)
             ypCorner(8) = yp(1)
             zpCorner(8) = zp(1)

             ! Compute Aux Variables (Normals/Area Ratios/Volume Ratio)

             ! \xi normal
             poly(:,1) = (/xpCorner(1),ypCorner(1),zpCorner(1)/)
             poly(:,2) = (/xpCorner(4),ypCorner(4),zpCorner(4)/)
             poly(:,3) = (/xpCorner(8),ypCorner(8),zpCorner(8)/)
             poly(:,4) = (/xpCorner(5),ypCorner(5),zpCorner(5)/)
             CALL unitNormal(poly(:,1),poly(:,2),poly(:,3),norm)
             aux(1,i,j,k) = norm(1) ! nx(i-1/2,j,k)
             aux(2,i,j,k) = norm(2) ! ny(i-1/2,j,k)
             aux(3,i,j,k) = norm(3) ! nz(i-1/2,j,k)
             CALL areaPolygon(poly,area)
             aux(4,i,j,k) = area/(dyc*dzc) ! gamma(i-1/2,j,k)
             
             ! \eta normal
             poly(:,1) = (/xpCorner(1),ypCorner(1),zpCorner(1)/)
             poly(:,2) = (/xpCorner(5),ypCorner(5),zpCorner(5)/)
             poly(:,3) = (/xpCorner(6),ypCorner(6),zpCorner(6)/)
             poly(:,4) = (/xpCorner(2),ypCorner(2),zpCorner(2)/)
             CALL unitNormal(poly(:,1),poly(:,2),poly(:,3),norm)
             aux(5,i,j,k) = norm(1) ! nx(i,j-1/2,k)
             aux(6,i,j,k) = norm(2) ! ny(i,j-1/2,k)
             aux(7,i,j,k) = norm(3) ! nz(i,j-1/2,k)
             CALL areaPolygon(poly,area)
             aux(8,i,j,k) = area/(dxc*dzc) ! gamma(i,j-1/2,k)
           
             ! \phi normal
             poly(:,1) = (/xpCorner(1),ypCorner(1),zpCorner(1)/)
             poly(:,2) = (/xpCorner(2),ypCorner(2),zpCorner(2)/)
             poly(:,3) = (/xpCorner(3),ypCorner(3),zpCorner(3)/)
             poly(:,4) = (/xpCorner(4),ypCorner(4),zpCorner(4)/)
             CALL unitNormal(poly(:,1),poly(:,2),poly(:,3),norm)
             aux(9,i,j,k) =  norm(1) ! nx(i,j,k-1/2)
             aux(10,i,j,k) = norm(2) ! ny(i,j,k-1/2)
             aux(11,i,j,k) = norm(3) ! nz(i,j,k-1/2)
             CALL areaPolygon(poly,area)
             aux(12,i,j,k) = area/(dxc*dyc) ! gamma(i,j,k-1/2)             

             ! Volume of hexahedron
             vHex = volumeOfHex(xpCorner,ypCorner,zpCorner)
             
             ! Capacity Function (Ratio of Physical:Computational volume)
             aux(13,i,j,k) = vHex/(dxc*dyc*dzc) ! kappa(i,j,k)
             
             ! Grid Delta in Direction of Gravity
             SELECT CASE(mapType)
             CASE("Identity","ZeroToOne","Rotate")
                aux(14,i,j,k) = ABS(zpCorner(5)-zpCorner(1)) ! dz
                aux(15,i,j,k) = (zpCorner(5) + zpCorner(1))/2.d0 ! z
             CASE("Spherical")
                aux(14,i,j,k) = &
                     ABS(SQRT(xpCorner(5)**2+ypCorner(5)**2+zpCorner(5)**2) &
                     - SQRT(xpCorner(1)**2+ypCorner(1)**2+zpCorner(1)**2)) ! dr
                aux(15,i,j,k) = (SQRT(xpCorner(5)**2+ypCorner(5)**2+zpCorner(5)**2) &
                     + SQRT(xpCorner(1)**2+ypCorner(1)**2+zpCorner(1)**2))/2.d0 &
                     - rEarth
             END SELECT
             
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE setAuxiliaryVariables

END MODULE euler3d_mappedgrid
