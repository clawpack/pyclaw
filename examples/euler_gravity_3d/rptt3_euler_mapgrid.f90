!-----------------------------------------------------------------------------------
! Name: rptt3(ixyz,icoor,imp,impt,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,
!            aux1,aux2,aux3,asdq,cmbsasdq,cpbsasdq)
!
! Description: Roe f-wave solver in the transverse direction for the 
!              3D Euler equations in general geometries
!
! Inputs: 
!         ixyz <INTEGER>   : direction to take slice x-direction if ixyz=1
!                                                    y-direction if ixyz=2.
!                                                    z-direction if ixyz=3.
!         icoor <INTEGER>  : direction in which the transverse solve should be
!                            performed.  icoor=2 (split in y-like direction)
!                                        icoor=3 (split in z-like direction)
!                            ixyz=1, icoor=3 means bsasdq=B*A*\Dq, split in z into:
!                                            cmbsasdq = C^-B*A*\Dq,
!                                            cpbsasdq = C^+B*A*\Dq.
!                            ixyz=2, icoor=3 means bsasdq=C*B*\Dq, split in x into: 
!                                            cmbsasdq = A^-C*B*\Dq,
!                                            cpbsasdq = A^+C*B*\Dq.
!         imp <INTEGER>    : index for aux arrays
!         impt <INTEGER>   : index for aux arrays
!         maxm <INTEGER>   : max number of grid cells (less the ghost cells)
!         meqn <INTEGER>   : number of equations in the system
!         mwaves <INTEGER> : nuber of waves in the system
!         maux <INTEGER>   : number of auxilary equations
!         mbc <INTEGER>    : number of ghost cells on either end
!         mx <INTEGER>     : number of elements
!         ql <REAL>        : state vector at left edge of each cell
!                            Note the i'th Riemann problem has left state qr(:,i-1)
!         qr <REAL>        : state vector at right edge of each cell
!                            Note the i'th Riemann problem has right state ql(:,i)
!         aux1 <REAL>      : 
!         aux2 <REAL>      :
!         aux3 <REAL>      :
!         asdq <REAL>      : array of flux differences (A*\Dq).  asdq(:,i) is the
!                            flux difference propagating away from the interface
!                            between cells i-1 and i.  Note: asdq represents
!                            B*\Dq if ixyz=2 or C*\Dq if ixyz=3.
!         
! Outputs: 
!          cmbsasdq <REAL> : left-going flux differences
!          cpbsasdq <REAL> : right-going flux differences
!
! Adapted from rptt3_euler.f90 in $CLAWHOME/riemann/src
!-----------------------------------------------------------------------------------
SUBROUTINE rptt3(ixyz,icoor,imp,impt,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,&
     aux1,aux2,aux3,bsasdq,cmbsasdq,cpbsasdq)

  IMPLICIT NONE

  ! Input
  INTEGER, INTENT(IN) :: ixyz,icoor,imp,impt,maxm,meqn,mwaves,maux,mbc,mx
  REAL(kind=8), DIMENSION(meqn,1-mbc:maxm+mbc), INTENT(IN) :: ql,qr,bsasdq
  REAL(kind=8), DIMENSION(maux,1-mbc:maxm+mbc,3), INTENT(IN) :: aux1,aux2,aux3

  ! Output
  REAL(kind=8), DIMENSION(meqn,1-mbc:maxm+mbc), INTENT(INOUT) :: cmbsasdq,cpbsasdq
  
  ! Local Storage
  REAL(kind=8), DIMENSION(meqn) :: beta
  REAL(kind=8), DIMENSION(-1:maxm+mbc) :: u2v2w2,u,v,w,enth,a,g1a2,euv
  REAL(kind=8) :: gamma1,asqrd,pl,pr,rhsq2,rhsqrtl,rhsqrtr,nx,ny,nz,vn,&
       areaRatioMinus,areaRatioPlus
  REAL(kind=8), DIMENSION(meqn,mwaves) :: waveb
  REAL(kind=8), DIMENSION(mwaves) :: sb
  INTEGER :: i,j,mu=2,mv=3,mw=4,mws,nu=1,nv=2,nw=3,iratio=4
  CHARACTER(20) :: nonsingularNormal

  ! Common block storage for ideal gas constant
  REAL(kind=8) :: gamma
  COMMON /cparam/ gamma

  ! Set (gamma-1)
  gamma1 = gamma - 1.d0

  ! Make sure the number of waves is 3 for the 3D Euler equations
  IF (mwaves /= 3) THEN
     WRITE(6,*) '*** Should have mwaves=3 for this Riemann solver'
     STOP
  ENDIF
    
  ! Set mu to point to  the component of the system that corresponds to momentum in 
  !  the direction of this slice, mv and mw to the orthogonal momentums.
  !
  !  ixyz indicates the direction of the original Riemann solve,
  !  called the x-like direction in the table below:
  ! 
  !          x-like direction   y-like direction   z-like direction
  !  ixyz=1:        x                  y                  z         
  !  ixyz=2:        y                  z                  x         
  !  ixyz=3:        z                  x                  y      
  IF(ixyz == 1)THEN ! x-like direction
     mu = 2
     mv = 3
     mw = 4
     ! Choose the coordinate direction to solve the Riemann problem
     !   Permute nu,nv,nw for the original ixyz incident direction (x-like)
     IF(icoor == 2)THEN ! y-like direction
        nu = 5
        nv = 6
        nw = 7
        iratio = 8
     ELSEIF(icoor == 3)THEN ! z-like direction
        nu = 9
        nv = 10
        nw = 11
        iratio = 12
     END IF
  ELSEIF(ixyz == 2)THEN ! y-like direction
     mu = 3
     mv = 4
     mw = 2
     ! Choose the coordinate direction to solve the Riemann problem
     !   Permute nu,nv,nw for the original ixyz incident direction (y-like)
     IF(icoor == 2)THEN ! z-like direction
        nu = 10
        nv = 11
        nw = 9
        iratio = 12
     ELSEIF(icoor == 3)THEN ! x-like direction
        nu = 2
        nv = 3
        nw = 1
        iratio = 4
     END IF
  ELSEIF(ixyz == 3)THEN ! z-like direction
     mu = 4
     mv = 2
     mw = 3
     ! Choose the coordinate direction to solve the Riemann problem
     !   Permute nu,nv,nw for the original ixyz incident direction (z-like)
     IF(icoor == 2)THEN ! x-like direction
        nu = 3
        nv = 1
        nw = 2
        iratio = 4
     ELSEIF(icoor == 3)THEN ! y-like direction
        nu = 7
        nv = 5
        nw = 6
        iratio = 8
     END IF
  ENDIF

  ! Note that notation for u,v, and w reflects assumption that the
  !   Riemann problems are in the x-direction with u in the normal
  !   direction and v and w in the orthogonal directions, but with the
  !   above definitions of mu, mv, and mw the routine also works with
  !   ixyz=2 and ixyz = 3
  !   and returns, for example, f0 as the Godunov flux g0 for the
  !   Riemann problems u_t + g(u)_y = 0 in the y-direction.
  !
  ! Compute the Roe-averaged variables needed in the Roe solver.

  ! Loop over grid cell interfaces
  DO i = 2-mbc, mx+mbc
     IF (qr(1,i-1) <= 0.d0 .OR. ql(1,i) <= 0.d0) THEN
        WRITE(*,*) i, mu, mv, mw
        WRITE(*,'(5e12.4)') (qr(j,i-1),j=1,5)
        WRITE(*,'(5e12.4)') (ql(j,i),j=1,5)
        IF (ixyz == 1) WRITE(6,*) '*** (rptt3) rho <= 0 in x-sweep at ',i
        IF (ixyz == 2) WRITE(6,*) '*** (rptt3) rho <= 0 in y-sweep at ',i
        IF (ixyz == 3) WRITE(6,*) '*** (rptt3) rho <= 0 in z-sweep at ',i
        WRITE(6,*) 'stopped with rho <= 0...'
        STOP
     ENDIF
     rhsqrtl = SQRT(qr(1,i-1))
     rhsqrtr = SQRT(ql(1,i))
     pl = gamma1*(qr(5,i-1) - 0.5d0*(qr(mu,i-1)**2 + &
          qr(mv,i-1)**2 + qr(mw,i-1)**2)/qr(1,i-1))
     pr = gamma1*(ql(5,i) - 0.5d0*(ql(mu,i)**2 + &
          ql(mv,i)**2 + ql(mw,i)**2)/ql(1,i))
     rhsq2 = rhsqrtl + rhsqrtr
     u(i) = (qr(mu,i-1)/rhsqrtl + ql(mu,i)/rhsqrtr) / rhsq2
     v(i) = (qr(mv,i-1)/rhsqrtl + ql(mv,i)/rhsqrtr) / rhsq2
     w(i) = (qr(mw,i-1)/rhsqrtl + ql(mw,i)/rhsqrtr) / rhsq2
     enth(i) = (((qr(5,i-1)+pl)/rhsqrtl &
          + (ql(5,i)+pr)/rhsqrtr)) / rhsq2
     u2v2w2(i) = u(i)**2 + v(i)**2 + w(i)**2
     asqrd = gamma1*(enth(i) - .5d0*u2v2w2(i))
     
     IF (i>=0 .AND. i<=mx .AND. asqrd <= 0.d0) THEN
        IF (ixyz == 1) WRITE(6,*) '*** (rptt3) a**2 <= 0 in x-sweep at ',i
        IF (ixyz == 2) WRITE(6,*) '*** (rptt3) a**2 <= 0 in y-sweep at ',i
        IF (ixyz == 3) WRITE(6,*) '*** (rptt3) a**2 <= 0 in z-sweep at ',i
        WRITE(6,*) 'stopped with a**2 < 0...'
        STOP
     ENDIF
     a(i) = SQRT(asqrd)
     g1a2(i) = gamma1 / asqrd
     euv(i) = enth(i) - u2v2w2(i)
  END DO

  ! Solve the Riemann Problem
  DO i = 2-mbc, mx+mbc

     nx = aux2(nu,i,2)
     ny = aux2(nv,i,2)
     nz = aux2(nw,i,2)
     vn = u(i)*nx + v(i)*ny + w(i)*nz
     IF(icoor == 2)THEN
        IF(impt==1)THEN
           areaRatioMinus = aux2(iratio,i-2+imp,1)
           areaRatioPlus  = aux3(iratio,i-2+imp,1)
        ELSEIF(impt==2)THEN
           areaRatioMinus = aux2(iratio,i-2+imp,3)
           areaRatioPlus  = aux3(iratio,i-2+imp,3)
        END IF
     ELSEIF(icoor == 3)THEN
        IF(impt == 1)THEN
           areaRatioMinus = aux1(iratio,i-2+imp,2)
           areaRatioPlus  = aux1(iratio,i-2+imp,3)
        ELSEIF(impt == 2)THEN
           areaRatioMinus = aux3(iratio,i-2+imp,2)
           areaRatioPlus  = aux3(iratio,i-2+imp,3)
        END IF
     END IF
     ! Pick the largest normal to avoid using a singular direction
     IF(ABS(nx) >= ABS(ny) .AND. ABS(nx) >= ABS(nz))THEN
        nonsingularNormal = "x"
     ELSEIF(ABS(ny) >= ABS(nx) .AND. ABS(ny) >= ABS(nz))THEN
        nonsingularNormal = "y"
     ELSEIF(ABS(nz) >= ABS(nx) .AND. ABS(nz) >= ABS(ny))THEN
        nonsingularNormal = "z"
     ELSE
        WRITE(*,*)"Invalid Normal Vector..."
        WRITE(*,*)"nx: ",nx,"ny: ",ny,"nz: ",nz
        STOP
     END IF
     
     beta(4) = g1a2(i)*(euv(i)*bsasdq(1,i) &
          + u(i)*bsasdq(mu,i) + v(i)*bsasdq(mv,i) + w(i)*bsasdq(mw,i) - bsasdq(5,i))
     beta(5) = ( (a(i)-vn)*bsasdq(1,i) - a(i)*beta(4) &
          + nx*bsasdq(mu,i) + ny*bsasdq(mv,i) + nz*bsasdq(mw,i) )/(2.d0*a(i))
     beta(1) = bsasdq(1,i) - beta(4) - beta(5)
     
     SELECT CASE(TRIM(nonsingularNormal))
     CASE("x")
        beta(2) = ( (v(i)-vn*ny)*bsasdq(1,i) + nx*ny*bsasdq(mu,i) &
             + (ny**2-1)*bsasdq(mv,i) + ny*nz*bsasdq(mw,i) )/nx
        beta(3) = (-(w(i)-vn*nz)*bsasdq(1,i) - nx*nz*bsasdq(mu,i) &
             - ny*nz*bsasdq(mv,i) - (nz**2-1)*bsasdq(mw,i) )/nx
     CASE("y")
        beta(2) = (-(u(i)-vn*nx)*bsasdq(1,i) - (nx**2-1)*bsasdq(mu,i) &
             - ny*nx*bsasdq(mv,i) - nx*nz*bsasdq(mw,i) )/ny
        beta(3) = ( (w(i)-vn*nz)*bsasdq(1,i) + nx*nz*bsasdq(mu,i) &
             + ny*nz*bsasdq(mv,i) + (nz**2-1)*bsasdq(mw,i) )/ny
     CASE("z")
        beta(2) = ( (u(i)-vn*nx)*bsasdq(1,i) + (nx**2-1)*bsasdq(mu,i) &
             + nx*ny*bsasdq(mv,i) + nz*nx*bsasdq(mw,i) )/nz
        beta(3) = (-(v(i)-vn*ny)*bsasdq(1,i) - nx*ny*bsasdq(mu,i) &
             - (ny**2-1)*bsasdq(mv,i) - nz*ny*bsasdq(mw,i) )/nz
     END SELECT
     
     waveb(1,1)  = beta(1)
     waveb(mu,1) = beta(1)*(u(i)-a(i)*nx)
     waveb(mv,1) = beta(1)*(v(i)-a(i)*ny)
     waveb(mw,1) = beta(1)*(w(i)-a(i)*nz)
     waveb(5,1)  = beta(1)*(enth(i) - vn*a(i))
     sb(1) = vn - a(i)
     
     waveb(1,2)  = beta(4)
     SELECT CASE(TRIM(nonsingularNormal))
     CASE("x")
        waveb(mu,2) = beta(4)*u(i) + beta(2)*ny - beta(3)*nz
        waveb(mv,2) = beta(4)*v(i) - beta(2)*nx
        waveb(mw,2) = beta(4)*w(i) + beta(3)*nx
        waveb(5,2)  = beta(4)*0.5d0*u2v2w2(i)  &
             + beta(2)*(u(i)*ny-v(i)*nx) + beta(3)*(w(i)*nx-u(i)*nz)
     CASE("y")
        waveb(mu,2) = beta(4)*u(i) + beta(2)*ny
        waveb(mv,2) = beta(4)*v(i) - beta(2)*nx + beta(3)*nz
        waveb(mw,2) = beta(4)*w(i) - beta(3)*ny
        waveb(5,2)  = beta(4)*0.5d0*u2v2w2(i)  &
             + beta(2)*(u(i)*ny-v(i)*nx) + beta(3)*(v(i)*nz-w(i)*ny)
     CASE("z")
        waveb(mu,2) = beta(4)*u(i) - beta(2)*nz
        waveb(mv,2) = beta(4)*v(i) + beta(3)*nz
        waveb(mw,2) = beta(4)*w(i) + beta(2)*nx - beta(3)*ny
        waveb(5,2)  = beta(4)*0.5d0*u2v2w2(i)  &
             + beta(2)*(w(i)*nx-u(i)*nz) + beta(3)*(v(i)*nz-w(i)*ny)
     END SELECT
     sb(2) = vn
     
     waveb(1,3)  = beta(5)
     waveb(mu,3) = beta(5)*(u(i) + a(i)*nx)
     waveb(mv,3) = beta(5)*(v(i) + a(i)*ny)
     waveb(mw,3) = beta(5)*(w(i) + a(i)*nz)
     waveb(5,3)  = beta(5)*(enth(i) + vn*a(i))
     sb(3) = vn + a(i)
     
     cmbsasdq(:,i) = 0.d0
     cpbsasdq(:,i) = 0.d0
     DO mws = 1,mwaves
        IF (sb(mws) < 0.d0)THEN
           waveb(:,mws) = waveb(:,mws)*areaRatioMinus
           cmbsasdq(:,i) = cmbsasdq(:,i) + sb(mws)*waveb(:,mws)
        ELSE
           waveb(:,mws) = waveb(:,mws)*areaRatioPlus
           cpbsasdq(:,i) = cpbsasdq(:,i) + sb(mws)*waveb(:,mws)
        END IF
     END DO
  END DO

END SUBROUTINE rptt3
