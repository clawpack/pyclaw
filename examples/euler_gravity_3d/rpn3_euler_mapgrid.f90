!-----------------------------------------------------------------------------------
! Name: rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)
!
! Description: Roe f-wave solver for the 3D Euler equations in general geometries
!
! Waves: 3
! Equations: 5
!
! Conserved Quantities:
!       1 density
!       2 x_momentum
!       3 y_momentum
!       4 z_momentum
!       5 energy
!
! Auxilliary Quantities:
!       1  nx xface (i-1/2,j,k)
!       2  ny xface (i-1/2,j,k)
!       3  nz xface (i-1/2,j,k)
!       4  gamma_x (area/(dyc*dzc))
!       5  nx yface (i,j-1/2,k)
!       6  ny yface (i,j-1/2,k)
!       7  nz yface (i,j-1/2,k)
!       8  gamma_x (area/(dxc*dzc))
!       9  nx zface (i,j,k-1/2)
!       10 ny zface (i,j,k-1/2)
!       11 nz zface (i,j,k-1/2)
!       12 gamma_z (area/(dxc*dyc))
!       13 kappa(i,j,k) volumeOfHex/(dxc*dyc*dzc)
!       14 dz (Cartesian), dr (Spherical)
!       15 z (Cartesian), r (Spherical)
!
! Inputs: 
!         ixyz <INTEGER>   : direction to take slice x-direction if ixyz=1
!                                                    y-direction if ixyz=2.
!                                                    z-direction if ixyz=3.
!         maxm <INTEGER>   : max number of grid cells (less the ghost cells)
!         meqn <INTEGER>   : number of equations in the system (=5)
!         mwaves <INTEGER> : nuber of waves in the system (=3)
!         maux <INTEGER>   : number of auxilary variables (=15)
!         mbc <INTEGER>    : number of ghost cells on either end
!         mx <INTEGER>     : number of elements
!         ql <REAL>        : state vector at left edge of each cell
!                            Note the i'th Riemann problem has left state qr(:,i-1)
!         qr <REAL>        : state vector at right edge of each cell
!                            Note the i'th Riemann problem has right state ql(:,i)
!         auxl <REAL>      : state of auxilary variable on left egde of cell 
!         auxr <REAL>      : state of auxilary variable on right egde of cell 
!         
! Outputs: 
!          fwave <REAL>    : f-wave vectors (eigenvectors of Roe matrix)
!          s <REAL>        : eigenvalues (eigenvalues of Roe matrix)
!          amdq <REAL>     : left-going fluctuations
!          apdq <REAL>     : right-going fluctuations
!
! Adapted from rpn3_euler.f90 in $CLAWHOME/riemann/src
!-----------------------------------------------------------------------------------
SUBROUTINE rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)

  IMPLICIT NONE

  ! Input
  INTEGER, INTENT(IN) :: ixyz,maxm,meqn,mwaves,maux,mbc,mx
  REAL(kind=8), DIMENSION(meqn,1-mbc:maxm+mbc), INTENT(IN) :: ql,qr
  REAL(kind=8), DIMENSION(maux,1-mbc:maxm+mbc), INTENT(IN) :: auxl,auxr

  ! Output
  REAL(kind=8), DIMENSION(meqn,mwaves,1-mbc:maxm+mbc), INTENT(INOUT) :: fwave
  REAL(kind=8), DIMENSION(mwaves,1-mbc:maxm+mbc), INTENT(INOUT) :: s
  REAL(kind=8), DIMENSION(meqn,1-mbc:maxm+mbc), INTENT(INOUT) :: amdq,apdq
  
  ! Local Storage
  REAL(kind=8), DIMENSION(meqn) :: delta,beta
  REAL(kind=8), DIMENSION(-1:maxm+mbc) :: u,v,w,enth,ek,a,a2g1,psi
  REAL(kind=8) :: gamma1,asqrd,pl,pr,vnl,vnr,rhsq2,rhsqrtl,rhsqrtr,nx,ny,nz,vn,&
       areaRatio,dr,finterp,qinterp,quinterp,qvinterp,qwinterp
  INTEGER :: i,j,mu=2,mv=3,mw=4,mws,nu,nv,nw,iratio=4,gFlux=0
  CHARACTER(20) :: nonsingularNormal

  ! Common block storage for ideal gas constant and gravity term
  REAL(kind=8) :: gamma ! Ideal gas constant
  REAL(kind=8) :: g_r  ! Gravitational Constant Magnitude
  LOGICAL :: gravity,gravityflux ! Turn Gravity Source Term On/Off
  COMMON /cparam/ gamma,g_r,gravity,gravityflux

  ! Set (gamma-1)
  gamma1 = gamma - 1.d0

  ! Gravity in Energy Flux On/Off
  gFlux = 0
  IF(gravityflux)gFlux=1

  ! Make sure the number of waves is 3 for this solver
  IF (mwaves /= 3) THEN
     WRITE(6,*) '*** Must set mwaves=3 for this Riemann solver'
     STOP
  ENDIF
    
  ! Set mu to point to  the component of the system that corresponds to momentum in 
  !  the direction of this slice, mv and mw to the orthogonal momentums.  Set the
  !  nu,nv,and nw indices for the normal vectors in a similar way.  The iratio is
  !  the aux array index for the gamma area scaling for mapped grids.
  IF(ixyz == 1)THEN
     mu = 2
     mv = 3
     mw = 4
     iratio = 4
  ELSE IF(ixyz == 2)THEN
     mu = 3
     mv = 4
     mw = 2
     iratio = 8
  ELSE IF (ixyz == 3)THEN
     mu = 4
     mv = 2
     mw = 3
     iratio = 12
  ELSE
     WRITE(*,*)"Something is wrong with ixyz..."
     STOP
  ENDIF
  nu = (ixyz-1)*4 + mu - 1
  nv = (ixyz-1)*4 + mv - 1
  nw = (ixyz-1)*4 + mw - 1

  ! Note that notation for u,v, and w reflects assumption that the
  !   Riemann problems are in the x-direction with u in the normal
  !   direction and v and w in the orthogonal directions, but with the
  !   above definitions of mu, mv, and mw the routine also works with
  !   ixyz=2 and ixyz = 3
  !   and returns, for example, f0 as the Godunov flux g0 for the
  !   Riemann problems u_t + g(u)_y = 0 in the y-direction.
  
  ! Initialize fwaves and speeds
  fwave = 0.d0
  s = 0.d0

  ! Loop over grid cell interfaces
  DO i = 2-mbc, mx+mbc
     IF (qr(1,i-1) <= 0.d0 .OR. ql(1,i) <= 0.d0) THEN
        WRITE(*,*) i, mu, mv, mw
        WRITE(*,'(5e12.4)') (qr(j,i-1),j=1,5)
        WRITE(*,'(5e12.4)') (ql(j,i),j=1,5)
        IF (ixyz == 1) WRITE(6,*) '*** (rpn3) rho <= 0 in x-sweep at ',i
        IF (ixyz == 2) WRITE(6,*) '*** (rpn3) rho <= 0 in y-sweep at ',i
        IF (ixyz == 3) WRITE(6,*) '*** (rpn3) rho <= 0 in z-sweep at ',i
        WRITE(6,*) 'stopped with rho <= 0...'
        STOP
     ENDIF
     rhsqrtl = SQRT(qr(1,i-1))
     rhsqrtr = SQRT(ql(1,i))
     pl = gamma1*(qr(5,i-1) - 0.5d0*(qr(mu,i-1)**2 + &
          qr(mv,i-1)**2 + qr(mw,i-1)**2)/qr(1,i-1) &
          - qr(1,i-1)*(g_r)*auxr(15,i-1)*gFlux )
     pr = gamma1*(ql(5,i) - 0.5d0*(ql(mu,i)**2 + &
          ql(mv,i)**2 + ql(mw,i)**2)/ql(1,i) &
          - ql(1,i)*(g_r)*auxl(15,i)*gFlux )
     rhsq2 = rhsqrtl + rhsqrtr
     u(i) = (qr(mu,i-1)/rhsqrtl + ql(mu,i)/rhsqrtr) / rhsq2
     v(i) = (qr(mv,i-1)/rhsqrtl + ql(mv,i)/rhsqrtr) / rhsq2
     w(i) = (qr(mw,i-1)/rhsqrtl + ql(mw,i)/rhsqrtr) / rhsq2
     enth(i) = ((qr(5,i-1)+pl)/rhsqrtl + (ql(5,i)+pr)/rhsqrtr) / rhsq2
     psi(i) = (qr(1,i-1)*(g_r)*auxr(15,i-1)/rhsqrtl &
          + ql(1,i)*(g_r)*auxl(15,i)/rhsqrtr)/rhsq2*gFlux
     ek(i) = 0.5d0*(u(i)**2 + v(i)**2 + w(i)**2)
     asqrd = gamma1*(enth(i) - ek(i))
     IF (i>=0 .AND. i<=mx .AND. asqrd <= 0.d0) THEN
        IF (ixyz == 1) WRITE(6,*) '*** (rpn3) a**2 <= 0 in x-sweep at ',i
        IF (ixyz == 2) WRITE(6,*) '*** (rpn3) a**2 <= 0 in y-sweep at ',i
        IF (ixyz == 3) WRITE(6,*) '*** (rpn3) a**2 <= 0 in z-sweep at ',i
        WRITE(6,*) 'stopped with a**2 < 0...'
        STOP
     ENDIF
     a(i) = SQRT(asqrd)
     a2g1(i)  = asqrd/gamma1
  END DO

  ! Now split the jump in q1d at each interface into f-waves
  !   find alpha(1) thru alpha(5), the coefficients of the 5 eigenvectors:
  DO i = 2-mbc, mx+mbc
     
     ! Normal Vectors, Normal Velocity, and Area Ratio
     nx = auxl(nu,i)
     ny = auxl(nv,i)
     nz = auxl(nw,i)
     vn = u(i)*nx + v(i)*ny + w(i)*nz
     areaRatio = auxl(iratio,i)

     ! Compute the jumps (delta F) in the flux functions
     pl = gamma1*(ql(5,i) &
          - 0.5d0*(ql(mu,i)**2+ql(mv,i)**2+ql(mw,i)**2)/ql(1,i) &
          - ql(1,i)*(g_r)*auxl(15,i)*gFlux)
     pr = gamma1*(qr(5,i-1) &
          - 0.5d0*(qr(mu,i-1)**2+qr(mv,i-1)**2+qr(mw,i-1)**2)/qr(1,i-1) &
          - qr(1,i-1)*(g_r)*auxr(15,i-1)*gFlux)
     vnl = (ql(mu,i)*nx + ql(mv,i)*ny + ql(mw,i)*nz)/ql(1,i)
     vnr = (qr(mu,i-1)*nx + qr(mv,i-1)*ny + qr(mw,i-1)*nz)/qr(1,i-1)
     delta(1) = ql(1,i)*vnl - qr(1,i-1)*vnr
     delta(2) = (vnl*ql(mu,i) + pl*nx) - (vnr*qr(mu,i-1) + pr*nx)
     delta(3) = (vnl*ql(mv,i) + pl*ny) - (vnr*qr(mv,i-1) + pr*ny)
     delta(4) = (vnl*ql(mw,i) + pl*nz) - (vnr*qr(mw,i-1) + pr*nz)
     delta(5) = ((ql(5,i) + pl)*vnl) - ((qr(5,i-1) + pr)*vnr)
     
     ! Modify delta(2/3/4/5) to account for gravitational term
     IF(gravity .AND. ixyz==3)THEN
        dr = (auxl(14,i) + auxl(14,i-1))*0.5d0
        finterp = auxl(14,i-1)/(auxl(14,i)+auxl(14,i-1))
        qinterp = ql(1,i)*finterp + qr(1,i-1)*(1.d0-finterp)
        quinterp = ql(mu,i)*finterp + qr(mu,i-1)*(1.d0-finterp)
        qvinterp = ql(mv,i)*finterp + qr(mv,i-1)*(1.d0-finterp)
        qwinterp = ql(mw,i)*finterp + qr(mw,i-1)*(1.d0-finterp)
        delta(2) = delta(2) + dr*g_r*nx*qinterp
        delta(3) = delta(3) + dr*g_r*ny*qinterp
        delta(4) = delta(4) + dr*g_r*nz*qinterp
        IF(.NOT.gravityflux)&
             delta(5) = delta(5) + dr*g_r*(quinterp*nx + qvinterp*ny + qwinterp*nz)
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

     beta(4) = ((enth(i) + psi(i) - 2.d0*ek(i))*delta(1) &
          + u(i)*delta(2) + v(i)*delta(3) + w(i)*delta(4) - delta(5))/a2g1(i)
     beta(5) = ( (a(i)-vn)*delta(1) - a(i)*beta(4) &
          + nx*delta(2) + ny*delta(3) + nz*delta(4) )/(2.d0*a(i))
     beta(1) = delta(1) - beta(4) - beta(5)

     SELECT CASE(TRIM(nonsingularNormal))
     CASE("x")
        beta(2) = ( (v(i)-vn*ny)*delta(1) + nx*ny*delta(2) &
             + (ny**2-1.d0)*delta(3) + ny*nz*delta(4) )/nx
        beta(3) = (-(w(i)-vn*nz)*delta(1) - nx*nz*delta(2) &
             - ny*nz*delta(3) - (nz**2-1.d0)*delta(4) )/nx
     CASE("y")
        beta(2) = (-(u(i)-vn*nx)*delta(1) - (nx**2-1.d0)*delta(2) &
             - ny*nx*delta(3) - nx*nz*delta(4) )/ny
        beta(3) = ( (w(i)-vn*nz)*delta(1) + nx*nz*delta(2) &
             + ny*nz*delta(3) + (nz**2-1.d0)*delta(4) )/ny
     CASE("z")
        beta(2) = ( (u(i)-vn*nx)*delta(1) + (nx**2-1.d0)*delta(2) &
             + nx*ny*delta(3) + nz*nx*delta(4) )/nz
        beta(3) = (-(v(i)-vn*ny)*delta(1) - nx*ny*delta(2) &
             - (ny**2-1.d0)*delta(3) - nz*ny*delta(4) )/nz
     END SELECT
     
     ! Compute the f-waves.
     !   Note that the 2-wave, 3-wave and 4-wave travel at the same speed
     !   and are lumped together in fwave(.,2,.).  The 5-wave is then stored
     !   in fwave(.,3,.).
     fwave(1,1,i)  = beta(1)
     fwave(mu,1,i) = beta(1)*(u(i)-a(i)*nx)
     fwave(mv,1,i) = beta(1)*(v(i)-a(i)*ny)
     fwave(mw,1,i) = beta(1)*(w(i)-a(i)*nz)
     fwave(5,1,i)  = beta(1)*(enth(i) + psi(i) - a(i)*vn)
     s(1,i) = (vn-a(i))*areaRatio
     
     fwave(1,2,i)  = beta(4)
     SELECT CASE(TRIM(nonsingularNormal))
     CASE("x")
        fwave(mu,2,i) = beta(4)*u(i) + beta(2)*ny - beta(3)*nz
        fwave(mv,2,i) = beta(4)*v(i) - beta(2)*nx
        fwave(mw,2,i) = beta(4)*w(i) + beta(3)*nx
        fwave(5,2,i)  = beta(4)*(ek(i) + psi(i))  &
             + beta(2)*(u(i)*ny-v(i)*nx) + beta(3)*(w(i)*nx-u(i)*nz)
     CASE("y")
        fwave(mu,2,i) = beta(4)*u(i) + beta(2)*ny
        fwave(mv,2,i) = beta(4)*v(i) - beta(2)*nx + beta(3)*nz
        fwave(mw,2,i) = beta(4)*w(i) - beta(3)*ny
        fwave(5,2,i)  = beta(4)*(ek(i) + psi(i))  &
             + beta(2)*(u(i)*ny-v(i)*nx) + beta(3)*(v(i)*nz-w(i)*ny)
     CASE("z")
        fwave(mu,2,i) = beta(4)*u(i) - beta(2)*nz
        fwave(mv,2,i) = beta(4)*v(i) + beta(3)*nz
        fwave(mw,2,i) = beta(4)*w(i) + beta(2)*nx - beta(3)*ny
        fwave(5,2,i)  = beta(4)*(ek(i) + psi(i))  &
             + beta(2)*(w(i)*nx-u(i)*nz) + beta(3)*(v(i)*nz-w(i)*ny)
     END SELECT
     s(2,i) = vn*areaRatio
     
     fwave(1,3,i)  = beta(5)
     fwave(mu,3,i) = beta(5)*(u(i)+a(i)*nx)
     fwave(mv,3,i) = beta(5)*(v(i)+a(i)*ny)
     fwave(mw,3,i) = beta(5)*(w(i)+a(i)*nz)
     fwave(5,3,i)  = beta(5)*(enth(i) + psi(i) + a(i)*vn)
     s(3,i) = (vn+a(i))*areaRatio

     ! Compute fluctuations amdq and apdq
     !  amdq = SUM s*fwave   over left-going waves
     !  apdq = SUM s*fwave   over right-going waves
     amdq(:,i) = 0.d0
     apdq(:,i) = 0.d0
     DO mws=1,mwaves
        fwave(:,mws,i) = fwave(:,mws,i)*areaRatio
        IF (s(mws,i) < 0.d0) THEN
           amdq(:,i) = amdq(:,i) + fwave(:,mws,i)
        ELSE
           apdq(:,i) = apdq(:,i) + fwave(:,mws,i)
        ENDIF
     END DO

  END DO
END SUBROUTINE rpn3
