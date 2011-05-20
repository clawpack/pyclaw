! ==========================================================
subroutine flux2(q,dq,q1d,dq1d,aux,dt,cfl,t,maux,meqn,mbc,maxnx,mx,my)
! ==========================================================

    ! Evaluate (delta t) *dq/dt
    ! On entry, q should give the initial data for this step
    ! dimension-by-dimension version
    ! (No genuinely multi-d terms)

    ! dq1d is used to return increments to q from flux1
    ! See the flux1 documentation for more information.

    use workspace
    use ClawParams
    implicit none

!f2py intent(in,out) dq  
!f2py intent(out) cfl  

    integer :: maux,meqn,mbc,maxnx,mx,my
    double precision, target, intent(in) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    double precision, intent(inout) :: dq(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    double precision, target, intent(in) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    double precision :: q1d(meqn,1-mbc:maxnx+mbc), dq1d(meqn,1-mbc:maxnx+mbc)
    double precision, intent(in) :: dt,t
    double precision, intent(out) :: cfl
    integer :: i,j,m
    double precision :: cfl1d
    double precision, pointer :: auxp(:,:),q1dp(:,:)
    
    cfl = 0.d0

    ! perform x-sweeps
    ! ==================

   
    do j = 0,my+1

        ! copy auxiliary data along a slice into 1d arrays:
        q1dp => q(:,:,j)
        if (maux .gt. 0)  then
            auxp => aux(:,:,j)
        endif


        ! compute modification dq1d along this slice:
        call flux1(q1dp,dq1d,auxp,dt,cfl1d,t,1,maux,meqn,mx,mbc,maxnx)
        cfl = dmax1(cfl,cfl1d)


        if (mcapa.eq.0) then
            ! no capa array.  Standard flux differencing:
            forall(i=1:mx, m=1:meqn)
                dq(m,i,j) = dq(m,i,j)+dq1d(m,i)
            end forall
        else
            ! with capa array.  Which is correct?
            forall(i=1:mx, m=1:meqn)
                dq(m,i,j) = dq(m,i,j)+dq1d(m,i)
            end forall
            ! dq(m,i,j) = dq(m,i,j)+dq1d(m,i)/aux(mcapa,i,j)
        endif
    enddo !end x sweeps


    ! perform y sweeps
    ! ==================

    do i = 0, mx+1
        ! copy auxiliary data along a slice into 1d arrays:
        q1dp => q(:,i,:)
        forall(j=1-mbc:my+mbc, m=1:meqn)
            q1d(m,j) = q(m,i,j)
        end forall

        if (maux .gt. 0)  then
            auxp => aux(:,i,:)
        endif

        call flux1(q1dp,dq1d,auxp,dt,cfl1d,t,2,maux,meqn,my,mbc,maxnx)
        cfl = dmax1(cfl,cfl1d)

        if (mcapa.eq.0) then
            ! no capa array.  Standard flux differencing:
            forall(j=1:my,m=1:meqn)
                dq(m,i,j) = dq(m,i,j)+dq1d(m,j)
            end forall
        else
            ! with capa array.  Which is correct?
            ! dq(m,i,j) = dq(m,i,j)+dq1d(m,j)/aux(mcapa,i,j)
            dq(:,i,:) = dq(:,i,:)+dq1d
        endif
    enddo !end y sweeps

end subroutine flux2
