! ==========================================================
subroutine flux3(q,dq,q1d,dq1d,aux,dt,cfl,t,num_aux,num_eqn,num_ghost,maxnx,mx,my,mz,rpn3)
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

    integer :: num_aux,num_eqn,num_ghost,maxnx,mx,my,mz
    double precision, target, intent(in) :: q(num_eqn, 1-num_ghost:mx+num_ghost, 1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
    double precision, intent(inout) :: dq(num_eqn, 1-num_ghost:mx+num_ghost, 1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
    double precision, target, intent(in) :: aux(num_aux,1-num_ghost:mx+num_ghost,1-num_ghost:my+num_ghost,1-num_ghost:mz+num_ghost)
    double precision :: q1d(num_eqn,1-num_ghost:maxnx+num_ghost), dq1d(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in) :: dt,t
    double precision, intent(out) :: cfl
    integer :: i,j,k,m
    double precision :: cfl1d
    double precision, pointer :: auxp(:,:),q1dp(:,:)
    external :: rpn3
    
!f2py intent(in,out) dq  
!f2py intent(out) cfl  
!f2py optional dq, q1d, dq1d

! Dummy interface just so f2py doesn't complain:
!f2py real(DP) x
!f2py x=rpn3(x)

    cfl = 0.d0

    ! perform x-sweeps
    ! ==================

   
    do j = 0,my+1
        do k = 0,mz+1
            ! copy data along a slice into 1d arrays:
            q1dp => q(:,:,j,k)
            if (num_aux.gt.0)  then
                auxp => aux(:,:,j,k)
            endif

            ! compute modification dq1d along this slice:
            call flux1(q1dp,dq1d,auxp,dt,cfl1d,t,1,num_aux,num_eqn,mx,num_ghost,maxnx,rpn3)
            cfl = dmax1(cfl,cfl1d)

            forall(i=1:mx, m=1:num_eqn)
                dq(m,i,j,k) = dq(m,i,j,k)+dq1d(m,i)
            end forall
        enddo
    enddo !end x sweeps


    ! perform y sweeps
    ! ==================

    do i = 0, mx+1
        do k = 0, mz+1
            ! copy data along a slice into 1d arrays:
            q1dp => q(:,i,:,k)
            if (num_aux.gt.0)  then
                auxp => aux(:,i,:,k)
            endif

            call flux1(q1dp,dq1d,auxp,dt,cfl1d,t,2,num_aux,num_eqn,my,num_ghost,maxnx,rpn3)
            cfl = dmax1(cfl,cfl1d)

            forall(j=1:my,m=1:num_eqn)
                dq(m,i,j,k) = dq(m,i,j,k)+dq1d(m,j)
            end forall
        enddo
    enddo !end y sweeps

    ! perform z sweeps
    ! ==================

    do i = 0, mx+1
        do j = 0, my+1
            ! copy data along a slice into 1d arrays:
            q1dp => q(:,i,j,:)
            if (num_aux.gt.0)  then
                auxp => aux(:,i,j,:)
            endif

            call flux1(q1dp,dq1d,auxp,dt,cfl1d,t,3,num_aux,num_eqn,mz,num_ghost,maxnx,rpn3)
            cfl = dmax1(cfl,cfl1d)

            forall(k=1:mz,m=1:num_eqn)
                dq(m,i,j,k) = dq(m,i,j,k)+dq1d(m,k)
            end forall
        enddo
    enddo !end y sweeps

end subroutine flux3
