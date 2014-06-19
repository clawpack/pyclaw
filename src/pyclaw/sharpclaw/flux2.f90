! ==========================================================
subroutine flux2(q,dq,aux,dt,cfl,t,num_aux,num_eqn,num_ghost,maxnx,mx,my,rpn2,tfluct2)
! ==========================================================
!    Evaluate (delta t) *dq/dt
!    On entry, q should give the initial data for this step
!    See the flux1 documentation for more information.
! ==========================================================

    use ClawParams
    implicit none
    
    ! Input (dummy) variables
    external :: rpn2, tfluct2
    integer, intent(in) :: num_aux, num_eqn, num_ghost, maxnx, mx, my
    double precision, intent(in) :: dt,t
    double precision, target, intent(in) :: q(num_eqn, 1-num_ghost:mx+num_ghost, 1-num_ghost:my+num_ghost)
    double precision, target, intent(in) :: aux(num_aux, 1-num_ghost:mx+num_ghost, 1-num_ghost:my+num_ghost)
    double precision, intent(inout) :: dq(num_eqn, 1-num_ghost:mx+num_ghost, 1-num_ghost:my+num_ghost)
    double precision, intent(out) :: cfl

!f2py intent(in,out) dq  
!f2py intent(out) cfl  
!f2py optional dq

! Dummy interface just so f2py doesn't complain:
!f2py real(DP) x
!f2py x=rpn2(x)
!f2py x=tfluct2(x)

    ! Local variables
    integer :: i,j,m
    double precision :: cfl1d
    double precision, pointer :: aux1d(:,:),q1d(:,:)
    double precision :: dq1d(num_eqn,1-num_ghost:maxnx+num_ghost)

! ==========================================================

    cfl = 0.d0

    ! perform x-sweeps
    ! ==================

    do j = 0,my+1
        ! copy data along a slice into 1d arrays:
        q1d => q(:,:,j)
        aux1d => aux(:,:,j)

        ! compute modification dq1d along this slice:
        call flux1(q1d,dq1d,aux1d,dt,cfl1d,t,1,num_aux,num_eqn,mx,num_ghost,maxnx,rpn2,tfluct2)
        cfl = dmax1(cfl,cfl1d)

        forall(i=1:mx, m=1:num_eqn)
            dq(m,i,j) = dq(m,i,j)+dq1d(m,i)
        end forall
    enddo !end x sweeps


    ! perform y sweeps
    ! ==================

    do i = 0, mx+1
        ! copy data along a slice into 1d arrays:
        q1d => q(:,i,:)
        aux1d => aux(:,i,:)

        call flux1(q1d,dq1d,aux1d,dt,cfl1d,t,2,num_aux,num_eqn,my,num_ghost,maxnx,rpn2,tfluct2)
        cfl = dmax1(cfl,cfl1d)

        forall(j=1:my, m=1:num_eqn)
            dq(m,i,j) = dq(m,i,j)+dq1d(m,j)
        end forall
    enddo !end y sweeps

end subroutine flux2
