! ==========================================================
subroutine flux2(q,dq,q1d,dq1d,aux,dt,cfl,t,num_aux,num_eqn,num_ghost,maxnx,mx,my,rpn2)
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

    integer :: num_aux,num_eqn,num_ghost,maxnx,mx,my
    double precision, target, intent(in) :: q(num_eqn, 1-num_ghost:mx+num_ghost, 1-num_ghost:my+num_ghost)
    double precision, intent(inout) :: dq(num_eqn, 1-num_ghost:mx+num_ghost, 1-num_ghost:my+num_ghost)
    double precision, target, intent(in) :: aux(num_aux, 1-num_ghost:mx+num_ghost, 1-num_ghost:my+num_ghost)
    double precision :: q1d(num_eqn,1-num_ghost:maxnx+num_ghost), dq1d(num_eqn,1-num_ghost:maxnx+num_ghost)
    double precision, intent(in) :: dt,t
    double precision, intent(out) :: cfl
    integer :: i,j,m
    double precision :: cfl1d
    double precision, pointer :: auxp(:,:),q1dp(:,:)
    external :: rpn2
    
!f2py intent(in,out) dq  
!f2py intent(out) cfl  
!f2py optional dq, q1d, dq1d

! Dummy interface just so f2py doesn't complain:
!f2py real(DP) x
!f2py x=rpn2(x)

    cfl = 0.d0

    ! perform x-sweeps
    ! ==================

   
    do j = 0,my+1

        ! copy auxiliary data along a slice into 1d arrays:
        q1dp => q(:,:,j)
        if (num_aux .gt. 0)  then
            auxp => aux(:,:,j)
        endif


        ! compute modification dq1d along this slice:
        call flux1(q1dp,dq1d,auxp,dt,cfl1d,t,1,num_aux,num_eqn,mx,num_ghost,maxnx,rpn2)
        cfl = dmax1(cfl,cfl1d)


        if (index_capa.eq.0) then
            ! no capa array.  Standard flux differencing:
            forall(i=1:mx, m=1:num_eqn)
                dq(m,i,j) = dq(m,i,j)+dq1d(m,i)
            end forall
        else
            ! with capa array.  Which is correct?
            forall(i=1:mx, m=1:num_eqn)
                dq(m,i,j) = dq(m,i,j)+dq1d(m,i)
            end forall
            ! dq(m,i,j) = dq(m,i,j)+dq1d(m,i)/aux(index_capa,i,j)
        endif
    enddo !end x sweeps


    ! perform y sweeps
    ! ==================

    do i = 0, mx+1
        ! copy auxiliary data along a slice into 1d arrays:
        q1dp => q(:,i,:)
        forall(j=1-num_ghost:my+num_ghost, m=1:num_eqn)
            q1d(m,j) = q(m,i,j)
        end forall

        if (num_aux .gt. 0)  then
            auxp => aux(:,i,:)
        endif

        call flux1(q1dp,dq1d,auxp,dt,cfl1d,t,2,num_aux,num_eqn,my,num_ghost,maxnx,rpn2)
        cfl = dmax1(cfl,cfl1d)

        if (index_capa.eq.0) then
            ! no capa array.  Standard flux differencing:
            forall(j=1:my,m=1:num_eqn)
                dq(m,i,j) = dq(m,i,j)+dq1d(m,j)
            end forall
        else
            ! with capa array.  Which is correct?
            ! dq(m,i,j) = dq(m,i,j)+dq1d(m,j)/aux(index_capa,i,j)
            dq(:,i,:) = dq(:,i,:)+dq1d
        endif
    enddo !end y sweeps

end subroutine flux2
