! ============================================================================
!  Program:     /Users/mandli/Documents/Research/notes/point_wise_rs/src
!  File:        main
!  Created:     2010-04-20
!  Author:      Kyle Mandli
! ============================================================================
!      Copyright (C) 2010-04-20 Kyle Mandli <mandli@amath.washington.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================
!  Stencils:                                   
!                         |          |          | 
!           ghost cells   |          |          | 
!    ----------+----------+----------+----------+-
!        1          2          3          4         cell centers  q,aux
!           1          2          3          4      rp locations  waves,s
! 
!      |          |          |          
!      |          |          |    ghost cells
!    --+----------+----------+----------+----------
!          mx-3      mx-2       mx-1        mx      cell centers  q,aux
!    mx-3       mx-2      mx-1        mx            rp locations  waves,s
!
!  First row are grid cell center indices
!  Second row are riemann problem indices
!  Riemann problems are solved in the indices at grid cell edges
!
module mod_rand

    implicit none
    
contains

    subroutine init_random_seed()
        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        call system_clock(count=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        call random_seed(put = seed)

        deallocate(seed)
    end subroutine

end module mod_rand

subroutine pw_rp(meqn,maux,mwaves,q_l,q_r,aux_l,aux_r,wave,s)

    implicit none
    integer, intent(in) :: meqn,maux,mwaves
    double precision, intent(in) :: q_l(1:meqn), q_r(1:meqn)
    double precision, intent(in) :: aux_l(1:maux), aux_r(1:maux)
    double precision, intent(out) :: wave(1:meqn,1:mwaves),s(1:mwaves)
    
    integer :: p
    
    do p=1,mwaves
        wave(:,p) = q_l(:) + q_r(:)
        s(p) = (aux_l(p) + aux_r(p)) / 2.d0
    enddo

end subroutine pw_rp

subroutine vec_rp(mx,meqn,maux,mwaves,q_l,q_r,aux_l,aux_r,wave,s)

    implicit none
    
    integer, intent(in) :: mx,meqn,maux,mwaves
    double precision, intent(in) :: q_l(1:meqn,1:mx-1), q_r(1:meqn,1:mx-1)
    double precision, intent(in) :: aux_l(1:maux,1:mx-1), aux_r(1:maux,1:mx-1)
    double precision, intent(out) :: wave(1:meqn,1:mwaves,1:mx),s(1:mwaves,1:mx)

    integer :: i,p
    
    do i=2,mx-1
        do p=1,mwaves
            wave(:,p,i) = q_l(:,i) + q_r(:,i)
            s(p,i) = (aux_l(p,i) + aux_r(p,i)) / 2.d0
        enddo
    enddo
            

end subroutine vec_rp

program clawpack_rs_test
    
    implicit none
    
    ! Parameters
    integer, parameter :: mx = 2**15, meqn = 2, maux = 2, mwaves = 2
    
    ! Arrays
    double precision, allocatable :: q(:,:), aux(:,:), wave(:,:,:), s(:,:)
    
    ! Locals
    integer :: i,n,m,count_rate,start,finish
    double precision :: pw_time, vec_time
    
    ! Allocate all arrays
    allocate(q(1:meqn,1:mx))
    allocate(aux(1:maux,1:mx))
    allocate(wave(1:meqn,1:mwaves,1:mx))
    allocate(s(1:mwaves,1:mx))
    
    ! Initialize arrays
    call random_number(q)
    call random_number(aux)
    wave = 0.d0
    s = 0.d0
    
    ! Point-wise test
    print *,"Point-wise test..."
    call system_clock(start,count_rate)
    do n=1,2**15
        !$OMP parallel do
        do i=2,mx-1
            call pw_rp(meqn,maux,mwaves,q(:,i-1),q(:,i),aux(:,i-1), &
                aux(:,i),wave(:,:,i),s(:,i))
        enddo
        !$OMP end parallel do
    enddo
    call system_clock(finish,count_rate)
    pw_time = float(finish - start) / float(count_rate)
    print *,"test done."
    
    ! Vectorized version
    print *,"Vectorized test..."
    wave = 0.d0
    s = 0.d0
    call system_clock(start,count_rate)
    do n=1,2**15
        call vec_rp(mx,meqn,maux,mwaves,q(:,1:mx-1),q(:,2:mx),aux(:,1:mx-1), &
            aux(:,2:mx),wave,s)
    enddo
    call system_clock(finish,count_rate)
    vec_time = float(finish - start) / float(count_rate)
    print *,"test done."

    print *, "Point-wise: ", pw_time
    print *, "Vectorized: ", vec_time
    print *, "Percent difference: ",100.d0 * abs(vec_time - pw_time) / vec_time,"%"
end program clawpack_rs_test