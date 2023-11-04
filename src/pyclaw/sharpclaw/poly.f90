module poly
contains

    ! ===================================================================
    subroutine poly4(q,ql,qr,num_eqn,maxnx,num_ghost)
    ! ===================================================================

        implicit none

        integer,          intent(in) :: num_eqn, maxnx, num_ghost
        double precision, intent(in) :: q(num_eqn,maxnx+2*num_ghost)
        double precision, intent(out) :: ql(num_eqn,maxnx+2*num_ghost),qr(num_eqn,maxnx+2*num_ghost)

        integer :: i,n
        
        do n = 1,num_eqn
            do i = num_ghost-1,maxnx+num_ghost+1
                ql(n,i) = (-5.d0*q(n,-2+i)+60.d0*q(n,-1+i)+90.d0*q(n,i)-20.d0*q(n,1+i)+3.d0*q(n,2+i))/128.d0
                qr(n,i) = (3.d0*q(n,-2+i)-20.d0*q(n,-1+i)+90.d0*q(n,i)+60.d0*q(n,1+i)-5.d0*q(n,2+i))/128.d0
            end do
        end do

    end subroutine poly4

    ! ===================================================================
    subroutine poly6(q,ql,qr,num_eqn,maxnx,num_ghost)
    ! ===================================================================

        implicit none

        integer,          intent(in) :: num_eqn, maxnx, num_ghost
        double precision, intent(in) :: q(num_eqn,maxnx+2*num_ghost)
        double precision, intent(out) :: ql(num_eqn,maxnx+2*num_ghost),qr(num_eqn,maxnx+2*num_ghost)

        integer :: i,n
        
        do n = 1,num_eqn
            do i = num_ghost-1,maxnx+num_ghost+1
                ql(n,i) = (7.d0*q(n,-3+i) - 70.d0*q(n,-2+i) + 525.d0*q(n,-1+i) + 700.d0*q(n,i) &
                           - 175.d0*q(n,1+i) + 42.d0*q(n,2+i) - 5.d0*q(n,3+i))/1024.d0
                qr(n,i) = (42.d0*q(n,-2+i) - 175.d0*q(n,-1+i) + 700.d0*q(n,i) &
                           + 525.d0*q(n,1+i) - 70.d0*q(n,2+i) + 7.d0*q(n,3+i) - 5.d0*q(n,-3+i))/1024.d0
            end do
        end do
    
    end subroutine poly6

    ! ===================================================================
    subroutine poly8(q,ql,qr,num_eqn,maxnx,num_ghost)
    ! ===================================================================

        implicit none

        integer,          intent(in) :: num_eqn, maxnx, num_ghost
        double precision, intent(in) :: q(num_eqn,maxnx+2*num_ghost)
        double precision, intent(out) :: ql(num_eqn,maxnx+2*num_ghost),qr(num_eqn,maxnx+2*num_ghost)

        integer :: i,n
        
        do n = 1,num_eqn
            do i = num_ghost-1,maxnx+num_ghost+1
                ql(n,i) = (-45.d0*q(n,-4+i) + 504.d0*q(n,-3+i) - 2940.d0*q(n,-2+i) + 17640.d0*q(n,-1+i) &
                           + 22050.d0*q(n,i) - 5880.d0*q(n,1+i) + 1764.d0*q(n,2+i) - 360.d0*q(n,3+i) + 35.d0*q(n,4+i))/32768.d0
                qr(n,i) = (35.d0*q(n,-4+i) - 360.d0*q(n,-3+i) + 1764.d0*q(n,-2+i) - 5880.d0*q(n,-1+i) + 22050.d0*q(n,i) &
                           + 17640.d0*q(n,1+i) - 2940.d0*q(n,2+i) + 504.d0*q(n,3+i) - 45.d0*q(n,4+i))/32768.d0
            end do
        end do
    
    end subroutine poly8

    ! ===================================================================
    subroutine poly10(q,ql,qr,num_eqn,maxnx,num_ghost)
    ! ===================================================================

        implicit none

        integer,          intent(in) :: num_eqn, maxnx, num_ghost
        double precision, intent(in) :: q(num_eqn,maxnx+2*num_ghost)
        double precision, intent(out) :: ql(num_eqn,maxnx+2*num_ghost),qr(num_eqn,maxnx+2*num_ghost)

        integer :: i,n

        do n = 1,num_eqn
            do i = num_ghost-1,maxnx+num_ghost+1
                ql(n,i) = (77.d0*q(n,-5+i) - 990.d0*q(n,-4+i) + 6237.d0*q(n,-3+i) - 27720.d0*q(n,-2+i) &
                           + 145530.d0*q(n,-1+i) + 174636.d0*q(n,i) - 48510.d0*q(n,1+i) + 16632.d0*q(n,2+i) - 4455.d0*q(n,3+i) &
                           + 770.d0*q(n,4+i) - 63.d0*q(n,5+i))/262144.d0
                qr(n,i) = (-63.d0*q(n,-5+i) + 770.d0*q(n,-4+i) - 4455.d0*q(n,-3+i) + 16632.d0*q(n,-2+i) & 
                           - 48510.d0*q(n,-1+i) + 174636.d0*q(n,i) + 145530.d0*q(n,1+i) - 27720.d0*q(n,2+i) + 6237.d0*q(n,3+i) &
                           - 990.d0*q(n,4+i) + 77.d0*q(n,5+i))/262144.d0
            end do
        end do
    end subroutine poly10

end module poly