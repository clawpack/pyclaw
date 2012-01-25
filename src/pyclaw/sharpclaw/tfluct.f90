!
!
! =========================================================
      subroutine tfluct(maxmx,num_eqn,num_waves,num_ghost,mx,ql,qr,auxl,auxr, &
                 s,adq)
! =========================================================
!
!
      implicit double precision (a-h,o-z)
!
!
      write(*,*) "Error: method(8)=1, but you have not defined"
      write(*,*) "a function tfluct1()."
      stop
!
      return
      end
