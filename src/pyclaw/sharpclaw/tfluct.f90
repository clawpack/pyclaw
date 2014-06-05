!
!
! =========================================================
      subroutine tfluct(maxnx,num_eqn,num_waves,num_aux,num_ghost,mx,&
                        ql,qr,auxl,auxr,amdq2)
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
