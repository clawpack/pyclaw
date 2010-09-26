      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /comaux/ rho(1000)
      common /comlim/ mylim,mrplim(2)
      common /combc/ omega
      common /comwall/ pi,t1,a1,tw1,t2,a2,tw2

c
c     # Set the material parameters for the acoustic equations
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)

      iunit = 8
      fname = 'rho.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)

      pi = 4.d0*datan(1.d0)

c     # parameters for wall
      read(7,*) t1,a1,tw1
      read(7,*) t2,a2,tw2

c     # Piecewise constant medium with interfaces at x = 1, 2, ...
c     # Random densities can be created in matlab, e.g.:
c     #       rho = rand(1000,1) + 1;
c     #       save setprob.rho rho -ascii

c     # for periodic:
c     #       rho = mod(1:1000,2)*2 + 1;

      do i=1,1000
         read(8,*) rho(i)
         enddo
c
c     # if mylim>0 then limiter is applied in rp1 rather than using
c     # default clawpack limiter.   
c     # mrplim(p) tells which limiter to use for the p'th wave.
c     # if mylim>0 then you should set mthlim(p)=0 in claw1ez.data!
c
      read(7,*) mylim
      read(7,*) mrplim(1), mrplim(2)

      return
      end
