      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /cqinit/ A1,beta1,x1,y1, A2,beta2,x2,y2
      common /cgrid/ igrid
c
c     # read data values for this problem
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                

c     # These parameters are used in qinit.f
      read(7,*) A1
      read(7,*) beta1
      read(7,*) x1
      read(7,*) y1
      read(7,*) A2
      read(7,*) beta2
      read(7,*) x2
      read(7,*) y2

      close(unit=7)

      return
      end
