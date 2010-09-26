      subroutine setprob
      implicit double precision (a-h,o-z)
      character*25 fname
      common /cparam/ rho,bulk,cc,zz

c
c     # Set the material parameters for the acoustic equations
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                
c
c     # Density and bulk modulus:

      read(7,*) rho
      read(7,*) bulk
c
c     # Compute sound speed and impendance:

      cc = dsqrt(bulk/rho)
      zz = rho*cc

      return
      end
