      program testrpn2
c
c     # Routine to test a Riemann solver rpn2 of the form used in 2d clawpack.
c     # You must first modify this program to set meqn and mwaves 
c     # appropriately for the solver being tested, and set
c     # any values that must be passed to rp in common blocks.
c
c     # This version is set to test 2d shallow water solvers.
c     # To use:
c          f77 testrpn2.f $CLAW/applications/shallow/2d/rp/rpn2sw.f
c          a.out
c
c     # Author: Randall J. LeVeque
c
c
      implicit double precision (a-h,o-z)
c
c     # The next two parameters need to be set appropriately for other rpn2:
      parameter (meqn = 3)
      parameter (mwaves = 3)

      parameter (maxmx = 2)
      parameter (mx = 2)
      parameter (mbc = 0)
      parameter (maux = 0)
c
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension wave(1-mbc:maxmx+mbc, meqn,mwaves)
      dimension s(1-mbc:maxmx+mbc, mwaves)
      dimension amdq(1-mbc:maxmx+mbc, meqn)
      dimension apdq(1-mbc:maxmx+mbc, meqn)

c     dimension  aux(1-mbc:maxmx+mbc, maux)

c     # gravitational constant needed in rp1sw.f:
      common /param/ grav
      grav = 1.d0

      write(6,*) ' '
      write(6,*) 'Riemann solver test routine for rp1'
      write(6,605) meqn,mwaves
  605 format('   meqn = ',i3,'    mwaves = ',i3)
      write(6,*) ' '


c     # read in data:
      write(6,*) 'Input the meqn components of ql'
      read(5,*) (q(1,m),m=1,meqn)
      write(6,*) 'Input the meqn components of qr'
      read(5,*) (q(2,m),m=1,meqn)


      write(6,600)
      do m=1,meqn
         write(6,604) q(1,m), q(2,m)
         enddo

      call rpn2(1,maxmx,meqn,mwaves,mbc,mx,q,q,aux,aux,
     &         wave,s,amdq,apdq)

      do mw=1,mwaves
         write(6,601) mw, s(2,mw)
         do m=1,meqn
            q(1,m) = q(1,m) + wave(2,m,mw)
            write(6,604) wave(2,m,mw),q(1,m)
            enddo
         enddo

      write(6,603)
      do m=1,meqn
         write(6,604) amdq(2,m), apdq(2,m)
         enddo
      write(6,*) ' '

  600 format(/,'           ql                    qr')
  601 format(/,i3,'-wave has speed s = ',d19.12,/,
     &   /,'           wave         q to right of this wave:')
  603 format(/,'           amdq                   apdq')
  604 format(2d22.12)

      stop 
      end
